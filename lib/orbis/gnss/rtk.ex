defmodule Orbis.GNSS.RTK do
  @moduledoc """
  RTK-facing carrier/code double-difference primitives.

  A base receiver and a rover receiver observing the same satellites have
  receiver-clock terms that differ by station but are common to every satellite.
  A *single difference* subtracts base from rover for the same satellite; a
  *double difference* subtracts a reference satellite's single difference:

      DD_s = (rover_s - base_s) - (rover_ref - base_ref)

  The receiver clocks cancel in the second subtraction. Satellite-clock,
  ephemeris, and short-baseline atmosphere errors that are common between base
  and rover also cancel in the receiver single difference. The remaining
  carrier-phase double differences are the measurement surface used by RTK
  baseline estimation and integer ambiguity fixing.

  `double_differences/3` returns the normalized measurements. The first solver
  layer, `solve_float_baseline_epochs/3`, estimates one static base-to-rover
  baseline from code and carrier-phase double differences, keeping one float
  carrier ambiguity per non-reference satellite across the arc.

  ## Example

      iex> base = [
      ...>   {"G01", 20_100.0, 20_103.0},
      ...>   {"G02", 21_105.0, 21_110.0}
      ...> ]
      iex> rover = [
      ...>   {"G01", 20_040.0, 20_044.0},
      ...>   {"G02", 21_060.0, 21_066.0}
      ...> ]
      iex> {:ok, result} = Orbis.GNSS.RTK.double_differences(base, rover, reference_satellite_id: "G01")
      iex> result.double_differences
      [%{satellite_id: "G02", reference_satellite_id: "G01", code_m: 15.0, phase_m: 17.0}]
  """

  import Orbis.GNSS.Core.LinearAlgebra,
    only: [
      correlated_normal_equations: 2,
      invert_matrix: 1,
      matmul: 2,
      matrix_sub: 2,
      solve_linear: 2,
      solve_matrix: 2,
      submatrix: 5,
      transpose: 1
    ]

  alias Orbis.GNSS.Core.Types

  @default_max_iterations 8
  @default_position_tolerance_m 1.0e-4
  @default_ambiguity_tolerance_m 1.0e-4
  @default_code_sigma_m 1.0
  @default_phase_sigma_m 0.02

  defmodule FloatBaselineSolution do
    @moduledoc """
    Float RTK baseline solution from code/carrier double differences.
    """

    @enforce_keys [
      :baseline_m,
      :rover_position_m,
      :reference_satellite_id,
      :used_sats,
      :ambiguities_m,
      :residuals_m,
      :metadata
    ]
    defstruct [
      :baseline_m,
      :rover_position_m,
      :reference_satellite_id,
      :used_sats,
      :ambiguities_m,
      :residuals_m,
      :metadata
    ]

    @type ecef :: %{x_m: float(), y_m: float(), z_m: float()}

    @type residual :: %{
            epoch: term(),
            satellite_id: String.t(),
            reference_satellite_id: String.t(),
            code_m: float(),
            phase_m: float()
          }

    @type t :: %__MODULE__{
            baseline_m: ecef(),
            rover_position_m: ecef(),
            reference_satellite_id: String.t(),
            used_sats: [String.t()],
            ambiguities_m: %{String.t() => float()},
            residuals_m: [residual()],
            metadata: %{
              iterations: pos_integer(),
              converged: boolean(),
              status: :state_tolerance | :max_iterations,
              ambiguity_float: %{
                order: [String.t()],
                covariance_m: [[float()]],
                covariance_inverse_m: [[float()]]
              },
              measurement_covariance: %{
                model: :double_difference,
                code_sigma_m: float(),
                phase_sigma_m: float()
              },
              code_rms_m: float(),
              phase_rms_m: float(),
              weighted_rms_m: float(),
              n_epochs: pos_integer(),
              n_observations: pos_integer(),
              dropped_sats: [String.t()]
            }
          }
  end

  @typedoc "Code and carrier-phase observation in metres."
  @type observation ::
          %{satellite_id: String.t(), code_m: number(), phase_m: number()}
          | {String.t(), number(), number()}

  @typedoc "ECEF position in metres."
  @type ecef_input ::
          {number(), number(), number()} | %{x_m: number(), y_m: number(), z_m: number()}

  @typedoc "Satellite ECEF position keyed by satellite id."
  @type satellite_positions :: %{required(String.t()) => ecef_input()}

  @typedoc """
  One RTK epoch carrying paired base/rover observations and satellite positions.

  `:epoch` is preserved in residual diagnostics; it is not interpreted by this
  first solver layer because satellite positions are supplied by the caller.
  """
  @type baseline_epoch :: %{
          required(:base_observations) => [observation()],
          required(:rover_observations) => [observation()],
          required(:satellite_positions_m) => satellite_positions(),
          optional(:epoch) => term()
        }

  @typedoc "One non-reference satellite's double-difference observation."
  @type double_difference :: %{
          satellite_id: String.t(),
          reference_satellite_id: String.t(),
          code_m: float(),
          phase_m: float()
        }

  @typedoc "Double-difference result with deterministic satellite ordering."
  @type result :: %{
          reference_satellite_id: String.t(),
          double_differences: [double_difference()],
          dropped_sats: [String.t()]
        }

  @doc """
  Solve a static float RTK baseline from multi-epoch double differences.

  `base_position` is the surveyed base ECEF position. Each epoch supplies base
  and rover code/carrier observations plus satellite ECEF positions at that
  receive epoch:

      epoch = %{
        epoch: ~N[2026-01-01 00:00:00],
        satellite_positions_m: %{"G01" => {21.0e6, 14.0e6, 20.0e6}, ...},
        base_observations: [%{satellite_id: "G01", code_m: p_base, phase_m: l_base}, ...],
        rover_observations: [%{satellite_id: "G01", code_m: p_rover, phase_m: l_rover}, ...]
      }

      {:ok, sol} = Orbis.GNSS.RTK.solve_float_baseline_epochs(base_position, [epoch])

  The model is

      DD = [rho_rover(s) - rho_base(s)] - [rho_rover(ref) - rho_base(ref)]

  for code, and the same geometry plus one float carrier ambiguity per
  non-reference satellite for phase. Receiver clocks and any satellite-common
  short-baseline errors cancel before the solve.

  The normal equations use the full double-difference covariance block for each
  epoch and measurement kind. Double-difference rows sharing the reference
  satellite are therefore correlated, and the returned
  `metadata.ambiguity_float` contains the resulting float ambiguity covariance
  and inverse covariance in metres.

  Options:

    * `:reference_satellite_id` - fixed reference satellite. When omitted, the
      highest-average-elevation satellite common to every epoch is used, with a
      lexicographic tie-break.
    * `:initial_baseline_m` - initial base-to-rover ECEF vector, default
      `{0.0, 0.0, 0.0}`.
    * `:code_sigma_m` / `:phase_sigma_m` - undifferenced receiver measurement
      sigmas in metres. The solver propagates them into the non-diagonal
      double-difference covariance where rows sharing the reference satellite
      are correlated. Defaults are `#{@default_code_sigma_m}` and
      `#{@default_phase_sigma_m}`.
    * `:max_iterations`, `:position_tolerance_m`,
      `:ambiguity_tolerance_m`.

  Returns `{:ok, %FloatBaselineSolution{}}` or a tagged error. The solver is a
  float-ambiguity RTK layer; integer double-difference fixing is a later step.
  """
  @spec solve_float_baseline_epochs(ecef_input(), [baseline_epoch()], keyword()) ::
          {:ok, FloatBaselineSolution.t()} | {:error, term()}
  def solve_float_baseline_epochs(base_position, epochs, opts \\ [])

  def solve_float_baseline_epochs(base_position, epochs, opts) when is_list(epochs) do
    with {:ok, base} <- Types.normalize_ecef(base_position, :invalid_base_position),
         :ok <- ensure_nonempty_epochs(epochs),
         {:ok, normalized_epochs} <- normalize_epochs(epochs),
         {:ok, common_sats, dropped_sats} <- common_epoch_sats(normalized_epochs),
         :ok <- ensure_baseline_satellites(common_sats),
         {:ok, reference_sat} <-
           baseline_reference_satellite(common_sats, opts, base, normalized_epochs),
         {:ok, solve_opts} <- baseline_solve_options(opts),
         {:ok, weights} <- baseline_weights(opts),
         {:ok, initial_baseline} <- initial_baseline(opts) do
      ambiguity_sats = Enum.reject(common_sats, &(&1 == reference_sat))

      if baseline_row_count(normalized_epochs, ambiguity_sats) <
           baseline_unknown_count(ambiguity_sats) do
        {:error,
         {:underdetermined, baseline_row_count(normalized_epochs, ambiguity_sats),
          baseline_unknown_count(ambiguity_sats)}}
      else
        state = %{
          baseline: initial_baseline,
          ambiguities: Map.new(ambiguity_sats, &{&1, 0.0})
        }

        iterate_baseline(
          base,
          normalized_epochs,
          reference_sat,
          ambiguity_sats,
          state,
          weights,
          solve_opts,
          dropped_sats,
          1
        )
      end
    end
  end

  def solve_float_baseline_epochs(_base_position, _epochs, _opts), do: {:error, :invalid_epochs}

  @doc """
  Build code and carrier-phase double differences from base and rover observations.

  Observations can be maps with `:satellite_id`, `:code_m`, and `:phase_m`, or
  `{satellite_id, code_m, phase_m}` tuples. Satellites are paired by id; any
  satellite not present at both receivers is reported in `:dropped_sats`.

  Options:

    * `:reference_satellite_id` - reference satellite for the second
      difference. When omitted, the lexicographically first common satellite is
      selected deterministically.

  Returns `{:ok, result}` or a tagged error. At least two common satellites are
  required so one non-reference double difference can be produced.
  """
  @spec double_differences([observation()], [observation()], keyword()) ::
          {:ok, result()} | {:error, term()}
  def double_differences(base_observations, rover_observations, opts \\ [])

  def double_differences(base_observations, rover_observations, opts)
      when is_list(base_observations) and is_list(rover_observations) do
    with {:ok, base} <- normalize_observations(base_observations, :invalid_base_observations),
         {:ok, rover} <- normalize_observations(rover_observations, :invalid_rover_observations),
         {:ok, common, dropped} <- common_observations(base, rover),
         {:ok, reference_sat} <- reference_satellite(common, opts) do
      ref_base = Map.fetch!(base, reference_sat)
      ref_rover = Map.fetch!(rover, reference_sat)
      ref_code_sd = ref_rover.code_m - ref_base.code_m
      ref_phase_sd = ref_rover.phase_m - ref_base.phase_m

      dds =
        common
        |> Enum.reject(&(&1 == reference_sat))
        |> Enum.map(fn sat ->
          base_obs = Map.fetch!(base, sat)
          rover_obs = Map.fetch!(rover, sat)

          %{
            satellite_id: sat,
            reference_satellite_id: reference_sat,
            code_m: rover_obs.code_m - base_obs.code_m - ref_code_sd,
            phase_m: rover_obs.phase_m - base_obs.phase_m - ref_phase_sd
          }
        end)

      {:ok,
       %{
         reference_satellite_id: reference_sat,
         double_differences: dds,
         dropped_sats: dropped
       }}
    end
  end

  def double_differences(_base_observations, _rover_observations, _opts),
    do: {:error, :invalid_observations}

  defp ensure_nonempty_epochs([]), do: {:error, :no_epochs}
  defp ensure_nonempty_epochs(_epochs), do: :ok

  defp normalize_epochs(epochs) do
    epochs
    |> Enum.with_index()
    |> Enum.reduce_while({:ok, []}, fn {epoch, idx}, {:ok, acc} ->
      case normalize_epoch(epoch, idx) do
        {:ok, normalized} -> {:cont, {:ok, [normalized | acc]}}
        {:error, _reason} = err -> {:halt, err}
      end
    end)
    |> case do
      {:ok, normalized} -> {:ok, Enum.reverse(normalized)}
      {:error, _reason} = err -> err
    end
  end

  defp normalize_epoch(
         %{
           base_observations: base_observations,
           rover_observations: rover_observations,
           satellite_positions_m: satellite_positions
         } = epoch,
         idx
       )
       when is_list(base_observations) and is_list(rover_observations) and
              is_map(satellite_positions) do
    with {:ok, base} <- normalize_observations(base_observations, :invalid_base_observations),
         {:ok, rover} <- normalize_observations(rover_observations, :invalid_rover_observations),
         {:ok, positions} <- normalize_satellite_positions(satellite_positions) do
      {:ok,
       %{
         idx: idx,
         epoch: Map.get(epoch, :epoch, idx),
         base: base,
         rover: rover,
         positions: positions
       }}
    end
  end

  defp normalize_epoch(_epoch, idx), do: {:error, {:invalid_epoch_observations, idx}}

  defp normalize_satellite_positions(positions) do
    positions
    |> Enum.reduce_while({:ok, %{}}, fn
      {sat, position}, {:ok, acc} when is_binary(sat) ->
        case Types.normalize_ecef(position, :invalid_satellite_position) do
          {:ok, ecef} -> {:cont, {:ok, Map.put(acc, sat, ecef)}}
          {:error, _reason} -> {:halt, {:error, {:invalid_satellite_position, sat}}}
        end

      {sat, _position}, {:ok, _acc} ->
        {:halt, {:error, {:invalid_satellite_position, sat}}}
    end)
  end

  defp common_epoch_sats(epochs) do
    per_epoch =
      Enum.map(epochs, fn epoch ->
        epoch.base
        |> Map.keys()
        |> MapSet.new()
        |> MapSet.intersection(epoch.rover |> Map.keys() |> MapSet.new())
        |> MapSet.intersection(epoch.positions |> Map.keys() |> MapSet.new())
      end)

    common =
      per_epoch
      |> Enum.reduce(fn sats, acc -> MapSet.intersection(acc, sats) end)
      |> MapSet.to_list()
      |> Enum.sort()

    all =
      epochs
      |> Enum.reduce(MapSet.new(), fn epoch, acc ->
        epoch.base
        |> Map.keys()
        |> MapSet.new()
        |> MapSet.union(epoch.rover |> Map.keys() |> MapSet.new())
        |> MapSet.union(epoch.positions |> Map.keys() |> MapSet.new())
        |> MapSet.union(acc)
      end)

    dropped =
      all
      |> MapSet.difference(MapSet.new(common))
      |> MapSet.to_list()
      |> Enum.sort()

    {:ok, common, dropped}
  end

  defp ensure_baseline_satellites(common_sats) do
    if length(common_sats) < 4,
      do: {:error, {:too_few_common_satellites, length(common_sats), 4}},
      else: :ok
  end

  defp baseline_solve_options(opts) do
    max_iterations = Keyword.get(opts, :max_iterations, @default_max_iterations)
    position_tolerance_m = Keyword.get(opts, :position_tolerance_m, @default_position_tolerance_m)

    ambiguity_tolerance_m =
      Keyword.get(opts, :ambiguity_tolerance_m, @default_ambiguity_tolerance_m)

    cond do
      not (is_integer(max_iterations) and max_iterations > 0) ->
        {:error, {:invalid_option, :max_iterations}}

      not (is_number(position_tolerance_m) and position_tolerance_m > 0.0) ->
        {:error, {:invalid_option, :position_tolerance_m}}

      not (is_number(ambiguity_tolerance_m) and ambiguity_tolerance_m > 0.0) ->
        {:error, {:invalid_option, :ambiguity_tolerance_m}}

      true ->
        {:ok,
         %{
           max_iterations: max_iterations,
           position_tolerance_m: position_tolerance_m / 1.0,
           ambiguity_tolerance_m: ambiguity_tolerance_m / 1.0
         }}
    end
  end

  defp baseline_weights(opts) do
    with {:ok, code_sigma_m} <- measurement_sigma(opts, :code_sigma_m, @default_code_sigma_m),
         {:ok, phase_sigma_m} <- measurement_sigma(opts, :phase_sigma_m, @default_phase_sigma_m) do
      {:ok,
       %{
         code_sigma_m: code_sigma_m,
         phase_sigma_m: phase_sigma_m,
         code_dd_weight: 1.0 / (2.0 * code_sigma_m),
         phase_dd_weight: 1.0 / (2.0 * phase_sigma_m)
       }}
    end
  end

  defp measurement_sigma(opts, key, default) do
    sigma = Keyword.get(opts, key, default)

    if is_number(sigma) and sigma > 0.0,
      do: {:ok, sigma / 1.0},
      else: {:error, {:invalid_sigma, key}}
  end

  defp initial_baseline(opts) do
    opts
    |> Keyword.get(:initial_baseline_m, {0.0, 0.0, 0.0})
    |> Types.normalize_ecef(:invalid_initial_baseline)
  end

  defp baseline_row_count(epochs, ambiguity_sats), do: 2 * length(epochs) * length(ambiguity_sats)

  defp baseline_unknown_count(ambiguity_sats), do: 3 + length(ambiguity_sats)

  defp iterate_baseline(
         base,
         epochs,
         reference_sat,
         ambiguity_sats,
         state,
         weights,
         opts,
         dropped_sats,
         iter
       ) do
    with {:ok, rows} <-
           build_baseline_rows(base, epochs, reference_sat, ambiguity_sats, state, weights),
         {:ok, dx} <-
           solve_baseline_normal_equations(rows, baseline_unknown_count(ambiguity_sats)) do
      next = apply_baseline_delta(state, ambiguity_sats, dx)
      {baseline_step, ambiguity_step} = baseline_step_norms(dx)

      cond do
        baseline_step <= opts.position_tolerance_m and
            ambiguity_step <= opts.ambiguity_tolerance_m ->
          finalize_baseline(
            base,
            epochs,
            reference_sat,
            ambiguity_sats,
            next,
            weights,
            dropped_sats,
            iter,
            true,
            :state_tolerance
          )

        iter >= opts.max_iterations ->
          finalize_baseline(
            base,
            epochs,
            reference_sat,
            ambiguity_sats,
            next,
            weights,
            dropped_sats,
            iter,
            false,
            :max_iterations
          )

        true ->
          iterate_baseline(
            base,
            epochs,
            reference_sat,
            ambiguity_sats,
            next,
            weights,
            opts,
            dropped_sats,
            iter + 1
          )
      end
    end
  end

  defp build_baseline_rows(base, epochs, reference_sat, ambiguity_sats, state, weights) do
    epochs
    |> Enum.reduce_while({:ok, []}, fn epoch, {:ok, acc} ->
      case build_epoch_baseline_rows(base, epoch, reference_sat, ambiguity_sats, state, weights) do
        {:ok, rows} -> {:cont, {:ok, rows ++ acc}}
        {:error, _reason} = err -> {:halt, err}
      end
    end)
    |> case do
      {:ok, rows} -> {:ok, Enum.reverse(rows)}
      {:error, _reason} = err -> err
    end
  end

  defp build_epoch_baseline_rows(base, epoch, reference_sat, ambiguity_sats, state, weights) do
    ref_dd = single_difference(epoch, reference_sat)
    ref_pos = Map.fetch!(epoch.positions, reference_sat)

    ambiguity_sats
    |> Enum.reduce({:ok, []}, fn sat, {:ok, acc} ->
      obs_dd = double_difference_measurement(epoch, sat, ref_dd)
      sat_pos = Map.fetch!(epoch.positions, sat)
      {geom_dd, deriv} = geometry_double_difference(base, state.baseline, sat_pos, ref_pos)
      ambiguity = Map.fetch!(state.ambiguities, sat)
      h_base = design_baseline_row(deriv, nil, ambiguity_sats)
      h_phase = design_baseline_row(deriv, sat, ambiguity_sats)

      code_row = %{
        kind: :code,
        epoch_idx: epoch.idx,
        epoch: epoch.epoch,
        sat: sat,
        h: h_base,
        y: obs_dd.code_m - geom_dd,
        sigma_m: weights.code_sigma_m,
        weight: weights.code_dd_weight
      }

      phase_row = %{
        kind: :phase,
        epoch_idx: epoch.idx,
        epoch: epoch.epoch,
        sat: sat,
        h: h_phase,
        y: obs_dd.phase_m - (geom_dd + ambiguity),
        sigma_m: weights.phase_sigma_m,
        weight: weights.phase_dd_weight
      }

      {:ok, [phase_row, code_row | acc]}
    end)
    |> case do
      {:ok, rows} -> {:ok, Enum.reverse(rows)}
      {:error, _reason} = err -> err
    end
  end

  defp single_difference(epoch, sat) do
    %{
      code_m: Map.fetch!(epoch.rover, sat).code_m - Map.fetch!(epoch.base, sat).code_m,
      phase_m: Map.fetch!(epoch.rover, sat).phase_m - Map.fetch!(epoch.base, sat).phase_m
    }
  end

  defp double_difference_measurement(epoch, sat, ref_dd) do
    sd = single_difference(epoch, sat)
    %{code_m: sd.code_m - ref_dd.code_m, phase_m: sd.phase_m - ref_dd.phase_m}
  end

  defp geometry_double_difference(base, baseline, sat_pos, ref_pos) do
    rover = add3(base, baseline)
    sat_sd = range(sat_pos, rover) - range(sat_pos, base)
    ref_sd = range(ref_pos, rover) - range(ref_pos, base)
    sat_deriv = range_derivative(rover, sat_pos)
    ref_deriv = range_derivative(rover, ref_pos)
    {sat_sd - ref_sd, sub3(sat_deriv, ref_deriv)}
  end

  defp design_baseline_row({dx, dy, dz}, ambiguity_sat, ambiguity_sats) do
    ambiguity_cols = Enum.map(ambiguity_sats, &if(&1 == ambiguity_sat, do: 1.0, else: 0.0))
    [dx, dy, dz | ambiguity_cols]
  end

  defp apply_baseline_delta(state, ambiguity_sats, dx) do
    {bx, by, bz} = state.baseline

    ambiguities =
      ambiguity_sats
      |> Enum.with_index()
      |> Map.new(fn {sat, idx} ->
        {sat, Map.fetch!(state.ambiguities, sat) + Enum.at(dx, 3 + idx)}
      end)

    %{
      state
      | baseline: {bx + Enum.at(dx, 0), by + Enum.at(dx, 1), bz + Enum.at(dx, 2)},
        ambiguities: ambiguities
    }
  end

  defp baseline_step_norms(dx) do
    baseline_step =
      dx
      |> Enum.take(3)
      |> Enum.map(&(&1 * &1))
      |> Enum.sum()
      |> :math.sqrt()

    ambiguity_step =
      dx
      |> Enum.drop(3)
      |> Enum.map(&abs/1)
      |> Enum.max(fn -> 0.0 end)

    {baseline_step, ambiguity_step}
  end

  defp finalize_baseline(
         base,
         epochs,
         reference_sat,
         ambiguity_sats,
         state,
         weights,
         dropped_sats,
         iterations,
         converged?,
         status
       ) do
    with {:ok, rows} <-
           build_baseline_rows(base, epochs, reference_sat, ambiguity_sats, state, weights),
         {:ok, covariance_m} <-
           baseline_ambiguity_covariance(
             rows,
             baseline_unknown_count(ambiguity_sats),
             3,
             length(ambiguity_sats)
           ),
         {:ok, covariance_inverse_m} <- invert_matrix(covariance_m) do
      residuals = baseline_residuals(rows, reference_sat)
      code_residuals = Enum.map(residuals, & &1.code_m)
      phase_residuals = Enum.map(residuals, & &1.phase_m)
      rover = add3(base, state.baseline)

      {:ok,
       %FloatBaselineSolution{
         baseline_m: ecef_map(state.baseline),
         rover_position_m: ecef_map(rover),
         reference_satellite_id: reference_sat,
         used_sats: ambiguity_sats,
         ambiguities_m: state.ambiguities,
         residuals_m: residuals,
         metadata: %{
           iterations: iterations,
           converged: converged?,
           status: status,
           ambiguity_float: %{
             order: ambiguity_sats,
             covariance_m: covariance_m,
             covariance_inverse_m: covariance_inverse_m
           },
           measurement_covariance: %{
             model: :double_difference,
             code_sigma_m: weights.code_sigma_m,
             phase_sigma_m: weights.phase_sigma_m
           },
           code_rms_m: rms(code_residuals),
           phase_rms_m: rms(phase_residuals),
           weighted_rms_m: weighted_rms(rows),
           n_epochs: length(epochs),
           n_observations: length(rows),
           dropped_sats: dropped_sats
         }
       }}
    end
  end

  defp baseline_residuals(rows, reference_sat) do
    rows
    |> Enum.group_by(&{&1.epoch, &1.sat})
    |> Enum.map(fn {{epoch, sat}, grouped} ->
      code = Enum.find(grouped, &(&1.kind == :code))
      phase = Enum.find(grouped, &(&1.kind == :phase))

      %{
        epoch: epoch,
        satellite_id: sat,
        reference_satellite_id: reference_sat,
        code_m: code.y,
        phase_m: phase.y
      }
    end)
    |> Enum.sort_by(&{inspect(&1.epoch), &1.satellite_id})
  end

  defp solve_baseline_normal_equations(rows, n) do
    {ata, aty} = baseline_normal_equations(rows, n)
    solve_linear(ata, aty)
  end

  defp baseline_normal_equations(rows, n) do
    rows
    |> baseline_covariance_blocks()
    |> correlated_normal_equations(n)
  end

  defp baseline_covariance_blocks(rows) do
    rows
    |> Enum.group_by(&{&1.epoch_idx, &1.kind})
    |> Enum.sort_by(fn {{epoch_idx, kind}, _rows} -> {epoch_idx, kind} end)
    |> Enum.map(fn {_key, block_rows} ->
      rows = Enum.sort_by(block_rows, & &1.sat)
      sigma_m = rows |> hd() |> Map.fetch!(:sigma_m)

      %{
        rows: rows,
        inverse_covariance: double_difference_inverse_covariance(length(rows), sigma_m)
      }
    end)
  end

  # For one epoch and one measurement kind, DD_i = SD_i - SD_ref. With
  # independent receiver observations of sigma^2 at each station, each single
  # difference has variance 2*sigma^2. Therefore the double-difference
  # covariance is v*(I + J), where v = 2*sigma^2: diagonal 4*sigma^2 and
  # off-diagonal 2*sigma^2 because every DD shares the same reference SD.
  defp double_difference_inverse_covariance(m, sigma_m) do
    sd_variance = 2.0 * sigma_m * sigma_m
    diagonal_scale = 1.0 / sd_variance * (1.0 - 1.0 / (m + 1.0))
    off_diagonal = -1.0 / (sd_variance * (m + 1.0))

    for i <- 0..(m - 1) do
      for j <- 0..(m - 1), do: if(i == j, do: diagonal_scale, else: off_diagonal)
    end
  end

  defp baseline_ambiguity_covariance(rows, n, start, n_ambiguities) do
    {normal, _rhs} = baseline_normal_equations(rows, n)
    ambiguity_covariance_from_normal(normal, start, n_ambiguities)
  end

  defp ambiguity_covariance_from_normal(normal, start, n_ambiguities) do
    a = submatrix(normal, 0, start, 0, start)
    b = submatrix(normal, 0, start, start, n_ambiguities)
    c = submatrix(normal, start, n_ambiguities, start, n_ambiguities)

    with {:ok, a_inv_b} <- solve_matrix(a, b) do
      bt_a_inv_b = matmul(transpose(b), a_inv_b)
      schur = matrix_sub(c, bt_a_inv_b)

      case invert_matrix(schur) do
        {:ok, covariance} -> {:ok, covariance}
        {:error, :singular_geometry} = err -> err
      end
    end
  end

  defp weighted_rms(rows) do
    values = Enum.map(rows, &(&1.y * &1.weight))
    rms(values)
  end

  defp rms([]), do: 0.0

  defp rms(values) do
    :math.sqrt(Enum.reduce(values, 0.0, fn v, acc -> acc + v * v end) / length(values))
  end

  defp range(sat_pos, receiver), do: norm(sub3(sat_pos, receiver))

  defp range_derivative(receiver, sat_pos) do
    rho = range(sat_pos, receiver)
    scale3(sub3(receiver, sat_pos), 1.0 / rho)
  end

  defp add3({ax, ay, az}, {bx, by, bz}), do: {ax + bx, ay + by, az + bz}
  defp sub3({ax, ay, az}, {bx, by, bz}), do: {ax - bx, ay - by, az - bz}
  defp scale3({x, y, z}, s), do: {x * s, y * s, z * s}
  defp norm({x, y, z}), do: :math.sqrt(x * x + y * y + z * z)
  defp ecef_map({x, y, z}), do: %{x_m: x, y_m: y, z_m: z}

  defp normalize_observations(observations, error_tag) do
    observations
    |> Enum.reduce_while({:ok, %{}}, fn observation, {:ok, acc} ->
      case normalize_observation(observation) do
        {:ok, %{satellite_id: sat} = obs} ->
          if Map.has_key?(acc, sat) do
            {:halt, {:error, {:duplicate_observation, sat}}}
          else
            {:cont, {:ok, Map.put(acc, sat, obs)}}
          end

        {:error, _} ->
          {:halt, {:error, {error_tag, observation}}}
      end
    end)
  end

  defp normalize_observation(%{satellite_id: sat, code_m: code, phase_m: phase})
       when is_binary(sat) and is_number(code) and is_number(phase) do
    {:ok, %{satellite_id: sat, code_m: code / 1.0, phase_m: phase / 1.0}}
  end

  defp normalize_observation({sat, code, phase})
       when is_binary(sat) and is_number(code) and is_number(phase) do
    {:ok, %{satellite_id: sat, code_m: code / 1.0, phase_m: phase / 1.0}}
  end

  defp normalize_observation(_observation), do: {:error, :invalid_observation}

  defp common_observations(base, rover) do
    base_sats = base |> Map.keys() |> MapSet.new()
    rover_sats = rover |> Map.keys() |> MapSet.new()

    common =
      base_sats
      |> MapSet.intersection(rover_sats)
      |> MapSet.to_list()
      |> Enum.sort()

    if length(common) < 2 do
      {:error, {:too_few_common_satellites, length(common), 2}}
    else
      dropped =
        base_sats
        |> MapSet.union(rover_sats)
        |> MapSet.difference(MapSet.new(common))
        |> MapSet.to_list()
        |> Enum.sort()

      {:ok, common, dropped}
    end
  end

  defp reference_satellite(common, opts) do
    case Keyword.get(opts, :reference_satellite_id) do
      nil ->
        {:ok, hd(common)}

      sat when is_binary(sat) ->
        if sat in common do
          {:ok, sat}
        else
          {:error, {:reference_satellite_missing, sat}}
        end

      _other ->
        {:error, {:invalid_option, :reference_satellite_id}}
    end
  end

  defp baseline_reference_satellite(common, opts, base, epochs) do
    case Keyword.get(opts, :reference_satellite_id) do
      nil ->
        {:ok, highest_elevation_reference(common, base, epochs)}

      sat when is_binary(sat) ->
        if sat in common do
          {:ok, sat}
        else
          {:error, {:reference_satellite_missing, sat}}
        end

      _other ->
        {:error, {:invalid_option, :reference_satellite_id}}
    end
  end

  defp highest_elevation_reference(common, base, epochs) do
    common
    |> Enum.map(fn sat -> {sat, average_elevation_score(base, epochs, sat)} end)
    |> Enum.sort_by(fn {sat, score} -> {-score, sat} end)
    |> hd()
    |> elem(0)
  end

  defp average_elevation_score(base, epochs, sat) do
    {:ok, up} = local_up(base)

    scores =
      Enum.map(epochs, fn epoch ->
        sat_pos = Map.fetch!(epoch.positions, sat)

        case unit3(sub3(sat_pos, base)) do
          {:ok, los} -> dot3(los, up)
          :zero -> -1.0
        end
      end)

    Enum.sum(scores) / length(scores)
  end

  defp local_up(base) do
    case norm(base) do
      n when n > 0.0 -> {:ok, scale3(base, 1.0 / n)}
      _zero -> {:ok, {0.0, 0.0, 1.0}}
    end
  end

  defp unit3(v) do
    case norm(v) do
      n when n > 0.0 -> {:ok, scale3(v, 1.0 / n)}
      _zero -> :zero
    end
  end

  defp dot3({ax, ay, az}, {bx, by, bz}), do: ax * bx + ay * by + az * bz
end
