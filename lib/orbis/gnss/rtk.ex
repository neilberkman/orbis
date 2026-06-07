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

  alias Orbis.GNSS.Core.Types

  @default_max_iterations 8
  @default_position_tolerance_m 1.0e-4
  @default_ambiguity_tolerance_m 1.0e-4
  @default_code_sigma_m 1.0
  @default_phase_sigma_m 0.02
  @pivot_epsilon 1.0e-12

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

  Options:

    * `:reference_satellite_id` - fixed reference satellite. When omitted, the
      lexicographically first satellite common to every epoch is used.
    * `:initial_baseline_m` - initial base-to-rover ECEF vector, default
      `{0.0, 0.0, 0.0}`.
    * `:code_sigma_m` / `:phase_sigma_m` - row sigmas in metres, defaults
      `#{@default_code_sigma_m}` and `#{@default_phase_sigma_m}`.
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
         {:ok, reference_sat} <- reference_satellite(common_sats, opts),
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
    with {:ok, code} <- sigma_weight(opts, :code_sigma_m, @default_code_sigma_m),
         {:ok, phase} <- sigma_weight(opts, :phase_sigma_m, @default_phase_sigma_m) do
      {:ok, %{code: code, phase: phase}}
    end
  end

  defp sigma_weight(opts, key, default) do
    sigma = Keyword.get(opts, key, default)

    if is_number(sigma) and sigma > 0.0,
      do: {:ok, 1.0 / sigma},
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
         {:ok, dx} <- solve_normal_equations(rows, baseline_unknown_count(ambiguity_sats)) do
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
        epoch: epoch.epoch,
        sat: sat,
        h: h_base,
        y: obs_dd.code_m - geom_dd,
        weight: weights.code
      }

      phase_row = %{
        kind: :phase,
        epoch: epoch.epoch,
        sat: sat,
        h: h_phase,
        y: obs_dd.phase_m - (geom_dd + ambiguity),
        weight: weights.phase
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
           build_baseline_rows(base, epochs, reference_sat, ambiguity_sats, state, weights) do
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

  defp solve_normal_equations(rows, n) do
    {ata, aty} = normal_equations(rows, n)
    solve_linear(ata, aty)
  end

  defp normal_equations(rows, n) do
    Enum.reduce(rows, {zero_matrix(n), zero_vector(n)}, fn row, {ata, aty} ->
      h = Enum.map(row.h, &(&1 * row.weight))
      y = row.y * row.weight
      {accumulate_ata(ata, h), accumulate_aty(aty, h, y)}
    end)
  end

  defp zero_matrix(n), do: for(_ <- 1..n, do: zero_vector(n))
  defp zero_vector(n), do: for(_ <- 1..n, do: 0.0)

  defp accumulate_ata(ata, h) do
    Enum.with_index(ata)
    |> Enum.map(fn {row, i} ->
      hi = Enum.at(h, i)

      row
      |> Enum.with_index()
      |> Enum.map(fn {aij, j} -> aij + hi * Enum.at(h, j) end)
    end)
  end

  defp accumulate_aty(aty, h, y) do
    aty
    |> Enum.zip(h)
    |> Enum.map(fn {acc, hi} -> acc + hi * y end)
  end

  defp solve_linear(a, b) do
    augmented = Enum.zip_with(a, b, fn row, bi -> row ++ [bi] end)

    case eliminate(augmented, 0, length(b)) do
      :singular -> {:error, :singular_geometry}
      upper -> {:ok, back_substitute(upper)}
    end
  end

  defp eliminate(rows, col, n) when col >= n, do: rows

  defp eliminate(rows, col, n) do
    {pivot_row, pivot_abs} =
      rows
      |> Enum.with_index()
      |> Enum.drop(col)
      |> Enum.map(fn {row, idx} -> {idx, abs(Enum.at(row, col))} end)
      |> Enum.max_by(fn {_idx, value} -> value end)

    if pivot_abs <= @pivot_epsilon do
      :singular
    else
      rows = swap_rows(rows, col, pivot_row)
      pivot = Enum.at(rows, col)
      pivot_value = Enum.at(pivot, col)

      rows =
        rows
        |> Enum.with_index()
        |> Enum.map(fn {row, idx} ->
          if idx <= col do
            row
          else
            factor = Enum.at(row, col) / pivot_value

            row
            |> Enum.zip(pivot)
            |> Enum.map(fn {rij, pij} -> rij - factor * pij end)
          end
        end)

      eliminate(rows, col + 1, n)
    end
  end

  defp swap_rows(rows, i, i), do: rows

  defp swap_rows(rows, i, j) do
    ri = Enum.at(rows, i)
    rj = Enum.at(rows, j)

    rows
    |> List.replace_at(i, rj)
    |> List.replace_at(j, ri)
  end

  defp back_substitute(rows) do
    n = length(rows)

    solved =
      (n - 1)..0//-1
      |> Enum.reduce(%{}, fn i, solved ->
        row = Enum.at(rows, i)

        known =
          if i == n - 1 do
            0.0
          else
            Enum.reduce((i + 1)..(n - 1), 0.0, fn j, acc ->
              acc + Enum.at(row, j) * Map.fetch!(solved, j)
            end)
          end

        xi = (Enum.at(row, n) - known) / Enum.at(row, i)
        Map.put(solved, i, xi)
      end)

    Enum.map(0..(n - 1), &Map.fetch!(solved, &1))
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
end
