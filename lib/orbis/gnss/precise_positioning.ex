defmodule Orbis.GNSS.PrecisePositioning do
  @moduledoc """
  Carrier-phase precise-positioning primitives.

  This is the first precise-positioning layer above the code and carrier-phase
  combinations in `Orbis.GNSS.IonosphereFree` / `Orbis.GNSS.CarrierPhase`. It
  solves one SP3-backed epoch from dual-frequency ionosphere-free code and phase
  observations:

      P_IF_i = rho_i(x) + b - c * dt_sat_i
      L_IF_i = rho_i(x) + b - c * dt_sat_i + N_i

  where `x` is the receiver ECEF position, `b` is the receiver clock in metres,
  and `N_i` is one float carrier-phase ambiguity per satellite, also in metres.
  The state is linearized and iterated over `[x, y, z, b, N_1, N_2, ...]`.

  `solve_float/4` solves one epoch. `solve_float_epochs/3` solves a static
  multi-epoch arc with one receiver position, one receiver clock per epoch, and
  one ambiguity per satellite held constant across the arc. That multi-epoch
  model is the first step where carrier phase can tighten position instead of
  being absorbed entirely by one ambiguity per epoch.

  `solve_fixed_epochs/3` starts from the same multi-epoch float model, searches
  integer ambiguity candidates on an explicit caller-supplied wavelength grid,
  then re-solves position and per-epoch clocks with those ambiguities held
  fixed. It is intentionally explicit about wavelengths because an
  ionosphere-free carrier phase in metres has no universal integer grid unless
  the caller supplies the effective combination wavelength.

  ## Observation shape

  Observations may be maps or tuples:

      %{satellite_id: "G05", code_m: 24_000_000.0, phase_m: 24_012_345.0}
      {"G05", 24_000_000.0, 24_012_345.0}

  `code_m` and `phase_m` should normally be ionosphere-free combinations. Use
  `Orbis.GNSS.IonosphereFree.iono_free/4` and
  `Orbis.GNSS.IonosphereFree.iono_free_phase_cycles/4` to form them from raw
  dual-frequency RINEX observations.
  """

  alias Orbis.GNSS.Core.Constants
  alias Orbis.GNSS.Observables
  alias Orbis.GNSS.Positioning
  alias Orbis.GNSS.SP3

  @default_max_iterations 8
  @default_position_tolerance_m 1.0e-4
  @default_clock_tolerance_m 1.0e-4
  @default_code_sigma_m 1.0
  @default_phase_sigma_m 0.01
  @default_integer_search_radius_cycles 1
  @default_integer_ratio_threshold 3.0
  @default_integer_candidate_limit 50_000
  @pivot_epsilon 1.0e-12

  defmodule Solution do
    @moduledoc """
    Float-ambiguity phase positioning solution for one epoch.
    """

    @enforce_keys [
      :position,
      :rx_clock_s,
      :rx_clock_m,
      :ambiguities_m,
      :residuals_m,
      :used_sats,
      :metadata
    ]
    defstruct [
      :position,
      :rx_clock_s,
      :rx_clock_m,
      :ambiguities_m,
      :residuals_m,
      :used_sats,
      :metadata
    ]

    @type position :: %{x_m: float(), y_m: float(), z_m: float()}
    @type residual :: %{code_m: float(), phase_m: float()}

    @type t :: %__MODULE__{
            position: position(),
            rx_clock_s: float(),
            rx_clock_m: float(),
            ambiguities_m: %{String.t() => float()},
            residuals_m: %{String.t() => residual()},
            used_sats: [String.t()],
            metadata: %{
              iterations: pos_integer(),
              converged: boolean(),
              status: :position_tolerance | :max_iterations,
              code_rms_m: float(),
              phase_rms_m: float(),
              weighted_rms_m: float()
            }
          }
  end

  defmodule MultiEpochSolution do
    @moduledoc """
    Static multi-epoch float-ambiguity phase positioning solution.
    """

    @enforce_keys [
      :position,
      :epoch_clocks,
      :ambiguities_m,
      :residuals_m,
      :used_sats,
      :epochs,
      :metadata
    ]
    defstruct [
      :position,
      :epoch_clocks,
      :ambiguities_m,
      :residuals_m,
      :used_sats,
      :epochs,
      :metadata
    ]

    @type position :: %{x_m: float(), y_m: float(), z_m: float()}

    @type epoch_clock :: %{
            epoch: NaiveDateTime.t(),
            rx_clock_s: float(),
            rx_clock_m: float()
          }

    @type residual :: %{
            epoch: NaiveDateTime.t(),
            satellite_id: String.t(),
            code_m: float(),
            phase_m: float()
          }

    @type t :: %__MODULE__{
            position: position(),
            epoch_clocks: [epoch_clock()],
            ambiguities_m: %{String.t() => float()},
            residuals_m: [residual()],
            used_sats: [String.t()],
            epochs: [NaiveDateTime.t()],
            metadata: %{
              iterations: pos_integer(),
              converged: boolean(),
              status: :state_tolerance | :max_iterations,
              n_epochs: pos_integer(),
              n_observations: pos_integer(),
              code_rms_m: float(),
              phase_rms_m: float(),
              weighted_rms_m: float()
            }
          }
  end

  defmodule FixedSolution do
    @moduledoc """
    Static multi-epoch integer-fixed carrier-phase solution.
    """

    @enforce_keys [
      :position,
      :epoch_clocks,
      :fixed_ambiguities_cycles,
      :fixed_ambiguities_m,
      :float_solution,
      :residuals_m,
      :used_sats,
      :epochs,
      :metadata
    ]
    defstruct [
      :position,
      :epoch_clocks,
      :fixed_ambiguities_cycles,
      :fixed_ambiguities_m,
      :float_solution,
      :residuals_m,
      :used_sats,
      :epochs,
      :metadata
    ]

    @type position :: %{x_m: float(), y_m: float(), z_m: float()}

    @type epoch_clock :: %{
            epoch: NaiveDateTime.t(),
            rx_clock_s: float(),
            rx_clock_m: float()
          }

    @type residual :: %{
            epoch: NaiveDateTime.t(),
            satellite_id: String.t(),
            code_m: float(),
            phase_m: float()
          }

    @type t :: %__MODULE__{
            position: position(),
            epoch_clocks: [epoch_clock()],
            fixed_ambiguities_cycles: %{String.t() => integer()},
            fixed_ambiguities_m: %{String.t() => float()},
            float_solution: MultiEpochSolution.t(),
            residuals_m: [residual()],
            used_sats: [String.t()],
            epochs: [NaiveDateTime.t()],
            metadata: %{
              iterations: pos_integer(),
              converged: boolean(),
              status: :state_tolerance | :max_iterations,
              n_epochs: pos_integer(),
              n_observations: pos_integer(),
              code_rms_m: float(),
              phase_rms_m: float(),
              weighted_rms_m: float(),
              integer_status: :fixed | :not_fixed,
              integer_ratio: float() | :infinity,
              integer_best_score: float(),
              integer_second_best_score: float() | nil,
              integer_candidates: pos_integer()
            }
          }
  end

  @typedoc "A dual-frequency ionosphere-free code/phase observation."
  @type observation ::
          %{satellite_id: String.t(), code_m: number(), phase_m: number()}
          | {String.t(), number(), number()}

  @typedoc "A receiver ECEF position in metres."
  @type receiver ::
          {number(), number(), number()} | %{x_m: number(), y_m: number(), z_m: number()}

  @typedoc "A set of code/phase observations for one epoch."
  @type epoch_observations ::
          %{epoch: NaiveDateTime.t(), observations: [observation()]}
          | {NaiveDateTime.t(), [observation()]}

  @doc """
  Solve a float-ambiguity carrier-phase position for one SP3-backed epoch.

  `source` is a loaded `Orbis.GNSS.SP3` product. `observations` is a list of
  ionosphere-free code/phase pairs for one epoch. `epoch` is interpreted in the
  SP3 product's time scale.

  ## Options

    * `:initial_guess` - `{x_m, y_m, z_m, clock_m}`. If omitted, the code
      observations are first passed through `Orbis.GNSS.Positioning.solve/4`
      with ionosphere/troposphere disabled, and that code-only solution seeds
      the float solve.
    * `:spp_initial_guess` - code-only SPP seed used only when `:initial_guess`
      is omitted (default `{0, 0, 0, 0}`).
    * `:code_sigma_m` - code row standard deviation (default `1.0` m).
    * `:phase_sigma_m` - phase row standard deviation (default `0.01` m).
    * `:max_iterations` - maximum nonlinear iterations (default `8`).
    * `:position_tolerance_m` - position-update convergence threshold
      (default `1.0e-4` m).
    * `:clock_tolerance_m` - receiver-clock update threshold (default
      `1.0e-4` m).

  Returns `{:ok, %Solution{}}` or `{:error, reason}`. Reasons include
  `:no_observations`, `{:too_few_satellites, used, 4}`,
  `{:duplicate_observation, sat}`, `{:invalid_observation, entry}`,
  `:invalid_initial_guess`, `{:invalid_sigma, key}`, `{:invalid_option, key}`,
  `{:code_seed_failed, reason}`, `{:no_ephemeris, sat, reason}`, and
  `:singular_geometry`. If the iteration limit is reached after a valid solve
  step, the function returns a solution with `metadata.converged == false` and
  `metadata.status == :max_iterations` so callers can inspect the residuals and
  decide whether to reject it.
  """
  @spec solve_float(SP3.t(), [observation()], NaiveDateTime.t(), keyword()) ::
          {:ok, Solution.t()} | {:error, term()}
  def solve_float(source, observations, epoch, opts \\ [])

  def solve_float(%SP3{} = sp3, observations, %NaiveDateTime{} = epoch, opts)
      when is_list(observations) do
    with :ok <- ensure_nonempty(observations),
         {:ok, obs} <- normalize_observations(observations),
         :ok <- ensure_enough(obs),
         {:ok, weights} <- weights(opts),
         {:ok, solve_opts} <- solve_options(opts),
         {:ok, state} <- initial_state(sp3, obs, epoch, opts) do
      iterate(sp3, obs, epoch, state, weights, solve_opts, 1)
    end
  end

  def solve_float(%SP3{}, observations, %NaiveDateTime{}, _opts) when not is_list(observations),
    do: {:error, :no_observations}

  @doc """
  Solve a static multi-epoch float-ambiguity carrier-phase position.

  `epoch_observations` is a list of `%{epoch: epoch, observations: obs}` maps or
  `{epoch, obs}` tuples. The receiver position is static across the whole arc,
  each epoch gets its own receiver clock, and each satellite gets one ambiguity
  held constant across every epoch where that satellite appears.

  This model is still float ambiguity only. It does not fix integer ambiguities
  or estimate a stochastic PPP process, but it lets changing geometry across the
  arc separate position from carrier ambiguities.

  Options are the same as `solve_float/4`, plus:

    * `:ambiguity_tolerance_m` - maximum ambiguity-update convergence threshold
      (default `1.0e-4` m).

  Returns `{:ok, %MultiEpochSolution{}}` or `{:error, reason}`. Reasons include
  `:no_epochs`, `{:too_few_epochs, used, 2}`, `{:duplicate_epoch, epoch}`,
  `{:too_few_epoch_observations, epoch, used, 4}`,
  `{:too_few_equations, equations, unknowns}`, and the same observation,
  option, ephemeris, seeding, and geometry errors as `solve_float/4`.
  """
  @spec solve_float_epochs(SP3.t(), [epoch_observations()], keyword()) ::
          {:ok, MultiEpochSolution.t()} | {:error, term()}
  def solve_float_epochs(source, epoch_observations, opts \\ [])

  def solve_float_epochs(%SP3{} = sp3, epoch_observations, opts)
      when is_list(epoch_observations) do
    with {:ok, epochs} <- normalize_epoch_observations(epoch_observations),
         :ok <- ensure_multi_enough(epochs),
         {:ok, weights} <- weights(opts),
         {:ok, solve_opts} <- solve_options(opts),
         {:ok, state} <- initial_multi_state(sp3, epochs, opts) do
      sat_ids = multi_satellite_ids(epochs)
      iterate_multi(sp3, epochs, sat_ids, state, weights, solve_opts, 1)
    end
  end

  def solve_float_epochs(%SP3{}, _epoch_observations, _opts), do: {:error, :no_epochs}

  @doc """
  Solve a static multi-epoch position with integer-fixed ambiguities.

  The function first solves the float multi-epoch model (`solve_float_epochs/3`),
  converts each float ambiguity from metres to cycles using the explicit
  `:ambiguity_wavelength_m` option, searches nearby integer candidates, and
  re-solves the receiver position and per-epoch clocks with the best integer
  ambiguities held fixed.

  ## Required option

    * `:ambiguity_wavelength_m` - either a positive scalar wavelength in metres
      for every satellite, or a map `%{"G05" => wavelength_m, ...}`.

  ## Additional options

    * `:integer_search_radius_cycles` - integer half-window around each rounded
      float ambiguity (default `1`).
    * `:integer_ratio_threshold` - minimum second-best / best weighted-score
      ratio for `metadata.integer_status == :fixed` (default `3.0`).
    * `:integer_candidate_limit` - maximum candidates to evaluate before
      returning `{:error, {:too_many_integer_candidates, count, limit}}`
      (default `50_000`).

  The fixed solution is returned even when the ratio test is not met; in that
  case `metadata.integer_status` is `:not_fixed` so callers can reject it.
  """
  @spec solve_fixed_epochs(SP3.t(), [epoch_observations()], keyword()) ::
          {:ok, FixedSolution.t()} | {:error, term()}
  def solve_fixed_epochs(source, epoch_observations, opts \\ [])

  def solve_fixed_epochs(%SP3{} = sp3, epoch_observations, opts)
      when is_list(epoch_observations) do
    with {:ok, epochs} <- normalize_epoch_observations(epoch_observations),
         :ok <- ensure_multi_enough(epochs),
         {:ok, weights} <- weights(opts),
         {:ok, solve_opts} <- solve_options(opts),
         {:ok, integer_opts} <- integer_options(opts),
         {:ok, float_sol} <- solve_float_epochs(sp3, epochs, opts),
         {:ok, wavelengths} <- ambiguity_wavelengths(float_sol.used_sats, opts),
         {:ok, fixed_cycles, fixed_meta} <-
           search_integer_ambiguities(float_sol, wavelengths, integer_opts),
         fixed_m = fixed_ambiguities_m(fixed_cycles, wavelengths),
         state = fixed_state_from_float(float_sol),
         {:ok, fixed_state, iterations, converged, status} <-
           iterate_fixed_multi(sp3, epochs, fixed_m, state, weights, solve_opts, 1) do
      finalize_fixed_multi(
        sp3,
        epochs,
        float_sol.used_sats,
        fixed_cycles,
        fixed_m,
        float_sol,
        fixed_state,
        weights,
        iterations,
        converged,
        status,
        fixed_meta
      )
    end
  end

  def solve_fixed_epochs(%SP3{}, _epoch_observations, _opts), do: {:error, :no_epochs}

  # --- input normalization -------------------------------------------------

  defp ensure_nonempty([]), do: {:error, :no_observations}
  defp ensure_nonempty(_), do: :ok

  defp normalize_observations(observations) do
    Enum.reduce_while(observations, {:ok, [], MapSet.new()}, fn entry, {:ok, acc, seen} ->
      case normalize_one(entry) do
        {:ok, %{satellite_id: sat} = obs} ->
          if MapSet.member?(seen, sat) do
            {:halt, {:error, {:duplicate_observation, sat}}}
          else
            {:cont, {:ok, [obs | acc], MapSet.put(seen, sat)}}
          end

        {:error, _} = err ->
          {:halt, err}
      end
    end)
    |> case do
      {:ok, acc, _seen} ->
        sorted = acc |> Enum.reverse() |> Enum.sort_by(& &1.satellite_id)
        {:ok, sorted}

      {:error, _} = err ->
        err
    end
  end

  defp normalize_one(%{satellite_id: sat, code_m: code, phase_m: phase} = obs)
       when is_binary(sat) and is_number(code) and is_number(phase) do
    {:ok, %{satellite_id: sat, code_m: code / 1.0, phase_m: phase / 1.0, raw: obs}}
  end

  defp normalize_one({sat, code, phase} = obs)
       when is_binary(sat) and is_number(code) and is_number(phase) do
    {:ok, %{satellite_id: sat, code_m: code / 1.0, phase_m: phase / 1.0, raw: obs}}
  end

  defp normalize_one(entry), do: {:error, {:invalid_observation, entry}}

  defp ensure_enough(obs) when length(obs) >= 4, do: :ok
  defp ensure_enough(obs), do: {:error, {:too_few_satellites, length(obs), 4}}

  defp normalize_epoch_observations([]), do: {:error, :no_epochs}

  defp normalize_epoch_observations(epoch_observations) do
    epoch_observations
    |> Enum.reduce_while({:ok, [], MapSet.new()}, fn entry, {:ok, acc, seen} ->
      case normalize_epoch_entry(entry) do
        {:ok, epoch, observations} ->
          if MapSet.member?(seen, epoch) do
            {:halt, {:error, {:duplicate_epoch, epoch}}}
          else
            with {:ok, obs} <- normalize_observations(observations),
                 :ok <- ensure_epoch_enough(epoch, obs) do
              {:cont, {:ok, [%{epoch: epoch, observations: obs} | acc], MapSet.put(seen, epoch)}}
            else
              {:error, _} = err -> {:halt, err}
            end
          end

        {:error, _} = err ->
          {:halt, err}
      end
    end)
    |> case do
      {:ok, acc, _seen} ->
        {:ok, Enum.sort_by(acc, &NaiveDateTime.to_iso8601(&1.epoch))}

      {:error, _} = err ->
        err
    end
  end

  defp normalize_epoch_entry(%{epoch: %NaiveDateTime{} = epoch, observations: observations})
       when is_list(observations), do: {:ok, epoch, observations}

  defp normalize_epoch_entry({%NaiveDateTime{} = epoch, observations}) when is_list(observations),
    do: {:ok, epoch, observations}

  defp normalize_epoch_entry(entry), do: {:error, {:invalid_epoch_observations, entry}}

  defp ensure_epoch_enough(_epoch, obs) when length(obs) >= 4, do: :ok

  defp ensure_epoch_enough(epoch, obs),
    do: {:error, {:too_few_epoch_observations, epoch, length(obs), 4}}

  defp ensure_multi_enough(epochs) when length(epochs) < 2,
    do: {:error, {:too_few_epochs, length(epochs), 2}}

  defp ensure_multi_enough(epochs) do
    n_epochs = length(epochs)
    n_sats = length(multi_satellite_ids(epochs))
    n_observations = multi_observation_count(epochs)
    equations = 2 * n_observations
    unknowns = 3 + n_epochs + n_sats

    cond do
      n_sats < 4 ->
        {:error, {:too_few_satellites, n_sats, 4}}

      equations < unknowns ->
        {:error, {:too_few_equations, equations, unknowns}}

      true ->
        :ok
    end
  end

  defp weights(opts) do
    code_sigma = Keyword.get(opts, :code_sigma_m, @default_code_sigma_m)
    phase_sigma = Keyword.get(opts, :phase_sigma_m, @default_phase_sigma_m)

    cond do
      not is_number(code_sigma) or code_sigma <= 0.0 ->
        {:error, {:invalid_sigma, :code_sigma_m}}

      not is_number(phase_sigma) or phase_sigma <= 0.0 ->
        {:error, {:invalid_sigma, :phase_sigma_m}}

      true ->
        {:ok, %{code: 1.0 / code_sigma, phase: 1.0 / phase_sigma}}
    end
  end

  defp solve_options(opts) do
    max_iterations = Keyword.get(opts, :max_iterations, @default_max_iterations)
    pos_tol = Keyword.get(opts, :position_tolerance_m, @default_position_tolerance_m)
    clock_tol = Keyword.get(opts, :clock_tolerance_m, @default_clock_tolerance_m)
    ambiguity_tol = Keyword.get(opts, :ambiguity_tolerance_m, @default_position_tolerance_m)

    cond do
      not is_integer(max_iterations) or max_iterations < 1 ->
        {:error, {:invalid_option, :max_iterations}}

      not is_number(pos_tol) or pos_tol < 0.0 ->
        {:error, {:invalid_option, :position_tolerance_m}}

      not is_number(clock_tol) or clock_tol < 0.0 ->
        {:error, {:invalid_option, :clock_tolerance_m}}

      not is_number(ambiguity_tol) or ambiguity_tol < 0.0 ->
        {:error, {:invalid_option, :ambiguity_tolerance_m}}

      true ->
        {:ok,
         %{
           max_iterations: max_iterations,
           position_tolerance_m: pos_tol / 1.0,
           clock_tolerance_m: clock_tol / 1.0,
           ambiguity_tolerance_m: ambiguity_tol / 1.0
         }}
    end
  end

  defp integer_options(opts) do
    radius =
      Keyword.get(opts, :integer_search_radius_cycles, @default_integer_search_radius_cycles)

    ratio = Keyword.get(opts, :integer_ratio_threshold, @default_integer_ratio_threshold)
    limit = Keyword.get(opts, :integer_candidate_limit, @default_integer_candidate_limit)

    cond do
      not is_integer(radius) or radius < 0 ->
        {:error, {:invalid_option, :integer_search_radius_cycles}}

      not is_number(ratio) or ratio < 0.0 ->
        {:error, {:invalid_option, :integer_ratio_threshold}}

      not is_integer(limit) or limit < 1 ->
        {:error, {:invalid_option, :integer_candidate_limit}}

      true ->
        {:ok,
         %{
           radius_cycles: radius,
           ratio_threshold: ratio / 1.0,
           candidate_limit: limit
         }}
    end
  end

  defp ambiguity_wavelengths(sat_ids, opts) do
    case Keyword.fetch(opts, :ambiguity_wavelength_m) do
      {:ok, wavelength} when is_number(wavelength) and wavelength > 0.0 ->
        {:ok, Map.new(sat_ids, &{&1, wavelength / 1.0})}

      {:ok, wavelength_by_sat} when is_map(wavelength_by_sat) ->
        sat_ids
        |> Enum.reduce_while({:ok, %{}}, fn sat, {:ok, acc} ->
          case Map.fetch(wavelength_by_sat, sat) do
            {:ok, value} when is_number(value) and value > 0.0 ->
              {:cont, {:ok, Map.put(acc, sat, value / 1.0)}}

            _ ->
              {:halt, {:error, {:invalid_ambiguity_wavelength, sat}}}
          end
        end)

      {:ok, _other} ->
        {:error, {:invalid_option, :ambiguity_wavelength_m}}

      :error ->
        {:error, :ambiguity_wavelength_required}
    end
  end

  # --- initialization ------------------------------------------------------

  defp initial_state(sp3, obs, epoch, opts) do
    case Keyword.fetch(opts, :initial_guess) do
      {:ok, guess} ->
        with {:ok, {x, y, z, clock_m}} <- normalize_guess(guess) do
          {:ok, state_from_guess(obs, {x, y, z}, clock_m)}
        end

      :error ->
        spp_seed(sp3, obs, epoch, opts)
    end
  end

  defp normalize_guess({x, y, z, clock_m})
       when is_number(x) and is_number(y) and is_number(z) and is_number(clock_m),
       do: {:ok, {x / 1.0, y / 1.0, z / 1.0, clock_m / 1.0}}

  defp normalize_guess(_guess), do: {:error, :invalid_initial_guess}

  defp state_from_guess(obs, position, clock_m) do
    ambiguities =
      Map.new(obs, fn o ->
        {o.satellite_id, o.phase_m - o.code_m}
      end)

    %{position: position, clock_m: clock_m, ambiguities: ambiguities}
  end

  defp spp_seed(sp3, obs, epoch, opts) do
    observations = Enum.map(obs, &{&1.satellite_id, &1.code_m})
    spp_initial = Keyword.get(opts, :spp_initial_guess, {0.0, 0.0, 0.0, 0.0})

    case Positioning.solve(sp3, observations, epoch,
           ionosphere: false,
           troposphere: false,
           initial_guess: spp_initial,
           with_geodetic: false
         ) do
      {:ok, sol} ->
        pos = {sol.position.x_m, sol.position.y_m, sol.position.z_m}
        state = state_from_guess(obs, pos, sol.rx_clock_s * Constants.speed_of_light_m_s())
        {:ok, state}

      {:error, reason} ->
        {:error, {:code_seed_failed, reason}}
    end
  end

  defp initial_multi_state(sp3, epochs, opts) do
    case Keyword.fetch(opts, :initial_guess) do
      {:ok, guess} ->
        with {:ok, {x, y, z, clock_m}} <- normalize_guess(guess) do
          {:ok,
           %{
             position: {x, y, z},
             clocks_m: List.duplicate(clock_m, length(epochs)),
             ambiguities: initial_ambiguities(epochs)
           }}
        end

      :error ->
        multi_spp_seed(sp3, epochs, opts)
    end
  end

  defp initial_ambiguities(epochs) do
    epochs
    |> Enum.flat_map(& &1.observations)
    |> Enum.reduce(%{}, fn obs, acc ->
      Map.put_new(acc, obs.satellite_id, obs.phase_m - obs.code_m)
    end)
  end

  defp multi_spp_seed(sp3, epochs, opts) do
    epochs
    |> Enum.reduce_while({:ok, [], []}, fn epoch_row, {:ok, positions, clocks} ->
      observations = Enum.map(epoch_row.observations, &{&1.satellite_id, &1.code_m})
      spp_initial = Keyword.get(opts, :spp_initial_guess, {0.0, 0.0, 0.0, 0.0})

      case Positioning.solve(sp3, observations, epoch_row.epoch,
             ionosphere: false,
             troposphere: false,
             initial_guess: spp_initial,
             with_geodetic: false
           ) do
        {:ok, sol} ->
          pos = {sol.position.x_m, sol.position.y_m, sol.position.z_m}
          clock_m = sol.rx_clock_s * Constants.speed_of_light_m_s()
          {:cont, {:ok, [pos | positions], [clock_m | clocks]}}

        {:error, reason} ->
          {:halt, {:error, {:code_seed_failed, epoch_row.epoch, reason}}}
      end
    end)
    |> case do
      {:ok, positions, clocks} ->
        {:ok,
         %{
           position: mean_position(positions),
           clocks_m: Enum.reverse(clocks),
           ambiguities: initial_ambiguities(epochs)
         }}

      {:error, _} = err ->
        err
    end
  end

  defp mean_position(positions) do
    {sx, sy, sz} =
      Enum.reduce(positions, {0.0, 0.0, 0.0}, fn {x, y, z}, {ax, ay, az} ->
        {ax + x, ay + y, az + z}
      end)

    n = length(positions)
    {sx / n, sy / n, sz / n}
  end

  # --- nonlinear solve -----------------------------------------------------

  defp iterate(sp3, obs, epoch, state, weights, opts, iter) do
    with {:ok, rows} <- build_rows(sp3, obs, epoch, state, weights),
         {:ok, dx} <- solve_normal_equations(rows, 4 + length(obs)) do
      next = apply_delta(state, obs, dx)
      {pos_step, clock_step} = step_norms(dx)

      cond do
        pos_step <= opts.position_tolerance_m and abs(clock_step) <= opts.clock_tolerance_m ->
          finalize(sp3, obs, epoch, next, weights, iter, true, :position_tolerance)

        iter >= opts.max_iterations ->
          finalize(sp3, obs, epoch, next, weights, iter, false, :max_iterations)

        true ->
          iterate(sp3, obs, epoch, next, weights, opts, iter + 1)
      end
    end
  end

  defp iterate_multi(sp3, epochs, sat_ids, state, weights, opts, iter) do
    with {:ok, rows} <- build_multi_rows(sp3, epochs, sat_ids, state, weights),
         {:ok, dx} <- solve_normal_equations(rows, 3 + length(epochs) + length(sat_ids)) do
      next = apply_multi_delta(state, epochs, sat_ids, dx)
      {pos_step, clock_step, ambiguity_step} = multi_step_norms(dx, length(epochs))

      cond do
        pos_step <= opts.position_tolerance_m and
          clock_step <= opts.clock_tolerance_m and
            ambiguity_step <= opts.ambiguity_tolerance_m ->
          finalize_multi(
            sp3,
            epochs,
            sat_ids,
            next,
            weights,
            iter,
            true,
            :state_tolerance
          )

        iter >= opts.max_iterations ->
          finalize_multi(sp3, epochs, sat_ids, next, weights, iter, false, :max_iterations)

        true ->
          iterate_multi(sp3, epochs, sat_ids, next, weights, opts, iter + 1)
      end
    end
  end

  defp iterate_fixed_multi(sp3, epochs, fixed_m, state, weights, opts, iter) do
    with {:ok, rows} <- build_fixed_multi_rows(sp3, epochs, fixed_m, state, weights),
         {:ok, dx} <- solve_normal_equations(rows, 3 + length(epochs)) do
      next = apply_fixed_multi_delta(state, dx, length(epochs))
      {pos_step, clock_step} = fixed_multi_step_norms(dx)

      cond do
        pos_step <= opts.position_tolerance_m and clock_step <= opts.clock_tolerance_m ->
          {:ok, next, iter, true, :state_tolerance}

        iter >= opts.max_iterations ->
          {:ok, next, iter, false, :max_iterations}

        true ->
          iterate_fixed_multi(sp3, epochs, fixed_m, next, weights, opts, iter + 1)
      end
    end
  end

  defp build_rows(sp3, obs, epoch, state, weights) do
    {rx, ry, rz} = state.position

    Enum.reduce_while(obs, {:ok, []}, fn o, {:ok, acc} ->
      case Observables.predict(sp3, o.satellite_id, {rx, ry, rz}, epoch) do
        {:ok, pred} ->
          {ex, ey, ez} = pred.los_unit
          sat_clock_m = Constants.speed_of_light_m_s() * (pred.sat_clock_s || 0.0)
          model_code = pred.geometric_range_m + state.clock_m - sat_clock_m
          ambiguity = Map.fetch!(state.ambiguities, o.satellite_id)

          code_row = %{
            kind: :code,
            sat: o.satellite_id,
            h: design_row({-ex, -ey, -ez, 1.0}, nil, obs),
            y: o.code_m - model_code,
            weight: weights.code
          }

          phase_row = %{
            kind: :phase,
            sat: o.satellite_id,
            h: design_row({-ex, -ey, -ez, 1.0}, o.satellite_id, obs),
            y: o.phase_m - (model_code + ambiguity),
            weight: weights.phase
          }

          {:cont, {:ok, [phase_row, code_row | acc]}}

        {:error, reason} ->
          {:halt, {:error, {:no_ephemeris, o.satellite_id, reason}}}
      end
    end)
    |> case do
      {:ok, rows} -> {:ok, Enum.reverse(rows)}
      {:error, _} = err -> err
    end
  end

  defp design_row({dx, dy, dz, dc}, ambiguity_sat, obs) do
    ambiguity_cols =
      Enum.map(obs, fn o ->
        if o.satellite_id == ambiguity_sat, do: 1.0, else: 0.0
      end)

    [dx, dy, dz, dc | ambiguity_cols]
  end

  defp build_multi_rows(sp3, epochs, sat_ids, state, weights) do
    {rx, ry, rz} = state.position

    epochs
    |> Enum.with_index()
    |> Enum.reduce_while({:ok, []}, fn {epoch_row, epoch_idx}, {:ok, acc} ->
      clock_m = Enum.at(state.clocks_m, epoch_idx)

      epoch_row.observations
      |> Enum.reduce_while({:ok, acc}, fn o, {:ok, rows_acc} ->
        case Observables.predict(sp3, o.satellite_id, {rx, ry, rz}, epoch_row.epoch) do
          {:ok, pred} ->
            {ex, ey, ez} = pred.los_unit
            sat_clock_m = Constants.speed_of_light_m_s() * (pred.sat_clock_s || 0.0)
            model_code = pred.geometric_range_m + clock_m - sat_clock_m
            ambiguity = Map.fetch!(state.ambiguities, o.satellite_id)
            base = {-ex, -ey, -ez}

            code_row = %{
              kind: :code,
              sat: o.satellite_id,
              epoch: epoch_row.epoch,
              h: multi_design_row(base, epoch_idx, nil, length(epochs), sat_ids),
              y: o.code_m - model_code,
              weight: weights.code
            }

            phase_row = %{
              kind: :phase,
              sat: o.satellite_id,
              epoch: epoch_row.epoch,
              h: multi_design_row(base, epoch_idx, o.satellite_id, length(epochs), sat_ids),
              y: o.phase_m - (model_code + ambiguity),
              weight: weights.phase
            }

            {:cont, {:ok, [phase_row, code_row | rows_acc]}}

          {:error, reason} ->
            {:halt, {:error, {:no_ephemeris, o.satellite_id, reason}}}
        end
      end)
      |> case do
        {:ok, rows_acc} -> {:cont, {:ok, rows_acc}}
        {:error, _} = err -> {:halt, err}
      end
    end)
    |> case do
      {:ok, rows} -> {:ok, Enum.reverse(rows)}
      {:error, _} = err -> err
    end
  end

  defp build_fixed_multi_rows(sp3, epochs, fixed_m, state, weights) do
    {rx, ry, rz} = state.position

    epochs
    |> Enum.with_index()
    |> Enum.reduce_while({:ok, []}, fn {epoch_row, epoch_idx}, {:ok, acc} ->
      clock_m = Enum.at(state.clocks_m, epoch_idx)

      epoch_row.observations
      |> Enum.reduce_while({:ok, acc}, fn o, {:ok, rows_acc} ->
        case Observables.predict(sp3, o.satellite_id, {rx, ry, rz}, epoch_row.epoch) do
          {:ok, pred} ->
            {ex, ey, ez} = pred.los_unit
            sat_clock_m = Constants.speed_of_light_m_s() * (pred.sat_clock_s || 0.0)
            model_code = pred.geometric_range_m + clock_m - sat_clock_m
            ambiguity = Map.fetch!(fixed_m, o.satellite_id)
            h = fixed_multi_design_row({-ex, -ey, -ez}, epoch_idx, length(epochs))

            code_row = %{
              kind: :code,
              sat: o.satellite_id,
              epoch: epoch_row.epoch,
              h: h,
              y: o.code_m - model_code,
              weight: weights.code
            }

            phase_row = %{
              kind: :phase,
              sat: o.satellite_id,
              epoch: epoch_row.epoch,
              h: h,
              y: o.phase_m - (model_code + ambiguity),
              weight: weights.phase
            }

            {:cont, {:ok, [phase_row, code_row | rows_acc]}}

          {:error, reason} ->
            {:halt, {:error, {:no_ephemeris, o.satellite_id, reason}}}
        end
      end)
      |> case do
        {:ok, rows_acc} -> {:cont, {:ok, rows_acc}}
        {:error, _} = err -> {:halt, err}
      end
    end)
    |> case do
      {:ok, rows} -> {:ok, Enum.reverse(rows)}
      {:error, _} = err -> err
    end
  end

  defp multi_design_row({dx, dy, dz}, epoch_idx, ambiguity_sat, n_epochs, sat_ids) do
    clock_cols = for idx <- 0..(n_epochs - 1), do: if(idx == epoch_idx, do: 1.0, else: 0.0)
    ambiguity_cols = Enum.map(sat_ids, &if(&1 == ambiguity_sat, do: 1.0, else: 0.0))
    [dx, dy, dz | clock_cols ++ ambiguity_cols]
  end

  defp fixed_multi_design_row({dx, dy, dz}, epoch_idx, n_epochs) do
    clock_cols = for idx <- 0..(n_epochs - 1), do: if(idx == epoch_idx, do: 1.0, else: 0.0)
    [dx, dy, dz | clock_cols]
  end

  defp apply_delta(state, obs, dx) do
    {rx, ry, rz} = state.position
    [dx0, dx1, dx2, dclock | ambiguity_deltas] = dx

    ambiguities =
      obs
      |> Enum.zip(ambiguity_deltas)
      |> Map.new(fn {o, d} ->
        {o.satellite_id, Map.fetch!(state.ambiguities, o.satellite_id) + d}
      end)

    %{
      position: {rx + dx0, ry + dx1, rz + dx2},
      clock_m: state.clock_m + dclock,
      ambiguities: ambiguities
    }
  end

  defp step_norms([dx, dy, dz, dclock | _]) do
    {:math.sqrt(dx * dx + dy * dy + dz * dz), dclock}
  end

  defp apply_multi_delta(state, epochs, sat_ids, dx) do
    {rx, ry, rz} = state.position
    [dx0, dx1, dx2 | rest] = dx
    {clock_deltas, ambiguity_deltas} = Enum.split(rest, length(epochs))

    clocks =
      state.clocks_m
      |> Enum.zip(clock_deltas)
      |> Enum.map(fn {clock, delta} -> clock + delta end)

    ambiguities =
      sat_ids
      |> Enum.zip(ambiguity_deltas)
      |> Map.new(fn {sat, delta} -> {sat, Map.fetch!(state.ambiguities, sat) + delta} end)

    %{
      position: {rx + dx0, ry + dx1, rz + dx2},
      clocks_m: clocks,
      ambiguities: ambiguities
    }
  end

  defp apply_fixed_multi_delta(state, dx, n_epochs) do
    {rx, ry, rz} = state.position
    [dx0, dx1, dx2 | clock_deltas] = dx

    clocks =
      state.clocks_m
      |> Enum.zip(Enum.take(clock_deltas, n_epochs))
      |> Enum.map(fn {clock, delta} -> clock + delta end)

    %{
      position: {rx + dx0, ry + dx1, rz + dx2},
      clocks_m: clocks
    }
  end

  defp multi_step_norms([dx, dy, dz | rest], n_epochs) do
    {clock_deltas, ambiguity_deltas} = Enum.split(rest, n_epochs)

    {
      :math.sqrt(dx * dx + dy * dy + dz * dz),
      max_abs(clock_deltas),
      max_abs(ambiguity_deltas)
    }
  end

  defp fixed_multi_step_norms([dx, dy, dz | clock_deltas]) do
    {:math.sqrt(dx * dx + dy * dy + dz * dz), max_abs(clock_deltas)}
  end

  defp max_abs([]), do: 0.0
  defp max_abs(xs), do: xs |> Enum.map(&abs/1) |> Enum.max()

  # --- integer ambiguity fixing -------------------------------------------

  defp search_integer_ambiguities(float_sol, wavelengths, opts) do
    sat_ids = float_sol.used_sats

    ranges =
      Enum.map(
        sat_ids,
        &candidate_range(&1, float_sol.ambiguities_m, wavelengths, opts.radius_cycles)
      )

    count = Enum.reduce(ranges, 1, fn range, acc -> acc * length(range) end)

    if count > opts.candidate_limit do
      {:error, {:too_many_integer_candidates, count, opts.candidate_limit}}
    else
      candidates =
        ranges
        |> cartesian_product()
        |> Enum.map(fn cycles ->
          fixed_cycles = Map.new(Enum.zip(sat_ids, cycles))
          score = ambiguity_score(float_sol.ambiguities_m, fixed_cycles, wavelengths, sat_ids)
          {score, fixed_cycles}
        end)
        |> Enum.sort_by(fn {score, fixed_cycles} ->
          {score, Enum.map(sat_ids, &Map.fetch!(fixed_cycles, &1))}
        end)

      [{best_score, fixed_cycles} | rest] = candidates

      second_best_score =
        case rest do
          [{score, _} | _] -> score
          [] -> nil
        end

      ratio = integer_ratio(best_score, second_best_score)
      status = if integer_ratio_pass?(ratio, opts.ratio_threshold), do: :fixed, else: :not_fixed

      {:ok, fixed_cycles,
       %{
         integer_status: status,
         integer_ratio: ratio,
         integer_best_score: best_score,
         integer_second_best_score: second_best_score,
         integer_candidates: count
       }}
    end
  end

  defp candidate_range(sat, ambiguities_m, wavelengths, radius) do
    wavelength = Map.fetch!(wavelengths, sat)
    center = round(Map.fetch!(ambiguities_m, sat) / wavelength)
    (center - radius)..(center + radius) |> Enum.to_list()
  end

  defp cartesian_product([]), do: [[]]

  defp cartesian_product([range | rest]) do
    for value <- range, tail <- cartesian_product(rest), do: [value | tail]
  end

  defp ambiguity_score(ambiguities_m, fixed_cycles, wavelengths, sat_ids) do
    Enum.reduce(sat_ids, 0.0, fn sat, acc ->
      wavelength = Map.fetch!(wavelengths, sat)
      delta = Map.fetch!(ambiguities_m, sat) - Map.fetch!(fixed_cycles, sat) * wavelength
      normalized = delta / wavelength
      acc + normalized * normalized
    end)
  end

  defp fixed_ambiguities_m(fixed_cycles, wavelengths) do
    Map.new(fixed_cycles, fn {sat, cycles} -> {sat, cycles * Map.fetch!(wavelengths, sat)} end)
  end

  defp fixed_state_from_float(float_sol) do
    %{
      position: {float_sol.position.x_m, float_sol.position.y_m, float_sol.position.z_m},
      clocks_m: Enum.map(float_sol.epoch_clocks, & &1.rx_clock_m)
    }
  end

  defp integer_ratio(_best_score, nil), do: :infinity

  defp integer_ratio(best_score, second_best_score) do
    cond do
      best_score == 0.0 and second_best_score > 0.0 -> :infinity
      best_score == 0.0 -> 0.0
      true -> second_best_score / best_score
    end
  end

  defp integer_ratio_pass?(:infinity, _threshold), do: true
  defp integer_ratio_pass?(ratio, threshold), do: ratio >= threshold

  # --- dynamic least squares ----------------------------------------------

  defp solve_normal_equations(rows, n) do
    {ata, aty} =
      Enum.reduce(rows, {zero_matrix(n), zero_vector(n)}, fn row, {ata, aty} ->
        h = Enum.map(row.h, &(&1 * row.weight))
        y = row.y * row.weight
        {accumulate_ata(ata, h), accumulate_aty(aty, h, y)}
      end)

    solve_linear(ata, aty)
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

  # Gaussian elimination with partial pivoting, returning :singular for
  # rank-deficient normal equations.
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

  # --- result assembly -----------------------------------------------------

  defp finalize(sp3, obs, epoch, state, weights, iterations, converged, status) do
    with {:ok, residual_rows} <- residual_rows(sp3, obs, epoch, state) do
      code = residual_rows |> Enum.map(& &1.code_m)
      phase = residual_rows |> Enum.map(& &1.phase_m)

      {x, y, z} = state.position

      {:ok,
       %Solution{
         position: %{x_m: x, y_m: y, z_m: z},
         rx_clock_s: state.clock_m / Constants.speed_of_light_m_s(),
         rx_clock_m: state.clock_m,
         ambiguities_m: state.ambiguities,
         residuals_m:
           Map.new(residual_rows, fn r -> {r.sat, %{code_m: r.code_m, phase_m: r.phase_m}} end),
         used_sats: Enum.map(obs, & &1.satellite_id),
         metadata: %{
           iterations: iterations,
           converged: converged,
           status: status,
           code_rms_m: rms(code),
           phase_rms_m: rms(phase),
           weighted_rms_m: weighted_rms(residual_rows, weights)
         }
       }}
    end
  end

  defp finalize_multi(sp3, epochs, sat_ids, state, weights, iterations, converged, status) do
    with {:ok, residual_rows} <- multi_residual_rows(sp3, epochs, state) do
      code = residual_rows |> Enum.map(& &1.code_m)
      phase = residual_rows |> Enum.map(& &1.phase_m)
      {x, y, z} = state.position

      {:ok,
       %MultiEpochSolution{
         position: %{x_m: x, y_m: y, z_m: z},
         epoch_clocks:
           epochs
           |> Enum.map(& &1.epoch)
           |> Enum.zip(state.clocks_m)
           |> Enum.map(fn {epoch, clock_m} ->
             %{
               epoch: epoch,
               rx_clock_s: clock_m / Constants.speed_of_light_m_s(),
               rx_clock_m: clock_m
             }
           end),
         ambiguities_m: state.ambiguities,
         residuals_m: residual_rows,
         used_sats: sat_ids,
         epochs: Enum.map(epochs, & &1.epoch),
         metadata: %{
           iterations: iterations,
           converged: converged,
           status: status,
           n_epochs: length(epochs),
           n_observations: multi_observation_count(epochs),
           code_rms_m: rms(code),
           phase_rms_m: rms(phase),
           weighted_rms_m: weighted_rms(residual_rows, weights)
         }
       }}
    end
  end

  defp finalize_fixed_multi(
         sp3,
         epochs,
         sat_ids,
         fixed_cycles,
         fixed_m,
         float_sol,
         state,
         weights,
         iterations,
         converged,
         status,
         fixed_meta
       ) do
    with {:ok, residual_rows} <- fixed_multi_residual_rows(sp3, epochs, fixed_m, state) do
      code = residual_rows |> Enum.map(& &1.code_m)
      phase = residual_rows |> Enum.map(& &1.phase_m)
      {x, y, z} = state.position

      {:ok,
       %FixedSolution{
         position: %{x_m: x, y_m: y, z_m: z},
         epoch_clocks:
           epochs
           |> Enum.map(& &1.epoch)
           |> Enum.zip(state.clocks_m)
           |> Enum.map(fn {epoch, clock_m} ->
             %{
               epoch: epoch,
               rx_clock_s: clock_m / Constants.speed_of_light_m_s(),
               rx_clock_m: clock_m
             }
           end),
         fixed_ambiguities_cycles: fixed_cycles,
         fixed_ambiguities_m: fixed_m,
         float_solution: float_sol,
         residuals_m: residual_rows,
         used_sats: sat_ids,
         epochs: Enum.map(epochs, & &1.epoch),
         metadata:
           Map.merge(fixed_meta, %{
             iterations: iterations,
             converged: converged,
             status: status,
             n_epochs: length(epochs),
             n_observations: multi_observation_count(epochs),
             code_rms_m: rms(code),
             phase_rms_m: rms(phase),
             weighted_rms_m: weighted_rms(residual_rows, weights)
           })
       }}
    end
  end

  defp residual_rows(sp3, obs, epoch, state) do
    {rx, ry, rz} = state.position

    Enum.reduce_while(obs, {:ok, []}, fn o, {:ok, acc} ->
      case Observables.predict(sp3, o.satellite_id, {rx, ry, rz}, epoch) do
        {:ok, pred} ->
          sat_clock_m = Constants.speed_of_light_m_s() * (pred.sat_clock_s || 0.0)
          model_code = pred.geometric_range_m + state.clock_m - sat_clock_m
          ambiguity = Map.fetch!(state.ambiguities, o.satellite_id)

          row = %{
            sat: o.satellite_id,
            code_m: o.code_m - model_code,
            phase_m: o.phase_m - (model_code + ambiguity)
          }

          {:cont, {:ok, [row | acc]}}

        {:error, reason} ->
          {:halt, {:error, {:no_ephemeris, o.satellite_id, reason}}}
      end
    end)
    |> case do
      {:ok, rows} -> {:ok, Enum.reverse(rows)}
      {:error, _} = err -> err
    end
  end

  defp multi_residual_rows(sp3, epochs, state) do
    {rx, ry, rz} = state.position

    epochs
    |> Enum.with_index()
    |> Enum.reduce_while({:ok, []}, fn {epoch_row, epoch_idx}, {:ok, acc} ->
      clock_m = Enum.at(state.clocks_m, epoch_idx)

      epoch_row.observations
      |> Enum.reduce_while({:ok, acc}, fn o, {:ok, rows_acc} ->
        case Observables.predict(sp3, o.satellite_id, {rx, ry, rz}, epoch_row.epoch) do
          {:ok, pred} ->
            sat_clock_m = Constants.speed_of_light_m_s() * (pred.sat_clock_s || 0.0)
            model_code = pred.geometric_range_m + clock_m - sat_clock_m
            ambiguity = Map.fetch!(state.ambiguities, o.satellite_id)

            row = %{
              epoch: epoch_row.epoch,
              satellite_id: o.satellite_id,
              code_m: o.code_m - model_code,
              phase_m: o.phase_m - (model_code + ambiguity)
            }

            {:cont, {:ok, [row | rows_acc]}}

          {:error, reason} ->
            {:halt, {:error, {:no_ephemeris, o.satellite_id, reason}}}
        end
      end)
      |> case do
        {:ok, rows_acc} -> {:cont, {:ok, rows_acc}}
        {:error, _} = err -> {:halt, err}
      end
    end)
    |> case do
      {:ok, rows} -> {:ok, Enum.reverse(rows)}
      {:error, _} = err -> err
    end
  end

  defp fixed_multi_residual_rows(sp3, epochs, fixed_m, state) do
    {rx, ry, rz} = state.position

    epochs
    |> Enum.with_index()
    |> Enum.reduce_while({:ok, []}, fn {epoch_row, epoch_idx}, {:ok, acc} ->
      clock_m = Enum.at(state.clocks_m, epoch_idx)

      epoch_row.observations
      |> Enum.reduce_while({:ok, acc}, fn o, {:ok, rows_acc} ->
        case Observables.predict(sp3, o.satellite_id, {rx, ry, rz}, epoch_row.epoch) do
          {:ok, pred} ->
            sat_clock_m = Constants.speed_of_light_m_s() * (pred.sat_clock_s || 0.0)
            model_code = pred.geometric_range_m + clock_m - sat_clock_m
            ambiguity = Map.fetch!(fixed_m, o.satellite_id)

            row = %{
              epoch: epoch_row.epoch,
              satellite_id: o.satellite_id,
              code_m: o.code_m - model_code,
              phase_m: o.phase_m - (model_code + ambiguity)
            }

            {:cont, {:ok, [row | rows_acc]}}

          {:error, reason} ->
            {:halt, {:error, {:no_ephemeris, o.satellite_id, reason}}}
        end
      end)
      |> case do
        {:ok, rows_acc} -> {:cont, {:ok, rows_acc}}
        {:error, _} = err -> {:halt, err}
      end
    end)
    |> case do
      {:ok, rows} -> {:ok, Enum.reverse(rows)}
      {:error, _} = err -> err
    end
  end

  defp multi_satellite_ids(epochs) do
    epochs
    |> Enum.flat_map(& &1.observations)
    |> Enum.map(& &1.satellite_id)
    |> Enum.uniq()
    |> Enum.sort()
  end

  defp multi_observation_count(epochs) do
    Enum.reduce(epochs, 0, fn epoch, acc -> acc + length(epoch.observations) end)
  end

  defp rms([]), do: 0.0

  defp rms(xs) do
    :math.sqrt(Enum.reduce(xs, 0.0, fn x, acc -> acc + x * x end) / length(xs))
  end

  defp weighted_rms(rows, weights) do
    values =
      Enum.flat_map(rows, fn r ->
        [r.code_m * weights.code, r.phase_m * weights.phase]
      end)

    rms(values)
  end
end
