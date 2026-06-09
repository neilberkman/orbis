defmodule Orbis.GNSS.PrecisePositioning do
  @moduledoc """
  Carrier-phase precise-positioning primitives.

  This is the first precise-positioning layer above the code and carrier-phase
  combinations in `Orbis.GNSS.IonosphereFree` / `Orbis.GNSS.CarrierPhase`. It
  solves one SP3-backed epoch from dual-frequency ionosphere-free code and phase
  observations:

      P_IF_i = rho_i(x) + b - c * dt_sat_i + T_i
      L_IF_i = rho_i(x) + b - c * dt_sat_i + T_i + N_i

  where `x` is the receiver ECEF position, `b` is the receiver clock in metres,
  `T_i` is the optional a-priori slant tropospheric delay plus any estimated
  residual zenith delay mapped to the line of sight, and `N_i` is one float
  carrier-phase ambiguity per satellite, also in metres. The single-epoch state
  is linearized and iterated over `[x, y, z, b, N_1, N_2, ...]`.

  `solve_float/4` solves one epoch. `solve_float_epochs/3` solves a static
  multi-epoch arc with one receiver position, one receiver clock per epoch, and
  one ambiguity per satellite held constant across the arc. That multi-epoch
  model is the first step where carrier phase can tighten position instead of
  being absorbed entirely by one ambiguity per epoch. Multi-epoch and fixed
  solves can also estimate one residual zenith troposphere delay over the arc
  (`estimate_ztd: true`) after the a-priori Saastamoinen/Niell correction.

  `solve_fixed_epochs/3` starts from the same multi-epoch float model, searches
  integer ambiguity candidates on an explicit caller-supplied wavelength grid,
  then re-solves position and per-epoch clocks with those ambiguities held
  fixed. `solve_widelane_fixed_epochs/3` is the dual-frequency convenience path:
  it fixes the Melbourne-Wubbena wide-lane integer first, subtracts that known
  contribution from the ionosphere-free phase ambiguity, then runs the bounded
  integer least-squares search on the remaining narrow-lane integer.

  ## Observation shape

  Observations may be maps or tuples:

      %{satellite_id: "G05", code_m: 24_000_000.0, phase_m: 24_012_345.0}
      {"G05", 24_000_000.0, 24_012_345.0}

  `code_m` and `phase_m` should normally be ionosphere-free combinations. Use
  `Orbis.GNSS.IonosphereFree.iono_free/4` and
  `Orbis.GNSS.IonosphereFree.iono_free_phase_cycles/4` to form them from raw
  dual-frequency RINEX observations.
  """

  import Orbis.GNSS.Core.LinearAlgebra,
    only: [
      invert_matrix: 1,
      matmul: 2,
      matrix_sub: 2,
      normal_equations: 2,
      solve_matrix: 2,
      solve_normal_equations: 2,
      submatrix: 5,
      transpose: 1
    ]

  alias Orbis.Coordinates
  alias Orbis.GNSS.CarrierPhase
  alias Orbis.GNSS.Core.Constants
  alias Orbis.GNSS.Core.IntegerLeastSquares
  alias Orbis.GNSS.IonosphereFree
  alias Orbis.GNSS.Observables
  alias Orbis.GNSS.Positioning
  alias Orbis.GNSS.SP3
  alias Orbis.GNSS.Troposphere

  @default_max_iterations 8
  @default_position_tolerance_m 1.0e-4
  @default_clock_tolerance_m 1.0e-4
  @default_code_sigma_m 1.0
  @default_phase_sigma_m 0.01
  @default_pressure_hpa 1013.25
  @default_temperature_k 288.15
  @default_relative_humidity 0.5
  @default_ztd_tolerance_m 1.0e-4
  @default_integer_search_radius_cycles 1
  @default_integer_ratio_threshold 3.0
  @default_integer_candidate_limit 200_000
  @default_cycle_slip_policy :error
  @min_elevation_weight_scale 1.0e-3

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
              weighted_rms_m: float(),
              troposphere_applied: boolean()
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
      :ztd_residual_m,
      :residuals_m,
      :used_sats,
      :epochs,
      :metadata
    ]
    defstruct [
      :position,
      :epoch_clocks,
      :ambiguities_m,
      :ztd_residual_m,
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
            required(:epoch) => NaiveDateTime.t(),
            required(:satellite_id) => String.t(),
            required(:code_m) => float(),
            required(:phase_m) => float(),
            optional(:code_weight) => float(),
            optional(:phase_weight) => float()
          }

    @type t :: %__MODULE__{
            position: position(),
            epoch_clocks: [epoch_clock()],
            ambiguities_m: %{String.t() => float()},
            ztd_residual_m: float() | nil,
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
              troposphere_applied: boolean(),
              ztd_estimated: boolean()
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
      :wide_lane_ambiguities_cycles,
      :ztd_residual_m,
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
      :wide_lane_ambiguities_cycles,
      :ztd_residual_m,
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
            required(:epoch) => NaiveDateTime.t(),
            required(:satellite_id) => String.t(),
            required(:code_m) => float(),
            required(:phase_m) => float(),
            optional(:code_weight) => float(),
            optional(:phase_weight) => float()
          }

    @type t :: %__MODULE__{
            position: position(),
            epoch_clocks: [epoch_clock()],
            fixed_ambiguities_cycles: %{String.t() => integer()},
            fixed_ambiguities_m: %{String.t() => float()},
            wide_lane_ambiguities_cycles: %{String.t() => integer()} | nil,
            ztd_residual_m: float() | nil,
            float_solution: MultiEpochSolution.t(),
            residuals_m: [residual()],
            used_sats: [String.t()],
            epochs: [NaiveDateTime.t()],
            metadata: %{
              required(:iterations) => pos_integer(),
              required(:converged) => boolean(),
              required(:status) => :state_tolerance | :max_iterations,
              required(:n_epochs) => pos_integer(),
              required(:n_observations) => pos_integer(),
              required(:code_rms_m) => float(),
              required(:phase_rms_m) => float(),
              required(:weighted_rms_m) => float(),
              required(:integer_status) => :fixed | :not_fixed,
              required(:integer_method) => :lambda | :widelane_narrowlane_lambda,
              required(:integer_ratio) => float() | :infinity,
              required(:integer_best_score) => float(),
              required(:integer_second_best_score) => float() | nil,
              required(:integer_candidates) => pos_integer(),
              required(:troposphere_applied) => boolean(),
              required(:ztd_estimated) => boolean(),
              optional(:wide_lane_fixed) => boolean(),
              optional(:dropped_cycle_slip_sats) => [String.t()],
              optional(:split_cycle_slip_arcs) => [map()],
              optional(:ambiguity_search) => %{
                required(:order) => [String.t()],
                required(:float_cycles) => %{String.t() => float()},
                required(:covariance_cycles) => [[float()]],
                required(:covariance_inverse_cycles) => [[float()]]
              }
            }
          }
  end

  @typedoc "A dual-frequency ionosphere-free code/phase observation."
  @type observation ::
          %{satellite_id: String.t(), code_m: number(), phase_m: number()}
          | {String.t(), number(), number()}

  @typedoc "Raw dual-frequency code/phase observation for wide-lane/narrow-lane fixing."
  @type dual_frequency_observation :: %{
          required(:satellite_id) => String.t(),
          required(:p1_m) => number(),
          required(:p2_m) => number(),
          required(:phi1_cyc) => number(),
          required(:phi2_cyc) => number(),
          required(:f1_hz) => number(),
          required(:f2_hz) => number(),
          optional(:lli1) => integer() | nil,
          optional(:lli2) => integer() | nil
        }

  @typedoc "A receiver ECEF position in metres."
  @type receiver ::
          {number(), number(), number()} | %{x_m: number(), y_m: number(), z_m: number()}

  @typedoc "A set of code/phase observations for one epoch."
  @type epoch_observations ::
          %{epoch: NaiveDateTime.t(), observations: [observation()]}
          | {NaiveDateTime.t(), [observation()]}

  @typedoc "A set of raw dual-frequency observations for one epoch."
  @type dual_frequency_epoch_observations ::
          %{epoch: NaiveDateTime.t(), observations: [dual_frequency_observation()]}
          | {NaiveDateTime.t(), [dual_frequency_observation()]}

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
    * `:elevation_weighting` - when `true`, scale both code and phase row
      standard deviations by `1 / sin(elevation)` so low-elevation
      observations contribute less to the float, fixed, and ambiguity-covariance
      solves (default `false`).
    * `:max_iterations` - maximum nonlinear iterations (default `8`).
    * `:position_tolerance_m` - position-update convergence threshold
      (default `1.0e-4` m).
    * `:clock_tolerance_m` - receiver-clock update threshold (default
      `1.0e-4` m).
    * `:troposphere` - apply an a-priori Saastamoinen/Niell slant
      tropospheric delay to both code and phase (default `false`).
    * `:pressure_hpa` - surface pressure in hPa when `:troposphere` is true
      (default `1013.25`).
    * `:temperature_k` - surface temperature in kelvin when `:troposphere` is
      true (default `288.15`).
    * `:relative_humidity` - relative humidity fraction when `:troposphere` is
      true (default `0.5`).
    * `:estimate_ztd` - on multi-epoch/fixed solves only, estimate one residual
      zenith troposphere delay in metres over the whole static arc, mapped with
      the Niell wet mapping factor. Requires `troposphere: true` (default
      `false`).
    * `:ztd_tolerance_m` - residual-ZTD update convergence threshold when
      `:estimate_ztd` is true (default `1.0e-4` m).

  Returns `{:ok, %Solution{}}` or `{:error, reason}`. Reasons include
  `:no_observations`, `{:too_few_satellites, used, 4}`,
  `{:duplicate_observation, sat}`, `{:invalid_observation, entry}`,
  `:invalid_initial_guess`, `{:invalid_sigma, key}`, `{:invalid_option, key}`,
  `{:code_seed_failed, reason}`, `{:no_ephemeris, sat, reason}`,
  `{:troposphere_failed, sat, reason}`, and `:singular_geometry`. If the
  iteration limit is reached after a valid solve step, the function returns a
  solution with `metadata.converged == false` and
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
         {:ok, tropo} <- troposphere_options(opts),
         :ok <- ensure_single_epoch_troposphere(tropo),
         {:ok, state} <- initial_state(sp3, obs, epoch, opts) do
      iterate(sp3, obs, epoch, state, weights, tropo, solve_opts, 1)
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
         {:ok, tropo} <- troposphere_options(opts),
         :ok <- ensure_multi_enough(epochs, tropo),
         {:ok, weights} <- weights(opts),
         {:ok, solve_opts} <- solve_options(opts),
         {:ok, state} <- initial_multi_state(sp3, epochs, opts) do
      sat_ids = multi_satellite_ids(epochs)
      state = state_with_ztd(state, tropo)
      iterate_multi(sp3, epochs, sat_ids, state, weights, tropo, solve_opts, 1)
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

    * `:integer_ratio_threshold` - minimum second-best / best weighted-score
      ratio for `metadata.integer_status == :fixed` (default `3.0`).
    * `:integer_search_radius_cycles` / `:integer_candidate_limit` - retained and
      still validated for backward compatibility, but no longer bound the search:
      integer resolution uses the LAMBDA method (decorrelation + reduction +
      mlambda search), which finds the true integer-least-squares optimum for any
      geometry with no search box, so it cannot return
      `{:error, {:too_many_integer_candidates, ...}}`.
    * `:ambiguity_offset_m` - optional scalar or `%{"G05" => offset_m, ...}` map
      subtracted from each float ambiguity before converting to cycles and added
      back after fixing (default `0.0`). This is mainly for affine carrier-phase
      combinations such as wide-lane/narrow-lane fixing.

  The fixed solution is returned even when the ratio test is not met; in that
  case `metadata.integer_status` is `:not_fixed` so callers can reject it.
  """
  @spec solve_fixed_epochs(SP3.t(), [epoch_observations()], keyword()) ::
          {:ok, FixedSolution.t()} | {:error, term()}
  def solve_fixed_epochs(source, epoch_observations, opts \\ [])

  def solve_fixed_epochs(%SP3{} = sp3, epoch_observations, opts)
      when is_list(epoch_observations) do
    with {:ok, epochs} <- normalize_epoch_observations(epoch_observations),
         {:ok, tropo} <- troposphere_options(opts),
         :ok <- ensure_multi_enough(epochs, tropo),
         {:ok, weights} <- weights(opts),
         {:ok, solve_opts} <- solve_options(opts),
         {:ok, integer_opts} <- integer_options(opts),
         {:ok, float_sol} <- solve_float_epochs(sp3, epochs, opts),
         {:ok, wavelengths} <- ambiguity_wavelengths(float_sol.used_sats, opts),
         {:ok, offsets} <- ambiguity_offsets(float_sol.used_sats, opts),
         {:ok, fixed_cycles, fixed_meta} <-
           search_integer_ambiguities(
             sp3,
             epochs,
             float_sol,
             wavelengths,
             offsets,
             weights,
             tropo,
             integer_opts
           ),
         fixed_m = fixed_ambiguities_m(fixed_cycles, wavelengths, offsets),
         state = fixed_state_from_float(float_sol),
         {:ok, fixed_state, iterations, converged, status} <-
           iterate_fixed_multi(sp3, epochs, fixed_m, state, weights, tropo, solve_opts, 1) do
      finalize_fixed_multi(
        sp3,
        epochs,
        float_sol.used_sats,
        fixed_cycles,
        fixed_m,
        float_sol,
        fixed_state,
        weights,
        tropo,
        iterations,
        converged,
        status,
        fixed_meta
      )
    end
  end

  def solve_fixed_epochs(%SP3{}, _epoch_observations, _opts), do: {:error, :no_epochs}

  @doc """
  Solve a static multi-epoch position from raw dual-frequency observations by
  fixing wide-lane then narrow-lane ambiguities.

  This is the real-data convenience layer above `solve_fixed_epochs/3`. Each
  observation must carry both code and carrier phase on two bands:

      %{
        satellite_id: "G05",
        p1_m: 24_000_000.0,
        p2_m: 24_000_004.0,
        phi1_cyc: 123_456_789.0,
        phi2_cyc: 96_123_456.0,
        f1_hz: 1_575_420_000.0,
        f2_hz: 1_227_600_000.0,
        lli1: 0,
        lli2: 0
      }

  For each satellite the function first estimates the Melbourne-Wubbena
  wide-lane integer `Nw = N1 - N2` over the arc. It then forms ionosphere-free
  code/phase observations and fixes the remaining band-1 narrow-lane integer
  with bounded integer least-squares using `lambda_NL = c / (f1 + f2)`. The returned
  `fixed_ambiguities_cycles` are those band-1 narrow-lane integers; the
  wide-lane integers are exposed as `wide_lane_ambiguities_cycles`.

  ## Options

  Accepts the same solve and integer-search options as `solve_fixed_epochs/3`,
  plus:

    * `:wide_lane_min_epochs` - minimum usable Melbourne-Wubbena epochs per
      satellite (default `2`).
    * `:wide_lane_tolerance_cycles` - maximum absolute distance between the
      averaged wide-lane float value and the nearest integer (default `0.5`
      cycles).
    * `:on_cycle_slip` - what to do when a satellite arc has a detected cycle
      slip: `:error` returns `{:error, {:cycle_slip_detected, sat, epoch,
      reasons}}` (default); `:drop_satellite` removes that satellite from the
      wide-lane and narrow-lane solve; `:split_arc` resets that satellite's
      ambiguity at each slip and keeps any resulting arc with at least
      `:wide_lane_min_epochs` usable epochs. Dropped satellites are reported in
      `metadata.dropped_cycle_slip_sats`; split fragments are reported in
      `metadata.split_cycle_slip_arcs`. Split fragments use ambiguity ids such
      as `"G21#2"` in `used_sats` and the ambiguity maps, while ephemeris lookup
      and residual rows continue to use the physical satellite id (`"G21"`).

  Cycle slips are detected with `Orbis.GNSS.CarrierPhase.detect_cycle_slips/2`;
  pass `:gf_threshold_m` / `:mw_threshold_cycles` to tune that detector.
  """
  @spec solve_widelane_fixed_epochs(SP3.t(), [dual_frequency_epoch_observations()], keyword()) ::
          {:ok, FixedSolution.t()} | {:error, term()}
  def solve_widelane_fixed_epochs(source, dual_epoch_observations, opts \\ [])

  def solve_widelane_fixed_epochs(%SP3{} = sp3, dual_epoch_observations, opts)
      when is_list(dual_epoch_observations) do
    with {:ok, dual_epochs} <- normalize_dual_epoch_observations(dual_epoch_observations),
         {:ok, prepared_dual_epochs, wide_lane_cycles, slip_meta} <-
           wide_lane_ambiguities(dual_epochs, opts),
         filtered_dual_epochs =
           filter_dual_epochs_by_wide_lanes(prepared_dual_epochs, wide_lane_cycles),
         {:ok, if_epochs, wavelengths, offsets} <-
           ionosphere_free_narrow_lane_epochs(filtered_dual_epochs, wide_lane_cycles),
         fixed_opts =
           opts
           |> Keyword.put(:ambiguity_wavelength_m, wavelengths)
           |> Keyword.put(:ambiguity_offset_m, offsets),
         {:ok, %FixedSolution{} = sol} <- solve_fixed_epochs(sp3, if_epochs, fixed_opts) do
      {:ok,
       %{
         sol
         | wide_lane_ambiguities_cycles: wide_lane_cycles,
           metadata:
             Map.merge(sol.metadata, %{
               integer_method: :widelane_narrowlane_lambda,
               wide_lane_fixed: true,
               dropped_cycle_slip_sats: slip_meta.dropped_sats,
               split_cycle_slip_arcs: slip_meta.split_arcs
             })
       }}
    end
  end

  def solve_widelane_fixed_epochs(%SP3{}, _dual_epoch_observations, _opts),
    do: {:error, :no_epochs}

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
    case normalize_ambiguity_id(obs, sat) do
      {:ok, ambiguity_id} ->
        {:ok,
         %{
           satellite_id: sat,
           ambiguity_id: ambiguity_id,
           code_m: code / 1.0,
           phase_m: phase / 1.0,
           raw: obs
         }}

      {:error, _} = err ->
        err
    end
  end

  defp normalize_one({sat, code, phase} = obs)
       when is_binary(sat) and is_number(code) and is_number(phase) do
    {:ok,
     %{satellite_id: sat, ambiguity_id: sat, code_m: code / 1.0, phase_m: phase / 1.0, raw: obs}}
  end

  defp normalize_one(entry), do: {:error, {:invalid_observation, entry}}

  defp normalize_ambiguity_id(obs, sat) do
    case Map.get(obs, :ambiguity_id, sat) do
      ambiguity_id when is_binary(ambiguity_id) -> {:ok, ambiguity_id}
      _other -> {:error, {:invalid_observation, obs}}
    end
  end

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

  defp normalize_dual_epoch_observations([]), do: {:error, :no_epochs}

  defp normalize_dual_epoch_observations(epoch_observations) do
    epoch_observations
    |> Enum.reduce_while({:ok, [], MapSet.new()}, fn entry, {:ok, acc, seen} ->
      case normalize_dual_epoch_entry(entry) do
        {:ok, epoch, observations} ->
          if MapSet.member?(seen, epoch) do
            {:halt, {:error, {:duplicate_epoch, epoch}}}
          else
            with {:ok, obs} <- normalize_dual_observations(observations),
                 :ok <- ensure_dual_epoch_enough(epoch, obs) do
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

  defp normalize_dual_epoch_entry(%{epoch: %NaiveDateTime{} = epoch, observations: observations})
       when is_list(observations), do: {:ok, epoch, observations}

  defp normalize_dual_epoch_entry({%NaiveDateTime{} = epoch, observations})
       when is_list(observations), do: {:ok, epoch, observations}

  defp normalize_dual_epoch_entry(entry), do: {:error, {:invalid_epoch_observations, entry}}

  defp normalize_dual_observations(observations) do
    observations
    |> Enum.reduce_while({:ok, [], MapSet.new()}, fn entry, {:ok, acc, seen} ->
      case normalize_dual_one(entry) do
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
        {:ok, acc |> Enum.reverse() |> Enum.sort_by(& &1.satellite_id)}

      {:error, _} = err ->
        err
    end
  end

  defp normalize_dual_one(
         %{
           satellite_id: sat,
           p1_m: p1,
           p2_m: p2,
           phi1_cyc: phi1,
           phi2_cyc: phi2,
           f1_hz: f1,
           f2_hz: f2
         } = obs
       )
       when is_binary(sat) and is_number(p1) and is_number(p2) and is_number(phi1) and
              is_number(phi2) and is_number(f1) and is_number(f2) and f1 > 0.0 and f2 > 0.0 do
    {:ok,
     %{
       satellite_id: sat,
       ambiguity_id: sat,
       p1_m: p1 / 1.0,
       p2_m: p2 / 1.0,
       phi1_cyc: phi1 / 1.0,
       phi2_cyc: phi2 / 1.0,
       f1_hz: f1 / 1.0,
       f2_hz: f2 / 1.0,
       lli1: Map.get(obs, :lli1),
       lli2: Map.get(obs, :lli2),
       raw: obs
     }}
  end

  defp normalize_dual_one(entry), do: {:error, {:invalid_dual_frequency_observation, entry}}

  defp ensure_dual_epoch_enough(_epoch, obs) when length(obs) >= 4, do: :ok

  defp ensure_dual_epoch_enough(epoch, obs),
    do: {:error, {:too_few_epoch_observations, epoch, length(obs), 4}}

  defp ensure_multi_enough(epochs, _tropo) when length(epochs) < 2,
    do: {:error, {:too_few_epochs, length(epochs), 2}}

  defp ensure_multi_enough(epochs, tropo) do
    n_epochs = length(epochs)
    n_sats = length(multi_satellite_ids(epochs))
    n_observations = multi_observation_count(epochs)
    equations = 2 * n_observations
    unknowns = 3 + n_epochs + ztd_unknown_count(tropo) + n_sats

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
    elevation_weighting = Keyword.get(opts, :elevation_weighting, false)

    cond do
      not is_number(code_sigma) or code_sigma <= 0.0 ->
        {:error, {:invalid_sigma, :code_sigma_m}}

      not is_number(phase_sigma) or phase_sigma <= 0.0 ->
        {:error, {:invalid_sigma, :phase_sigma_m}}

      elevation_weighting not in [true, false] ->
        {:error, {:invalid_option, :elevation_weighting}}

      true ->
        {:ok,
         %{
           code: 1.0 / code_sigma,
           phase: 1.0 / phase_sigma,
           elevation_weighting?: elevation_weighting
         }}
    end
  end

  defp measurement_weight(weights, kind, elevation_deg) do
    base = Map.fetch!(weights, kind)

    if weights.elevation_weighting? do
      base * elevation_weight_scale(elevation_deg)
    else
      base
    end
  end

  defp elevation_weight_scale(elevation_deg) do
    sin_el = :math.sin(elevation_deg * :math.pi() / 180.0)

    cond do
      not is_number(sin_el) -> @min_elevation_weight_scale
      sin_el < @min_elevation_weight_scale -> @min_elevation_weight_scale
      true -> sin_el
    end
  end

  defp solve_options(opts) do
    max_iterations = Keyword.get(opts, :max_iterations, @default_max_iterations)
    pos_tol = Keyword.get(opts, :position_tolerance_m, @default_position_tolerance_m)
    clock_tol = Keyword.get(opts, :clock_tolerance_m, @default_clock_tolerance_m)
    ambiguity_tol = Keyword.get(opts, :ambiguity_tolerance_m, @default_position_tolerance_m)
    ztd_tol = Keyword.get(opts, :ztd_tolerance_m, @default_ztd_tolerance_m)

    cond do
      not is_integer(max_iterations) or max_iterations < 1 ->
        {:error, {:invalid_option, :max_iterations}}

      not is_number(pos_tol) or pos_tol < 0.0 ->
        {:error, {:invalid_option, :position_tolerance_m}}

      not is_number(clock_tol) or clock_tol < 0.0 ->
        {:error, {:invalid_option, :clock_tolerance_m}}

      not is_number(ambiguity_tol) or ambiguity_tol < 0.0 ->
        {:error, {:invalid_option, :ambiguity_tolerance_m}}

      not is_number(ztd_tol) or ztd_tol < 0.0 ->
        {:error, {:invalid_option, :ztd_tolerance_m}}

      true ->
        {:ok,
         %{
           max_iterations: max_iterations,
           position_tolerance_m: pos_tol / 1.0,
           clock_tolerance_m: clock_tol / 1.0,
           ambiguity_tolerance_m: ambiguity_tol / 1.0,
           ztd_tolerance_m: ztd_tol / 1.0
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

  defp troposphere_options(opts) do
    estimate_ztd = Keyword.get(opts, :estimate_ztd, false)

    case Keyword.get(opts, :troposphere, false) do
      false ->
        if estimate_ztd == false do
          {:ok, %{enabled?: false, met: nil, estimate_ztd?: false}}
        else
          {:error, {:invalid_option, :estimate_ztd}}
        end

      true ->
        pressure = Keyword.get(opts, :pressure_hpa, @default_pressure_hpa)
        temperature = Keyword.get(opts, :temperature_k, @default_temperature_k)
        humidity = Keyword.get(opts, :relative_humidity, @default_relative_humidity)

        cond do
          not is_number(pressure) or pressure <= 0.0 ->
            {:error, {:invalid_option, :pressure_hpa}}

          not is_number(temperature) or temperature <= 0.0 ->
            {:error, {:invalid_option, :temperature_k}}

          not is_number(humidity) or humidity < 0.0 or humidity > 1.0 ->
            {:error, {:invalid_option, :relative_humidity}}

          estimate_ztd not in [true, false] ->
            {:error, {:invalid_option, :estimate_ztd}}

          true ->
            {:ok,
             %{
               enabled?: true,
               estimate_ztd?: estimate_ztd,
               met: %{
                 pressure_hpa: pressure / 1.0,
                 temperature_k: temperature / 1.0,
                 relative_humidity: humidity / 1.0
               }
             }}
        end

      _other ->
        {:error, {:invalid_option, :troposphere}}
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

  defp ambiguity_offsets(sat_ids, opts) do
    case Keyword.fetch(opts, :ambiguity_offset_m) do
      {:ok, offset} when is_number(offset) ->
        {:ok, Map.new(sat_ids, &{&1, offset / 1.0})}

      {:ok, offset_by_sat} when is_map(offset_by_sat) ->
        sat_ids
        |> Enum.reduce_while({:ok, %{}}, fn sat, {:ok, acc} ->
          case Map.fetch(offset_by_sat, sat) do
            {:ok, value} when is_number(value) ->
              {:cont, {:ok, Map.put(acc, sat, value / 1.0)}}

            _ ->
              {:halt, {:error, {:invalid_ambiguity_offset, sat}}}
          end
        end)

      {:ok, _other} ->
        {:error, {:invalid_option, :ambiguity_offset_m}}

      :error ->
        {:ok, Map.new(sat_ids, &{&1, 0.0})}
    end
  end

  defp wide_lane_ambiguities(epochs, opts) do
    with {:ok, wl_opts} <- wide_lane_options(opts),
         {:ok, slip_policy} <- cycle_slip_policy(opts) do
      arcs =
        epochs
        |> Enum.flat_map(fn epoch_row ->
          Enum.map(epoch_row.observations, &{&1.satellite_id, epoch_row.epoch, &1})
        end)
        |> Enum.group_by(fn {sat, _epoch, _obs} -> sat end)

      arcs
      |> Enum.sort_by(fn {sat, _arc} -> sat end)
      |> Enum.reduce_while({:ok, [], %{}, [], []}, fn {sat, arc},
                                                      {:ok, entries, cycles, dropped, split_arcs} ->
        arc = Enum.sort_by(arc, fn {_sat, epoch, _obs} -> NaiveDateTime.to_iso8601(epoch) end)

        case prepare_wide_lane_arc(sat, arc, opts, wl_opts, slip_policy) do
          {:ok, arc_entries, arc_cycles, arc_dropped, arc_splits} ->
            {:cont,
             {:ok, arc_entries ++ entries, Map.merge(cycles, arc_cycles), arc_dropped ++ dropped,
              arc_splits ++ split_arcs}}

          {:error, _} = err ->
            {:halt, err}
        end
      end)
      |> case do
        {:ok, entries, cycles, dropped, split_arcs} ->
          {:ok, dual_epochs_from_entries(entries), cycles,
           %{
             dropped_sats: dropped |> Enum.uniq() |> Enum.sort(),
             split_arcs: split_arcs |> Enum.sort_by(&{&1.satellite_id, &1.ambiguity_id})
           }}

        {:error, _} = err ->
          err
      end
    end
  end

  defp wide_lane_options(opts) do
    min_epochs = Keyword.get(opts, :wide_lane_min_epochs, 2)
    tolerance = Keyword.get(opts, :wide_lane_tolerance_cycles, 0.5)

    cond do
      not is_integer(min_epochs) or min_epochs < 1 ->
        {:error, {:invalid_option, :wide_lane_min_epochs}}

      not is_number(tolerance) or tolerance < 0.0 ->
        {:error, {:invalid_option, :wide_lane_tolerance_cycles}}

      true ->
        {:ok, %{min_epochs: min_epochs, tolerance_cycles: tolerance / 1.0}}
    end
  end

  defp cycle_slip_policy(opts) do
    case Keyword.get(opts, :on_cycle_slip, @default_cycle_slip_policy) do
      :error -> {:ok, :error}
      :drop_satellite -> {:ok, :drop_satellite}
      :split_arc -> {:ok, :split_arc}
      _other -> {:error, {:invalid_option, :on_cycle_slip}}
    end
  end

  defp prepare_wide_lane_arc(sat, arc, opts, wl_opts, :split_arc) do
    case cycle_slips_for_arc(arc, opts) do
      [] ->
        estimate_tagged_wide_lane_arc(sat, sat, arc, wl_opts, [])

      slips ->
        segments = split_wide_lane_arc(arc, MapSet.new(Enum.map(slips, & &1.epoch)))

        segments
        |> Enum.reduce_while({:ok, [], %{}, [], []}, fn {segment_idx, segment},
                                                        {:ok, entries, cycles, dropped,
                                                         split_arcs} ->
          if length(segment) < wl_opts.min_epochs do
            {:cont, {:ok, entries, cycles, dropped, split_arcs}}
          else
            ambiguity_id = split_ambiguity_id(sat, segment_idx)

            case estimate_tagged_wide_lane_arc(ambiguity_id, ambiguity_id, segment, wl_opts, [
                   split_arc_metadata(sat, ambiguity_id, segment)
                 ]) do
              {:ok, arc_entries, arc_cycles, _arc_dropped, arc_splits} ->
                {:cont,
                 {:ok, arc_entries ++ entries, Map.merge(cycles, arc_cycles), dropped,
                  arc_splits ++ split_arcs}}

              {:error, _} = err ->
                {:halt, err}
            end
          end
        end)
        |> case do
          {:ok, entries, cycles, dropped, split_arcs} when map_size(cycles) > 0 ->
            {:ok, entries, cycles, dropped, split_arcs}

          {:ok, _entries, _cycles, dropped, split_arcs} ->
            {:ok, [], %{}, [sat | dropped], split_arcs}

          {:error, _} = err ->
            err
        end
    end
  end

  defp prepare_wide_lane_arc(sat, arc, opts, wl_opts, slip_policy) do
    case cycle_slips_for_arc(arc, opts) do
      [] ->
        estimate_tagged_wide_lane_arc(sat, sat, arc, wl_opts, [])

      [_slip | _] when slip_policy == :drop_satellite ->
        {:ok, [], %{}, [sat], []}

      [slip | _] ->
        {:error, {:cycle_slip_detected, sat, slip.epoch, slip.reasons}}
    end
  end

  defp estimate_tagged_wide_lane_arc(error_id, ambiguity_id, arc, wl_opts, split_arcs) do
    case estimate_wide_lane_integer(error_id, arc, wl_opts) do
      {:ok, fixed} ->
        entries =
          Enum.map(arc, fn {_sat, epoch, obs} ->
            {epoch, Map.put(obs, :ambiguity_id, ambiguity_id)}
          end)

        {:ok, entries, %{ambiguity_id => fixed}, [], split_arcs}

      {:error, _} = err ->
        err
    end
  end

  defp cycle_slips_for_arc(arc, opts) do
    arc
    |> carrier_phase_arc()
    |> CarrierPhase.detect_cycle_slips(opts)
    |> Enum.filter(& &1.slip)
  end

  defp carrier_phase_arc(arc) do
    Enum.map(arc, fn {_sat, epoch, obs} ->
      %{
        epoch: epoch,
        phi1: obs.phi1_cyc,
        phi2: obs.phi2_cyc,
        p1: obs.p1_m,
        p2: obs.p2_m,
        f1: obs.f1_hz,
        f2: obs.f2_hz,
        lli1: obs.lli1,
        lli2: obs.lli2
      }
    end)
  end

  defp split_wide_lane_arc(arc, slip_epochs) do
    {segments, current, current_idx} =
      Enum.reduce(arc, {[], [], 1}, fn {_sat, epoch, _obs} = row, {segments, current, idx} ->
        if MapSet.member?(slip_epochs, epoch) do
          segments =
            if current == [], do: segments, else: [{idx, Enum.reverse(current)} | segments]

          {segments, [row], idx + 1}
        else
          {segments, [row | current], idx}
        end
      end)

    segments =
      if current == [], do: segments, else: [{current_idx, Enum.reverse(current)} | segments]

    Enum.reverse(segments)
  end

  defp split_ambiguity_id(sat, segment_idx), do: "#{sat}##{segment_idx}"

  defp split_arc_metadata(sat, ambiguity_id, segment) do
    epochs = Enum.map(segment, fn {_sat, epoch, _obs} -> epoch end)

    %{
      satellite_id: sat,
      ambiguity_id: ambiguity_id,
      start_epoch: List.first(epochs),
      end_epoch: List.last(epochs),
      n_epochs: length(epochs)
    }
  end

  defp dual_epochs_from_entries(entries) do
    entries
    |> Enum.group_by(fn {epoch, _obs} -> epoch end, fn {_epoch, obs} -> obs end)
    |> Enum.map(fn {epoch, observations} ->
      %{
        epoch: epoch,
        observations: Enum.sort_by(observations, &{&1.satellite_id, ambiguity_id(&1)})
      }
    end)
    |> Enum.sort_by(&NaiveDateTime.to_iso8601(&1.epoch))
  end

  defp ambiguity_id(obs), do: Map.get(obs, :ambiguity_id, obs.satellite_id)

  defp estimate_wide_lane_integer(sat, arc, opts) do
    arc
    |> Enum.reduce_while({:ok, []}, fn {_sat, _epoch, obs}, {:ok, acc} ->
      case wide_lane_cycles(obs) do
        {:ok, value} -> {:cont, {:ok, [value | acc]}}
        {:error, reason} -> {:halt, {:error, {:wide_lane_failed, sat, reason}}}
      end
    end)
    |> case do
      {:ok, cycles} ->
        if length(cycles) < opts.min_epochs do
          {:error, {:too_few_wide_lane_epochs, sat, length(cycles), opts.min_epochs}}
        else
          mean = Enum.sum(cycles) / length(cycles)
          fixed = round(mean)

          if abs(mean - fixed) <= opts.tolerance_cycles do
            {:ok, fixed}
          else
            {:error, {:wide_lane_not_integer, sat, mean, fixed}}
          end
        end

      {:error, _} = err ->
        err
    end
  end

  defp wide_lane_cycles(obs) do
    with {:ok, mw_m} <-
           CarrierPhase.melbourne_wubbena(
             obs.phi1_cyc,
             obs.phi2_cyc,
             obs.p1_m,
             obs.p2_m,
             obs.f1_hz,
             obs.f2_hz
           ),
         {:ok, lambda_wl} <- CarrierPhase.wide_lane_wavelength(obs.f1_hz, obs.f2_hz) do
      {:ok, mw_m / lambda_wl}
    end
  end

  defp filter_dual_epochs_by_wide_lanes(dual_epochs, wide_lane_cycles) do
    keep = MapSet.new(Map.keys(wide_lane_cycles))

    dual_epochs
    |> Enum.map(fn epoch_row ->
      observations =
        Enum.filter(epoch_row.observations, &MapSet.member?(keep, ambiguity_id(&1)))

      %{epoch_row | observations: observations}
    end)
    |> Enum.reject(&(&1.observations == []))
  end

  defp ionosphere_free_narrow_lane_epochs(dual_epochs, wide_lane_cycles) do
    with {:ok, params} <- narrow_lane_params(dual_epochs, wide_lane_cycles),
         {:ok, if_epochs} <- ionosphere_free_epochs(dual_epochs) do
      wavelengths = Map.new(params, fn {sat, p} -> {sat, p.wavelength_m} end)
      offsets = Map.new(params, fn {sat, p} -> {sat, p.offset_m} end)
      {:ok, if_epochs, wavelengths, offsets}
    end
  end

  defp narrow_lane_params(dual_epochs, wide_lane_cycles) do
    dual_epochs
    |> Enum.flat_map(& &1.observations)
    |> Enum.reduce_while({:ok, %{}}, fn obs, {:ok, acc} ->
      sat = ambiguity_id(obs)

      with {:ok, wide_lane} <- fetch_wide_lane(wide_lane_cycles, sat),
           {:ok, params} <- narrow_lane_param(obs.f1_hz, obs.f2_hz, wide_lane),
           :ok <- ensure_consistent_narrow_lane_params(sat, params, Map.get(acc, sat)) do
        {:cont, {:ok, Map.put_new(acc, sat, params)}}
      else
        {:error, _} = err -> {:halt, err}
      end
    end)
  end

  defp fetch_wide_lane(wide_lane_cycles, sat) do
    case Map.fetch(wide_lane_cycles, sat) do
      {:ok, wide_lane} -> {:ok, wide_lane}
      :error -> {:error, {:missing_wide_lane_ambiguity, sat}}
    end
  end

  defp narrow_lane_param(f1, f2, wide_lane_cycles) do
    with {:ok, gamma} <- IonosphereFree.gamma(f1, f2) do
      c = Constants.speed_of_light_m_s()
      beta = gamma - 1.0
      lambda2 = c / f2

      {:ok,
       %{
         wavelength_m: c / (f1 + f2),
         offset_m: beta * lambda2 * wide_lane_cycles,
         f1_hz: f1,
         f2_hz: f2
       }}
    end
  end

  defp ensure_consistent_narrow_lane_params(_sat, _params, nil), do: :ok

  defp ensure_consistent_narrow_lane_params(sat, params, prev) do
    if same_frequency?(params.f1_hz, prev.f1_hz) and same_frequency?(params.f2_hz, prev.f2_hz) do
      :ok
    else
      {:error, {:inconsistent_frequencies, sat}}
    end
  end

  defp same_frequency?(a, b), do: abs(a - b) <= 1.0e-6

  defp ionosphere_free_epochs(dual_epochs) do
    dual_epochs
    |> Enum.reduce_while({:ok, []}, fn epoch_row, {:ok, acc} ->
      case ionosphere_free_observations(epoch_row.observations) do
        {:ok, observations} ->
          {:cont, {:ok, [%{epoch: epoch_row.epoch, observations: observations} | acc]}}

        {:error, _} = err ->
          {:halt, err}
      end
    end)
    |> case do
      {:ok, rows} -> {:ok, Enum.reverse(rows)}
      {:error, _} = err -> err
    end
  end

  defp ionosphere_free_observations(observations) do
    observations
    |> Enum.reduce_while({:ok, []}, fn obs, {:ok, acc} ->
      with {:ok, code_m} <- IonosphereFree.iono_free(obs.p1_m, obs.p2_m, obs.f1_hz, obs.f2_hz),
           {:ok, phase_m} <-
             IonosphereFree.iono_free_phase_cycles(
               obs.phi1_cyc,
               obs.phi2_cyc,
               obs.f1_hz,
               obs.f2_hz
             ) do
        {:cont,
         {:ok,
          [
            %{
              satellite_id: obs.satellite_id,
              ambiguity_id: ambiguity_id(obs),
              code_m: code_m,
              phase_m: phase_m
            }
            | acc
          ]}}
      else
        {:error, reason} -> {:halt, {:error, {:ionosphere_free_failed, obs.satellite_id, reason}}}
      end
    end)
    |> case do
      {:ok, obs} -> {:ok, Enum.reverse(obs)}
      {:error, _} = err -> err
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
        {ambiguity_id(o), o.phase_m - o.code_m}
      end)

    %{position: position, clock_m: clock_m, ambiguities: ambiguities}
  end

  defp spp_seed(sp3, obs, epoch, opts) do
    observations = Enum.map(obs, &{&1.satellite_id, &1.code_m})
    spp_initial = Keyword.get(opts, :spp_initial_guess, {0.0, 0.0, 0.0, 0.0})

    case Positioning.solve(sp3, observations, epoch, spp_seed_options(opts, spp_initial)) do
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
      Map.put_new(acc, ambiguity_id(obs), obs.phase_m - obs.code_m)
    end)
  end

  defp multi_spp_seed(sp3, epochs, opts) do
    epochs
    |> Enum.reduce_while({:ok, [], []}, fn epoch_row, {:ok, positions, clocks} ->
      observations = Enum.map(epoch_row.observations, &{&1.satellite_id, &1.code_m})
      spp_initial = Keyword.get(opts, :spp_initial_guess, {0.0, 0.0, 0.0, 0.0})

      case Positioning.solve(
             sp3,
             observations,
             epoch_row.epoch,
             spp_seed_options(opts, spp_initial)
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

  defp spp_seed_options(opts, initial_guess) do
    [
      ionosphere: false,
      troposphere: Keyword.get(opts, :troposphere, false),
      pressure_hpa: Keyword.get(opts, :pressure_hpa, @default_pressure_hpa),
      temperature_k: Keyword.get(opts, :temperature_k, @default_temperature_k),
      relative_humidity: Keyword.get(opts, :relative_humidity, @default_relative_humidity),
      initial_guess: initial_guess,
      with_geodetic: false
    ]
  end

  defp mean_position(positions) do
    {sx, sy, sz} =
      Enum.reduce(positions, {0.0, 0.0, 0.0}, fn {x, y, z}, {ax, ay, az} ->
        {ax + x, ay + y, az + z}
      end)

    n = length(positions)
    {sx / n, sy / n, sz / n}
  end

  defp state_with_ztd(state, %{estimate_ztd?: true}), do: Map.put_new(state, :ztd_m, 0.0)
  defp state_with_ztd(state, _tropo), do: state

  defp ztd_unknown_count(%{estimate_ztd?: true}), do: 1
  defp ztd_unknown_count(_tropo), do: 0

  defp ensure_single_epoch_troposphere(%{estimate_ztd?: true}),
    do: {:error, {:invalid_option, :estimate_ztd}}

  defp ensure_single_epoch_troposphere(_tropo), do: :ok

  # --- nonlinear solve -----------------------------------------------------

  defp iterate(sp3, obs, epoch, state, weights, tropo, opts, iter) do
    with {:ok, rows} <- build_rows(sp3, obs, epoch, state, weights, tropo),
         {:ok, dx} <- solve_normal_equations(rows, 4 + length(obs)) do
      next = apply_delta(state, obs, dx)
      {pos_step, clock_step} = step_norms(dx)

      cond do
        pos_step <= opts.position_tolerance_m and abs(clock_step) <= opts.clock_tolerance_m ->
          finalize(sp3, obs, epoch, next, weights, tropo, iter, true, :position_tolerance)

        iter >= opts.max_iterations ->
          finalize(sp3, obs, epoch, next, weights, tropo, iter, false, :max_iterations)

        true ->
          iterate(sp3, obs, epoch, next, weights, tropo, opts, iter + 1)
      end
    end
  end

  defp iterate_multi(sp3, epochs, sat_ids, state, weights, tropo, opts, iter) do
    with {:ok, rows} <- build_multi_rows(sp3, epochs, sat_ids, state, weights, tropo),
         {:ok, dx} <-
           solve_normal_equations(
             rows,
             3 + length(epochs) + ztd_unknown_count(tropo) + length(sat_ids)
           ) do
      next = apply_multi_delta(state, epochs, sat_ids, dx, tropo)

      {pos_step, clock_step, ztd_step, ambiguity_step} =
        multi_step_norms(dx, length(epochs), tropo)

      cond do
        pos_step <= opts.position_tolerance_m and
          clock_step <= opts.clock_tolerance_m and
          ztd_step <= opts.ztd_tolerance_m and
            ambiguity_step <= opts.ambiguity_tolerance_m ->
          finalize_multi(
            sp3,
            epochs,
            sat_ids,
            next,
            weights,
            tropo,
            iter,
            true,
            :state_tolerance
          )

        iter >= opts.max_iterations ->
          finalize_multi(sp3, epochs, sat_ids, next, weights, tropo, iter, false, :max_iterations)

        true ->
          iterate_multi(sp3, epochs, sat_ids, next, weights, tropo, opts, iter + 1)
      end
    end
  end

  defp iterate_fixed_multi(sp3, epochs, fixed_m, state, weights, tropo, opts, iter) do
    with {:ok, rows} <- build_fixed_multi_rows(sp3, epochs, fixed_m, state, weights, tropo),
         {:ok, dx} <- solve_normal_equations(rows, 3 + length(epochs) + ztd_unknown_count(tropo)) do
      next = apply_fixed_multi_delta(state, dx, length(epochs), tropo)
      {pos_step, clock_step, ztd_step} = fixed_multi_step_norms(dx, tropo)

      cond do
        pos_step <= opts.position_tolerance_m and
          clock_step <= opts.clock_tolerance_m and
            ztd_step <= opts.ztd_tolerance_m ->
          {:ok, next, iter, true, :state_tolerance}

        iter >= opts.max_iterations ->
          {:ok, next, iter, false, :max_iterations}

        true ->
          iterate_fixed_multi(sp3, epochs, fixed_m, next, weights, tropo, opts, iter + 1)
      end
    end
  end

  defp build_rows(sp3, obs, epoch, state, weights, tropo) do
    {rx, ry, rz} = state.position

    Enum.reduce_while(obs, {:ok, []}, fn o, {:ok, acc} ->
      case Observables.predict(sp3, o.satellite_id, {rx, ry, rz}, epoch) do
        {:ok, pred} ->
          case model_troposphere(pred, {rx, ry, rz}, epoch, tropo) do
            {:ok, tropo_model} ->
              {ex, ey, ez} = pred.los_unit
              sat_clock_m = Constants.speed_of_light_m_s() * (pred.sat_clock_s || 0.0)
              tropo_m = applied_troposphere_m(tropo_model, state)
              model_code = pred.geometric_range_m + state.clock_m - sat_clock_m + tropo_m
              arc_id = ambiguity_id(o)
              ambiguity = Map.fetch!(state.ambiguities, arc_id)

              code_row = %{
                kind: :code,
                sat: o.satellite_id,
                h: design_row({-ex, -ey, -ez, 1.0}, nil, obs),
                y: o.code_m - model_code,
                weight: measurement_weight(weights, :code, pred.elevation_deg)
              }

              phase_row = %{
                kind: :phase,
                sat: o.satellite_id,
                h: design_row({-ex, -ey, -ez, 1.0}, arc_id, obs),
                y: o.phase_m - (model_code + ambiguity),
                weight: measurement_weight(weights, :phase, pred.elevation_deg)
              }

              {:cont, {:ok, [phase_row, code_row | acc]}}

            {:error, reason} ->
              {:halt, {:error, {:troposphere_failed, o.satellite_id, reason}}}
          end

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
        if ambiguity_id(o) == ambiguity_sat, do: 1.0, else: 0.0
      end)

    [dx, dy, dz, dc | ambiguity_cols]
  end

  defp build_multi_rows(sp3, epochs, sat_ids, state, weights, tropo) do
    {rx, ry, rz} = state.position

    epochs
    |> Enum.with_index()
    |> Enum.reduce_while({:ok, []}, fn {epoch_row, epoch_idx}, {:ok, acc} ->
      clock_m = Enum.at(state.clocks_m, epoch_idx)

      epoch_row.observations
      |> Enum.reduce_while({:ok, acc}, fn o, {:ok, rows_acc} ->
        case Observables.predict(sp3, o.satellite_id, {rx, ry, rz}, epoch_row.epoch) do
          {:ok, pred} ->
            case model_troposphere(pred, {rx, ry, rz}, epoch_row.epoch, tropo) do
              {:ok, tropo_model} ->
                {ex, ey, ez} = pred.los_unit
                sat_clock_m = Constants.speed_of_light_m_s() * (pred.sat_clock_s || 0.0)
                tropo_m = applied_troposphere_m(tropo_model, state)
                model_code = pred.geometric_range_m + clock_m - sat_clock_m + tropo_m
                arc_id = ambiguity_id(o)
                ambiguity = Map.fetch!(state.ambiguities, arc_id)
                base = {-ex, -ey, -ez}

                code_row = %{
                  kind: :code,
                  sat: o.satellite_id,
                  epoch: epoch_row.epoch,
                  h:
                    multi_design_row(
                      base,
                      epoch_idx,
                      nil,
                      length(epochs),
                      sat_ids,
                      tropo_model.ztd_mapping,
                      tropo
                    ),
                  y: o.code_m - model_code,
                  weight: measurement_weight(weights, :code, pred.elevation_deg)
                }

                phase_row = %{
                  kind: :phase,
                  sat: o.satellite_id,
                  epoch: epoch_row.epoch,
                  h:
                    multi_design_row(
                      base,
                      epoch_idx,
                      arc_id,
                      length(epochs),
                      sat_ids,
                      tropo_model.ztd_mapping,
                      tropo
                    ),
                  y: o.phase_m - (model_code + ambiguity),
                  weight: measurement_weight(weights, :phase, pred.elevation_deg)
                }

                {:cont, {:ok, [phase_row, code_row | rows_acc]}}

              {:error, reason} ->
                {:halt, {:error, {:troposphere_failed, o.satellite_id, reason}}}
            end

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

  defp build_fixed_multi_rows(sp3, epochs, fixed_m, state, weights, tropo) do
    {rx, ry, rz} = state.position

    epochs
    |> Enum.with_index()
    |> Enum.reduce_while({:ok, []}, fn {epoch_row, epoch_idx}, {:ok, acc} ->
      clock_m = Enum.at(state.clocks_m, epoch_idx)

      epoch_row.observations
      |> Enum.reduce_while({:ok, acc}, fn o, {:ok, rows_acc} ->
        case Observables.predict(sp3, o.satellite_id, {rx, ry, rz}, epoch_row.epoch) do
          {:ok, pred} ->
            case model_troposphere(pred, {rx, ry, rz}, epoch_row.epoch, tropo) do
              {:ok, tropo_model} ->
                {ex, ey, ez} = pred.los_unit
                sat_clock_m = Constants.speed_of_light_m_s() * (pred.sat_clock_s || 0.0)
                tropo_m = applied_troposphere_m(tropo_model, state)
                model_code = pred.geometric_range_m + clock_m - sat_clock_m + tropo_m
                ambiguity = Map.fetch!(fixed_m, ambiguity_id(o))

                h =
                  fixed_multi_design_row(
                    {-ex, -ey, -ez},
                    epoch_idx,
                    length(epochs),
                    tropo_model.ztd_mapping,
                    tropo
                  )

                code_row = %{
                  kind: :code,
                  sat: o.satellite_id,
                  epoch: epoch_row.epoch,
                  h: h,
                  y: o.code_m - model_code,
                  weight: measurement_weight(weights, :code, pred.elevation_deg)
                }

                phase_row = %{
                  kind: :phase,
                  sat: o.satellite_id,
                  epoch: epoch_row.epoch,
                  h: h,
                  y: o.phase_m - (model_code + ambiguity),
                  weight: measurement_weight(weights, :phase, pred.elevation_deg)
                }

                {:cont, {:ok, [phase_row, code_row | rows_acc]}}

              {:error, reason} ->
                {:halt, {:error, {:troposphere_failed, o.satellite_id, reason}}}
            end

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

  defp multi_design_row(
         {dx, dy, dz},
         epoch_idx,
         ambiguity_sat,
         n_epochs,
         sat_ids,
         ztd_mapping,
         tropo
       ) do
    clock_cols = for idx <- 0..(n_epochs - 1), do: if(idx == epoch_idx, do: 1.0, else: 0.0)
    ztd_cols = ztd_design_cols(ztd_mapping, tropo)
    ambiguity_cols = Enum.map(sat_ids, &if(&1 == ambiguity_sat, do: 1.0, else: 0.0))
    [dx, dy, dz | clock_cols ++ ztd_cols ++ ambiguity_cols]
  end

  defp fixed_multi_design_row({dx, dy, dz}, epoch_idx, n_epochs, ztd_mapping, tropo) do
    clock_cols = for idx <- 0..(n_epochs - 1), do: if(idx == epoch_idx, do: 1.0, else: 0.0)
    [dx, dy, dz | clock_cols ++ ztd_design_cols(ztd_mapping, tropo)]
  end

  defp ztd_design_cols(ztd_mapping, %{estimate_ztd?: true}), do: [ztd_mapping]
  defp ztd_design_cols(_ztd_mapping, _tropo), do: []

  defp model_troposphere(_pred, _receiver_m, _epoch, %{enabled?: false}),
    do: {:ok, %{slant_m: 0.0, ztd_mapping: 0.0}}

  defp model_troposphere(pred, {rx, ry, rz}, epoch, %{enabled?: true, met: met} = tropo) do
    geo = Coordinates.to_geodetic({rx / 1000.0, ry / 1000.0, rz / 1000.0})

    height_m = geo.altitude_km * 1000.0

    with {:ok, slant_m} <-
           Troposphere.slant_delay(
             pred.elevation_deg,
             geo.latitude,
             geo.longitude,
             height_m,
             met,
             epoch
           ),
         {:ok, ztd_mapping} <-
           ztd_mapping(pred.elevation_deg, geo.latitude, height_m, epoch, tropo) do
      {:ok, %{slant_m: slant_m, ztd_mapping: ztd_mapping}}
    end
  end

  defp ztd_mapping(elevation_deg, latitude_deg, height_m, epoch, %{estimate_ztd?: true}) do
    with {:ok, %{wet: wet_mapping}} <-
           Troposphere.mapping(elevation_deg, latitude_deg, height_m, epoch) do
      {:ok, wet_mapping}
    end
  end

  defp ztd_mapping(_elevation_deg, _latitude_deg, _height_m, _epoch, _tropo), do: {:ok, 0.0}

  defp applied_troposphere_m(tropo_model, state) do
    tropo_model.slant_m + Map.get(state, :ztd_m, 0.0) * tropo_model.ztd_mapping
  end

  defp apply_delta(state, obs, dx) do
    {rx, ry, rz} = state.position
    [dx0, dx1, dx2, dclock | ambiguity_deltas] = dx

    ambiguities =
      obs
      |> Enum.zip(ambiguity_deltas)
      |> Map.new(fn {o, d} ->
        arc_id = ambiguity_id(o)
        {arc_id, Map.fetch!(state.ambiguities, arc_id) + d}
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

  defp apply_multi_delta(state, epochs, sat_ids, dx, tropo) do
    {rx, ry, rz} = state.position
    [dx0, dx1, dx2 | rest] = dx
    {clock_deltas, rest} = Enum.split(rest, length(epochs))
    {ztd_deltas, ambiguity_deltas} = Enum.split(rest, ztd_unknown_count(tropo))

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
      ztd_m: state_ztd_m(state) + ztd_delta(ztd_deltas),
      ambiguities: ambiguities
    }
  end

  defp apply_fixed_multi_delta(state, dx, n_epochs, tropo) do
    {rx, ry, rz} = state.position
    [dx0, dx1, dx2 | rest] = dx
    {clock_deltas, ztd_deltas} = Enum.split(rest, n_epochs)

    clocks =
      state.clocks_m
      |> Enum.zip(Enum.take(clock_deltas, n_epochs))
      |> Enum.map(fn {clock, delta} -> clock + delta end)

    %{
      position: {rx + dx0, ry + dx1, rz + dx2},
      clocks_m: clocks,
      ztd_m: state_ztd_m(state) + ztd_delta(Enum.take(ztd_deltas, ztd_unknown_count(tropo)))
    }
  end

  defp multi_step_norms([dx, dy, dz | rest], n_epochs, tropo) do
    {clock_deltas, rest} = Enum.split(rest, n_epochs)
    {ztd_deltas, ambiguity_deltas} = Enum.split(rest, ztd_unknown_count(tropo))

    {
      :math.sqrt(dx * dx + dy * dy + dz * dz),
      max_abs(clock_deltas),
      max_abs(ztd_deltas),
      max_abs(ambiguity_deltas)
    }
  end

  defp fixed_multi_step_norms([dx, dy, dz | rest], tropo) do
    # Fixed arcs have `[x, y, z, clocks..., ztd?]`.
    n_ztd = ztd_unknown_count(tropo)
    n_clocks = length(rest) - n_ztd
    {clock_deltas, ztd_deltas} = Enum.split(rest, n_clocks)

    {:math.sqrt(dx * dx + dy * dy + dz * dz), max_abs(clock_deltas), max_abs(ztd_deltas)}
  end

  defp max_abs([]), do: 0.0
  defp max_abs(xs), do: xs |> Enum.map(&abs/1) |> Enum.max()
  defp state_ztd_m(state), do: Map.get(state, :ztd_m, 0.0)
  defp ztd_delta([delta]), do: delta
  defp ztd_delta(_deltas), do: 0.0

  # --- integer ambiguity fixing -------------------------------------------

  defp search_integer_ambiguities(
         sp3,
         epochs,
         float_sol,
         wavelengths,
         offsets,
         weights,
         tropo,
         opts
       ) do
    sat_ids = float_sol.used_sats

    case ambiguity_covariance_cycles(sp3, epochs, sat_ids, float_sol, wavelengths, weights, tropo) do
      {:ok, q_cycles} ->
        float_cycles_by_sat = float_ambiguities_cycles(float_sol, wavelengths, offsets)
        IntegerLeastSquares.search(float_cycles_by_sat, q_cycles, opts)

      {:error, _} = err ->
        err
    end
  end

  defp fixed_ambiguities_m(fixed_cycles, wavelengths, offsets) do
    Map.new(fixed_cycles, fn {sat, cycles} ->
      {sat, Map.fetch!(offsets, sat) + cycles * Map.fetch!(wavelengths, sat)}
    end)
  end

  defp fixed_state_from_float(float_sol) do
    %{
      position: {float_sol.position.x_m, float_sol.position.y_m, float_sol.position.z_m},
      clocks_m: Enum.map(float_sol.epoch_clocks, & &1.rx_clock_m),
      ztd_m: float_sol.ztd_residual_m || 0.0
    }
  end

  defp float_ambiguities_cycles(float_sol, wavelengths, offsets) do
    Map.new(float_sol.used_sats, fn sat ->
      {sat,
       (Map.fetch!(float_sol.ambiguities_m, sat) - Map.fetch!(offsets, sat)) /
         Map.fetch!(wavelengths, sat)}
    end)
  end

  defp ambiguity_covariance_cycles(sp3, epochs, sat_ids, float_sol, wavelengths, weights, tropo) do
    state = %{
      position: {float_sol.position.x_m, float_sol.position.y_m, float_sol.position.z_m},
      clocks_m: Enum.map(float_sol.epoch_clocks, & &1.rx_clock_m),
      ztd_m: float_sol.ztd_residual_m || 0.0,
      ambiguities: float_sol.ambiguities_m
    }

    n = 3 + length(epochs) + ztd_unknown_count(tropo) + length(sat_ids)
    start = 3 + length(epochs) + ztd_unknown_count(tropo)

    with {:ok, rows} <- build_multi_rows(sp3, epochs, sat_ids, state, weights, tropo),
         {normal, _rhs} = normal_equations(rows, n),
         {:ok, covariance_m} <- ambiguity_covariance_from_normal(normal, start, length(sat_ids)) do
      covariance_cycles =
        for i <- 0..(length(sat_ids) - 1) do
          sat_i = Enum.at(sat_ids, i)
          lambda_i = Map.fetch!(wavelengths, sat_i)

          for j <- 0..(length(sat_ids) - 1) do
            sat_j = Enum.at(sat_ids, j)
            lambda_j = Map.fetch!(wavelengths, sat_j)
            (covariance_m |> Enum.at(i) |> Enum.at(j)) / (lambda_i * lambda_j)
          end
        end

      {:ok, covariance_cycles}
    end
  end

  defp ambiguity_covariance_from_normal(normal, start, n_ambiguities) do
    a = submatrix(normal, 0, start, 0, start)
    b = submatrix(normal, 0, start, start, n_ambiguities)
    c = submatrix(normal, start, n_ambiguities, start, n_ambiguities)

    with {:ok, a_inv_b} <- solve_matrix(a, b) do
      bt_a_inv_b = matmul(transpose(b), a_inv_b)
      schur = matrix_sub(c, bt_a_inv_b)

      case invert_matrix(schur) do
        {:ok, covariance} ->
          {:ok, covariance}

        {:error, :singular_geometry} = err ->
          err
      end
    end
  end

  # --- result assembly -----------------------------------------------------

  defp finalize(sp3, obs, epoch, state, weights, tropo, iterations, converged, status) do
    with {:ok, residual_rows} <- residual_rows(sp3, obs, epoch, state, tropo, weights) do
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
           weighted_rms_m: weighted_rms(residual_rows, weights),
           troposphere_applied: tropo.enabled?
         }
       }}
    end
  end

  defp finalize_multi(sp3, epochs, sat_ids, state, weights, tropo, iterations, converged, status) do
    with {:ok, residual_rows} <- multi_residual_rows(sp3, epochs, state, tropo, weights) do
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
         ztd_residual_m: solution_ztd_residual_m(state, tropo),
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
           weighted_rms_m: weighted_rms(residual_rows, weights),
           troposphere_applied: tropo.enabled?,
           ztd_estimated: tropo.estimate_ztd?
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
         tropo,
         iterations,
         converged,
         status,
         fixed_meta
       ) do
    with {:ok, residual_rows} <-
           fixed_multi_residual_rows(sp3, epochs, fixed_m, state, tropo, weights) do
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
         wide_lane_ambiguities_cycles: nil,
         ztd_residual_m: solution_ztd_residual_m(state, tropo),
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
             weighted_rms_m: weighted_rms(residual_rows, weights),
             troposphere_applied: tropo.enabled?,
             ztd_estimated: tropo.estimate_ztd?
           })
       }}
    end
  end

  defp solution_ztd_residual_m(state, %{estimate_ztd?: true}), do: state_ztd_m(state)
  defp solution_ztd_residual_m(_state, _tropo), do: nil

  defp residual_rows(sp3, obs, epoch, state, tropo, weights) do
    {rx, ry, rz} = state.position

    Enum.reduce_while(obs, {:ok, []}, fn o, {:ok, acc} ->
      case Observables.predict(sp3, o.satellite_id, {rx, ry, rz}, epoch) do
        {:ok, pred} ->
          case model_troposphere(pred, {rx, ry, rz}, epoch, tropo) do
            {:ok, tropo_model} ->
              sat_clock_m = Constants.speed_of_light_m_s() * (pred.sat_clock_s || 0.0)
              tropo_m = applied_troposphere_m(tropo_model, state)
              model_code = pred.geometric_range_m + state.clock_m - sat_clock_m + tropo_m
              ambiguity = Map.fetch!(state.ambiguities, ambiguity_id(o))

              row = %{
                sat: o.satellite_id,
                code_m: o.code_m - model_code,
                phase_m: o.phase_m - (model_code + ambiguity),
                code_weight: measurement_weight(weights, :code, pred.elevation_deg),
                phase_weight: measurement_weight(weights, :phase, pred.elevation_deg)
              }

              {:cont, {:ok, [row | acc]}}

            {:error, reason} ->
              {:halt, {:error, {:troposphere_failed, o.satellite_id, reason}}}
          end

        {:error, reason} ->
          {:halt, {:error, {:no_ephemeris, o.satellite_id, reason}}}
      end
    end)
    |> case do
      {:ok, rows} -> {:ok, Enum.reverse(rows)}
      {:error, _} = err -> err
    end
  end

  defp multi_residual_rows(sp3, epochs, state, tropo, weights) do
    {rx, ry, rz} = state.position

    epochs
    |> Enum.with_index()
    |> Enum.reduce_while({:ok, []}, fn {epoch_row, epoch_idx}, {:ok, acc} ->
      clock_m = Enum.at(state.clocks_m, epoch_idx)

      epoch_row.observations
      |> Enum.reduce_while({:ok, acc}, fn o, {:ok, rows_acc} ->
        case Observables.predict(sp3, o.satellite_id, {rx, ry, rz}, epoch_row.epoch) do
          {:ok, pred} ->
            case model_troposphere(pred, {rx, ry, rz}, epoch_row.epoch, tropo) do
              {:ok, tropo_model} ->
                sat_clock_m = Constants.speed_of_light_m_s() * (pred.sat_clock_s || 0.0)
                tropo_m = applied_troposphere_m(tropo_model, state)
                model_code = pred.geometric_range_m + clock_m - sat_clock_m + tropo_m
                ambiguity = Map.fetch!(state.ambiguities, ambiguity_id(o))

                row = %{
                  epoch: epoch_row.epoch,
                  satellite_id: o.satellite_id,
                  code_m: o.code_m - model_code,
                  phase_m: o.phase_m - (model_code + ambiguity),
                  code_weight: measurement_weight(weights, :code, pred.elevation_deg),
                  phase_weight: measurement_weight(weights, :phase, pred.elevation_deg)
                }

                {:cont, {:ok, [row | rows_acc]}}

              {:error, reason} ->
                {:halt, {:error, {:troposphere_failed, o.satellite_id, reason}}}
            end

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

  defp fixed_multi_residual_rows(sp3, epochs, fixed_m, state, tropo, weights) do
    {rx, ry, rz} = state.position

    epochs
    |> Enum.with_index()
    |> Enum.reduce_while({:ok, []}, fn {epoch_row, epoch_idx}, {:ok, acc} ->
      clock_m = Enum.at(state.clocks_m, epoch_idx)

      epoch_row.observations
      |> Enum.reduce_while({:ok, acc}, fn o, {:ok, rows_acc} ->
        case Observables.predict(sp3, o.satellite_id, {rx, ry, rz}, epoch_row.epoch) do
          {:ok, pred} ->
            case model_troposphere(pred, {rx, ry, rz}, epoch_row.epoch, tropo) do
              {:ok, tropo_model} ->
                sat_clock_m = Constants.speed_of_light_m_s() * (pred.sat_clock_s || 0.0)
                tropo_m = applied_troposphere_m(tropo_model, state)
                model_code = pred.geometric_range_m + clock_m - sat_clock_m + tropo_m
                ambiguity = Map.fetch!(fixed_m, ambiguity_id(o))

                row = %{
                  epoch: epoch_row.epoch,
                  satellite_id: o.satellite_id,
                  code_m: o.code_m - model_code,
                  phase_m: o.phase_m - (model_code + ambiguity),
                  code_weight: measurement_weight(weights, :code, pred.elevation_deg),
                  phase_weight: measurement_weight(weights, :phase, pred.elevation_deg)
                }

                {:cont, {:ok, [row | rows_acc]}}

              {:error, reason} ->
                {:halt, {:error, {:troposphere_failed, o.satellite_id, reason}}}
            end

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
    |> Enum.map(&ambiguity_id/1)
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
        [
          r.code_m * Map.get(r, :code_weight, weights.code),
          r.phase_m * Map.get(r, :phase_weight, weights.phase)
        ]
      end)

    rms(values)
  end
end
