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

  `double_differences/3` returns the normalized measurements. The float solver,
  `solve_float_baseline_epochs/3`, estimates one static base-to-rover baseline
  from code and carrier-phase double differences, keeping one float carrier
  ambiguity per non-reference double-difference arc across the data. Clean arcs
  use the physical satellite id (for example `"G05"`) as the ambiguity id; split
  arcs use explicit ids so a cycle slip resets the ambiguity without pretending
  the satellite disappeared.
  `solve_fixed_baseline_epochs/3` adds bounded integer least-squares ambiguity
  fixing on top of the same correlated double-difference covariance and
  re-solves the baseline with the selected integers held fixed.

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
      [%{satellite_id: "G02", reference_satellite_id: "G01", ambiguity_id: "G02", code_m: 15.0, phase_m: 17.0}]
  """

  import Bitwise, only: [band: 2]

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

  alias Orbis.GNSS.{CarrierPhase, IonosphereFree}
  alias Orbis.GNSS.Core.{Constants, Types}
  alias Orbis.GNSS.Core.IntegerLeastSquares

  @default_max_iterations 8
  @default_position_tolerance_m 1.0e-4
  @default_ambiguity_tolerance_m 1.0e-4
  @default_code_sigma_m 1.0
  @default_phase_sigma_m 0.02
  @default_integer_search_radius_cycles 1
  @default_integer_ratio_threshold 3.0
  @default_integer_candidate_limit 200_000
  @default_partial_min_ambiguities 4
  @default_cycle_slip_policy :error
  @default_hatch_window_cap 100
  @min_elevation_sin 0.05

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
            ambiguity_id: String.t(),
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
              physical_sats: [String.t()],
              ambiguity_satellites: %{String.t() => String.t()},
              ambiguity_float: %{
                order: [String.t()],
                covariance_m: [[float()]],
                covariance_inverse_m: [[float()]]
              },
              measurement_covariance: %{
                model: :double_difference,
                code_sigma_m: float(),
                phase_sigma_m: float(),
                elevation_weighting: boolean(),
                min_elevation_sin: float()
              },
              code_rms_m: float(),
              phase_rms_m: float(),
              weighted_rms_m: float(),
              n_epochs: pos_integer(),
              n_observations: pos_integer(),
              dropped_sats: [String.t()],
              dropped_cycle_slip_sats: [String.t()],
              split_cycle_slip_arcs: [map()]
            }
          }
  end

  defmodule FixedBaselineSolution do
    @moduledoc """
    Integer-fixed RTK baseline solution from code/carrier double differences.
    """

    @enforce_keys [
      :baseline_m,
      :rover_position_m,
      :reference_satellite_id,
      :used_sats,
      :fixed_ambiguities_cycles,
      :fixed_ambiguities_m,
      :wide_lane_ambiguities_cycles,
      :float_solution,
      :residuals_m,
      :metadata
    ]
    defstruct [
      :baseline_m,
      :rover_position_m,
      :reference_satellite_id,
      :used_sats,
      :fixed_ambiguities_cycles,
      :fixed_ambiguities_m,
      :wide_lane_ambiguities_cycles,
      :float_solution,
      :residuals_m,
      :metadata
    ]

    @type ecef :: %{x_m: float(), y_m: float(), z_m: float()}

    @type t :: %__MODULE__{
            baseline_m: ecef(),
            rover_position_m: ecef(),
            reference_satellite_id: String.t(),
            used_sats: [String.t()],
            fixed_ambiguities_cycles: %{String.t() => integer()},
            fixed_ambiguities_m: %{String.t() => float()},
            wide_lane_ambiguities_cycles: %{String.t() => integer()} | nil,
            float_solution: FloatBaselineSolution.t(),
            residuals_m: [FloatBaselineSolution.residual()],
            metadata: %{
              required(:iterations) => pos_integer(),
              required(:converged) => boolean(),
              required(:status) => :state_tolerance | :max_iterations,
              required(:integer_status) => :fixed | :not_fixed,
              required(:integer_method) => :lambda | :widelane_narrowlane_lambda,
              required(:integer_ratio) => float() | :infinity,
              required(:integer_best_score) => float(),
              required(:integer_second_best_score) => float() | nil,
              required(:integer_candidates) => pos_integer(),
              required(:code_rms_m) => float(),
              required(:phase_rms_m) => float(),
              required(:weighted_rms_m) => float(),
              required(:n_epochs) => pos_integer(),
              required(:n_observations) => pos_integer(),
              required(:measurement_covariance) => %{
                model: :double_difference,
                code_sigma_m: float(),
                phase_sigma_m: float(),
                elevation_weighting: boolean(),
                min_elevation_sin: float()
              },
              required(:ambiguity_search) => %{
                order: [String.t()],
                float_cycles: %{String.t() => float()},
                covariance_cycles: [[float()]],
                covariance_inverse_cycles: [[float()]]
              },
              required(:ambiguity_offsets_m) => %{String.t() => float()},
              optional(:wide_lane_fixed) => boolean(),
              optional(:wide_lane_ambiguities_cycles) => %{String.t() => integer()},
              optional(:physical_sats) => [String.t()],
              optional(:ambiguity_satellites) => %{String.t() => String.t()},
              optional(:partial_ambiguity_resolution) => boolean(),
              optional(:partial_fixed) => boolean(),
              optional(:partial_fixed_ambiguities) => [String.t()],
              optional(:partial_free_ambiguities) => [String.t()],
              optional(:partial_full_set) => map(),
              optional(:dropped_cycle_slip_sats) => [String.t()],
              optional(:split_cycle_slip_arcs) => [map()]
            }
          }
  end

  @typedoc """
  Code and carrier-phase observation in metres.

  Map observations may optionally carry `:ambiguity_id` to identify a carrier
  arc and `:lli` (or `:loss_of_lock_indicator`) for single-frequency
  loss-of-lock handling. Tuple observations use the satellite id as the
  ambiguity id and have no LLI.
  """
  @type observation ::
          %{
            required(:satellite_id) => String.t(),
            required(:code_m) => number(),
            required(:phase_m) => number(),
            optional(:ambiguity_id) => String.t(),
            optional(:lli) => integer() | nil,
            optional(:loss_of_lock_indicator) => integer() | nil
          }
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

  @typedoc """
  Raw dual-frequency code/carrier observation for wide-lane/narrow-lane RTK.

  `p1_m` / `p2_m` are code pseudoranges in metres, `phi1_cyc` / `phi2_cyc` are
  carrier phases in cycles, and `f1_hz` / `f2_hz` are the corresponding carrier
  frequencies. `:ambiguity_id` is normally omitted; the wide-lane solver sets it
  internally when `:on_cycle_slip` is `:split_arc`.
  """
  @type dual_frequency_observation :: %{
          required(:satellite_id) => String.t(),
          required(:p1_m) => number(),
          required(:p2_m) => number(),
          required(:phi1_cyc) => number(),
          required(:phi2_cyc) => number(),
          required(:f1_hz) => number(),
          required(:f2_hz) => number(),
          optional(:ambiguity_id) => String.t(),
          optional(:lli1) => integer() | nil,
          optional(:lli2) => integer() | nil
        }

  @typedoc "One RTK epoch carrying raw dual-frequency base/rover observations."
  @type dual_frequency_baseline_epoch :: %{
          required(:base_observations) => [dual_frequency_observation()],
          required(:rover_observations) => [dual_frequency_observation()],
          required(:satellite_positions_m) => satellite_positions(),
          optional(:epoch) => term()
        }

  @typedoc "One non-reference satellite's double-difference observation."
  @type double_difference :: %{
          satellite_id: String.t(),
          reference_satellite_id: String.t(),
          ambiguity_id: String.t(),
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

  The reference satellite must be available in every epoch. Other satellites may
  appear in only part of the arc; each available non-reference observation
  contributes a row, and split cycle-slip fragments become independent ambiguity
  ids.

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
    * `:on_cycle_slip` - what to do when a base or rover observation carries an
      LLI loss-of-lock bit: `:error` returns
      `{:error, {:cycle_slip_detected, receiver, sat, epoch, [:lli]}}`
      (default); `:drop_satellite` removes that satellite from the arc;
      `:split_arc` starts a new ambiguity arc at the slipped epoch.
    * `:elevation_weighting` - when `true`, scales each undifferenced
      measurement sigma by `1 / max(sin(elevation), #{@min_elevation_sin})`
      before propagating the double-difference covariance. Default `false`
      preserves the constant-sigma, transcendental-free solve path.
    * `:code_smoothing` - when `true`, applies per-receiver/per-ambiguity-arc
      Hatch carrier smoothing to code observations before forming double
      differences. Default `false`.
    * `:hatch_window_cap` - maximum Hatch smoothing window when
      `:code_smoothing` is enabled (default `#{@default_hatch_window_cap}`).
    * `:max_iterations`, `:position_tolerance_m`,
      `:ambiguity_tolerance_m`.

  Returns `{:ok, %FloatBaselineSolution{}}` or a tagged error.
  """
  @spec solve_float_baseline_epochs(ecef_input(), [baseline_epoch()], keyword()) ::
          {:ok, FloatBaselineSolution.t()} | {:error, term()}
  def solve_float_baseline_epochs(base_position, epochs, opts \\ [])

  def solve_float_baseline_epochs(base_position, epochs, opts) when is_list(epochs) do
    with {:ok, base} <- Types.normalize_ecef(base_position, :invalid_base_position),
         :ok <- ensure_nonempty_epochs(epochs),
         {:ok, normalized_epochs} <- normalize_epochs(epochs),
         {:ok, normalized_epochs, slip_meta} <-
           prepare_epochs_for_cycle_slips(normalized_epochs, opts),
         {:ok, normalized_epochs, smoothing_meta} <-
           prepare_epochs_for_code_smoothing(normalized_epochs, opts),
         {:ok, common_sats, _not_common_sats} <- common_epoch_sats(normalized_epochs),
         {:ok, all_sats} <- all_epoch_sats(normalized_epochs),
         :ok <- ensure_baseline_satellites(all_sats),
         {:ok, reference_sat} <-
           baseline_reference_satellite(common_sats, opts, base, normalized_epochs),
         {:ok, solve_opts} <- baseline_solve_options(opts),
         {:ok, weights} <- baseline_weights(opts),
         {:ok, initial_baseline} <- initial_baseline(opts),
         {:ok, ambiguity_ids, ambiguity_satellites} <-
           baseline_ambiguity_index(normalized_epochs, all_sats, reference_sat) do
      physical_sats = Enum.reject(all_sats, &(&1 == reference_sat))

      if baseline_row_count(normalized_epochs, reference_sat) <
           baseline_unknown_count(ambiguity_ids) do
        {:error,
         {:underdetermined, baseline_row_count(normalized_epochs, reference_sat),
          baseline_unknown_count(ambiguity_ids)}}
      else
        state = %{
          baseline: initial_baseline,
          ambiguities: Map.new(ambiguity_ids, &{&1, 0.0})
        }

        iterate_baseline(
          base,
          normalized_epochs,
          reference_sat,
          physical_sats,
          ambiguity_ids,
          ambiguity_satellites,
          state,
          weights,
          solve_opts,
          [],
          Map.merge(slip_meta, smoothing_meta),
          1
        )
      end
    end
  end

  def solve_float_baseline_epochs(_base_position, _epochs, _opts), do: {:error, :invalid_epochs}

  @doc """
  Solve a static RTK baseline with integer-fixed double-difference ambiguities.

  The function first runs `solve_float_baseline_epochs/3`, converts the float
  double-difference ambiguities from metres to cycles using
  `:ambiguity_wavelength_m`, runs the shared bounded integer least-squares
  search with the correlated float ambiguity covariance, and then re-solves the
  baseline with the selected integer ambiguities held fixed.

  Required option:

    * `:ambiguity_wavelength_m` - either a positive scalar wavelength in metres
      for every non-reference satellite, or a map `%{"G05" => wavelength_m, ...}`.
    * `:ambiguity_offset_m` - optional fixed ambiguity offset in metres, either
      a scalar or a map keyed by ambiguity id / physical satellite id. The fixed
      carrier ambiguity model is `offset_m + integer * wavelength_m`. Defaults
      to zero and is useful for dual-frequency wide-lane/narrow-lane workflows
      where the wide-lane integer contributes a known ionosphere-free offset.

  Integer search options mirror `Orbis.GNSS.PrecisePositioning`:

    * `:integer_search_radius_cycles` - default
      `#{@default_integer_search_radius_cycles}`.
    * `:integer_ratio_threshold` - default `#{@default_integer_ratio_threshold}`.
    * `:integer_candidate_limit` - default `#{@default_integer_candidate_limit}`.
    * `:partial_ambiguity_resolution` - when `true`, a rejected full-set
      integer fix is followed by confidence-ranked subset searches. A passing
      subset is held fixed while the remaining ambiguities stay in the re-solve
      as float states (default `false`).
    * `:partial_min_ambiguities` - minimum subset size for partial ambiguity
      resolution (default `#{@default_partial_min_ambiguities}`).

  The fixed solution is returned even when the ratio test fails; in that case
  `metadata.integer_status` is `:not_fixed`.
  """
  @spec solve_fixed_baseline_epochs(ecef_input(), [baseline_epoch()], keyword()) ::
          {:ok, FixedBaselineSolution.t()} | {:error, term()}
  def solve_fixed_baseline_epochs(base_position, epochs, opts \\ [])

  def solve_fixed_baseline_epochs(base_position, epochs, opts) when is_list(epochs) do
    with {:ok, float_sol} <- solve_float_baseline_epochs(base_position, epochs, opts),
         {:ok, wavelengths} <-
           ambiguity_wavelengths(
             float_sol.used_sats,
             Map.fetch!(float_sol.metadata, :ambiguity_satellites),
             opts
           ),
         {:ok, offsets} <-
           ambiguity_offsets(
             float_sol.used_sats,
             Map.fetch!(float_sol.metadata, :ambiguity_satellites),
             opts
           ),
         {:ok, integer_opts} <- integer_options(opts),
         {:ok, fixed_cycles, fixed_meta} <-
           search_baseline_ambiguities(float_sol, wavelengths, offsets, integer_opts),
         fixed_m = fixed_ambiguities_m(fixed_cycles, wavelengths, offsets),
         free_ambiguity_ids = free_ambiguity_ids(float_sol.used_sats, fixed_cycles),
         {:ok, base} <- Types.normalize_ecef(base_position, :invalid_base_position),
         {:ok, normalized_epochs} <- normalize_epochs(epochs),
         {:ok, normalized_epochs, _slip_meta} <-
           prepare_epochs_for_cycle_slips(normalized_epochs, opts),
         {:ok, weights} <- baseline_weights(opts),
         {:ok, solve_opts} <- baseline_solve_options(opts) do
      state = %{
        baseline: {float_sol.baseline_m.x_m, float_sol.baseline_m.y_m, float_sol.baseline_m.z_m},
        ambiguities: Map.take(float_sol.ambiguities_m, free_ambiguity_ids)
      }

      iterate_fixed_baseline(
        base,
        normalized_epochs,
        float_sol.reference_satellite_id,
        Map.fetch!(float_sol.metadata, :physical_sats),
        float_sol.used_sats,
        free_ambiguity_ids,
        fixed_m,
        state,
        weights,
        solve_opts,
        float_sol,
        fixed_cycles,
        fixed_meta,
        1
      )
    end
  end

  def solve_fixed_baseline_epochs(_base_position, _epochs, _opts), do: {:error, :invalid_epochs}

  @doc """
  Solve a static RTK baseline from raw dual-frequency observations by fixing
  wide-lane then narrow-lane double-difference ambiguities.

  This is the dual-frequency convenience layer above
  `solve_fixed_baseline_epochs/3`. Each base and rover observation must carry
  two code and phase measurements:

      %{
        satellite_id: "G05",
        p1_m: 20_200_000.0,
        p2_m: 20_200_004.0,
        phi1_cyc: 106_000_000.0,
        phi2_cyc: 82_000_000.0,
        f1_hz: 1_575_420_000.0,
        f2_hz: 1_227_600_000.0,
        lli1: 0,
        lli2: 0
      }

  For every non-reference double-difference arc the function estimates the
  Melbourne-Wubbena wide-lane integer first. It then forms ionosphere-free code
  and phase double differences and fixes the remaining narrow-lane integer with
  bounded integer least-squares. The returned `fixed_ambiguities_cycles` are the
  narrow-lane integers; `wide_lane_ambiguities_cycles` reports the fixed
  wide-lane integers.

  Options are the same as `solve_fixed_baseline_epochs/3`, except
  `:ambiguity_wavelength_m` and `:ambiguity_offset_m` are derived internally.
  Additional wide-lane options:

    * `:wide_lane_min_epochs` - minimum Melbourne-Wubbena epochs per
      double-difference arc (default `2`).
    * `:wide_lane_tolerance_cycles` - maximum absolute distance between the
      averaged wide-lane float value and the nearest integer (default `0.5`
      cycles).
    * `:on_cycle_slip` - `:error` (default), `:drop_satellite`, or `:split_arc`.
      Split arcs get fresh ambiguity ids and are fixed independently.
    * `:partial_ambiguity_resolution` - when `true`, a rejected full narrow-lane
      integer fix is followed by confidence-ranked subset searches (and, when
      the greedy ranking finds nothing, a bounded largest-first exhaustive
      combination search). Holding the wide-lane integers fixed collapses the
      per-ambiguity bias for most satellites, so the dual-frequency partial fix
      can safely cover a larger subset than the single-frequency partial. The
      ratio threshold is never weakened (default `false`).
    * `:partial_min_ambiguities` - minimum subset size for partial ambiguity
      resolution (default `#{@default_partial_min_ambiguities}`).

  Returns `{:ok, %FixedBaselineSolution{}}` or a tagged error.
  """
  @spec solve_widelane_fixed_baseline_epochs(
          ecef_input(),
          [dual_frequency_baseline_epoch()],
          keyword()
        ) :: {:ok, FixedBaselineSolution.t()} | {:error, term()}
  def solve_widelane_fixed_baseline_epochs(base_position, dual_epochs, opts \\ [])

  def solve_widelane_fixed_baseline_epochs(base_position, dual_epochs, opts)
      when is_list(dual_epochs) do
    with {:ok, base} <- Types.normalize_ecef(base_position, :invalid_base_position),
         :ok <- ensure_nonempty_epochs(dual_epochs),
         {:ok, normalized_dual_epochs} <- normalize_dual_baseline_epochs(dual_epochs),
         {:ok, prepared_dual_epochs, slip_meta} <-
           prepare_dual_baseline_cycle_slips(normalized_dual_epochs, opts),
         {:ok, common_sats, _dropped_sats} <- common_epoch_sats(prepared_dual_epochs),
         :ok <- ensure_baseline_satellites(common_sats),
         {:ok, reference_sat} <-
           baseline_reference_satellite(common_sats, opts, base, prepared_dual_epochs),
         {:ok, wide_lane_cycles} <-
           estimate_dual_baseline_wide_lanes(prepared_dual_epochs, reference_sat, opts),
         {:ok, if_epochs, wavelengths, offsets} <-
           ionosphere_free_baseline_epochs(
             prepared_dual_epochs,
             reference_sat,
             wide_lane_cycles
           ),
         fixed_opts =
           opts
           |> Keyword.put(:reference_satellite_id, reference_sat)
           |> Keyword.put(:ambiguity_wavelength_m, wavelengths)
           |> Keyword.put(:ambiguity_offset_m, offsets),
         {:ok, %FixedBaselineSolution{} = sol} <-
           solve_fixed_baseline_epochs(base_position, if_epochs, fixed_opts) do
      used_wide_lane_cycles = Map.take(wide_lane_cycles, sol.used_sats)

      {:ok,
       %{
         sol
         | wide_lane_ambiguities_cycles: used_wide_lane_cycles,
           metadata:
             Map.merge(sol.metadata, %{
               integer_method: :widelane_narrowlane_lambda,
               wide_lane_fixed: true,
               wide_lane_ambiguities_cycles: used_wide_lane_cycles,
               dropped_cycle_slip_sats: slip_meta.dropped_sats,
               split_cycle_slip_arcs: slip_meta.split_arcs
             })
       }}
    end
  end

  def solve_widelane_fixed_baseline_epochs(_base_position, _dual_epochs, _opts),
    do: {:error, :invalid_epochs}

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
      ref_sd_id = single_difference_ambiguity_id(reference_sat, ref_base, ref_rover)
      ref_dd = %{satellite_id: reference_sat, ambiguity_id: ref_sd_id}

      dds =
        common
        |> Enum.reject(&(&1 == reference_sat))
        |> Enum.map(fn sat ->
          base_obs = Map.fetch!(base, sat)
          rover_obs = Map.fetch!(rover, sat)
          sat_sd_id = single_difference_ambiguity_id(sat, base_obs, rover_obs)

          %{
            satellite_id: sat,
            reference_satellite_id: reference_sat,
            ambiguity_id: double_difference_ambiguity_id(sat, sat_sd_id, ref_dd),
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

  defp normalize_dual_baseline_epochs(epochs) do
    epochs
    |> Enum.with_index()
    |> Enum.reduce_while({:ok, []}, fn {epoch, idx}, {:ok, acc} ->
      case normalize_dual_baseline_epoch(epoch, idx) do
        {:ok, normalized} -> {:cont, {:ok, [normalized | acc]}}
        {:error, _reason} = err -> {:halt, err}
      end
    end)
    |> case do
      {:ok, normalized} -> {:ok, Enum.reverse(normalized)}
      {:error, _reason} = err -> err
    end
  end

  defp normalize_dual_baseline_epoch(
         %{
           base_observations: base_observations,
           rover_observations: rover_observations,
           satellite_positions_m: satellite_positions
         } = epoch,
         idx
       )
       when is_list(base_observations) and is_list(rover_observations) and
              is_map(satellite_positions) do
    with {:ok, base} <-
           normalize_dual_observations(base_observations, :invalid_base_observations),
         {:ok, rover} <-
           normalize_dual_observations(rover_observations, :invalid_rover_observations),
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

  defp normalize_dual_baseline_epoch(_epoch, idx),
    do: {:error, {:invalid_epoch_observations, idx}}

  defp normalize_dual_observations(observations, error_tag) do
    observations
    |> Enum.reduce_while({:ok, %{}}, fn observation, {:ok, acc} ->
      case normalize_dual_observation(observation) do
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

  defp normalize_dual_observation(
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
    with {:ok, ambiguity_id} <- normalize_observation_ambiguity_id(obs, sat),
         {:ok, lli1} <- normalize_dual_lli(obs, :lli1),
         {:ok, lli2} <- normalize_dual_lli(obs, :lli2) do
      {:ok,
       %{
         satellite_id: sat,
         ambiguity_id: ambiguity_id,
         p1_m: p1 / 1.0,
         p2_m: p2 / 1.0,
         phi1_cyc: phi1 / 1.0,
         phi2_cyc: phi2 / 1.0,
         f1_hz: f1 / 1.0,
         f2_hz: f2 / 1.0,
         lli1: lli1,
         lli2: lli2
       }}
    end
  end

  defp normalize_dual_observation(_observation), do: {:error, :invalid_observation}

  defp normalize_dual_lli(obs, key) do
    case Map.get(obs, key) do
      nil -> {:ok, nil}
      value when is_integer(value) -> {:ok, value}
      _other -> {:error, :invalid_observation}
    end
  end

  defp prepare_dual_baseline_cycle_slips(epochs, opts) do
    with {:ok, policy} <- rtk_cycle_slip_policy(opts) do
      slips = dual_cycle_slip_events(epochs, opts)

      case {policy, slips} do
        {_policy, []} ->
          {:ok, epochs, empty_cycle_slip_meta()}

        {:error, [slip | _]} ->
          {:error,
           {:cycle_slip_detected, slip.receiver, slip.satellite_id, slip.epoch, slip.reasons}}

        {:drop_satellite, slips} ->
          dropped = slips |> Enum.map(& &1.satellite_id) |> Enum.uniq() |> Enum.sort()

          {:ok, drop_dual_cycle_slip_sats(epochs, dropped),
           %{dropped_sats: dropped, split_arcs: []}}

        {:split_arc, slips} ->
          split_epochs = split_dual_cycle_slip_arcs(epochs, slips)

          {:ok, split_epochs,
           %{dropped_sats: [], split_arcs: dual_cycle_slip_split_metadata(split_epochs, slips)}}
      end
    end
  end

  defp dual_cycle_slip_events(epochs, opts) do
    dual_cycle_slip_events_for_receiver(:base, epochs, opts) ++
      dual_cycle_slip_events_for_receiver(:rover, epochs, opts)
  end

  defp dual_cycle_slip_events_for_receiver(receiver, epochs, opts) do
    epochs
    |> Enum.flat_map(fn epoch ->
      epoch
      |> dual_receiver_observations(receiver)
      |> Enum.map(fn {sat, obs} -> {sat, epoch.epoch, obs} end)
    end)
    |> Enum.group_by(fn {sat, _epoch, _obs} -> sat end)
    |> Enum.sort_by(fn {sat, _arc} -> sat end)
    |> Enum.flat_map(fn {sat, arc} ->
      arc =
        Enum.sort_by(arc, fn {_sat, epoch, _obs} -> inspect(epoch) end)

      arc
      |> dual_carrier_phase_arc()
      |> CarrierPhase.detect_cycle_slips(opts)
      |> Enum.filter(& &1.slip)
      |> Enum.map(fn slip ->
        %{
          receiver: receiver,
          satellite_id: sat,
          epoch: slip.epoch,
          reasons: slip.reasons
        }
      end)
    end)
  end

  defp dual_carrier_phase_arc(arc) do
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

  defp drop_dual_cycle_slip_sats(epochs, dropped) do
    dropped = MapSet.new(dropped)

    Enum.map(epochs, fn epoch ->
      %{
        epoch
        | base: Map.drop(epoch.base, MapSet.to_list(dropped)),
          rover: Map.drop(epoch.rover, MapSet.to_list(dropped)),
          positions: Map.drop(epoch.positions, MapSet.to_list(dropped))
      }
    end)
  end

  defp split_dual_cycle_slip_arcs(epochs, slips) do
    split_sides = MapSet.new(slips, &{&1.receiver, &1.satellite_id})
    slip_epochs = MapSet.new(slips, &{&1.receiver, &1.satellite_id, &1.epoch})

    {split_epochs, _segments} =
      Enum.map_reduce(epochs, %{}, fn epoch, segments ->
        {base, segments} =
          split_dual_receiver_cycle_slip_arcs(
            :base,
            epoch.epoch,
            epoch.base,
            split_sides,
            slip_epochs,
            segments
          )

        {rover, segments} =
          split_dual_receiver_cycle_slip_arcs(
            :rover,
            epoch.epoch,
            epoch.rover,
            split_sides,
            slip_epochs,
            segments
          )

        {%{epoch | base: base, rover: rover}, segments}
      end)

    split_epochs
  end

  defp split_dual_receiver_cycle_slip_arcs(
         receiver,
         epoch,
         observations,
         split_sides,
         slip_epochs,
         segments
       ) do
    observations
    |> Enum.sort_by(fn {sat, _obs} -> sat end)
    |> Enum.reduce({%{}, segments}, fn {sat, obs}, {acc, segments} ->
      key = {receiver, sat}

      if MapSet.member?(split_sides, key) do
        current_segment = Map.get(segments, key, 1)

        segment =
          if MapSet.member?(slip_epochs, {receiver, sat, epoch}),
            do: current_segment + 1,
            else: current_segment

        ambiguity_id = split_side_ambiguity_id(sat, receiver, segment)
        obs = %{obs | ambiguity_id: ambiguity_id}

        {Map.put(acc, sat, obs), Map.put(segments, key, segment)}
      else
        {Map.put(acc, sat, obs), segments}
      end
    end)
  end

  defp dual_cycle_slip_split_metadata(epochs, slips) do
    split_sides = MapSet.new(slips, &{&1.receiver, &1.satellite_id})

    epochs
    |> Enum.flat_map(fn epoch ->
      dual_split_metadata_entries(:base, epoch.epoch, epoch.base, split_sides) ++
        dual_split_metadata_entries(:rover, epoch.epoch, epoch.rover, split_sides)
    end)
    |> Enum.group_by(fn entry -> {entry.receiver, entry.satellite_id, entry.ambiguity_id} end)
    |> Enum.map(fn {{receiver, sat, ambiguity_id}, entries} ->
      sorted = Enum.sort_by(entries, &inspect(&1.epoch))

      %{
        receiver: receiver,
        satellite_id: sat,
        ambiguity_id: ambiguity_id,
        start_epoch: hd(sorted).epoch,
        end_epoch: List.last(sorted).epoch,
        n_epochs: length(sorted)
      }
    end)
    |> Enum.sort_by(&{&1.satellite_id, &1.receiver, &1.ambiguity_id})
  end

  defp dual_split_metadata_entries(receiver, epoch, observations, split_sides) do
    observations
    |> Enum.sort_by(fn {sat, _obs} -> sat end)
    |> Enum.flat_map(fn {sat, obs} ->
      if MapSet.member?(split_sides, {receiver, sat}) do
        [
          %{
            receiver: receiver,
            satellite_id: sat,
            ambiguity_id: obs.ambiguity_id,
            epoch: epoch
          }
        ]
      else
        []
      end
    end)
  end

  defp estimate_dual_baseline_wide_lanes(epochs, reference_sat, opts) do
    with {:ok, wl_opts} <- dual_wide_lane_options(opts),
         {:ok, samples} <- dual_wide_lane_samples(epochs, reference_sat) do
      split_arc? = Keyword.get(opts, :on_cycle_slip, @default_cycle_slip_policy) == :split_arc

      samples
      |> Enum.group_by(& &1.ambiguity_id)
      |> Enum.sort_by(fn {ambiguity_id, _samples} -> ambiguity_id end)
      |> Enum.reduce_while({:ok, %{}}, fn {ambiguity_id, grouped}, {:ok, acc} ->
        cycles = Enum.map(grouped, & &1.cycles)

        case estimate_dual_wide_lane_integer(ambiguity_id, cycles, wl_opts) do
          {:ok, fixed} ->
            {:cont, {:ok, Map.put(acc, ambiguity_id, fixed)}}

          {:error, {:too_few_wide_lane_epochs, _id, _got, _required}} when split_arc? ->
            {:cont, {:ok, acc}}

          {:error, _} = err ->
            {:halt, err}
        end
      end)
    end
  end

  defp dual_wide_lane_options(opts) do
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

  defp dual_wide_lane_samples(epochs, reference_sat) do
    epochs
    |> Enum.reduce_while({:ok, []}, fn epoch, {:ok, acc} ->
      common = dual_epoch_common_sats(epoch)

      if reference_sat in common do
        ref_sd = dual_single_difference(epoch, reference_sat)

        common
        |> Enum.reject(&(&1 == reference_sat))
        |> Enum.reduce_while({:ok, acc}, fn sat, {:ok, samples} ->
          case dual_wide_lane_double_difference(epoch, sat, ref_sd) do
            {:ok, sample} -> {:cont, {:ok, [sample | samples]}}
            {:error, _} = err -> {:halt, err}
          end
        end)
        |> case do
          {:ok, samples} -> {:cont, {:ok, samples}}
          {:error, _} = err -> {:halt, err}
        end
      else
        {:halt, {:error, {:reference_satellite_missing, reference_sat}}}
      end
    end)
    |> case do
      {:ok, samples} -> {:ok, Enum.reverse(samples)}
      {:error, _} = err -> err
    end
  end

  defp dual_epoch_common_sats(epoch) do
    epoch.base
    |> Map.keys()
    |> MapSet.new()
    |> MapSet.intersection(epoch.rover |> Map.keys() |> MapSet.new())
    |> MapSet.intersection(epoch.positions |> Map.keys() |> MapSet.new())
    |> MapSet.to_list()
    |> Enum.sort()
  end

  defp dual_single_difference(epoch, sat) do
    base_obs = Map.fetch!(epoch.base, sat)
    rover_obs = Map.fetch!(epoch.rover, sat)

    %{
      satellite_id: sat,
      ambiguity_id: single_difference_ambiguity_id(sat, base_obs, rover_obs),
      wide_lane_cycles:
        dual_observation_wide_lane_cycles(rover_obs) - dual_observation_wide_lane_cycles(base_obs)
    }
  end

  defp dual_wide_lane_double_difference(epoch, sat, ref_sd) do
    sd = dual_single_difference(epoch, sat)

    {:ok,
     %{
       satellite_id: sat,
       reference_satellite_id: ref_sd.satellite_id,
       ambiguity_id: double_difference_ambiguity_id(sat, sd.ambiguity_id, ref_sd),
       cycles: sd.wide_lane_cycles - ref_sd.wide_lane_cycles
     }}
  rescue
    MatchError -> {:error, {:wide_lane_failed, sat, :equal_frequencies}}
    ArithmeticError -> {:error, {:wide_lane_failed, sat, :equal_frequencies}}
  end

  defp dual_observation_wide_lane_cycles(obs) do
    {:ok, mw_m} =
      CarrierPhase.melbourne_wubbena(
        obs.phi1_cyc,
        obs.phi2_cyc,
        obs.p1_m,
        obs.p2_m,
        obs.f1_hz,
        obs.f2_hz
      )

    {:ok, lambda_wl} = CarrierPhase.wide_lane_wavelength(obs.f1_hz, obs.f2_hz)
    mw_m / lambda_wl
  end

  defp estimate_dual_wide_lane_integer(ambiguity_id, cycles, opts) do
    if length(cycles) < opts.min_epochs do
      {:error, {:too_few_wide_lane_epochs, ambiguity_id, length(cycles), opts.min_epochs}}
    else
      mean = Enum.sum(cycles) / length(cycles)
      fixed = round(mean)

      if abs(mean - fixed) <= opts.tolerance_cycles do
        {:ok, fixed}
      else
        {:error, {:wide_lane_not_integer, ambiguity_id, mean, fixed}}
      end
    end
  end

  defp ionosphere_free_baseline_epochs(dual_epochs, reference_sat, wide_lane_cycles) do
    with {:ok, params} <- dual_narrow_lane_params(dual_epochs, reference_sat, wide_lane_cycles),
         {:ok, if_epochs} <-
           dual_ionosphere_free_baseline_epochs(dual_epochs, reference_sat, wide_lane_cycles) do
      wavelengths = Map.new(params, fn {id, p} -> {id, p.wavelength_m} end)
      offsets = Map.new(params, fn {id, p} -> {id, p.offset_m} end)
      {:ok, if_epochs, wavelengths, offsets}
    end
  end

  defp dual_narrow_lane_params(dual_epochs, reference_sat, wide_lane_cycles) do
    dual_epochs
    |> Enum.reduce_while({:ok, %{}}, fn epoch, {:ok, acc} ->
      ref_sd = dual_single_difference(epoch, reference_sat)

      epoch
      |> dual_epoch_common_sats()
      |> Enum.reject(&(&1 == reference_sat))
      |> Enum.reduce_while({:ok, acc}, fn sat, {:ok, params_acc} ->
        with {:ok, sample} <- dual_wide_lane_double_difference(epoch, sat, ref_sd),
             {:ok, wide_lane} <- fetch_dual_wide_lane(wide_lane_cycles, sample.ambiguity_id),
             {:ok, params} <-
               dual_narrow_lane_param_from_epoch(
                 epoch,
                 sat,
                 reference_sat,
                 sample.ambiguity_id,
                 wide_lane
               ),
             :ok <-
               ensure_consistent_dual_narrow_lane_params(
                 sample.ambiguity_id,
                 params,
                 Map.get(params_acc, sample.ambiguity_id)
               ) do
          {:cont, {:ok, Map.put_new(params_acc, sample.ambiguity_id, params)}}
        else
          {:error, {:missing_wide_lane_ambiguity, _id}} ->
            {:cont, {:ok, params_acc}}

          {:error, _} = err ->
            {:halt, err}
        end
      end)
      |> case do
        {:ok, params_acc} -> {:cont, {:ok, params_acc}}
        {:error, _} = err -> {:halt, err}
      end
    end)
  end

  defp fetch_dual_wide_lane(wide_lane_cycles, ambiguity_id) do
    case Map.fetch(wide_lane_cycles, ambiguity_id) do
      {:ok, wide_lane} -> {:ok, wide_lane}
      :error -> {:error, {:missing_wide_lane_ambiguity, ambiguity_id}}
    end
  end

  defp dual_narrow_lane_param_from_epoch(
         epoch,
         sat,
         reference_sat,
         ambiguity_id,
         wide_lane_cycles
       ) do
    sat_base = Map.fetch!(epoch.base, sat)
    sat_rover = Map.fetch!(epoch.rover, sat)
    ref_base = Map.fetch!(epoch.base, reference_sat)
    ref_rover = Map.fetch!(epoch.rover, reference_sat)

    with :ok <-
           ensure_same_dual_frequencies(ambiguity_id, [sat_base, sat_rover, ref_base, ref_rover]) do
      dual_narrow_lane_param(sat_base.f1_hz, sat_base.f2_hz, wide_lane_cycles)
    end
  end

  defp ensure_same_dual_frequencies(ambiguity_id, [first | rest]) do
    if Enum.all?(
         rest,
         &(same_frequency?(&1.f1_hz, first.f1_hz) and same_frequency?(&1.f2_hz, first.f2_hz))
       ) do
      :ok
    else
      {:error, {:inconsistent_frequencies, ambiguity_id}}
    end
  end

  defp dual_narrow_lane_param(f1, f2, wide_lane_cycles) do
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

  defp ensure_consistent_dual_narrow_lane_params(_id, _params, nil), do: :ok

  defp ensure_consistent_dual_narrow_lane_params(ambiguity_id, params, prev) do
    if same_frequency?(params.f1_hz, prev.f1_hz) and same_frequency?(params.f2_hz, prev.f2_hz) do
      :ok
    else
      {:error, {:inconsistent_frequencies, ambiguity_id}}
    end
  end

  defp dual_ionosphere_free_baseline_epochs(dual_epochs, reference_sat, wide_lane_cycles) do
    dual_epochs
    |> Enum.reduce_while({:ok, []}, fn epoch, {:ok, acc} ->
      keep_sats = dual_ionosphere_free_keep_sats(epoch, reference_sat, wide_lane_cycles)

      if length(keep_sats) < 2 do
        {:cont, {:ok, acc}}
      else
        with {:ok, base_obs} <- dual_ionosphere_free_observations(epoch.base, keep_sats),
             {:ok, rover_obs} <- dual_ionosphere_free_observations(epoch.rover, keep_sats) do
          {:cont,
           {:ok,
            [
              %{
                epoch: epoch.epoch,
                base_observations: base_obs,
                rover_observations: rover_obs,
                satellite_positions_m: Map.take(epoch.positions, keep_sats)
              }
              | acc
            ]}}
        else
          {:error, _} = err -> {:halt, err}
        end
      end
    end)
    |> case do
      {:ok, []} -> {:error, :no_epochs}
      {:ok, rows} -> {:ok, Enum.reverse(rows)}
      {:error, _} = err -> err
    end
  end

  defp dual_ionosphere_free_keep_sats(epoch, reference_sat, wide_lane_cycles) do
    common = dual_epoch_common_sats(epoch)

    if reference_sat in common do
      ref_sd = dual_single_difference(epoch, reference_sat)

      kept_nonrefs =
        common
        |> Enum.reject(&(&1 == reference_sat))
        |> Enum.filter(fn sat ->
          case dual_wide_lane_double_difference(epoch, sat, ref_sd) do
            {:ok, sample} -> Map.has_key?(wide_lane_cycles, sample.ambiguity_id)
            {:error, _} -> false
          end
        end)

      if kept_nonrefs == [], do: [], else: [reference_sat | kept_nonrefs]
    else
      []
    end
  end

  defp dual_ionosphere_free_observations(observations, keep_sats) do
    keep = MapSet.new(keep_sats)

    observations
    |> Enum.sort_by(fn {sat, _obs} -> sat end)
    |> Enum.reduce_while({:ok, []}, fn {sat, obs}, {:ok, acc} ->
      if MapSet.member?(keep, sat) do
        case dual_ionosphere_free_observation(obs) do
          {:ok, converted} -> {:cont, {:ok, [converted | acc]}}
          {:error, _} = err -> {:halt, err}
        end
      else
        {:cont, {:ok, acc}}
      end
    end)
    |> case do
      {:ok, obs} -> {:ok, Enum.reverse(obs)}
      {:error, _} = err -> err
    end
  end

  defp dual_ionosphere_free_observation(obs) do
    with {:ok, code_m} <- IonosphereFree.iono_free(obs.p1_m, obs.p2_m, obs.f1_hz, obs.f2_hz),
         {:ok, phase_m} <-
           IonosphereFree.iono_free_phase_cycles(
             obs.phi1_cyc,
             obs.phi2_cyc,
             obs.f1_hz,
             obs.f2_hz
           ) do
      {:ok,
       %{
         satellite_id: obs.satellite_id,
         ambiguity_id: obs.ambiguity_id,
         code_m: code_m,
         phase_m: phase_m
       }}
    else
      {:error, reason} -> {:error, {:ionosphere_free_failed, obs.satellite_id, reason}}
    end
  end

  defp dual_receiver_observations(epoch, :base), do: epoch.base
  defp dual_receiver_observations(epoch, :rover), do: epoch.rover

  defp same_frequency?(a, b), do: abs(a - b) <= 1.0e-6

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

  defp prepare_epochs_for_cycle_slips(epochs, opts) do
    with {:ok, policy} <- rtk_cycle_slip_policy(opts) do
      slips = cycle_slip_events(epochs)

      case {policy, slips} do
        {_policy, []} ->
          {:ok, epochs, empty_cycle_slip_meta()}

        {:error, [slip | _]} ->
          {:error,
           {:cycle_slip_detected, slip.receiver, slip.satellite_id, slip.epoch, slip.reasons}}

        {:drop_satellite, slips} ->
          dropped = slips |> Enum.map(& &1.satellite_id) |> Enum.uniq() |> Enum.sort()
          {:ok, drop_cycle_slip_sats(epochs, dropped), %{dropped_sats: dropped, split_arcs: []}}

        {:split_arc, slips} ->
          split_keys = MapSet.new(slips, &{&1.receiver, &1.satellite_id})
          split_epochs = split_cycle_slip_arcs(epochs, split_keys)

          {:ok, split_epochs,
           %{dropped_sats: [], split_arcs: cycle_slip_split_metadata(split_epochs, split_keys)}}
      end
    end
  end

  defp empty_cycle_slip_meta, do: %{dropped_sats: [], split_arcs: []}

  defp prepare_epochs_for_code_smoothing(epochs, opts) do
    with {:ok, smoothing} <- rtk_code_smoothing(opts) do
      case smoothing do
        :none ->
          {:ok, epochs, %{code_smoothing: false, code_smoothing_window_cap: nil}}

        {:hatch, cap} ->
          epochs =
            epochs
            |> smooth_receiver_codes(:base, cap)
            |> smooth_receiver_codes(:rover, cap)

          {:ok, epochs, %{code_smoothing: true, code_smoothing_window_cap: cap}}
      end
    end
  end

  defp rtk_code_smoothing(opts) do
    case Keyword.get(opts, :code_smoothing, false) do
      value when value in [false, nil] ->
        {:ok, :none}

      value when value in [true, :hatch] ->
        cap = Keyword.get(opts, :hatch_window_cap, @default_hatch_window_cap)

        if is_integer(cap) and cap >= 1 do
          {:ok, {:hatch, cap}}
        else
          {:error, {:invalid_option, :hatch_window_cap}}
        end

      _other ->
        {:error, {:invalid_option, :code_smoothing}}
    end
  end

  defp smooth_receiver_codes(epochs, receiver, cap) do
    {smoothed, _states} =
      Enum.map_reduce(epochs, %{}, fn epoch, states ->
        {observations, states} =
          epoch
          |> receiver_observations(receiver)
          |> Enum.sort_by(fn {sat, _obs} -> sat end)
          |> Enum.map_reduce(states, fn {sat, obs}, state_acc ->
            key = {receiver, obs.ambiguity_id}
            state = Map.get(state_acc, key)
            {obs, next_state} = smooth_observation_code(obs, state, cap)

            state_acc =
              if is_nil(next_state),
                do: Map.delete(state_acc, key),
                else: Map.put(state_acc, key, next_state)

            {{sat, obs}, state_acc}
          end)

        {put_receiver_observations(epoch, receiver, Map.new(observations)), states}
      end)

    smoothed
  end

  defp receiver_observations(epoch, :base), do: epoch.base
  defp receiver_observations(epoch, :rover), do: epoch.rover

  defp put_receiver_observations(epoch, :base, observations), do: %{epoch | base: observations}
  defp put_receiver_observations(epoch, :rover, observations), do: %{epoch | rover: observations}

  defp smooth_observation_code(%{code_m: code_m, phase_m: phase_m} = obs, state, cap)
       when is_number(code_m) and is_number(phase_m) do
    if is_nil(state) or lli_set?(obs.lli) do
      {%{obs | code_m: code_m}, %{p_smooth: code_m, phase_m: phase_m, window: 1}}
    else
      %{p_smooth: prev_smooth, phase_m: prev_phase, window: prev_window} = state
      n = min(prev_window + 1, cap)
      p_smooth = code_m / n + (n - 1) / n * (prev_smooth + (phase_m - prev_phase))
      {%{obs | code_m: p_smooth}, %{p_smooth: p_smooth, phase_m: phase_m, window: n}}
    end
  end

  defp smooth_observation_code(obs, _state, _cap), do: {obs, nil}

  defp rtk_cycle_slip_policy(opts) do
    case Keyword.get(opts, :on_cycle_slip, @default_cycle_slip_policy) do
      :error -> {:ok, :error}
      :drop_satellite -> {:ok, :drop_satellite}
      :split_arc -> {:ok, :split_arc}
      _other -> {:error, {:invalid_option, :on_cycle_slip}}
    end
  end

  defp cycle_slip_events(epochs) do
    Enum.flat_map(epochs, fn epoch ->
      cycle_slip_events_for_receiver(:base, epoch.epoch, epoch.base) ++
        cycle_slip_events_for_receiver(:rover, epoch.epoch, epoch.rover)
    end)
  end

  defp cycle_slip_events_for_receiver(receiver, epoch, observations) do
    observations
    |> Enum.sort_by(fn {sat, _obs} -> sat end)
    |> Enum.flat_map(fn {sat, obs} ->
      if lli_set?(obs.lli) do
        [%{receiver: receiver, satellite_id: sat, epoch: epoch, reasons: [:lli]}]
      else
        []
      end
    end)
  end

  defp drop_cycle_slip_sats(epochs, dropped) do
    dropped = MapSet.new(dropped)

    Enum.map(epochs, fn epoch ->
      %{
        epoch
        | base: Map.drop(epoch.base, MapSet.to_list(dropped)),
          rover: Map.drop(epoch.rover, MapSet.to_list(dropped)),
          positions: Map.drop(epoch.positions, MapSet.to_list(dropped))
      }
    end)
  end

  defp split_cycle_slip_arcs(epochs, split_keys) do
    {split_epochs, _segments} =
      Enum.map_reduce(epochs, %{}, fn epoch, segments ->
        {base, segments} =
          split_receiver_cycle_slip_arcs(:base, epoch.epoch, epoch.base, split_keys, segments)

        {rover, segments} =
          split_receiver_cycle_slip_arcs(:rover, epoch.epoch, epoch.rover, split_keys, segments)

        {%{epoch | base: base, rover: rover}, segments}
      end)

    split_epochs
  end

  defp split_receiver_cycle_slip_arcs(receiver, _epoch, observations, split_keys, segments) do
    observations
    |> Enum.sort_by(fn {sat, _obs} -> sat end)
    |> Enum.reduce({%{}, segments}, fn {sat, obs}, {acc, segments} ->
      key = {receiver, sat}

      if MapSet.member?(split_keys, key) do
        current_segment = Map.get(segments, key, 1)
        segment = if lli_set?(obs.lli), do: current_segment + 1, else: current_segment
        ambiguity_id = split_side_ambiguity_id(sat, receiver, segment)
        obs = %{obs | ambiguity_id: ambiguity_id}

        {Map.put(acc, sat, obs), Map.put(segments, key, segment)}
      else
        {Map.put(acc, sat, obs), segments}
      end
    end)
  end

  defp split_side_ambiguity_id(sat, receiver, segment), do: "#{sat}@#{receiver}##{segment}"

  defp cycle_slip_split_metadata(epochs, split_keys) do
    epochs
    |> Enum.flat_map(fn epoch ->
      split_metadata_entries(:base, epoch.epoch, epoch.base, split_keys) ++
        split_metadata_entries(:rover, epoch.epoch, epoch.rover, split_keys)
    end)
    |> Enum.group_by(fn entry -> {entry.receiver, entry.satellite_id, entry.ambiguity_id} end)
    |> Enum.map(fn {{receiver, sat, ambiguity_id}, entries} ->
      sorted = Enum.sort_by(entries, &inspect(&1.epoch))

      %{
        receiver: receiver,
        satellite_id: sat,
        ambiguity_id: ambiguity_id,
        start_epoch: hd(sorted).epoch,
        end_epoch: List.last(sorted).epoch,
        n_epochs: length(sorted)
      }
    end)
    |> Enum.sort_by(&{&1.satellite_id, &1.receiver, &1.ambiguity_id})
  end

  defp split_metadata_entries(receiver, epoch, observations, split_keys) do
    observations
    |> Enum.sort_by(fn {sat, _obs} -> sat end)
    |> Enum.flat_map(fn {sat, obs} ->
      if MapSet.member?(split_keys, {receiver, sat}) do
        [
          %{
            receiver: receiver,
            satellite_id: sat,
            ambiguity_id: obs.ambiguity_id,
            epoch: epoch
          }
        ]
      else
        []
      end
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

  defp all_epoch_sats(epochs) do
    all =
      epochs
      |> Enum.reduce(MapSet.new(), fn epoch, acc ->
        epoch.base
        |> Map.keys()
        |> MapSet.new()
        |> MapSet.intersection(epoch.rover |> Map.keys() |> MapSet.new())
        |> MapSet.intersection(epoch.positions |> Map.keys() |> MapSet.new())
        |> MapSet.union(acc)
      end)
      |> MapSet.to_list()
      |> Enum.sort()

    {:ok, all}
  end

  defp epoch_available_nonrefs(epoch, reference_sat, physical_sats \\ nil) do
    available =
      epoch.base
      |> Map.keys()
      |> MapSet.new()
      |> MapSet.intersection(epoch.rover |> Map.keys() |> MapSet.new())
      |> MapSet.intersection(epoch.positions |> Map.keys() |> MapSet.new())
      |> MapSet.delete(reference_sat)

    case physical_sats do
      nil -> available
      sats -> MapSet.intersection(available, MapSet.new(sats))
    end
    |> MapSet.to_list()
    |> Enum.sort()
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

  defp integer_options(opts) do
    radius =
      Keyword.get(opts, :integer_search_radius_cycles, @default_integer_search_radius_cycles)

    ratio = Keyword.get(opts, :integer_ratio_threshold, @default_integer_ratio_threshold)
    limit = Keyword.get(opts, :integer_candidate_limit, @default_integer_candidate_limit)
    partial? = Keyword.get(opts, :partial_ambiguity_resolution, false)
    partial_min = Keyword.get(opts, :partial_min_ambiguities, @default_partial_min_ambiguities)

    cond do
      not is_integer(radius) or radius < 0 ->
        {:error, {:invalid_option, :integer_search_radius_cycles}}

      not is_number(ratio) or ratio < 0.0 ->
        {:error, {:invalid_option, :integer_ratio_threshold}}

      not is_integer(limit) or limit < 1 ->
        {:error, {:invalid_option, :integer_candidate_limit}}

      not is_boolean(partial?) ->
        {:error, {:invalid_option, :partial_ambiguity_resolution}}

      not is_integer(partial_min) or partial_min < 1 ->
        {:error, {:invalid_option, :partial_min_ambiguities}}

      true ->
        {:ok,
         %{
           radius_cycles: radius,
           ratio_threshold: ratio / 1.0,
           candidate_limit: limit,
           partial_ambiguity_resolution?: partial?,
           partial_min_ambiguities: partial_min
         }}
    end
  end

  defp baseline_weights(opts) do
    with {:ok, code_sigma_m} <- measurement_sigma(opts, :code_sigma_m, @default_code_sigma_m),
         {:ok, phase_sigma_m} <- measurement_sigma(opts, :phase_sigma_m, @default_phase_sigma_m),
         {:ok, elevation_weighting?} <- elevation_weighting(opts) do
      {:ok,
       %{
         code_sigma_m: code_sigma_m,
         phase_sigma_m: phase_sigma_m,
         elevation_weighting?: elevation_weighting?
       }}
    end
  end

  defp measurement_sigma(opts, key, default) do
    sigma = Keyword.get(opts, key, default)

    if is_number(sigma) and sigma > 0.0,
      do: {:ok, sigma / 1.0},
      else: {:error, {:invalid_sigma, key}}
  end

  defp elevation_weighting(opts) do
    case Keyword.get(opts, :elevation_weighting, false) do
      value when is_boolean(value) -> {:ok, value}
      _other -> {:error, {:invalid_option, :elevation_weighting}}
    end
  end

  defp initial_baseline(opts) do
    opts
    |> Keyword.get(:initial_baseline_m, {0.0, 0.0, 0.0})
    |> Types.normalize_ecef(:invalid_initial_baseline)
  end

  defp baseline_ambiguity_index(epochs, common_sats, reference_sat) do
    physical_sats = Enum.reject(common_sats, &(&1 == reference_sat))

    epochs
    |> Enum.reduce_while({:ok, %{}}, fn epoch, {:ok, acc} ->
      ref_dd = single_difference(epoch, reference_sat)

      epoch
      |> epoch_available_nonrefs(reference_sat, physical_sats)
      |> Enum.reduce_while({:ok, acc}, fn sat, {:ok, acc} ->
        ambiguity_id = double_difference_measurement(epoch, sat, ref_dd).ambiguity_id

        case Map.fetch(acc, ambiguity_id) do
          {:ok, ^sat} ->
            {:cont, {:ok, acc}}

          {:ok, other_sat} ->
            {:halt, {:error, {:duplicate_ambiguity_id, ambiguity_id, other_sat, sat}}}

          :error ->
            {:cont, {:ok, Map.put(acc, ambiguity_id, sat)}}
        end
      end)
      |> case do
        {:ok, acc} -> {:cont, {:ok, acc}}
        {:error, _reason} = err -> {:halt, err}
      end
    end)
    |> case do
      {:ok, ambiguity_satellites} ->
        ambiguity_ids =
          ambiguity_satellites
          |> Enum.sort_by(fn {ambiguity_id, sat} -> {sat, ambiguity_id} end)
          |> Enum.map(&elem(&1, 0))

        {:ok, ambiguity_ids, Map.new(ambiguity_satellites)}

      {:error, _reason} = err ->
        err
    end
  end

  defp ambiguity_wavelengths(ambiguity_ids, ambiguity_satellites, opts) do
    case Keyword.fetch(opts, :ambiguity_wavelength_m) do
      {:ok, wavelength} when is_number(wavelength) and wavelength > 0.0 ->
        {:ok, Map.new(ambiguity_ids, &{&1, wavelength / 1.0})}

      {:ok, wavelengths} when is_map(wavelengths) ->
        ambiguity_ids
        |> Enum.reduce_while({:ok, %{}}, fn sat, {:ok, acc} ->
          physical_sat = Map.fetch!(ambiguity_satellites, sat)

          case ambiguity_wavelength(wavelengths, sat, physical_sat) do
            {:ok, wavelength} when is_number(wavelength) and wavelength > 0.0 ->
              {:cont, {:ok, Map.put(acc, sat, wavelength / 1.0)}}

            _other ->
              {:halt, {:error, {:invalid_ambiguity_wavelength, sat}}}
          end
        end)

      {:ok, _other} ->
        {:error, {:invalid_option, :ambiguity_wavelength_m}}

      :error ->
        {:error, :ambiguity_wavelength_required}
    end
  end

  defp ambiguity_offsets(ambiguity_ids, ambiguity_satellites, opts) do
    case Keyword.fetch(opts, :ambiguity_offset_m) do
      {:ok, offset} when is_number(offset) ->
        {:ok, Map.new(ambiguity_ids, &{&1, offset / 1.0})}

      {:ok, offsets} when is_map(offsets) ->
        ambiguity_ids
        |> Enum.reduce_while({:ok, %{}}, fn sat, {:ok, acc} ->
          physical_sat = Map.fetch!(ambiguity_satellites, sat)

          case ambiguity_value(offsets, sat, physical_sat) do
            {:ok, offset} when is_number(offset) ->
              {:cont, {:ok, Map.put(acc, sat, offset / 1.0)}}

            _other ->
              {:halt, {:error, {:invalid_ambiguity_offset, sat}}}
          end
        end)

      {:ok, _other} ->
        {:error, {:invalid_option, :ambiguity_offset_m}}

      :error ->
        {:ok, Map.new(ambiguity_ids, &{&1, 0.0})}
    end
  end

  defp ambiguity_wavelength(wavelengths, ambiguity_id, physical_sat) do
    case ambiguity_value(wavelengths, ambiguity_id, physical_sat) do
      {:ok, _wavelength} = ok -> ok
      :error -> :error
    end
  end

  defp ambiguity_value(values, ambiguity_id, physical_sat) do
    case Map.fetch(values, ambiguity_id) do
      {:ok, _value} = ok -> ok
      :error -> Map.fetch(values, physical_sat)
    end
  end

  defp baseline_row_count(epochs, reference_sat) do
    epochs
    |> Enum.map(fn epoch -> length(epoch_available_nonrefs(epoch, reference_sat)) end)
    |> Enum.sum()
    |> Kernel.*(2)
  end

  defp baseline_unknown_count(ambiguity_ids), do: 3 + length(ambiguity_ids)

  defp iterate_baseline(
         base,
         epochs,
         reference_sat,
         physical_sats,
         ambiguity_ids,
         ambiguity_satellites,
         state,
         weights,
         opts,
         dropped_sats,
         slip_meta,
         iter
       ) do
    with {:ok, rows} <-
           build_baseline_rows(
             base,
             epochs,
             reference_sat,
             physical_sats,
             ambiguity_ids,
             state,
             weights
           ),
         {:ok, dx} <-
           solve_baseline_normal_equations(rows, baseline_unknown_count(ambiguity_ids)) do
      next = apply_baseline_delta(state, ambiguity_ids, dx)
      {baseline_step, ambiguity_step} = baseline_step_norms(dx)

      cond do
        baseline_step <= opts.position_tolerance_m and
            ambiguity_step <= opts.ambiguity_tolerance_m ->
          finalize_baseline(
            base,
            epochs,
            reference_sat,
            physical_sats,
            ambiguity_ids,
            ambiguity_satellites,
            next,
            weights,
            dropped_sats,
            slip_meta,
            iter,
            true,
            :state_tolerance
          )

        iter >= opts.max_iterations ->
          finalize_baseline(
            base,
            epochs,
            reference_sat,
            physical_sats,
            ambiguity_ids,
            ambiguity_satellites,
            next,
            weights,
            dropped_sats,
            slip_meta,
            iter,
            false,
            :max_iterations
          )

        true ->
          iterate_baseline(
            base,
            epochs,
            reference_sat,
            physical_sats,
            ambiguity_ids,
            ambiguity_satellites,
            next,
            weights,
            opts,
            dropped_sats,
            slip_meta,
            iter + 1
          )
      end
    end
  end

  defp build_baseline_rows(
         base,
         epochs,
         reference_sat,
         physical_sats,
         ambiguity_ids,
         state,
         weights
       ) do
    epochs
    |> Enum.reduce_while({:ok, []}, fn epoch, {:ok, acc} ->
      case build_epoch_baseline_rows(
             base,
             epoch,
             reference_sat,
             physical_sats,
             ambiguity_ids,
             state,
             weights
           ) do
        {:ok, rows} -> {:cont, {:ok, rows ++ acc}}
        {:error, _reason} = err -> {:halt, err}
      end
    end)
    |> case do
      {:ok, rows} -> {:ok, Enum.reverse(rows)}
      {:error, _reason} = err -> err
    end
  end

  defp build_epoch_baseline_rows(
         base,
         epoch,
         reference_sat,
         physical_sats,
         ambiguity_ids,
         state,
         weights
       ) do
    ref_dd = single_difference(epoch, reference_sat)
    ref_pos = Map.fetch!(epoch.positions, reference_sat)

    epoch
    |> epoch_available_nonrefs(reference_sat, physical_sats)
    |> Enum.reduce({:ok, []}, fn sat, {:ok, acc} ->
      obs_dd = double_difference_measurement(epoch, sat, ref_dd)
      sat_pos = Map.fetch!(epoch.positions, sat)
      {geom_dd, deriv} = geometry_double_difference(base, state.baseline, sat_pos, ref_pos)
      ambiguity_id = obs_dd.ambiguity_id
      ambiguity = Map.fetch!(state.ambiguities, ambiguity_id)
      h_base = design_baseline_row(deriv, nil, ambiguity_ids)
      h_phase = design_baseline_row(deriv, ambiguity_id, ambiguity_ids)
      code_variance = row_variance(weights, :code, base, sat_pos, ref_pos)
      phase_variance = row_variance(weights, :phase, base, sat_pos, ref_pos)

      code_row = %{
        kind: :code,
        epoch_idx: epoch.idx,
        epoch: epoch.epoch,
        sat: sat,
        ambiguity_id: ambiguity_id,
        h: h_base,
        y: obs_dd.code_m - geom_dd,
        sd_variance_m2: code_variance.sd_variance_m2,
        ref_sd_variance_m2: code_variance.ref_sd_variance_m2,
        weight: code_variance.weight
      }

      phase_row = %{
        kind: :phase,
        epoch_idx: epoch.idx,
        epoch: epoch.epoch,
        sat: sat,
        ambiguity_id: ambiguity_id,
        h: h_phase,
        y: obs_dd.phase_m - (geom_dd + ambiguity),
        sd_variance_m2: phase_variance.sd_variance_m2,
        ref_sd_variance_m2: phase_variance.ref_sd_variance_m2,
        weight: phase_variance.weight
      }

      {:ok, [phase_row, code_row | acc]}
    end)
    |> case do
      {:ok, rows} -> {:ok, Enum.reverse(rows)}
      {:error, _reason} = err -> err
    end
  end

  defp single_difference(epoch, sat) do
    base_obs = Map.fetch!(epoch.base, sat)
    rover_obs = Map.fetch!(epoch.rover, sat)

    %{
      satellite_id: sat,
      code_m: rover_obs.code_m - base_obs.code_m,
      phase_m: rover_obs.phase_m - base_obs.phase_m,
      ambiguity_id: single_difference_ambiguity_id(sat, base_obs, rover_obs)
    }
  end

  defp double_difference_measurement(epoch, sat, ref_dd) do
    sd = single_difference(epoch, sat)

    %{
      code_m: sd.code_m - ref_dd.code_m,
      phase_m: sd.phase_m - ref_dd.phase_m,
      ambiguity_id: double_difference_ambiguity_id(sat, sd.ambiguity_id, ref_dd)
    }
  end

  defp single_difference_ambiguity_id(sat, base_obs, rover_obs) do
    case {base_obs.ambiguity_id, rover_obs.ambiguity_id} do
      {^sat, ^sat} -> sat
      {^sat, rover_id} -> rover_id
      {base_id, ^sat} -> base_id
      {same_id, same_id} -> same_id
      {base_id, rover_id} -> "#{sat}:base=#{base_id},rover=#{rover_id}"
    end
  end

  defp double_difference_ambiguity_id(sat, sat_sd_id, ref_dd) do
    if sat_sd_id == sat and ref_dd.ambiguity_id == ref_dd.satellite_id,
      do: sat,
      else: "#{sat_sd_id}|ref=#{ref_dd.ambiguity_id}"
  end

  defp geometry_double_difference(base, baseline, sat_pos, ref_pos) do
    rover = add3(base, baseline)
    sat_sd = range(sat_pos, rover) - range(sat_pos, base)
    ref_sd = range(ref_pos, rover) - range(ref_pos, base)
    sat_deriv = range_derivative(rover, sat_pos)
    ref_deriv = range_derivative(rover, ref_pos)
    {sat_sd - ref_sd, sub3(sat_deriv, ref_deriv)}
  end

  defp row_variance(weights, kind, base, sat_pos, ref_pos) do
    sigma =
      case kind do
        :code -> weights.code_sigma_m
        :phase -> weights.phase_sigma_m
      end

    sat_sd_variance =
      single_difference_variance(sigma, weights.elevation_weighting?, base, sat_pos)

    ref_sd_variance =
      single_difference_variance(sigma, weights.elevation_weighting?, base, ref_pos)

    dd_variance = sat_sd_variance + ref_sd_variance

    %{
      sd_variance_m2: sat_sd_variance,
      ref_sd_variance_m2: ref_sd_variance,
      weight: 1.0 / :math.sqrt(dd_variance)
    }
  end

  defp single_difference_variance(sigma_m, false, _base, _sat_pos) do
    2.0 * sigma_m * sigma_m
  end

  defp single_difference_variance(sigma_m, true, base, sat_pos) do
    sin_el = elevation_sin(base, sat_pos) |> max(@min_elevation_sin)
    scaled_sigma = sigma_m / sin_el
    2.0 * scaled_sigma * scaled_sigma
  end

  defp design_baseline_row({dx, dy, dz}, ambiguity_id, ambiguity_ids) do
    ambiguity_cols = Enum.map(ambiguity_ids, &if(&1 == ambiguity_id, do: 1.0, else: 0.0))
    [dx, dy, dz | ambiguity_cols]
  end

  defp apply_baseline_delta(state, ambiguity_ids, dx) do
    {bx, by, bz} = state.baseline

    ambiguities =
      ambiguity_ids
      |> Enum.with_index()
      |> Map.new(fn {ambiguity_id, idx} ->
        {ambiguity_id, Map.fetch!(state.ambiguities, ambiguity_id) + Enum.at(dx, 3 + idx)}
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
         physical_sats,
         ambiguity_ids,
         ambiguity_satellites,
         state,
         weights,
         dropped_sats,
         slip_meta,
         iterations,
         converged?,
         status
       ) do
    with {:ok, rows} <-
           build_baseline_rows(
             base,
             epochs,
             reference_sat,
             physical_sats,
             ambiguity_ids,
             state,
             weights
           ),
         {:ok, covariance_m} <-
           baseline_ambiguity_covariance(
             rows,
             baseline_unknown_count(ambiguity_ids),
             3,
             length(ambiguity_ids)
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
         used_sats: ambiguity_ids,
         ambiguities_m: state.ambiguities,
         residuals_m: residuals,
         metadata: %{
           iterations: iterations,
           converged: converged?,
           status: status,
           physical_sats: physical_sats,
           ambiguity_satellites: ambiguity_satellites,
           ambiguity_float: %{
             order: ambiguity_ids,
             covariance_m: covariance_m,
             covariance_inverse_m: covariance_inverse_m
           },
           measurement_covariance: %{
             model: :double_difference,
             code_sigma_m: weights.code_sigma_m,
             phase_sigma_m: weights.phase_sigma_m,
             elevation_weighting: weights.elevation_weighting?,
             min_elevation_sin: @min_elevation_sin
           },
           code_rms_m: rms(code_residuals),
           phase_rms_m: rms(phase_residuals),
           weighted_rms_m: weighted_rms(rows),
           n_epochs: length(epochs),
           n_observations: length(rows),
           dropped_sats: Enum.uniq(dropped_sats ++ slip_meta.dropped_sats) |> Enum.sort(),
           dropped_cycle_slip_sats: slip_meta.dropped_sats,
           split_cycle_slip_arcs: slip_meta.split_arcs,
           code_smoothing: Map.get(slip_meta, :code_smoothing, false),
           code_smoothing_window_cap: Map.get(slip_meta, :code_smoothing_window_cap)
         }
       }}
    end
  end

  defp search_baseline_ambiguities(float_sol, wavelengths, offsets, integer_opts) do
    covariance_cycles =
      covariance_m_to_cycles(
        float_sol.metadata.ambiguity_float.covariance_m,
        float_sol.used_sats,
        wavelengths
      )

    float_cycles =
      Map.new(float_sol.used_sats, fn sat ->
        ambiguity_m = Map.fetch!(float_sol.ambiguities_m, sat)
        offset_m = Map.fetch!(offsets, sat)
        {sat, (ambiguity_m - offset_m) / Map.fetch!(wavelengths, sat)}
      end)

    with {:ok, fixed_cycles, meta} <-
           IntegerLeastSquares.search(float_cycles, covariance_cycles, integer_opts) do
      meta = Map.put(meta, :ambiguity_offsets_m, offsets)

      if meta.integer_status == :fixed or not integer_opts.partial_ambiguity_resolution? do
        {:ok, fixed_cycles, Map.merge(meta, partial_meta(false, false, fixed_cycles, []))}
      else
        search_partial_baseline_ambiguities(
          float_sol.used_sats,
          float_cycles,
          covariance_cycles,
          offsets,
          integer_opts,
          fixed_cycles,
          meta
        )
      end
    end
  end

  defp partial_meta(enabled?, fixed?, fixed_cycles, free_ambiguities, extra \\ %{}) do
    Map.merge(
      %{
        partial_ambiguity_resolution: enabled?,
        partial_fixed: fixed?,
        partial_fixed_ambiguities: fixed_cycles |> Map.keys() |> Enum.sort(),
        partial_free_ambiguities: Enum.sort(free_ambiguities)
      },
      extra
    )
  end

  defp search_partial_baseline_ambiguities(
         ambiguity_ids,
         float_cycles,
         covariance_cycles,
         offsets,
         integer_opts,
         full_fixed_cycles,
         full_meta
       ) do
    ranked_ids =
      ambiguity_ids_ranked_by_integer_confidence(ambiguity_ids, float_cycles, covariance_cycles)

    min_size = min(integer_opts.partial_min_ambiguities, length(ambiguity_ids) - 1)

    if min_size < 1 do
      {:ok, full_fixed_cycles,
       Map.merge(
         full_meta,
         partial_meta(true, false, full_fixed_cycles, [], %{
           ambiguity_offsets_m: offsets,
           partial_full_set: full_set_integer_summary(full_meta)
         })
       )}
    else
      (length(ambiguity_ids) - 1)..min_size//-1
      |> Enum.reduce_while(:not_fixed, fn subset_size, :not_fixed ->
        subset_ids = ranked_ids |> Enum.take(subset_size) |> Enum.sort()

        case search_ambiguity_subset(
               subset_ids,
               ambiguity_ids,
               float_cycles,
               covariance_cycles,
               offsets,
               integer_opts,
               full_meta
             ) do
          {:ok, _fixed_cycles, %{integer_status: :fixed}} = ok -> {:halt, ok}
          {:ok, _fixed_cycles, _meta} -> {:cont, :not_fixed}
          {:error, _reason} = err -> {:halt, err}
        end
      end)
      |> case do
        {:ok, _fixed_cycles, _meta} = ok ->
          ok

        :not_fixed ->
          # The confidence-ranked nested-take is a greedy strategy: it keeps the
          # highest-margin ambiguities and only ever drops the worst-ranked
          # tail. On a dual-frequency narrow-lane arc the highest-margin
          # ambiguity can be the very one whose inclusion destroys the ratio
          # test, so the greedy ranking can step over a safe larger subset that
          # excludes it. When the greedy ranking finds nothing, fall back to a
          # bounded largest-first exhaustive combination search and accept the
          # highest-ratio subset of the largest size that passes the (unchanged)
          # ratio threshold. Single-frequency arcs reach a passing greedy subset
          # first, so they never enter this fallback and keep their historical
          # behavior.
          search_partial_baseline_ambiguities_exhaustive(
            ambiguity_ids,
            float_cycles,
            covariance_cycles,
            offsets,
            integer_opts,
            full_fixed_cycles,
            full_meta,
            min_size
          )
      end
    end
  end

  # Maximum ambiguity count for which the largest-first exhaustive subset
  # fallback runs. Above this the combination count grows too large to be worth
  # the search, and the historical not-fixed result is preserved.
  @partial_exhaustive_max_ambiguities 14

  # Hard global ceiling on the number of candidate subsets the largest-first
  # exhaustive fallback evaluates, so worst-case work is bounded directly rather
  # than only indirectly via the ambiguity-count cap. Covers a full N=14 search
  # (~16k subsets); if the budget is exhausted before a passing subset is found,
  # the safe historical not-fixed result stands. The actual count evaluated is
  # reported in `metadata.partial_exhaustive_subsets_evaluated`.
  @partial_exhaustive_max_subsets 20_000

  defp search_partial_baseline_ambiguities_exhaustive(
         ambiguity_ids,
         float_cycles,
         covariance_cycles,
         offsets,
         integer_opts,
         full_fixed_cycles,
         full_meta,
         min_size
       ) do
    not_fixed = fn evaluated ->
      {:ok, full_fixed_cycles,
       Map.merge(
         full_meta,
         partial_meta(true, false, full_fixed_cycles, [], %{
           ambiguity_offsets_m: offsets,
           partial_full_set: full_set_integer_summary(full_meta),
           partial_exhaustive_subsets_evaluated: evaluated
         })
       )}
    end

    if length(ambiguity_ids) > @partial_exhaustive_max_ambiguities do
      not_fixed.(0)
    else
      result =
        (length(ambiguity_ids) - 1)..min_size//-1
        |> Enum.reduce_while({:not_fixed, 0}, fn subset_size, {:not_fixed, evaluated} ->
          combos = ambiguity_subset_combinations(ambiguity_ids, subset_size)

          if evaluated + length(combos) > @partial_exhaustive_max_subsets do
            # Evaluating this size would exceed the global subset-work budget;
            # stop and keep the safe historical not-fixed result.
            {:halt, {:not_fixed, evaluated}}
          else
            best_passing =
              Enum.reduce(combos, nil, fn subset_ids, best ->
                subset_ids = Enum.sort(subset_ids)

                case search_ambiguity_subset(
                       subset_ids,
                       ambiguity_ids,
                       float_cycles,
                       covariance_cycles,
                       offsets,
                       integer_opts,
                       full_meta
                     ) do
                  {:ok, _fixed_cycles, %{integer_status: :fixed, integer_ratio: ratio}} = ok ->
                    case best do
                      {best_ratio, _} when best_ratio >= ratio -> best
                      _ -> {ratio, ok}
                    end

                  _ ->
                    best
                end
              end)

            evaluated = evaluated + length(combos)

            case best_passing do
              {_ratio, {:ok, fixed_cycles, meta}} ->
                meta = Map.put(meta, :partial_exhaustive_subsets_evaluated, evaluated)
                {:halt, {:fixed, {:ok, fixed_cycles, meta}}}

              nil ->
                {:cont, {:not_fixed, evaluated}}
            end
          end
        end)

      case result do
        {:fixed, ok} -> ok
        {:not_fixed, evaluated} -> not_fixed.(evaluated)
      end
    end
  end

  defp ambiguity_subset_combinations(_ids, 0), do: [[]]
  defp ambiguity_subset_combinations([], _size), do: []

  defp ambiguity_subset_combinations([head | tail], size) do
    with_head =
      tail
      |> ambiguity_subset_combinations(size - 1)
      |> Enum.map(&[head | &1])

    with_head ++ ambiguity_subset_combinations(tail, size)
  end

  defp search_ambiguity_subset(
         subset_ids,
         all_ids,
         float_cycles,
         covariance_cycles,
         offsets,
         integer_opts,
         full_meta
       ) do
    subset_cycles = Map.take(float_cycles, subset_ids)
    subset_covariance = covariance_submatrix(all_ids, covariance_cycles, subset_ids)

    case IntegerLeastSquares.search(subset_cycles, subset_covariance, integer_opts) do
      {:ok, fixed_cycles, meta} ->
        free_ids = all_ids -- subset_ids

        meta =
          meta
          |> Map.put(:ambiguity_offsets_m, Map.take(offsets, subset_ids))
          |> Map.merge(
            partial_meta(true, meta.integer_status == :fixed, fixed_cycles, free_ids, %{
              partial_full_set: full_set_integer_summary(full_meta)
            })
          )

        {:ok, fixed_cycles, meta}

      {:error, _reason} = err ->
        err
    end
  end

  defp ambiguity_ids_ranked_by_integer_confidence(ambiguity_ids, float_cycles, covariance_cycles) do
    ambiguity_ids
    |> Enum.with_index()
    |> Enum.sort_by(fn {id, idx} ->
      variance = covariance_cycles |> Enum.at(idx) |> Enum.at(idx)

      distance_to_integer =
        abs(Map.fetch!(float_cycles, id) - round(Map.fetch!(float_cycles, id)))

      sigma = :math.sqrt(max(variance, 0.0))
      # Larger margin-to-boundary per sigma is more reliable; sort ascending.
      normalized_margin = (0.5 - distance_to_integer) / sigma
      {-normalized_margin, variance, id}
    end)
    |> Enum.map(&elem(&1, 0))
  end

  defp covariance_submatrix(all_ids, covariance, subset_ids) do
    indices = Enum.map(subset_ids, fn id -> Enum.find_index(all_ids, &(&1 == id)) end)

    for i <- indices do
      for j <- indices, do: covariance |> Enum.at(i) |> Enum.at(j)
    end
  end

  defp full_set_integer_summary(meta) do
    %{
      integer_status: meta.integer_status,
      integer_ratio: meta.integer_ratio,
      integer_best_score: meta.integer_best_score,
      integer_second_best_score: meta.integer_second_best_score,
      integer_candidates: meta.integer_candidates,
      order: meta.ambiguity_search.order
    }
  end

  defp free_ambiguity_ids(ambiguity_ids, fixed_cycles) do
    fixed = fixed_cycles |> Map.keys() |> MapSet.new()
    Enum.reject(ambiguity_ids, &MapSet.member?(fixed, &1))
  end

  defp covariance_m_to_cycles(covariance_m, sat_ids, wavelengths) do
    for i <- 0..(length(sat_ids) - 1) do
      sat_i = Enum.at(sat_ids, i)
      lambda_i = Map.fetch!(wavelengths, sat_i)

      for j <- 0..(length(sat_ids) - 1) do
        sat_j = Enum.at(sat_ids, j)
        lambda_j = Map.fetch!(wavelengths, sat_j)
        (covariance_m |> Enum.at(i) |> Enum.at(j)) / (lambda_i * lambda_j)
      end
    end
  end

  defp fixed_ambiguities_m(fixed_cycles, wavelengths, offsets) do
    Map.new(fixed_cycles, fn {sat, cycles} ->
      {sat, Map.fetch!(offsets, sat) + cycles * Map.fetch!(wavelengths, sat)}
    end)
  end

  defp iterate_fixed_baseline(
         base,
         epochs,
         reference_sat,
         physical_sats,
         ambiguity_ids,
         free_ambiguity_ids,
         fixed_m,
         state,
         weights,
         opts,
         float_sol,
         fixed_cycles,
         fixed_meta,
         iter
       ) do
    with {:ok, rows} <-
           build_fixed_baseline_rows(
             base,
             epochs,
             reference_sat,
             physical_sats,
             ambiguity_ids,
             free_ambiguity_ids,
             fixed_m,
             state,
             weights
           ),
         {:ok, dx} <-
           solve_baseline_normal_equations(rows, 3 + length(free_ambiguity_ids)) do
      next = apply_fixed_baseline_delta(state, free_ambiguity_ids, dx)
      {baseline_step, ambiguity_step} = fixed_baseline_step_norms(dx)

      cond do
        baseline_step <= opts.position_tolerance_m and
            ambiguity_step <= opts.ambiguity_tolerance_m ->
          finalize_fixed_baseline(
            base,
            epochs,
            reference_sat,
            physical_sats,
            ambiguity_ids,
            free_ambiguity_ids,
            fixed_m,
            next,
            weights,
            float_sol,
            fixed_cycles,
            fixed_meta,
            iter,
            true,
            :state_tolerance
          )

        iter >= opts.max_iterations ->
          finalize_fixed_baseline(
            base,
            epochs,
            reference_sat,
            physical_sats,
            ambiguity_ids,
            free_ambiguity_ids,
            fixed_m,
            next,
            weights,
            float_sol,
            fixed_cycles,
            fixed_meta,
            iter,
            false,
            :max_iterations
          )

        true ->
          iterate_fixed_baseline(
            base,
            epochs,
            reference_sat,
            physical_sats,
            ambiguity_ids,
            free_ambiguity_ids,
            fixed_m,
            next,
            weights,
            opts,
            float_sol,
            fixed_cycles,
            fixed_meta,
            iter + 1
          )
      end
    end
  end

  defp build_fixed_baseline_rows(
         base,
         epochs,
         reference_sat,
         physical_sats,
         ambiguity_ids,
         free_ambiguity_ids,
         fixed_m,
         state,
         weights
       ) do
    epochs
    |> Enum.reduce_while({:ok, []}, fn epoch, {:ok, acc} ->
      case build_fixed_epoch_baseline_rows(
             base,
             epoch,
             reference_sat,
             physical_sats,
             ambiguity_ids,
             free_ambiguity_ids,
             fixed_m,
             state,
             weights
           ) do
        {:ok, rows} -> {:cont, {:ok, rows ++ acc}}
        {:error, _reason} = err -> {:halt, err}
      end
    end)
    |> case do
      {:ok, rows} -> {:ok, Enum.reverse(rows)}
      {:error, _reason} = err -> err
    end
  end

  defp build_fixed_epoch_baseline_rows(
         base,
         epoch,
         reference_sat,
         physical_sats,
         _ambiguity_ids,
         free_ambiguity_ids,
         fixed_m,
         state,
         weights
       ) do
    ref_dd = single_difference(epoch, reference_sat)
    ref_pos = Map.fetch!(epoch.positions, reference_sat)

    epoch
    |> epoch_available_nonrefs(reference_sat, physical_sats)
    |> Enum.reduce({:ok, []}, fn sat, {:ok, acc} ->
      obs_dd = double_difference_measurement(epoch, sat, ref_dd)
      sat_pos = Map.fetch!(epoch.positions, sat)
      {geom_dd, deriv} = geometry_double_difference(base, state.baseline, sat_pos, ref_pos)
      ambiguity_id = obs_dd.ambiguity_id

      {phase_ambiguity, phase_h} =
        case Map.fetch(fixed_m, ambiguity_id) do
          {:ok, fixed} ->
            {fixed, fixed_design_baseline_row(deriv, nil, free_ambiguity_ids)}

          :error ->
            {Map.fetch!(state.ambiguities, ambiguity_id),
             fixed_design_baseline_row(deriv, ambiguity_id, free_ambiguity_ids)}
        end

      code_h = fixed_design_baseline_row(deriv, nil, free_ambiguity_ids)
      code_variance = row_variance(weights, :code, base, sat_pos, ref_pos)
      phase_variance = row_variance(weights, :phase, base, sat_pos, ref_pos)

      code_row = %{
        kind: :code,
        epoch_idx: epoch.idx,
        epoch: epoch.epoch,
        sat: sat,
        ambiguity_id: ambiguity_id,
        h: code_h,
        y: obs_dd.code_m - geom_dd,
        sd_variance_m2: code_variance.sd_variance_m2,
        ref_sd_variance_m2: code_variance.ref_sd_variance_m2,
        weight: code_variance.weight
      }

      phase_row = %{
        kind: :phase,
        epoch_idx: epoch.idx,
        epoch: epoch.epoch,
        sat: sat,
        ambiguity_id: ambiguity_id,
        h: phase_h,
        y: obs_dd.phase_m - (geom_dd + phase_ambiguity),
        sd_variance_m2: phase_variance.sd_variance_m2,
        ref_sd_variance_m2: phase_variance.ref_sd_variance_m2,
        weight: phase_variance.weight
      }

      {:ok, [phase_row, code_row | acc]}
    end)
    |> case do
      {:ok, rows} -> {:ok, Enum.reverse(rows)}
      {:error, _reason} = err -> err
    end
  end

  defp fixed_design_baseline_row({dx, dy, dz}, ambiguity_id, free_ambiguity_ids) do
    ambiguity_cols =
      Enum.map(free_ambiguity_ids, &if(&1 == ambiguity_id, do: 1.0, else: 0.0))

    [dx, dy, dz | ambiguity_cols]
  end

  defp apply_fixed_baseline_delta(state, free_ambiguity_ids, [dx, dy, dz | ambiguity_deltas]) do
    {bx, by, bz} = state.baseline

    ambiguities =
      free_ambiguity_ids
      |> Enum.with_index()
      |> Map.new(fn {ambiguity_id, idx} ->
        {ambiguity_id,
         Map.fetch!(state.ambiguities, ambiguity_id) + Enum.at(ambiguity_deltas, idx)}
      end)

    %{state | baseline: {bx + dx, by + dy, bz + dz}, ambiguities: ambiguities}
  end

  defp fixed_baseline_step_norms([dx, dy, dz | ambiguity_deltas]) do
    baseline_step = :math.sqrt(dx * dx + dy * dy + dz * dz)
    ambiguity_step = Enum.map(ambiguity_deltas, &abs/1) |> Enum.max(fn -> 0.0 end)
    {baseline_step, ambiguity_step}
  end

  defp finalize_fixed_baseline(
         base,
         epochs,
         reference_sat,
         physical_sats,
         ambiguity_ids,
         free_ambiguity_ids,
         fixed_m,
         state,
         weights,
         float_sol,
         fixed_cycles,
         fixed_meta,
         iterations,
         converged?,
         status
       ) do
    with {:ok, rows} <-
           build_fixed_baseline_rows(
             base,
             epochs,
             reference_sat,
             physical_sats,
             ambiguity_ids,
             free_ambiguity_ids,
             fixed_m,
             state,
             weights
           ) do
      residuals = baseline_residuals(rows, reference_sat)
      code_residuals = Enum.map(residuals, & &1.code_m)
      phase_residuals = Enum.map(residuals, & &1.phase_m)
      rover = add3(base, state.baseline)

      {:ok,
       %FixedBaselineSolution{
         baseline_m: ecef_map(state.baseline),
         rover_position_m: ecef_map(rover),
         reference_satellite_id: reference_sat,
         used_sats: ambiguity_ids,
         fixed_ambiguities_cycles: fixed_cycles,
         fixed_ambiguities_m: fixed_m,
         wide_lane_ambiguities_cycles: Map.get(fixed_meta, :wide_lane_ambiguities_cycles),
         float_solution: float_sol,
         residuals_m: residuals,
         metadata:
           Map.merge(fixed_meta, %{
             iterations: iterations,
             converged: converged?,
             status: status,
             code_rms_m: rms(code_residuals),
             phase_rms_m: rms(phase_residuals),
             weighted_rms_m: weighted_rms(rows),
             n_epochs: length(epochs),
             n_observations: length(rows),
             physical_sats: physical_sats,
             ambiguity_satellites: Map.fetch!(float_sol.metadata, :ambiguity_satellites),
             dropped_cycle_slip_sats: Map.get(float_sol.metadata, :dropped_cycle_slip_sats, []),
             split_cycle_slip_arcs: Map.get(float_sol.metadata, :split_cycle_slip_arcs, []),
             measurement_covariance: %{
               model: :double_difference,
               code_sigma_m: weights.code_sigma_m,
               phase_sigma_m: weights.phase_sigma_m,
               elevation_weighting: weights.elevation_weighting?,
               min_elevation_sin: @min_elevation_sin
             }
           })
       }}
    end
  end

  defp baseline_residuals(rows, reference_sat) do
    rows
    |> Enum.group_by(&{&1.epoch, &1.sat, &1.ambiguity_id})
    |> Enum.map(fn {{epoch, sat, ambiguity_id}, grouped} ->
      code = Enum.find(grouped, &(&1.kind == :code))
      phase = Enum.find(grouped, &(&1.kind == :phase))

      %{
        epoch: epoch,
        satellite_id: sat,
        reference_satellite_id: reference_sat,
        ambiguity_id: ambiguity_id,
        code_m: code.y,
        phase_m: phase.y
      }
    end)
    |> Enum.sort_by(&{inspect(&1.epoch), &1.satellite_id, &1.ambiguity_id})
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

      %{
        rows: rows,
        inverse_covariance: double_difference_inverse_covariance(rows)
      }
    end)
  end

  # For one epoch and one measurement kind, DD_i = SD_i - SD_ref. Independent
  # receiver observations give each satellite single-difference variance
  # `row.sd_variance_m2`; all DD rows in the block share the same reference
  # single difference with variance `row.ref_sd_variance_m2`. With constant
  # sigmas this is the familiar v*(I + J) covariance; with elevation weighting,
  # the diagonal terms become satellite-specific but the reference covariance
  # remains shared by every row.
  defp double_difference_inverse_covariance(rows) do
    if constant_double_difference_variance?(rows) do
      equal_double_difference_inverse_covariance(length(rows), hd(rows).sd_variance_m2)
    else
      covariance = double_difference_covariance(rows)
      {:ok, inverse} = invert_matrix(covariance)
      inverse
    end
  end

  defp constant_double_difference_variance?(rows) do
    first = hd(rows).sd_variance_m2

    Enum.all?(rows, fn row ->
      row.sd_variance_m2 == first and row.ref_sd_variance_m2 == first
    end)
  end

  defp equal_double_difference_inverse_covariance(m, sd_variance_m2) do
    diagonal_scale = 1.0 / sd_variance_m2 * (1.0 - 1.0 / (m + 1.0))
    off_diagonal = -1.0 / (sd_variance_m2 * (m + 1.0))

    for i <- 0..(m - 1) do
      for j <- 0..(m - 1), do: if(i == j, do: diagonal_scale, else: off_diagonal)
    end
  end

  defp double_difference_covariance(rows) do
    ref_variance = rows |> hd() |> Map.fetch!(:ref_sd_variance_m2)

    rows
    |> Enum.with_index()
    |> Enum.map(fn {row, i} ->
      rows
      |> Enum.with_index()
      |> Enum.map(fn {_other, j} ->
        if i == j, do: row.sd_variance_m2 + ref_variance, else: ref_variance
      end)
    end)
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

  defp normalize_observation(%{satellite_id: sat, code_m: code, phase_m: phase} = obs)
       when is_binary(sat) and is_number(code) and is_number(phase) do
    with {:ok, ambiguity_id} <- normalize_observation_ambiguity_id(obs, sat),
         {:ok, lli} <- normalize_observation_lli(obs) do
      {:ok,
       %{
         satellite_id: sat,
         ambiguity_id: ambiguity_id,
         code_m: code / 1.0,
         phase_m: phase / 1.0,
         lli: lli
       }}
    end
  end

  defp normalize_observation({sat, code, phase})
       when is_binary(sat) and is_number(code) and is_number(phase) do
    {:ok,
     %{satellite_id: sat, ambiguity_id: sat, code_m: code / 1.0, phase_m: phase / 1.0, lli: nil}}
  end

  defp normalize_observation(_observation), do: {:error, :invalid_observation}

  defp normalize_observation_ambiguity_id(obs, sat) do
    case Map.get(obs, :ambiguity_id, sat) do
      ambiguity_id when is_binary(ambiguity_id) -> {:ok, ambiguity_id}
      _other -> {:error, :invalid_observation}
    end
  end

  defp normalize_observation_lli(obs) do
    lli = Map.get(obs, :lli, Map.get(obs, :loss_of_lock_indicator))

    case lli do
      nil -> {:ok, nil}
      value when is_integer(value) -> {:ok, value}
      _other -> {:error, :invalid_observation}
    end
  end

  defp lli_set?(lli) when is_integer(lli), do: band(lli, 1) == 1
  defp lli_set?(_lli), do: false

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

  defp elevation_sin(base, sat_pos) do
    {:ok, up} = local_up(base)

    case unit3(sub3(sat_pos, base)) do
      {:ok, los} -> dot3(los, up)
      :zero -> -1.0
    end
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
