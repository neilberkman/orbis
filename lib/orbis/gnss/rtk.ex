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
      matvec: 2,
      matrix_add: 2,
      matrix_sub: 2,
      solve_linear: 2,
      solve_matrix: 2,
      submatrix: 5,
      transpose: 1,
      zero_matrix: 1,
      zero_vector: 1
    ]

  alias Orbis.GNSS.{CarrierPhase, IonosphereFree}
  alias Orbis.GNSS.Core.{Constants, Types}
  alias Orbis.GNSS.Core.IntegerLeastSquares
  alias Orbis.NIF

  @default_max_iterations 8
  @default_position_tolerance_m 1.0e-4
  @default_ambiguity_tolerance_m 1.0e-4
  @default_code_sigma_m 1.0
  @default_phase_sigma_m 0.02
  @default_stochastic_model :simple
  @default_integer_search_radius_cycles 1
  @default_integer_ratio_threshold 3.0
  @default_integer_candidate_limit 200_000
  @default_partial_min_ambiguities 4
  @default_max_residual_exclusions 1
  @default_cycle_slip_policy :error
  @default_sagnac true
  @default_hatch_window_cap 100
  @default_filter_baseline_prior_sigma_m 100.0
  @default_filter_ambiguity_prior_sigma_m 1_000.0
  @default_filter_hold_sigma_m 1.0e-4
  @default_filter_process_noise_baseline_sigma_m 0.0
  @default_filter_dynamics_model :constant_position
  @default_filter_innovation_screen_min_rows 8
  @rtk_filter_state_version 3
  @min_elevation_sin 0.05
  @double_difference_options [:reference_satellite_id]
  @float_baseline_options [
    :reference_satellite_id,
    :initial_baseline_m,
    :code_sigma_m,
    :phase_sigma_m,
    :stochastic_model,
    :on_cycle_slip,
    :elevation_weighting,
    :sagnac,
    :elevation_mask_deg,
    :code_smoothing,
    :hatch_window_cap,
    :max_iterations,
    :position_tolerance_m,
    :ambiguity_tolerance_m
  ]
  @integer_baseline_options [
    :ambiguity_wavelength_m,
    :ambiguity_offset_m,
    :integer_search_radius_cycles,
    :integer_ratio_threshold,
    :integer_candidate_limit,
    :partial_ambiguity_resolution,
    :partial_min_ambiguities,
    :float_only_systems
  ]
  @residual_validation_options [:residual_threshold_sigma, :max_residual_exclusions]
  @fixed_baseline_options @float_baseline_options ++
                            @integer_baseline_options ++ @residual_validation_options
  @filter_baseline_options @fixed_baseline_options ++
                             [
                               :baseline_prior_sigma_m,
                               :ambiguity_prior_sigma_m,
                               :hold_sigma_m,
                               :process_noise_baseline_sigma_m,
                               :dynamics_model,
                               :innovation_screen_sigma,
                               :innovation_screen_min_rows,
                               :filter_kernel
                             ]
  @dual_wide_lane_options [
    :wide_lane_min_epochs,
    :wide_lane_tolerance_cycles,
    :gf_threshold_m,
    :mw_threshold_cycles
  ]
  @widelane_baseline_options (@fixed_baseline_options --
                                [:ambiguity_wavelength_m, :ambiguity_offset_m]) ++
                               @dual_wide_lane_options
  @widelane_delegate_drop_options @dual_wide_lane_options ++
                                    [:ambiguity_wavelength_m, :ambiguity_offset_m]

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
            phase_m: float(),
            code_sigma_m: float(),
            phase_sigma_m: float(),
            code_normalized: float(),
            phase_normalized: float()
          }

    @type t :: %__MODULE__{
            baseline_m: ecef(),
            rover_position_m: ecef(),
            reference_satellite_id: String.t() | %{String.t() => String.t()},
            used_sats: [String.t()],
            ambiguities_m: %{String.t() => float()},
            residuals_m: [residual()],
            metadata: %{
              iterations: pos_integer(),
              converged: boolean(),
              status: :state_tolerance | :max_iterations,
              physical_sats: [String.t()],
              reference_satellites: %{String.t() => String.t()},
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
                stochastic_model: :simple | :rtklib,
                elevation_weighting: boolean(),
                sagnac: boolean(),
                min_elevation_sin: float()
              },
              code_rms_m: float(),
              phase_rms_m: float(),
              weighted_rms_m: float(),
              n_epochs: pos_integer(),
              n_observations: pos_integer(),
              dropped_sats: [String.t()],
              dropped_cycle_slip_sats: [String.t()],
              elevation_mask_deg: float() | nil,
              elevation_masked_sats: [String.t()],
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
            reference_satellite_id: String.t() | %{String.t() => String.t()},
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
                stochastic_model: :simple | :rtklib,
                elevation_weighting: boolean(),
                sagnac: boolean(),
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
              optional(:reference_satellites) => %{String.t() => String.t()},
              optional(:ambiguity_satellites) => %{String.t() => String.t()},
              optional(:partial_ambiguity_resolution) => boolean(),
              optional(:partial_fixed) => boolean(),
              optional(:partial_fixed_ambiguities) => [String.t()],
              optional(:partial_free_ambiguities) => [String.t()],
              optional(:partial_full_set) => map(),
              optional(:dropped_cycle_slip_sats) => [String.t()],
              optional(:elevation_mask_deg) => float() | nil,
              optional(:elevation_masked_sats) => [String.t()],
              optional(:split_cycle_slip_arcs) => [map()]
            }
          }
  end

  defmodule FilterBaselineSolution do
    @moduledoc """
    Sequential RTK baseline-filter result.
    """

    @enforce_keys [
      :baseline_m,
      :rover_position_m,
      :reference_satellite_id,
      :fixed_ambiguities_cycles,
      :epochs,
      :metadata
    ]
    defstruct [
      :baseline_m,
      :rover_position_m,
      :reference_satellite_id,
      :fixed_ambiguities_cycles,
      :epochs,
      :metadata
    ]

    @type ecef :: %{x_m: float(), y_m: float(), z_m: float()}

    @type epoch_result :: %{
            epoch: term(),
            index: non_neg_integer(),
            baseline_m: ecef(),
            integer_status: :fixed | :not_fixed,
            integer_ratio: float() | :infinity | nil,
            integer_best_score: float() | nil,
            integer_second_best_score: float() | nil,
            integer_candidates: non_neg_integer() | nil,
            ambiguity_search: map() | nil,
            residuals_m: [FloatBaselineSolution.residual()],
            newly_fixed_ambiguities: [String.t()],
            fixed_ambiguities: [String.t()]
          }

    @type t :: %__MODULE__{
            baseline_m: ecef(),
            rover_position_m: ecef(),
            reference_satellite_id: String.t() | %{String.t() => String.t()},
            fixed_ambiguities_cycles: %{String.t() => integer()},
            epochs: [epoch_result()],
            metadata: map()
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
  `:satellite_positions_m` is used for satellite selection and elevation
  weighting. When the caller has receiver-specific transmit-time positions, it
  may also provide `:base_satellite_positions_m` and
  `:rover_satellite_positions_m`; otherwise both default to
  `:satellite_positions_m`.
  """
  @type baseline_epoch :: %{
          required(:base_observations) => [observation()],
          required(:rover_observations) => [observation()],
          required(:satellite_positions_m) => satellite_positions(),
          optional(:base_satellite_positions_m) => satellite_positions(),
          optional(:rover_satellite_positions_m) => satellite_positions(),
          optional(:velocity_mps) => ecef_input(),
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
          optional(:base_satellite_positions_m) => satellite_positions(),
          optional(:rover_satellite_positions_m) => satellite_positions(),
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
        # Optional when base/rover transmit-time satellite positions differ:
        base_satellite_positions_m: %{"G01" => {21.0e6, 14.0e6, 20.0e6}, ...},
        rover_satellite_positions_m: %{"G01" => {21.0e6, 14.0e6, 20.0e6}, ...},
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

    * `:reference_satellite_id` - fixed double-difference reference. Accepts a
      satellite id binary (single-system data only) or a per-system map such as
      `%{"G" => "G04", "E" => "E11"}` covering every observed system. When
      omitted, each system uses its highest-average-elevation satellite common
      to every epoch in which the system appears, with a lexicographic
      tie-break. Every non-reference satellite forms its double difference
      against its own system's reference; there are no cross-system double
      differences.
    * `:initial_baseline_m` - initial base-to-rover ECEF vector, default
      `{0.0, 0.0, 0.0}`.
    * `:code_sigma_m` / `:phase_sigma_m` - undifferenced receiver measurement
      sigmas in metres. The solver propagates them into the non-diagonal
      double-difference covariance where rows sharing the reference satellite
      are correlated. Defaults are `#{@default_code_sigma_m}` and
      `#{@default_phase_sigma_m}`.
    * `:stochastic_model` - `:simple` (default) uses constant sigmas, optionally
      scaled by `:elevation_weighting`. `:rtklib` uses RTKLIB's floor-plus-
      elevation single-difference variance shape, treating `:code_sigma_m` and
      `:phase_sigma_m` as the model's constant/elevation coefficients in
      metres.
    * `:on_cycle_slip` - what to do when a base or rover observation carries an
      LLI loss-of-lock bit: `:error` returns
      `{:error, {:cycle_slip_detected, receiver, sat, epoch, [:lli]}}`
      (default); `:drop_satellite` removes that satellite from the arc;
      `:split_arc` starts a new ambiguity arc at the slipped epoch.
    * `:elevation_weighting` - when `true`, scales each undifferenced
      measurement sigma by `1 / max(sin(elevation), #{@min_elevation_sin})`
      before propagating the double-difference covariance. Default `false`
      preserves the constant-sigma, transcendental-free solve path.
    * `:sagnac` - when `true` (default), applies the standard first-order
      Earth-rotation correction to each receiver-satellite range before
      forming double differences. Set `false` only for synthetic fixtures whose
      observations were generated from plain Euclidean range.
    * `:elevation_mask_deg` - optional elevation mask in degrees. Satellites
      below the mask at the base station are removed before reference selection
      and ambiguity construction.
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
    with :ok <- validate_options(opts, @float_baseline_options),
         {:ok, base} <- Types.normalize_ecef(base_position, :invalid_base_position),
         {:ok, normalized_epochs, prep_meta} <- prepare_baseline_epochs(base, epochs, opts),
         {:ok, all_sats} <- all_epoch_sats(normalized_epochs),
         :ok <- ensure_baseline_satellites(all_sats),
         {:ok, refs} <-
           baseline_reference_satellites(opts, base, normalized_epochs, all_sats),
         {:ok, solve_opts} <- baseline_solve_options(opts),
         {:ok, weights} <- baseline_weights(opts),
         {:ok, initial_baseline} <- initial_baseline(opts),
         {:ok, ambiguity_ids, ambiguity_satellites} <-
           baseline_ambiguity_index(normalized_epochs, all_sats, refs) do
      physical_sats = nonreference_sats(all_sats, refs)

      if baseline_row_count(normalized_epochs, refs) <
           baseline_unknown_count(ambiguity_ids) do
        {:error,
         {:underdetermined, baseline_row_count(normalized_epochs, refs),
          baseline_unknown_count(ambiguity_ids)}}
      else
        state = %{
          baseline: initial_baseline,
          ambiguities: Map.new(ambiguity_ids, &{&1, 0.0})
        }

        iterate_baseline(
          base,
          normalized_epochs,
          refs,
          physical_sats,
          ambiguity_ids,
          ambiguity_satellites,
          state,
          weights,
          solve_opts,
          [],
          prep_meta,
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
    * `:float_only_systems` - list of constellation letters (for example
      `["R"]`) whose double-difference ambiguities are never entered into the
      integer search; they contribute float measurement rows only. GLONASS is
      the canonical use: FDMA inter-channel biases break the clean DD integer
      assumption (default `[]`).
    * `:residual_threshold_sigma` - optional normalized-residual gate. When set,
      the float solve is checked before integer search; the worst offending
      satellite is excluded and the solve retried up to `:max_residual_exclusions`.
    * `:max_residual_exclusions` - maximum satellites the residual gate may
      exclude (default `#{@default_max_residual_exclusions}` when the residual
      gate is enabled).

  The fixed solution is returned even when the ratio test fails; in that case
  `metadata.integer_status` is `:not_fixed`.
  """
  @spec solve_fixed_baseline_epochs(ecef_input(), [baseline_epoch()], keyword()) ::
          {:ok, FixedBaselineSolution.t()} | {:error, term()}
  def solve_fixed_baseline_epochs(base_position, epochs, opts \\ [])

  def solve_fixed_baseline_epochs(base_position, epochs, opts) when is_list(epochs) do
    float_opts = Keyword.take(opts, @float_baseline_options)

    with :ok <- validate_options(opts, @fixed_baseline_options),
         {:ok, _float_only_systems} <- float_only_systems(opts),
         {:ok, residual_opts} <- residual_validation_options(opts) do
      solve_fixed_baseline_epochs_attempt(
        base_position,
        epochs,
        opts,
        float_opts,
        residual_opts,
        []
      )
    end
  end

  def solve_fixed_baseline_epochs(_base_position, _epochs, _opts), do: {:error, :invalid_epochs}

  @doc """
  Run a sequential static RTK baseline filter with per-epoch ambiguity fixing.

  This is the RTKLIB-style real-time path: it carries one static baseline and
  one single-difference ambiguity state per satellite arc across epochs,
  performs a correlated double-difference measurement update at each epoch,
  attempts integer fixing from the posterior covariance of the corresponding
  double-difference ambiguity combinations, and holds accepted integers as tight
  pseudo-measurements on those combinations in later epochs.

  Options are the fixed-baseline options plus the filter parameters below.
  `:partial_ambiguity_resolution` is deliberately rejected for this entry
  point: the sequential filter only holds a full-set fix until partial
  sequential AR has post-fix validation against real data.

    * `:baseline_prior_sigma_m` - initial baseline prior sigma in metres
      (default `#{@default_filter_baseline_prior_sigma_m}`).
    * `:ambiguity_prior_sigma_m` - initial ambiguity prior sigma in metres
      (default `#{@default_filter_ambiguity_prior_sigma_m}`).
    * `:hold_sigma_m` - pseudo-measurement sigma for fixed ambiguity holds
      (default `#{@default_filter_hold_sigma_m}`).
    * `:dynamics_model` - `:constant_position` (default) keeps the carried
      baseline mean fixed between epochs. `:velocity_propagated` advances the
      prediction mean by each epoch's optional `:velocity_mps` times elapsed
      seconds; process-noise meaning is unchanged.
    * `:innovation_screen_sigma` - optional predicted-residual screen in the
      Rust kernel. When set, epoch rows with `abs(innovation * weight)` above
      this value are rejected before the measurement update.
    * `:innovation_screen_min_rows` - minimum accepted row count for the
      innovation screen. If fewer rows survive, the epoch coasts on the
      predicted state (default `#{@default_filter_innovation_screen_min_rows}`).
    * `:filter_kernel` - `:rust` (default) or `:elixir` — the Elixir reference
      implementation, bit-identical to the kernel; every kernel capability is
      gated by === trace tests against it. The kernel carries the per-system
      references and honors `:float_only_systems`.

  Returns `{:ok, %FilterBaselineSolution{}}` or a tagged error.
  """
  @spec solve_filter_baseline_epochs(ecef_input(), [baseline_epoch()], keyword()) ::
          {:ok, FilterBaselineSolution.t()} | {:error, term()}
  def solve_filter_baseline_epochs(base_position, epochs, opts \\ [])

  def solve_filter_baseline_epochs(base_position, epochs, opts) when is_list(epochs) do
    float_opts = Keyword.take(opts, @float_baseline_options)

    with :ok <- validate_options(opts, @filter_baseline_options),
         {:ok, base} <- Types.normalize_ecef(base_position, :invalid_base_position),
         {:ok, filter_opts} <- sequential_filter_options(opts),
         {:ok, filter_kernel} <- filter_kernel(opts),
         {:ok, float_only_systems} <- float_only_systems(opts),
         {:ok, normalized_epochs, prep_meta} <- prepare_baseline_epochs(base, epochs, float_opts),
         {:ok, all_sats} <- all_epoch_sats(normalized_epochs),
         :ok <- ensure_baseline_satellites(all_sats),
         {:ok, refs} <-
           baseline_reference_satellites(opts, base, normalized_epochs, all_sats),
         {:ok, weights} <- baseline_weights(float_opts),
         {:ok, solve_opts} <- baseline_solve_options(float_opts),
         {:ok, initial_baseline} <- initial_baseline(float_opts),
         {:ok, sd_ambiguity_ids, sd_ambiguity_satellites} <-
           single_difference_ambiguity_index(normalized_epochs, all_sats),
         {:ok, dd_ambiguity_ids, dd_ambiguity_satellites, dd_ambiguity_pairs} <-
           sequential_dd_ambiguity_index(normalized_epochs, all_sats, refs),
         {:ok, wavelengths} <-
           ambiguity_wavelengths(dd_ambiguity_ids, dd_ambiguity_satellites, opts),
         {:ok, offsets} <- ambiguity_offsets(dd_ambiguity_ids, dd_ambiguity_satellites, opts),
         {:ok, integer_opts} <- integer_options(opts),
         :ok <- validate_filter_integer_options(integer_opts) do
      physical_sats = nonreference_sats(all_sats, refs)

      run_sequential_baseline_filter(
        base,
        normalized_epochs,
        refs,
        physical_sats,
        sd_ambiguity_ids,
        sd_ambiguity_satellites,
        dd_ambiguity_ids,
        dd_ambiguity_satellites,
        dd_ambiguity_pairs,
        float_only_systems,
        wavelengths,
        offsets,
        weights,
        solve_opts,
        integer_opts,
        filter_opts,
        initial_baseline,
        prep_meta,
        filter_kernel
      )
    end
  end

  def solve_filter_baseline_epochs(_base_position, _epochs, _opts), do: {:error, :invalid_epochs}

  defp validate_filter_integer_options(%{partial_ambiguity_resolution?: true}),
    do: {:error, {:unsupported_option, :partial_ambiguity_resolution}}

  defp validate_filter_integer_options(_integer_opts), do: :ok

  # The satellite system is the constellation letter, the first grapheme of the
  # RINEX satellite id ("G01" -> "G", "R12" -> "R").
  defp satellite_system(satellite_id), do: String.first(satellite_id)

  defp satellite_systems(sats) do
    sats
    |> Enum.map(&satellite_system/1)
    |> Enum.uniq()
    |> Enum.sort()
  end

  defp ensure_single_widelane_system(epochs) do
    epochs
    |> Enum.reduce_while(MapSet.new(), fn epoch, systems ->
      systems =
        epoch
        |> dual_observation_sats()
        |> Enum.reduce(systems, fn sat, acc -> MapSet.put(acc, satellite_system(sat)) end)

      if MapSet.size(systems) > 1 do
        {:halt, {:error, {:unsupported_widelane, :multi_gnss}}}
      else
        {:cont, systems}
      end
    end)
    |> case do
      {:error, _reason} = err -> err
      _systems -> :ok
    end
  end

  defp dual_observation_sats(epoch) do
    epoch.base
    |> Map.keys()
    |> MapSet.new()
    |> MapSet.union(epoch.rover |> Map.keys() |> MapSet.new())
    |> MapSet.to_list()
  end

  defp reference_satellite_set(refs), do: refs |> Map.values() |> MapSet.new()

  defp nonreference_sats(all_sats, refs) do
    ref_set = reference_satellite_set(refs)
    Enum.reject(all_sats, &MapSet.member?(ref_set, &1))
  end

  # Reported reference shape: a single-system solve keeps today's bare satellite
  # id; a multi-system solve reports the per-system reference map.
  defp reference_satellite_report(refs) when map_size(refs) == 1, do: refs |> Map.values() |> hd()

  defp reference_satellite_report(refs), do: refs

  defp float_only_systems(opts) do
    case Keyword.get(opts, :float_only_systems, []) do
      systems when is_list(systems) ->
        if Enum.all?(systems, &system_letter?/1),
          do: {:ok, systems},
          else: {:error, {:invalid_option, :float_only_systems}}

      _other ->
        {:error, {:invalid_option, :float_only_systems}}
    end
  end

  defp system_letter?(<<letter>>) when letter in ?A..?Z, do: true
  defp system_letter?(_other), do: false

  defp float_only_ambiguity_ids(ambiguity_satellites, []) when is_map(ambiguity_satellites),
    do: MapSet.new()

  defp float_only_ambiguity_ids(ambiguity_satellites, float_only_systems)
       when is_map(ambiguity_satellites) do
    ambiguity_satellites
    |> Enum.filter(fn {_ambiguity_id, sat} -> satellite_system(sat) in float_only_systems end)
    |> MapSet.new(fn {ambiguity_id, _sat} -> ambiguity_id end)
  end

  # Per-system common satellites: for each system, the satellites common to
  # every epoch in which that system appears. This is the single-reference
  # invariant ("the reference must be available in every epoch") applied per
  # system.
  defp per_system_common_sats(epochs) do
    epochs
    |> Enum.reduce(%{}, fn epoch, acc ->
      epoch
      |> epoch_sats()
      |> Enum.group_by(&satellite_system/1)
      |> Enum.reduce(acc, fn {system, sats}, acc ->
        sats = MapSet.new(sats)
        Map.update(acc, system, sats, &MapSet.intersection(&1, sats))
      end)
    end)
    |> Map.new(fn {system, sats} -> {system, sats |> MapSet.to_list() |> Enum.sort()} end)
  end

  defp per_system_epochs(epochs, system) do
    Enum.filter(epochs, fn epoch ->
      epoch |> epoch_sats() |> Enum.any?(&(satellite_system(&1) == system))
    end)
  end

  # Per-system double-difference reference selection.
  #
  #   * no option       — per system, the highest-average-elevation satellite of
  #     that system's common set (lexicographic tie-break), scored over the
  #     epochs in which the system appears;
  #   * binary (legacy) — valid only when a single system is observed;
  #   * map             — explicit per-system references covering every observed
  #     system.
  #
  # Returns `{:ok, %{"G" => "G04", ...}}`. With a single system this degenerates
  # to exactly today's selection (same candidate list, same scoring epochs).
  defp baseline_reference_satellites(opts, base, epochs, all_sats) do
    systems = satellite_systems(all_sats)
    common_by_system = per_system_common_sats(epochs)

    case Keyword.get(opts, :reference_satellite_id) do
      nil ->
        Enum.reduce_while(systems, {:ok, %{}}, fn system, {:ok, refs} ->
          case Map.get(common_by_system, system, []) do
            [] ->
              {:halt, {:error, {:no_common_reference_satellite, system}}}

            common ->
              reference =
                highest_elevation_reference(common, base, per_system_epochs(epochs, system))

              {:cont, {:ok, Map.put(refs, system, reference)}}
          end
        end)

      sat when is_binary(sat) ->
        case systems do
          [system] ->
            if sat in Map.get(common_by_system, system, []) do
              {:ok, %{system => sat}}
            else
              {:error, {:reference_satellite_missing, sat}}
            end

          _multiple ->
            {:error, {:reference_satellite_single_system, sat}}
        end

      refs when is_map(refs) ->
        Enum.reduce_while(systems, {:ok, %{}}, fn system, {:ok, acc} ->
          case Map.fetch(refs, system) do
            :error ->
              {:halt, {:error, {:reference_satellite_missing_system, system}}}

            {:ok, sat} when is_binary(sat) ->
              if sat in Map.get(common_by_system, system, []) do
                {:cont, {:ok, Map.put(acc, system, sat)}}
              else
                {:halt, {:error, {:reference_satellite_missing, sat}}}
              end

            {:ok, _other} ->
              {:halt, {:error, {:invalid_option, :reference_satellite_id}}}
          end
        end)

      _other ->
        {:error, {:invalid_option, :reference_satellite_id}}
    end
  end

  defp solve_fixed_baseline_epochs_attempt(
         base_position,
         epochs,
         opts,
         float_opts,
         residual_opts,
         exclusions
       ) do
    with {:ok, float_sol} <- solve_float_baseline_epochs(base_position, epochs, float_opts) do
      case residual_validation_outlier(float_sol, residual_opts) do
        nil ->
          finish_fixed_baseline_epochs(
            base_position,
            epochs,
            opts,
            float_opts,
            float_sol,
            residual_validation_meta(residual_opts, exclusions)
          )

        outlier ->
          if length(exclusions) < residual_opts.max_exclusions do
            epochs = drop_baseline_satellite_from_epochs(epochs, outlier.satellite_id)

            solve_fixed_baseline_epochs_attempt(
              base_position,
              epochs,
              opts,
              float_opts,
              residual_opts,
              [outlier | exclusions]
            )
          else
            {:error, {:residual_validation_failed, outlier, Enum.reverse(exclusions)}}
          end
      end
    end
  end

  defp finish_fixed_baseline_epochs(
         base_position,
         epochs,
         opts,
         float_opts,
         float_sol,
         residual_validation_meta
       ) do
    with {:ok, wavelengths} <-
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
         {:ok, float_only_systems} <- float_only_systems(opts),
         {:ok, fixed_cycles, fixed_meta} <-
           search_baseline_ambiguities(
             float_sol,
             wavelengths,
             offsets,
             integer_opts,
             float_only_systems
           ),
         fixed_m = fixed_ambiguities_m(fixed_cycles, wavelengths, offsets),
         free_ambiguity_ids = free_ambiguity_ids(float_sol.used_sats, fixed_cycles),
         {:ok, base} <- Types.normalize_ecef(base_position, :invalid_base_position),
         {:ok, normalized_epochs, _prep_meta} <- prepare_baseline_epochs(base, epochs, float_opts),
         {:ok, weights} <- baseline_weights(float_opts),
         {:ok, solve_opts} <- baseline_solve_options(float_opts) do
      state = %{
        baseline: {float_sol.baseline_m.x_m, float_sol.baseline_m.y_m, float_sol.baseline_m.z_m},
        ambiguities: Map.take(float_sol.ambiguities_m, free_ambiguity_ids)
      }

      case iterate_fixed_baseline(
             base,
             normalized_epochs,
             Map.fetch!(float_sol.metadata, :reference_satellites),
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
           ) do
        {:ok, sol} -> {:ok, attach_residual_validation_metadata(sol, residual_validation_meta)}
        {:error, _reason} = err -> err
      end
    end
  end

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

  This path is intentionally limited to one constellation at a time. If the
  normalized dual-frequency observations contain multiple constellation
  letters, it returns `{:error, {:unsupported_widelane, :multi_gnss}}` before
  wide-lane estimation.

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
    with :ok <- validate_options(opts, @widelane_baseline_options),
         {:ok, base} <- Types.normalize_ecef(base_position, :invalid_base_position),
         :ok <- ensure_nonempty_epochs(dual_epochs),
         {:ok, normalized_dual_epochs} <- normalize_dual_baseline_epochs(dual_epochs),
         :ok <- ensure_single_widelane_system(normalized_dual_epochs),
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
           |> Keyword.drop(@widelane_delegate_drop_options)
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
      difference: a satellite id binary (single-system data only) or a
      per-system map covering every observed system. When omitted, each
      system's lexicographically first common satellite is selected
      deterministically. Non-reference satellites difference against their own
      system's reference.

  Returns `{:ok, result}` or a tagged error. At least two common satellites are
  required so one non-reference double difference can be produced.
  """
  @spec double_differences([observation()], [observation()], keyword()) ::
          {:ok, result()} | {:error, term()}
  def double_differences(base_observations, rover_observations, opts \\ [])

  def double_differences(base_observations, rover_observations, opts)
      when is_list(base_observations) and is_list(rover_observations) do
    with :ok <- validate_options(opts, @double_difference_options),
         {:ok, base} <- normalize_observations(base_observations, :invalid_base_observations),
         {:ok, rover} <- normalize_observations(rover_observations, :invalid_rover_observations),
         {:ok, common, dropped} <- common_observations(base, rover),
         {:ok, refs} <- reference_satellites(common, opts) do
      ref_set = reference_satellite_set(refs)

      ref_data =
        Map.new(refs, fn {system, reference_sat} ->
          ref_base = Map.fetch!(base, reference_sat)
          ref_rover = Map.fetch!(rover, reference_sat)

          {system,
           %{
             satellite_id: reference_sat,
             ambiguity_id: single_difference_ambiguity_id(reference_sat, ref_base, ref_rover),
             code_sd: ref_rover.code_m - ref_base.code_m,
             phase_sd: ref_rover.phase_m - ref_base.phase_m
           }}
        end)

      dds =
        common
        |> Enum.reject(&MapSet.member?(ref_set, &1))
        |> Enum.map(fn sat ->
          ref = Map.fetch!(ref_data, satellite_system(sat))
          base_obs = Map.fetch!(base, sat)
          rover_obs = Map.fetch!(rover, sat)
          sat_sd_id = single_difference_ambiguity_id(sat, base_obs, rover_obs)

          %{
            satellite_id: sat,
            reference_satellite_id: ref.satellite_id,
            ambiguity_id: double_difference_ambiguity_id(sat, sat_sd_id, ref),
            code_m: rover_obs.code_m - base_obs.code_m - ref.code_sd,
            phase_m: rover_obs.phase_m - base_obs.phase_m - ref.phase_sd
          }
        end)

      {:ok,
       %{
         reference_satellite_id: reference_satellite_report(refs),
         double_differences: dds,
         dropped_sats: dropped
       }}
    end
  end

  def double_differences(_base_observations, _rover_observations, _opts),
    do: {:error, :invalid_observations}

  defp ensure_nonempty_epochs([]), do: {:error, :no_epochs}
  defp ensure_nonempty_epochs(_epochs), do: :ok

  defp validate_options(opts, allowed) when is_list(opts) do
    if Keyword.keyword?(opts) do
      allowed = MapSet.new(allowed)

      case Enum.find(Keyword.keys(opts), &(not MapSet.member?(allowed, &1))) do
        nil -> :ok
        key -> {:error, {:invalid_option, key}}
      end
    else
      {:error, {:invalid_option, :opts}}
    end
  end

  defp validate_options(_opts, _allowed), do: {:error, {:invalid_option, :opts}}

  defp residual_validation_options(opts) do
    threshold = Keyword.get(opts, :residual_threshold_sigma)
    max_exclusions = Keyword.get(opts, :max_residual_exclusions, @default_max_residual_exclusions)

    cond do
      not (is_integer(max_exclusions) and max_exclusions >= 0) ->
        {:error, {:invalid_option, :max_residual_exclusions}}

      is_nil(threshold) ->
        {:ok, %{enabled?: false, threshold_sigma: nil, max_exclusions: max_exclusions}}

      not (is_number(threshold) and threshold > 0.0) ->
        {:error, {:invalid_option, :residual_threshold_sigma}}

      true ->
        {:ok,
         %{
           enabled?: true,
           threshold_sigma: threshold / 1.0,
           max_exclusions: max_exclusions
         }}
    end
  end

  defp residual_validation_outlier(_float_sol, %{enabled?: false}), do: nil

  defp residual_validation_outlier(float_sol, %{threshold_sigma: threshold}) do
    float_sol.residuals_m
    |> Enum.flat_map(&residual_components/1)
    |> Enum.max_by(&abs(&1.normalized_residual), fn -> nil end)
    |> case do
      nil ->
        nil

      %{normalized_residual: normalized} = outlier when abs(normalized) > threshold ->
        Map.put(outlier, :threshold_sigma, threshold)

      _within_threshold ->
        nil
    end
  end

  defp residual_components(residual) do
    [
      %{
        epoch: residual.epoch,
        satellite_id: residual.satellite_id,
        reference_satellite_id: residual.reference_satellite_id,
        ambiguity_id: residual.ambiguity_id,
        kind: :code,
        residual_m: residual.code_m,
        sigma_m: residual.code_sigma_m,
        normalized_residual: residual.code_normalized
      },
      %{
        epoch: residual.epoch,
        satellite_id: residual.satellite_id,
        reference_satellite_id: residual.reference_satellite_id,
        ambiguity_id: residual.ambiguity_id,
        kind: :phase,
        residual_m: residual.phase_m,
        sigma_m: residual.phase_sigma_m,
        normalized_residual: residual.phase_normalized
      }
    ]
  end

  defp residual_validation_meta(%{enabled?: false}, _exclusions), do: nil

  defp residual_validation_meta(opts, exclusions) do
    exclusions = Enum.reverse(exclusions)

    %{
      threshold_sigma: opts.threshold_sigma,
      max_exclusions: opts.max_exclusions,
      excluded_sats: exclusions |> Enum.map(& &1.satellite_id) |> Enum.uniq() |> Enum.sort(),
      exclusions: exclusions
    }
  end

  defp attach_residual_validation_metadata(sol, nil), do: sol

  defp attach_residual_validation_metadata(%FixedBaselineSolution{} = sol, meta) do
    %{sol | metadata: Map.put(sol.metadata, :residual_validation, meta)}
  end

  defp drop_baseline_satellite_from_epochs(epochs, sat) do
    Enum.map(epochs, &drop_baseline_satellite_from_epoch(&1, sat))
  end

  defp drop_baseline_satellite_from_epoch(epoch, sat) when is_map(epoch) do
    epoch
    |> Map.update!(:base_observations, &drop_observations_satellite(&1, sat))
    |> Map.update!(:rover_observations, &drop_observations_satellite(&1, sat))
    |> Map.update!(:satellite_positions_m, &Map.delete(&1, sat))
    |> drop_optional_satellite_positions(:base_satellite_positions_m, sat)
    |> drop_optional_satellite_positions(:rover_satellite_positions_m, sat)
  end

  defp drop_optional_satellite_positions(epoch, key, sat) do
    if Map.has_key?(epoch, key), do: Map.update!(epoch, key, &Map.delete(&1, sat)), else: epoch
  end

  defp drop_observations_satellite(observations, sat) do
    Enum.reject(observations, &(observation_satellite_id(&1) == sat))
  end

  defp observation_satellite_id({sat, _code, _phase}), do: sat
  defp observation_satellite_id(%{satellite_id: sat}), do: sat
  defp observation_satellite_id(_obs), do: nil

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
      {:ok, normalized} -> {:ok, normalized |> Enum.reverse() |> attach_prediction_dts()}
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
    base_satellite_positions = Map.get(epoch, :base_satellite_positions_m, satellite_positions)
    rover_satellite_positions = Map.get(epoch, :rover_satellite_positions_m, satellite_positions)

    with {:ok, base} <- normalize_observations(base_observations, :invalid_base_observations),
         {:ok, rover} <- normalize_observations(rover_observations, :invalid_rover_observations),
         {:ok, positions} <- normalize_satellite_positions(satellite_positions),
         {:ok, base_positions} <- normalize_satellite_positions(base_satellite_positions),
         {:ok, rover_positions} <- normalize_satellite_positions(rover_satellite_positions),
         {:ok, velocity_mps} <- normalize_epoch_velocity(epoch, idx) do
      {:ok,
       %{
         idx: idx,
         epoch: Map.get(epoch, :epoch, idx),
         prediction_time: Map.get(epoch, :epoch),
         base: base,
         rover: rover,
         positions: positions,
         base_positions: base_positions,
         rover_positions: rover_positions,
         velocity_mps: velocity_mps
       }}
    end
  end

  defp normalize_epoch(_epoch, idx), do: {:error, {:invalid_epoch_observations, idx}}

  defp normalize_epoch_velocity(epoch, idx) do
    case Map.get(epoch, :velocity_mps) do
      nil ->
        {:ok, nil}

      velocity ->
        case Types.normalize_ecef(velocity, {:invalid_velocity_mps, idx}) do
          {:ok, {vx, vy, vz} = normalized} ->
            if finite_number?(vx) and finite_number?(vy) and finite_number?(vz),
              do: {:ok, normalized},
              else: {:error, {:invalid_velocity_mps, idx}}

          {:error, _reason} = err ->
            err
        end
    end
  end

  defp finite_number?(value) when is_number(value), do: value - value == 0.0
  defp finite_number?(_value), do: false

  defp attach_prediction_dts(epochs) do
    {epochs, _prev_time} =
      Enum.map_reduce(epochs, nil, fn epoch, prev_time ->
        time = Map.get(epoch, :prediction_time)
        dt_s = prediction_dt_seconds(prev_time, time)

        epoch =
          epoch
          |> Map.put(:prediction_dt_s, dt_s)
          |> Map.delete(:prediction_time)

        {epoch, time}
      end)

    epochs
  end

  defp prediction_dt_seconds(nil, _time), do: 0.0
  defp prediction_dt_seconds(_prev_time, nil), do: 0.0

  defp prediction_dt_seconds(prev, time) when is_number(prev) and is_number(time),
    do: (time - prev) / 1.0

  defp prediction_dt_seconds(%DateTime{} = prev, %DateTime{} = time),
    do: DateTime.diff(time, prev, :microsecond) / 1_000_000.0

  defp prediction_dt_seconds(%NaiveDateTime{} = prev, %NaiveDateTime{} = time),
    do: NaiveDateTime.diff(time, prev, :microsecond) / 1_000_000.0

  defp prediction_dt_seconds(%Date{} = prev, %Date{} = time), do: Date.diff(time, prev) * 86_400.0

  defp prediction_dt_seconds(_prev_time, _time), do: 0.0

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
    base_satellite_positions = Map.get(epoch, :base_satellite_positions_m, satellite_positions)
    rover_satellite_positions = Map.get(epoch, :rover_satellite_positions_m, satellite_positions)

    with {:ok, base} <-
           normalize_dual_observations(base_observations, :invalid_base_observations),
         {:ok, rover} <-
           normalize_dual_observations(rover_observations, :invalid_rover_observations),
         {:ok, positions} <- normalize_satellite_positions(satellite_positions),
         {:ok, base_positions} <- normalize_satellite_positions(base_satellite_positions),
         {:ok, rover_positions} <- normalize_satellite_positions(rover_satellite_positions) do
      {:ok,
       %{
         idx: idx,
         epoch: Map.get(epoch, :epoch, idx),
         base: base,
         rover: rover,
         positions: positions,
         base_positions: base_positions,
         rover_positions: rover_positions
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

      result =
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

      # As in the single-frequency prep (see `prepare_epochs_for_cycle_slips/2`),
      # a satellite that disappears and later reappears must start a FRESH
      # ambiguity arc: lock was lost during the outage so the integer can differ.
      # This must happen HERE, before wide-lane estimation, so the wide-lane
      # integers, narrow-lane wavelengths, and offsets are all keyed by the same
      # segmented `~raN` ambiguity ids that the delegated fixed-solve will see.
      # Dual epochs share the `%{base/rover => %{sat => %{ambiguity_id: ...}}}`
      # shape, so `segment_reacquired_arcs/1` applies unchanged. The delegated
      # `solve_fixed_baseline_epochs/3` re-segments the derived ionosphere-free
      # epochs, but `segment_reacquired_arcs/1` is idempotent so the ids match.
      case result do
        {:ok, prepared, meta} -> {:ok, segment_reacquired_arcs(prepared), meta}
        {:error, _reason} = err -> err
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
          positions: Map.drop(epoch.positions, MapSet.to_list(dropped)),
          base_positions: Map.drop(epoch.base_positions, MapSet.to_list(dropped)),
          rover_positions: Map.drop(epoch.rover_positions, MapSet.to_list(dropped))
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
    epoch_sats(epoch)
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
                satellite_positions_m: Map.take(epoch.positions, keep_sats),
                base_satellite_positions_m: Map.take(epoch.base_positions, keep_sats),
                rover_satellite_positions_m: Map.take(epoch.rover_positions, keep_sats)
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

      result =
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

      # A satellite that disappears (sets below the horizon, or loses lock with no
      # LLI flag) and later reappears must NOT reuse its previous carrier-phase
      # ambiguity: lock was lost during the outage, so the integer can differ.
      # This is ALWAYS a fresh ambiguity arc, independent of the LLI cycle-slip
      # policy above — reusing the pre-outage ambiguity corrupts the solution
      # (the sequential filter would hold a stale integer). Continuous arcs are
      # left unchanged.
      #
      # The wide-lane path also segments dual epochs up-front (see
      # `prepare_dual_baseline_cycle_slips/2`) and then delegates here, so this
      # may run a second time on already-segmented ids. That is safe:
      # `segment_reacquired_arcs/1` is idempotent (see `reacquired_ambiguity_id/2`).
      case result do
        {:ok, prepared, meta} -> {:ok, segment_reacquired_arcs(prepared), meta}
        {:error, _reason} = err -> err
      end
    end
  end

  defp empty_cycle_slip_meta, do: %{dropped_sats: [], split_arcs: []}

  defp segment_reacquired_arcs(epochs) do
    epochs
    |> segment_receiver_reacquisitions(:base)
    |> segment_receiver_reacquisitions(:rover)
  end

  # Walk one receiver's observation stream in order; each time a satellite
  # reappears after being absent, bump its arc counter and rewrite its
  # `ambiguity_id` so the post-outage arc is a distinct ambiguity. A satellite's
  # first arc keeps its original id, so continuous tracking is unaffected.
  defp segment_receiver_reacquisitions(epochs, receiver) do
    {epochs, _state} =
      Enum.map_reduce(epochs, %{present_last: MapSet.new(), arcs: %{}}, fn epoch, state ->
        observations = Map.fetch!(epoch, receiver)
        present = observations |> Map.keys() |> MapSet.new()

        {observations, arcs} =
          Enum.reduce(observations, {observations, state.arcs}, fn {sat, obs},
                                                                   {observations, arcs} ->
            reacquired? =
              Map.has_key?(arcs, sat) and not MapSet.member?(state.present_last, sat)

            arc = Map.get(arcs, sat, 0) + if(reacquired?, do: 1, else: 0)
            arcs = Map.put(arcs, sat, arc)

            observations =
              if arc > 0 do
                Map.put(observations, sat, %{
                  obs
                  | ambiguity_id: reacquired_ambiguity_id(obs.ambiguity_id, arc)
                })
              else
                observations
              end

            {observations, arcs}
          end)

        {Map.put(epoch, receiver, observations), %{present_last: present, arcs: arcs}}
      end)

    epochs
  end

  # Idempotent: strip any trailing `~raN` before re-appending, so segmenting the
  # same gap pattern twice (the wide-lane path segments dual epochs, then the
  # delegated fixed-solve segments the derived ionosphere-free epochs) yields
  # identical ids and the offset/wavelength maps stay keyed consistently.
  defp reacquired_ambiguity_id(ambiguity_id, arc) do
    base = String.replace(ambiguity_id, ~r/~ra\d+$/, "")
    "#{base}~ra#{arc}"
  end

  defp prepare_baseline_epochs(base, epochs, opts) do
    with :ok <- ensure_nonempty_epochs(epochs),
         {:ok, normalized_epochs} <- normalize_epochs(epochs),
         {:ok, normalized_epochs, slip_meta} <-
           prepare_epochs_for_cycle_slips(normalized_epochs, opts),
         {:ok, normalized_epochs, smoothing_meta} <-
           prepare_epochs_for_code_smoothing(normalized_epochs, opts),
         {:ok, normalized_epochs, mask_meta} <-
           apply_elevation_mask(base, normalized_epochs, opts) do
      {:ok, normalized_epochs, Map.merge(Map.merge(slip_meta, smoothing_meta), mask_meta)}
    end
  end

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

  defp apply_elevation_mask(base, epochs, opts) do
    with {:ok, mask_deg} <- elevation_mask_deg(opts) do
      case mask_deg do
        nil ->
          {:ok, epochs, %{elevation_mask_deg: nil, elevation_masked_sats: []}}

        mask_deg ->
          min_sin = :math.sin(mask_deg * :math.pi() / 180.0)

          {epochs, masked} =
            Enum.map_reduce(epochs, MapSet.new(), fn epoch, masked ->
              kept =
                epoch.positions
                |> Enum.filter(fn {_sat, sat_pos} -> elevation_sin(base, sat_pos) >= min_sin end)
                |> Enum.map(fn {sat, _sat_pos} -> sat end)

              kept_set = MapSet.new(kept)

              masked =
                epoch.positions
                |> Map.keys()
                |> MapSet.new()
                |> MapSet.difference(kept_set)
                |> MapSet.union(masked)

              {%{
                 epoch
                 | base: Map.take(epoch.base, kept),
                   rover: Map.take(epoch.rover, kept),
                   positions: Map.take(epoch.positions, kept),
                   base_positions: Map.take(epoch.base_positions, kept),
                   rover_positions: Map.take(epoch.rover_positions, kept)
               }, masked}
            end)

          {:ok, epochs,
           %{
             elevation_mask_deg: mask_deg,
             elevation_masked_sats: masked |> MapSet.to_list() |> Enum.sort()
           }}
      end
    end
  end

  defp elevation_mask_deg(opts) do
    case Keyword.get(opts, :elevation_mask_deg) do
      nil ->
        {:ok, nil}

      deg when is_number(deg) and deg >= 0.0 and deg < 90.0 ->
        {:ok, deg / 1.0}

      _other ->
        {:error, {:invalid_option, :elevation_mask_deg}}
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
          positions: Map.drop(epoch.positions, MapSet.to_list(dropped)),
          base_positions: Map.drop(epoch.base_positions, MapSet.to_list(dropped)),
          rover_positions: Map.drop(epoch.rover_positions, MapSet.to_list(dropped))
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
        epoch_sats(epoch)
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
        |> MapSet.union(epoch.base_positions |> Map.keys() |> MapSet.new())
        |> MapSet.union(epoch.rover_positions |> Map.keys() |> MapSet.new())
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
        epoch_sats(epoch)
        |> MapSet.union(acc)
      end)
      |> MapSet.to_list()
      |> Enum.sort()

    {:ok, all}
  end

  defp epoch_available_nonrefs(epoch, refs, physical_sats \\ nil) do
    available =
      epoch
      |> epoch_sats()
      |> MapSet.difference(reference_satellite_set(refs))

    case physical_sats do
      nil -> available
      sats -> MapSet.intersection(available, MapSet.new(sats))
    end
    |> MapSet.to_list()
    |> Enum.sort()
  end

  defp epoch_available_sats(epoch, physical_sats) do
    epoch
    |> epoch_sats()
    |> MapSet.intersection(MapSet.new(physical_sats))
    |> MapSet.to_list()
    |> Enum.sort()
  end

  defp epoch_sats(epoch) do
    epoch.base
    |> Map.keys()
    |> MapSet.new()
    |> MapSet.intersection(epoch.rover |> Map.keys() |> MapSet.new())
    |> MapSet.intersection(epoch.positions |> Map.keys() |> MapSet.new())
    |> MapSet.intersection(epoch.base_positions |> Map.keys() |> MapSet.new())
    |> MapSet.intersection(epoch.rover_positions |> Map.keys() |> MapSet.new())
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

  defp sequential_filter_options(opts) do
    with {:ok, baseline_prior_sigma_m} <-
           positive_option(
             opts,
             :baseline_prior_sigma_m,
             @default_filter_baseline_prior_sigma_m
           ),
         {:ok, ambiguity_prior_sigma_m} <-
           positive_option(
             opts,
             :ambiguity_prior_sigma_m,
             @default_filter_ambiguity_prior_sigma_m
           ),
         {:ok, hold_sigma_m} <-
           positive_option(opts, :hold_sigma_m, @default_filter_hold_sigma_m),
         {:ok, process_noise_baseline_sigma_m} <-
           nonnegative_option(
             opts,
             :process_noise_baseline_sigma_m,
             @default_filter_process_noise_baseline_sigma_m
           ),
         {:ok, dynamics_model} <- dynamics_model(opts),
         {:ok, innovation_screen} <- innovation_screen_options(opts) do
      {:ok,
       %{
         baseline_prior_sigma_m: baseline_prior_sigma_m,
         ambiguity_prior_sigma_m: ambiguity_prior_sigma_m,
         hold_sigma_m: hold_sigma_m,
         process_noise_baseline_sigma_m: process_noise_baseline_sigma_m,
         dynamics_model: dynamics_model,
         innovation_screen: innovation_screen
       }}
    end
  end

  defp dynamics_model(opts) do
    case Keyword.get(opts, :dynamics_model, @default_filter_dynamics_model) do
      value when value in [:constant_position, :velocity_propagated] -> {:ok, value}
      _other -> {:error, {:invalid_option, :dynamics_model}}
    end
  end

  defp innovation_screen_options(opts) do
    threshold = Keyword.get(opts, :innovation_screen_sigma)

    min_rows =
      Keyword.get(opts, :innovation_screen_min_rows, @default_filter_innovation_screen_min_rows)

    cond do
      not is_nil(threshold) and (not is_number(threshold) or threshold <= 0.0) ->
        {:error, {:invalid_option, :innovation_screen_sigma}}

      not is_integer(min_rows) or min_rows < 1 ->
        {:error, {:invalid_option, :innovation_screen_min_rows}}

      is_nil(threshold) ->
        {:ok, %{enabled?: false, threshold_sigma: nil, min_rows: min_rows}}

      true ->
        {:ok, %{enabled?: true, threshold_sigma: threshold / 1.0, min_rows: min_rows}}
    end
  end

  defp filter_kernel(opts) do
    case Keyword.get(opts, :filter_kernel, :rust) do
      value when value in [:elixir, :rust] -> {:ok, value}
      _other -> {:error, {:invalid_option, :filter_kernel}}
    end
  end

  defp positive_option(opts, key, default) do
    value = Keyword.get(opts, key, default)

    if is_number(value) and value > 0.0,
      do: {:ok, value / 1.0},
      else: {:error, {:invalid_option, key}}
  end

  defp nonnegative_option(opts, key, default) do
    value = Keyword.get(opts, key, default)

    if is_number(value) and value >= 0.0,
      do: {:ok, value / 1.0},
      else: {:error, {:invalid_option, key}}
  end

  defp baseline_weights(opts) do
    with {:ok, code_sigma_m} <- measurement_sigma(opts, :code_sigma_m, @default_code_sigma_m),
         {:ok, phase_sigma_m} <- measurement_sigma(opts, :phase_sigma_m, @default_phase_sigma_m),
         {:ok, stochastic_model} <- stochastic_model(opts),
         {:ok, elevation_weighting?} <- elevation_weighting(opts),
         {:ok, sagnac?} <- sagnac(opts) do
      {:ok,
       %{
         code_sigma_m: code_sigma_m,
         phase_sigma_m: phase_sigma_m,
         stochastic_model: stochastic_model,
         elevation_weighting?: elevation_weighting?,
         sagnac?: sagnac?
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

  defp sagnac(opts) do
    case Keyword.get(opts, :sagnac, @default_sagnac) do
      value when is_boolean(value) -> {:ok, value}
      _other -> {:error, {:invalid_option, :sagnac}}
    end
  end

  defp stochastic_model(opts) do
    case Keyword.get(opts, :stochastic_model, @default_stochastic_model) do
      value when value in [:simple, :rtklib] -> {:ok, value}
      _other -> {:error, {:invalid_option, :stochastic_model}}
    end
  end

  defp initial_baseline(opts) do
    opts
    |> Keyword.get(:initial_baseline_m, {0.0, 0.0, 0.0})
    |> Types.normalize_ecef(:invalid_initial_baseline)
  end

  defp baseline_ambiguity_index(epochs, common_sats, refs) do
    physical_sats = nonreference_sats(common_sats, refs)

    epochs
    |> Enum.reduce_while({:ok, %{}}, fn epoch, {:ok, acc} ->
      ref_sds = epoch_reference_sds(epoch, refs)

      epoch
      |> epoch_available_nonrefs(refs, physical_sats)
      |> Enum.reduce_while({:ok, acc}, fn sat, {:ok, acc} ->
        ref_dd = Map.fetch!(ref_sds, satellite_system(sat))
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

  defp single_difference_ambiguity_index(epochs, all_sats) do
    epochs
    |> Enum.reduce_while({:ok, %{}}, fn epoch, {:ok, acc} ->
      epoch
      |> epoch_available_sats(all_sats)
      |> Enum.reduce_while({:ok, acc}, fn sat, {:ok, acc} ->
        ambiguity_id = single_difference(epoch, sat).ambiguity_id

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

  defp sequential_dd_ambiguity_index(epochs, all_sats, refs) do
    physical_sats = nonreference_sats(all_sats, refs)

    epochs
    |> Enum.reduce_while({:ok, {%{}, %{}}}, fn epoch, {:ok, {satellites, pairs}} ->
      ref_sds = epoch_reference_sds(epoch, refs)

      epoch
      |> epoch_available_nonrefs(refs, physical_sats)
      |> Enum.reduce_while({:ok, {satellites, pairs}}, fn sat, {:ok, {satellites, pairs}} ->
        ref_dd = Map.fetch!(ref_sds, satellite_system(sat))
        sd = single_difference(epoch, sat)
        dd_id = double_difference_ambiguity_id(sat, sd.ambiguity_id, ref_dd)
        pair = %{sat_sd_id: sd.ambiguity_id, ref_sd_id: ref_dd.ambiguity_id}

        case Map.fetch(satellites, dd_id) do
          {:ok, ^sat} ->
            {:cont, {:ok, {satellites, Map.put_new(pairs, dd_id, pair)}}}

          {:ok, other_sat} ->
            {:halt, {:error, {:duplicate_ambiguity_id, dd_id, other_sat, sat}}}

          :error ->
            {:cont, {:ok, {Map.put(satellites, dd_id, sat), Map.put(pairs, dd_id, pair)}}}
        end
      end)
      |> case do
        {:ok, acc} -> {:cont, {:ok, acc}}
        {:error, _reason} = err -> {:halt, err}
      end
    end)
    |> case do
      {:ok, {ambiguity_satellites, pairs}} ->
        ambiguity_ids =
          ambiguity_satellites
          |> Enum.sort_by(fn {ambiguity_id, sat} -> {sat, ambiguity_id} end)
          |> Enum.map(&elem(&1, 0))

        {:ok, ambiguity_ids, Map.new(ambiguity_satellites), pairs}

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

  defp baseline_row_count(epochs, refs) do
    epochs
    |> Enum.map(fn epoch -> length(epoch_available_nonrefs(epoch, refs)) end)
    |> Enum.sum()
    |> Kernel.*(2)
  end

  defp baseline_unknown_count(ambiguity_ids), do: 3 + length(ambiguity_ids)

  defp iterate_baseline(
         base,
         epochs,
         refs,
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
             refs,
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
            refs,
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
            refs,
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
            refs,
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

  defp build_baseline_rows(base, epochs, refs, physical_sats, ambiguity_ids, state, weights) do
    epochs
    |> Enum.reduce_while({:ok, []}, fn epoch, {:ok, acc} ->
      case build_epoch_baseline_rows(
             base,
             epoch,
             refs,
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

  defp build_epoch_baseline_rows(base, epoch, refs, physical_sats, ambiguity_ids, state, weights) do
    ref_data = epoch_reference_data(epoch, refs)

    epoch
    |> epoch_available_nonrefs(refs, physical_sats)
    |> Enum.reduce({:ok, []}, fn sat, {:ok, acc} ->
      %{sd: ref_dd, pos: ref_pos, base_pos: ref_base_pos, rover_pos: ref_rover_pos} =
        Map.fetch!(ref_data, satellite_system(sat))

      obs_dd = double_difference_measurement(epoch, sat, ref_dd)
      sat_pos = Map.fetch!(epoch.positions, sat)
      sat_base_pos = Map.fetch!(epoch.base_positions, sat)
      sat_rover_pos = Map.fetch!(epoch.rover_positions, sat)

      {geom_dd, deriv} =
        geometry_double_difference(
          base,
          state.baseline,
          sat_base_pos,
          sat_rover_pos,
          ref_base_pos,
          ref_rover_pos,
          weights
        )

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
        ref_sat: ref_dd.satellite_id,
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
        ref_sat: ref_dd.satellite_id,
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

  defp build_epoch_sequential_baseline_rows(
         base,
         epoch,
         refs,
         physical_sats,
         sd_ambiguity_ids,
         dd_ambiguity_pairs,
         state,
         weights
       ) do
    ref_data = epoch_reference_data(epoch, refs)

    epoch
    |> epoch_available_nonrefs(refs, physical_sats)
    |> Enum.reduce({:ok, []}, fn sat, {:ok, acc} ->
      %{sd: ref_dd, pos: ref_pos, base_pos: ref_base_pos, rover_pos: ref_rover_pos} =
        Map.fetch!(ref_data, satellite_system(sat))

      obs_dd = double_difference_measurement(epoch, sat, ref_dd)
      sat_pos = Map.fetch!(epoch.positions, sat)
      sat_base_pos = Map.fetch!(epoch.base_positions, sat)
      sat_rover_pos = Map.fetch!(epoch.rover_positions, sat)

      {geom_dd, deriv} =
        geometry_double_difference(
          base,
          state.baseline,
          sat_base_pos,
          sat_rover_pos,
          ref_base_pos,
          ref_rover_pos,
          weights
        )

      %{sat_sd_id: sat_sd_id, ref_sd_id: ref_sd_id} =
        Map.fetch!(dd_ambiguity_pairs, obs_dd.ambiguity_id)

      ambiguity = dd_state_ambiguity_m(state, sat_sd_id, ref_sd_id)
      h_base = design_baseline_row(deriv, nil, sd_ambiguity_ids)

      h_phase =
        design_single_difference_baseline_row(deriv, sat_sd_id, ref_sd_id, sd_ambiguity_ids)

      code_variance = row_variance(weights, :code, base, sat_pos, ref_pos)
      phase_variance = row_variance(weights, :phase, base, sat_pos, ref_pos)

      code_row = %{
        kind: :code,
        epoch_idx: epoch.idx,
        epoch: epoch.epoch,
        sat: sat,
        ref_sat: ref_dd.satellite_id,
        ambiguity_id: obs_dd.ambiguity_id,
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
        ref_sat: ref_dd.satellite_id,
        ambiguity_id: obs_dd.ambiguity_id,
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

  # Per-epoch reference-satellite data, one entry per system whose reference is
  # observed in this epoch. The per-system common invariant guarantees the
  # reference is present in every epoch in which its system appears.
  defp epoch_reference_data(epoch, refs) do
    available = epoch_sats(epoch)

    refs
    |> Enum.filter(fn {_system, ref} -> MapSet.member?(available, ref) end)
    |> Map.new(fn {system, ref} ->
      {system,
       %{
         sd: single_difference(epoch, ref),
         pos: Map.fetch!(epoch.positions, ref),
         base_pos: Map.fetch!(epoch.base_positions, ref),
         rover_pos: Map.fetch!(epoch.rover_positions, ref)
       }}
    end)
  end

  defp epoch_reference_sds(epoch, refs) do
    available = epoch_sats(epoch)

    refs
    |> Enum.filter(fn {_system, ref} -> MapSet.member?(available, ref) end)
    |> Map.new(fn {system, ref} -> {system, single_difference(epoch, ref)} end)
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

  defp dd_state_ambiguity_m(state, sat_sd_id, ref_sd_id) do
    Map.fetch!(state.ambiguities, sat_sd_id) - Map.fetch!(state.ambiguities, ref_sd_id)
  end

  defp geometry_double_difference(
         base,
         baseline,
         sat_base_pos,
         sat_rover_pos,
         ref_base_pos,
         ref_rover_pos,
         weights
       ) do
    rover = add3(base, baseline)

    sat_sd =
      geometric_range(sat_rover_pos, rover, weights) -
        geometric_range(sat_base_pos, base, weights)

    ref_sd =
      geometric_range(ref_rover_pos, rover, weights) -
        geometric_range(ref_base_pos, base, weights)

    sat_deriv = geometric_range_derivative(rover, sat_rover_pos, weights)
    ref_deriv = geometric_range_derivative(rover, ref_rover_pos, weights)
    {sat_sd - ref_sd, sub3(sat_deriv, ref_deriv)}
  end

  defp row_variance(weights, kind, base, sat_pos, ref_pos) do
    sigma =
      case kind do
        :code -> weights.code_sigma_m
        :phase -> weights.phase_sigma_m
      end

    sat_sd_variance = single_difference_variance(sigma, weights, base, sat_pos)
    ref_sd_variance = single_difference_variance(sigma, weights, base, ref_pos)

    dd_variance = sat_sd_variance + ref_sd_variance

    %{
      sd_variance_m2: sat_sd_variance,
      ref_sd_variance_m2: ref_sd_variance,
      weight: 1.0 / :math.sqrt(dd_variance)
    }
  end

  defp single_difference_variance(
         sigma_m,
         %{stochastic_model: :simple, elevation_weighting?: false},
         _base,
         _sat_pos
       ) do
    2.0 * sigma_m * sigma_m
  end

  defp single_difference_variance(
         sigma_m,
         %{stochastic_model: :simple, elevation_weighting?: true},
         base,
         sat_pos
       ) do
    sin_el = elevation_sin(base, sat_pos) |> max(@min_elevation_sin)
    scaled_sigma = sigma_m / sin_el
    2.0 * scaled_sigma * scaled_sigma
  end

  defp single_difference_variance(sigma_m, %{stochastic_model: :rtklib}, base, sat_pos) do
    sin_el = elevation_sin(base, sat_pos) |> max(@min_elevation_sin)
    2.0 * (sigma_m * sigma_m + sigma_m * sigma_m / (sin_el * sin_el))
  end

  defp design_baseline_row({dx, dy, dz}, ambiguity_id, ambiguity_ids) do
    ambiguity_cols = Enum.map(ambiguity_ids, &if(&1 == ambiguity_id, do: 1.0, else: 0.0))
    [dx, dy, dz | ambiguity_cols]
  end

  defp design_single_difference_baseline_row({dx, dy, dz}, sat_sd_id, ref_sd_id, sd_ambiguity_ids) do
    ambiguity_cols =
      Enum.map(sd_ambiguity_ids, fn id ->
        cond do
          id == sat_sd_id -> 1.0
          id == ref_sd_id -> -1.0
          true -> 0.0
        end
      end)

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
         refs,
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
             refs,
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
         {:ok, covariance_inverse_m} <- invert_matrix(covariance_m),
         {:ok, residuals} <- baseline_residuals(rows) do
      code_residuals = Enum.map(residuals, & &1.code_m)
      phase_residuals = Enum.map(residuals, & &1.phase_m)
      rover = add3(base, state.baseline)

      {:ok,
       %FloatBaselineSolution{
         baseline_m: ecef_map(state.baseline),
         rover_position_m: ecef_map(rover),
         reference_satellite_id: reference_satellite_report(refs),
         used_sats: ambiguity_ids,
         ambiguities_m: state.ambiguities,
         residuals_m: residuals,
         metadata: %{
           iterations: iterations,
           converged: converged?,
           status: status,
           physical_sats: physical_sats,
           reference_satellites: refs,
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
             stochastic_model: weights.stochastic_model,
             elevation_weighting: weights.elevation_weighting?,
             sagnac: weights.sagnac?,
             min_elevation_sin: @min_elevation_sin
           },
           code_rms_m: rms(code_residuals),
           phase_rms_m: rms(phase_residuals),
           weighted_rms_m: weighted_rms(rows),
           n_epochs: length(epochs),
           n_observations: length(rows),
           dropped_sats:
             Enum.uniq(
               dropped_sats ++
                 slip_meta.dropped_sats ++
                 Map.get(slip_meta, :elevation_masked_sats, [])
             )
             |> Enum.sort(),
           dropped_cycle_slip_sats: slip_meta.dropped_sats,
           elevation_mask_deg: Map.get(slip_meta, :elevation_mask_deg),
           elevation_masked_sats: Map.get(slip_meta, :elevation_masked_sats, []),
           split_cycle_slip_arcs: slip_meta.split_arcs,
           code_smoothing: Map.get(slip_meta, :code_smoothing, false),
           code_smoothing_window_cap: Map.get(slip_meta, :code_smoothing_window_cap)
         }
       }}
    end
  end

  defp search_baseline_ambiguities(
         float_sol,
         wavelengths,
         offsets,
         integer_opts,
         float_only_systems
       ) do
    ambiguity_satellites = Map.fetch!(float_sol.metadata, :ambiguity_satellites)
    float_only_ids = float_only_ambiguity_ids(ambiguity_satellites, float_only_systems)
    search_sats = Enum.reject(float_sol.used_sats, &MapSet.member?(float_only_ids, &1))

    if search_sats == [] do
      {:ok, %{},
       Map.merge(empty_integer_search_meta(offsets), partial_meta(false, false, %{}, []))}
    else
      covariance_cycles =
        covariance_m_to_cycles(
          float_sol.metadata.ambiguity_float.covariance_m,
          float_sol.used_sats,
          wavelengths
        )

      # Float-only systems' ambiguities never enter the integer search set
      # (full-set or partial). With none configured this is the identity
      # selection, leaving the historical arithmetic untouched.
      search_covariance =
        if search_sats == float_sol.used_sats,
          do: covariance_cycles,
          else: covariance_submatrix(float_sol.used_sats, covariance_cycles, search_sats)

      float_cycles =
        Map.new(search_sats, fn sat ->
          ambiguity_m = Map.fetch!(float_sol.ambiguities_m, sat)
          offset_m = Map.fetch!(offsets, sat)
          {sat, (ambiguity_m - offset_m) / Map.fetch!(wavelengths, sat)}
        end)

      with {:ok, fixed_cycles, meta} <-
             IntegerLeastSquares.search(float_cycles, search_covariance, integer_opts) do
        meta = Map.put(meta, :ambiguity_offsets_m, offsets)

        if meta.integer_status == :fixed or not integer_opts.partial_ambiguity_resolution? do
          {:ok, fixed_cycles, Map.merge(meta, partial_meta(false, false, fixed_cycles, []))}
        else
          search_partial_baseline_ambiguities(
            search_sats,
            float_cycles,
            search_covariance,
            offsets,
            integer_opts,
            fixed_cycles,
            meta
          )
        end
      end
    end
  end

  defp empty_integer_search_meta(offsets) do
    %{
      integer_status: :not_fixed,
      integer_method: :lambda,
      integer_ratio: nil,
      integer_best_score: nil,
      integer_second_best_score: nil,
      integer_candidates: 0,
      ambiguity_search: %{
        order: [],
        float_cycles: %{},
        covariance_cycles: [],
        covariance_inverse_cycles: []
      },
      ambiguity_offsets_m: offsets
    }
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

  defp dd_covariance_m_to_cycles(sd_covariance_m, sd_ids, dd_ids, dd_pairs, wavelengths) do
    sd_index = sd_ids |> Enum.with_index() |> Map.new()

    dd_ids
    |> Enum.map(fn dd_i ->
      lambda_i = Map.fetch!(wavelengths, dd_i)
      ti = dd_transform_indices(Map.fetch!(dd_pairs, dd_i), sd_index)

      dd_ids
      |> Enum.map(fn dd_j ->
        lambda_j = Map.fetch!(wavelengths, dd_j)
        tj = dd_transform_indices(Map.fetch!(dd_pairs, dd_j), sd_index)

        dd_covariance_m(sd_covariance_m, ti, tj) / (lambda_i * lambda_j)
      end)
    end)
  end

  defp dd_transform_indices(%{sat_sd_id: sat_sd_id, ref_sd_id: ref_sd_id}, sd_index) do
    [{Map.fetch!(sd_index, sat_sd_id), 1.0}, {Map.fetch!(sd_index, ref_sd_id), -1.0}]
  end

  defp dd_covariance_m(covariance_m, ti, tj) do
    Enum.reduce(ti, 0.0, fn {i, wi}, acc ->
      acc +
        Enum.reduce(tj, 0.0, fn {j, wj}, inner ->
          inner + wi * wj * (covariance_m |> Enum.at(i) |> Enum.at(j))
        end)
    end)
  end

  defp fixed_ambiguities_m(fixed_cycles, wavelengths, offsets) do
    Map.new(fixed_cycles, fn {sat, cycles} ->
      {sat, Map.fetch!(offsets, sat) + cycles * Map.fetch!(wavelengths, sat)}
    end)
  end

  defp run_sequential_baseline_filter(
         base,
         epochs,
         refs,
         physical_sats,
         sd_ambiguity_ids,
         sd_ambiguity_satellites,
         dd_ambiguity_ids,
         dd_ambiguity_satellites,
         dd_ambiguity_pairs,
         float_only_systems,
         wavelengths,
         offsets,
         weights,
         solve_opts,
         integer_opts,
         filter_opts,
         initial_baseline,
         prep_meta,
         filter_kernel
       ) do
    n = baseline_unknown_count(sd_ambiguity_ids)

    float_only_dd_ids =
      float_only_ambiguity_ids(dd_ambiguity_satellites, float_only_systems)

    {initial_ambiguities, initial_ambiguity_count} =
      initial_sequential_ambiguities(
        epochs,
        Map.values(refs) ++ physical_sats,
        sd_ambiguity_ids
      )

    initial = %{
      state: %{
        baseline: initial_baseline,
        ambiguities: initial_ambiguities
      },
      information:
        sequential_initial_information(
          n,
          filter_opts.baseline_prior_sigma_m,
          filter_opts.ambiguity_prior_sigma_m
        ),
      fixed_cycles: %{},
      fixed_m: %{},
      rust_state: nil,
      epochs: []
    }

    case filter_kernel do
      :elixir ->
        epochs
        |> Enum.reduce_while({:ok, initial}, fn epoch, {:ok, acc} ->
          case sequential_filter_epoch(
                 base,
                 epoch,
                 refs,
                 physical_sats,
                 sd_ambiguity_ids,
                 dd_ambiguity_ids,
                 dd_ambiguity_pairs,
                 float_only_dd_ids,
                 wavelengths,
                 offsets,
                 weights,
                 solve_opts,
                 integer_opts,
                 filter_opts,
                 acc
               ) do
            {:ok, next} -> {:cont, {:ok, next}}
            {:error, _reason} = err -> {:halt, err}
          end
        end)

      :rust ->
        initial = %{
          initial
          | rust_state:
              rust_initial_filter_state(
                epochs,
                refs,
                initial_baseline,
                sd_ambiguity_ids,
                initial_ambiguities,
                filter_opts.baseline_prior_sigma_m,
                filter_opts.ambiguity_prior_sigma_m
              )
        }

        sequential_filter_epochs_rust(
          base,
          epochs,
          refs,
          physical_sats,
          dd_ambiguity_pairs,
          dd_ambiguity_satellites,
          float_only_dd_ids,
          float_only_systems,
          wavelengths,
          offsets,
          weights,
          solve_opts,
          integer_opts,
          filter_opts,
          initial
        )
    end
    |> case do
      {:ok, acc} ->
        baseline = acc.state.baseline
        rover = add3(base, baseline)
        epoch_results = Enum.reverse(acc.epochs)
        first_fixed = Enum.find(epoch_results, &(&1.integer_status == :fixed))
        fixed_epochs = Enum.count(epoch_results, &(&1.integer_status == :fixed))

        {:ok,
         %FilterBaselineSolution{
           baseline_m: ecef_map(baseline),
           rover_position_m: ecef_map(rover),
           reference_satellite_id: reference_satellite_report(refs),
           fixed_ambiguities_cycles: acc.fixed_cycles,
           epochs: epoch_results,
           metadata: %{
             integer_method: :sequential_lambda,
             ambiguity_state: :single_difference,
             first_fixed_epoch: first_fixed && first_fixed.epoch,
             first_fixed_index: first_fixed && first_fixed.index,
             fixed_epoch_count: fixed_epochs,
             n_epochs: length(epochs),
             physical_sats: physical_sats,
             reference_satellites: refs,
             float_only_systems: float_only_systems,
             ambiguity_satellites: dd_ambiguity_satellites,
             single_difference_ambiguity_satellites: sd_ambiguity_satellites,
             single_difference_ambiguity_count: length(sd_ambiguity_ids),
             measurement_covariance: %{
               model: :double_difference,
               code_sigma_m: weights.code_sigma_m,
               phase_sigma_m: weights.phase_sigma_m,
               stochastic_model: weights.stochastic_model,
               elevation_weighting: weights.elevation_weighting?,
               sagnac: weights.sagnac?,
               min_elevation_sin: @min_elevation_sin
             },
             dropped_sats:
               Enum.uniq(prep_meta.dropped_sats ++ Map.get(prep_meta, :elevation_masked_sats, []))
               |> Enum.sort(),
             dropped_cycle_slip_sats: prep_meta.dropped_sats,
             elevation_mask_deg: Map.get(prep_meta, :elevation_mask_deg),
             elevation_masked_sats: Map.get(prep_meta, :elevation_masked_sats, []),
             split_cycle_slip_arcs: prep_meta.split_arcs,
             hold_sigma_m: filter_opts.hold_sigma_m,
             baseline_prior_sigma_m: filter_opts.baseline_prior_sigma_m,
             ambiguity_prior_sigma_m: filter_opts.ambiguity_prior_sigma_m,
             dynamics_model: filter_opts.dynamics_model,
             ambiguity_initialization: :phase_code,
             initialized_ambiguity_count: initial_ambiguity_count,
             filter_kernel: filter_kernel
           }
         }}

      {:error, _reason} = err ->
        err
    end
  end

  defp initial_sequential_ambiguities(epochs, physical_sats, sd_ambiguity_ids) do
    zero = Map.new(sd_ambiguity_ids, &{&1, 0.0})
    ambiguity_set = MapSet.new(sd_ambiguity_ids)

    seeded =
      Enum.reduce_while(epochs, %{}, fn epoch, acc ->
        if map_size(acc) == length(sd_ambiguity_ids) do
          {:halt, acc}
        else
          acc =
            epoch
            |> epoch_available_sats(physical_sats)
            |> Enum.reduce(acc, fn sat, sat_acc ->
              sd = single_difference(epoch, sat)

              if MapSet.member?(ambiguity_set, sd.ambiguity_id) and
                   not Map.has_key?(sat_acc, sd.ambiguity_id) do
                # RTK filters conventionally seed carrier ambiguities from the
                # phase-code difference; the prior remains intentionally broad,
                # but starting near the code-pinned level avoids spending early
                # epochs pulling kilometre-scale zero-state ambiguities into place.
                Map.put(sat_acc, sd.ambiguity_id, sd.phase_m - sd.code_m)
              else
                sat_acc
              end
            end)

          {:cont, acc}
        end
      end)

    {Map.merge(zero, seeded), map_size(seeded)}
  end

  defp sequential_filter_epoch(
         base,
         epoch,
         refs,
         physical_sats,
         sd_ambiguity_ids,
         dd_ambiguity_ids,
         dd_ambiguity_pairs,
         float_only_dd_ids,
         wavelengths,
         offsets,
         weights,
         solve_opts,
         integer_opts,
         filter_opts,
         acc
       ) do
    {predicted_acc, fallback_acc} = time_update_information(acc, epoch, filter_opts)

    case do_sequential_filter_epoch(
           base,
           epoch,
           refs,
           physical_sats,
           sd_ambiguity_ids,
           dd_ambiguity_ids,
           dd_ambiguity_pairs,
           float_only_dd_ids,
           wavelengths,
           offsets,
           weights,
           solve_opts,
           integer_opts,
           filter_opts,
           predicted_acc
         ) do
      {:error, :singular_geometry} = err ->
        if fallback_acc do
          do_sequential_filter_epoch(
            base,
            epoch,
            refs,
            physical_sats,
            sd_ambiguity_ids,
            dd_ambiguity_ids,
            dd_ambiguity_pairs,
            float_only_dd_ids,
            wavelengths,
            offsets,
            weights,
            solve_opts,
            integer_opts,
            filter_opts,
            fallback_acc
          )
        else
          err
        end

      result ->
        result
    end
  end

  defp do_sequential_filter_epoch(
         base,
         epoch,
         refs,
         physical_sats,
         sd_ambiguity_ids,
         dd_ambiguity_ids,
         dd_ambiguity_pairs,
         float_only_dd_ids,
         wavelengths,
         offsets,
         weights,
         solve_opts,
         integer_opts,
         filter_opts,
         acc
       ) do
    with {:ok, screen_meta, row_mask} <-
           sequential_innovation_screen(
             base,
             epoch,
             refs,
             physical_sats,
             sd_ambiguity_ids,
             dd_ambiguity_pairs,
             acc.state,
             weights,
             filter_opts.innovation_screen
           ) do
      if screen_meta && screen_meta.coasted? do
        sequential_coasted_epoch(
          base,
          epoch,
          refs,
          physical_sats,
          sd_ambiguity_ids,
          dd_ambiguity_pairs,
          float_only_dd_ids,
          weights,
          screen_meta,
          acc
        )
      else
        do_sequential_filter_epoch_update(
          base,
          epoch,
          refs,
          physical_sats,
          sd_ambiguity_ids,
          dd_ambiguity_ids,
          dd_ambiguity_pairs,
          float_only_dd_ids,
          wavelengths,
          offsets,
          weights,
          solve_opts,
          integer_opts,
          filter_opts,
          screen_meta,
          row_mask,
          acc
        )
      end
    end
  end

  # Coasted epoch (kernel `coasted_update`): the screen rejected too many rows,
  # so the measurement update is skipped entirely. The carried state keeps the
  # PREDICTED prior (time update already applied by the caller); the reported
  # baseline is that prior; nothing is searched or newly fixed. Residuals are
  # reported from the full (unmasked) rows at the prior, mirroring the rust
  # path's post-hoc residual construction.
  defp sequential_coasted_epoch(
         base,
         epoch,
         refs,
         physical_sats,
         sd_ambiguity_ids,
         dd_ambiguity_pairs,
         float_only_dd_ids,
         weights,
         screen_meta,
         acc
       ) do
    with {:ok, residual_rows} <-
           build_epoch_sequential_baseline_rows(
             base,
             epoch,
             refs,
             physical_sats,
             sd_ambiguity_ids,
             dd_ambiguity_pairs,
             acc.state,
             weights
           ),
         {:ok, residuals} <- baseline_residuals(residual_rows) do
      all_fixed = acc.fixed_cycles |> Map.keys() |> Enum.sort()

      epoch_result = %{
        epoch: epoch.epoch,
        index: epoch.idx,
        baseline_m: ecef_map(acc.state.baseline),
        integer_status: :coasted,
        integer_ratio: rust_public_integer_ratio(0.0, residual_rows, acc, float_only_dd_ids),
        integer_best_score: nil,
        integer_second_best_score: nil,
        integer_candidates: nil,
        ambiguity_search: nil,
        residuals_m: residuals,
        newly_fixed_ambiguities: [],
        fixed_ambiguities: all_fixed,
        innovation_screen: screen_meta
      }

      {:ok, %{acc | epochs: [epoch_result | acc.epochs]}}
    end
  end

  defp do_sequential_filter_epoch_update(
         base,
         epoch,
         refs,
         physical_sats,
         sd_ambiguity_ids,
         dd_ambiguity_ids,
         dd_ambiguity_pairs,
         float_only_dd_ids,
         wavelengths,
         offsets,
         weights,
         solve_opts,
         integer_opts,
         filter_opts,
         screen_meta,
         row_mask,
         acc
       ) do
    with {:ok, state, information, rows} <-
           iterate_sequential_filter_epoch(
             base,
             epoch,
             refs,
             physical_sats,
             sd_ambiguity_ids,
             dd_ambiguity_pairs,
             acc.state,
             acc.state,
             acc.information,
             acc.fixed_m,
             weights,
             solve_opts,
             filter_opts,
             row_mask,
             1
           ),
         {:ok, covariance} <- invert_matrix(information),
         {:ok, fixed_cycles, fixed_m, search_meta} <-
           sequential_search_and_hold(
             state,
             covariance,
             rows,
             sd_ambiguity_ids,
             dd_ambiguity_ids,
             dd_ambiguity_pairs,
             float_only_dd_ids,
             acc.fixed_cycles,
             wavelengths,
             offsets,
             integer_opts
           ),
         newly_fixed =
           fixed_cycles |> Map.keys() |> Kernel.--(Map.keys(acc.fixed_cycles)) |> Enum.sort(),
         {:ok, report_state} <-
           conditioned_fixed_state(
             state,
             newly_fixed,
             base,
             epoch,
             refs,
             physical_sats,
             sd_ambiguity_ids,
             dd_ambiguity_pairs,
             acc,
             fixed_m,
             weights,
             solve_opts,
             filter_opts,
             row_mask
           ),
         {:ok, residual_rows} <-
           build_epoch_sequential_baseline_rows(
             base,
             epoch,
             refs,
             physical_sats,
             sd_ambiguity_ids,
             dd_ambiguity_pairs,
             report_state,
             weights
           ),
         {:ok, residuals} <- baseline_residuals(residual_rows) do
      all_fixed = fixed_cycles |> Map.keys() |> Enum.sort()
      fixed? = all_fixed != []

      epoch_result = %{
        epoch: epoch.epoch,
        index: epoch.idx,
        baseline_m: ecef_map(report_state.baseline),
        integer_status: if(fixed?, do: :fixed, else: :not_fixed),
        integer_ratio: Map.get(search_meta, :integer_ratio),
        integer_best_score: Map.get(search_meta, :integer_best_score),
        integer_second_best_score: Map.get(search_meta, :integer_second_best_score),
        integer_candidates: Map.get(search_meta, :integer_candidates),
        ambiguity_search: Map.get(search_meta, :ambiguity_search),
        residuals_m: residuals,
        newly_fixed_ambiguities: newly_fixed,
        fixed_ambiguities: all_fixed,
        innovation_screen: screen_meta
      }

      {:ok,
       %{
         acc
         | state: state,
           information: information,
           fixed_cycles: fixed_cycles,
           fixed_m: fixed_m,
           epochs: [epoch_result | acc.epochs]
       }}
    end
  end

  defp sequential_filter_epochs_rust(
         base,
         epochs,
         refs,
         physical_sats,
         dd_ambiguity_pairs,
         dd_ambiguity_satellites,
         float_only_dd_ids,
         float_only_systems,
         wavelengths,
         offsets,
         weights,
         solve_opts,
         integer_opts,
         filter_opts,
         acc
       ) do
    rust_wavelengths = rust_sd_keyed_values(dd_ambiguity_pairs, wavelengths)
    rust_offsets = rust_sd_keyed_values(dd_ambiguity_pairs, offsets)
    rust_epochs = Enum.map(epochs, &rust_epoch_term(&1, refs, physical_sats))

    case NIF.rtk_filter_update_epochs(
           acc.rust_state,
           rust_epochs,
           base,
           rust_model_term(weights),
           rust_wavelengths,
           rust_offsets,
           rust_update_opts_term(filter_opts, solve_opts, integer_opts, float_only_systems)
         ) do
      {:ok, updates} ->
        epochs
        |> Enum.zip(updates)
        |> Enum.reduce_while({:ok, acc}, fn {epoch, update}, {:ok, acc} ->
          case apply_rust_filter_update(
                 update,
                 base,
                 epoch,
                 refs,
                 physical_sats,
                 dd_ambiguity_pairs,
                 dd_ambiguity_satellites,
                 float_only_dd_ids,
                 weights,
                 acc
               ) do
            {:ok, next} -> {:cont, {:ok, next}}
            {:error, _reason} = err -> {:halt, err}
          end
        end)

      # The batch NIF tags a mid-arc failure with its epoch index
      # ({:error, epoch_index, reason}); carry both to the caller.
      {:error, epoch_index, reason} ->
        {:error, {reason, epoch_index: epoch_index}}

      {:error, _reason} = err ->
        err
    end
  end

  defp apply_rust_filter_update(
         {rust_state, reported_baseline, ratio, fixed?, newly_fixed_sd, fixed_sd_ids,
          innovation_screen},
         base,
         epoch,
         refs,
         physical_sats,
         dd_ambiguity_pairs,
         dd_ambiguity_satellites,
         float_only_dd_ids,
         weights,
         acc
       ) do
    with {:ok, state, information, header_refs, sd_fixed_cycles, sd_fixed_m, sd_ambiguity_ids} <-
           rust_filter_state_to_elixir(rust_state),
         # The carried state keeps the float baseline; the reported solution is
         # the kernel's ambiguity-conditioned baseline (matches the Elixir path).
         report_state = %{state | baseline: reported_baseline},
         {:ok, fixed_cycles} <-
           rust_sd_fixed_map_to_dd_ids(
             sd_fixed_cycles,
             dd_ambiguity_pairs,
             dd_ambiguity_satellites,
             refs,
             header_refs
           ),
         {:ok, fixed_m} <-
           rust_sd_fixed_map_to_dd_ids(
             sd_fixed_m,
             dd_ambiguity_pairs,
             dd_ambiguity_satellites,
             refs,
             header_refs
           ),
         {:ok, newly_fixed} <-
           rust_sd_ids_to_dd_ids(
             newly_fixed_sd,
             dd_ambiguity_pairs,
             dd_ambiguity_satellites,
             refs,
             header_refs
           ),
         {:ok, all_fixed} <-
           rust_sd_ids_to_dd_ids(
             fixed_sd_ids,
             dd_ambiguity_pairs,
             dd_ambiguity_satellites,
             refs,
             header_refs
           ),
         {:ok, residual_rows} <-
           build_epoch_sequential_baseline_rows(
             base,
             epoch,
             refs,
             physical_sats,
             sd_ambiguity_ids,
             dd_ambiguity_pairs,
             report_state,
             weights
           ),
         integer_ratio = rust_public_integer_ratio(ratio, residual_rows, acc, float_only_dd_ids),
         {:ok, residuals} <- baseline_residuals(residual_rows) do
      screen_meta = rust_innovation_screen_meta(innovation_screen)

      epoch_result = %{
        epoch: epoch.epoch,
        index: epoch.idx,
        baseline_m: ecef_map(report_state.baseline),
        integer_status:
          if(screen_meta && screen_meta.coasted?,
            do: :coasted,
            else: if(fixed?, do: :fixed, else: :not_fixed)
          ),
        integer_ratio: integer_ratio,
        integer_best_score: nil,
        integer_second_best_score: nil,
        integer_candidates: nil,
        ambiguity_search: nil,
        residuals_m: residuals,
        newly_fixed_ambiguities: newly_fixed,
        fixed_ambiguities: all_fixed,
        innovation_screen: screen_meta
      }

      {:ok,
       %{
         acc
         | state: state,
           information: information,
           fixed_cycles: fixed_cycles,
           fixed_m: fixed_m,
           rust_state: rust_state,
           epochs: [epoch_result | acc.epochs]
       }}
    end
  end

  defp rust_public_integer_ratio(ratio, residual_rows, acc, float_only_dd_ids) do
    fixed_set = acc.fixed_cycles |> Map.keys() |> MapSet.new()

    search_ids =
      residual_rows
      |> Enum.map(& &1.ambiguity_id)
      |> Enum.uniq()
      |> Enum.sort()
      |> Enum.reject(&(MapSet.member?(fixed_set, &1) or MapSet.member?(float_only_dd_ids, &1)))

    if search_ids != [], do: ratio
  end

  defp rust_initial_filter_state(
         epochs,
         refs,
         initial_baseline,
         sd_ambiguity_ids,
         initial_ambiguities,
         baseline_prior_sigma_m,
         ambiguity_prior_sigma_m
       ) do
    # Per-system reference SD ambiguity ids, each taken from the first epoch in
    # which that system's reference is observed (a reference is present in every
    # epoch where its system appears, but a system may join the arc late).
    header_refs =
      refs
      |> Enum.sort()
      |> Enum.map(fn {system, reference_sat} ->
        epoch = Enum.find(epochs, &MapSet.member?(epoch_sats(&1), reference_sat))
        {system, single_difference(epoch, reference_sat).ambiguity_id}
      end)

    # Pre-size the kernel state with the reference's globally-sorted column order
    # (sd_ambiguity_ids) and seeds, so the kernel's information matrix is
    # column-identical to the Elixir reference. Otherwise the kernel grows columns
    # by first-sighting insertion (reference first) — a permutation of the sorted
    # order that takes different partial pivots in the solve, breaking 0-ULP.
    # `ensure_ambiguity` is idempotent, so the per-epoch seeding is a no-op.
    n = baseline_unknown_count(sd_ambiguity_ids)

    information =
      sequential_initial_information(n, baseline_prior_sigma_m, ambiguity_prior_sigma_m)

    ambiguities_m = Enum.map(sd_ambiguity_ids, &Map.fetch!(initial_ambiguities, &1))

    {
      {@rtk_filter_state_version, header_refs, sd_ambiguity_ids, ambiguity_prior_sigma_m, 0},
      initial_baseline,
      ambiguities_m,
      List.flatten(information),
      [],
      []
    }
  end

  defp rust_epoch_term(epoch, refs, physical_sats) do
    available = epoch_sats(epoch)

    references =
      refs
      |> Enum.sort()
      |> Enum.filter(fn {_system, reference_sat} -> MapSet.member?(available, reference_sat) end)
      |> Enum.map(fn {_system, reference_sat} -> rust_sat_term(epoch, reference_sat) end)

    nonref =
      epoch
      |> epoch_available_nonrefs(refs, physical_sats)
      |> Enum.map(&rust_sat_term(epoch, &1))

    {references, nonref, Map.get(epoch, :velocity_mps), Map.get(epoch, :prediction_dt_s, 0.0)}
  end

  defp rust_sat_term(epoch, sat) do
    base = Map.fetch!(epoch.base, sat)
    rover = Map.fetch!(epoch.rover, sat)

    {
      {sat, single_difference_ambiguity_id(sat, base, rover)},
      {base.code_m, base.phase_m, rover.code_m, rover.phase_m},
      {
        Map.fetch!(epoch.base_positions, sat),
        Map.fetch!(epoch.rover_positions, sat),
        Map.fetch!(epoch.positions, sat)
      }
    }
  end

  defp rust_model_term(weights) do
    {
      weights.code_sigma_m,
      weights.phase_sigma_m,
      Atom.to_string(weights.stochastic_model),
      weights.elevation_weighting?,
      weights.sagnac?
    }
  end

  defp rust_update_opts_term(filter_opts, solve_opts, integer_opts, float_only_systems) do
    {
      filter_opts.hold_sigma_m,
      solve_opts.position_tolerance_m,
      solve_opts.ambiguity_tolerance_m,
      solve_opts.max_iterations,
      filter_opts.process_noise_baseline_sigma_m,
      integer_opts.ratio_threshold,
      {
        Atom.to_string(filter_opts.dynamics_model),
        float_only_systems,
        rust_innovation_screen_sigma(filter_opts.innovation_screen),
        filter_opts.innovation_screen.min_rows
      }
    }
  end

  defp rust_innovation_screen_sigma(%{enabled?: true, threshold_sigma: threshold}), do: threshold
  defp rust_innovation_screen_sigma(_screen), do: 0.0

  defp rust_innovation_screen_meta(nil), do: nil

  defp rust_innovation_screen_meta(
         {threshold, min_rows, input_rows, accepted_rows, rejected_rows, rejected_code_rows,
          {rejected_phase_rows, max_normalized, max_rejected_normalized, coasted?}}
       ) do
    %{
      threshold_sigma: threshold,
      min_rows: min_rows,
      input_rows: input_rows,
      accepted_rows: accepted_rows,
      rejected_rows: rejected_rows,
      rejected_code_rows: rejected_code_rows,
      rejected_phase_rows: rejected_phase_rows,
      max_abs_normalized_innovation: max_normalized,
      max_rejected_abs_normalized_innovation: max_rejected_normalized,
      coasted?: coasted?
    }
  end

  defp rust_sd_keyed_values(dd_ambiguity_pairs, values) do
    dd_ambiguity_pairs
    |> Map.new(fn {dd_id, %{sat_sd_id: sat_sd_id}} ->
      {sat_sd_id, Map.fetch!(values, dd_id)}
    end)
    |> Enum.sort()
  end

  defp rust_filter_state_to_elixir(
         {{@rtk_filter_state_version, header_refs, sd_ids, _ambiguity_prior_sigma_m,
           _epoch_count}, baseline, sd_ambiguities, information, fixed_cycles, fixed_m}
       )
       when length(sd_ids) == length(sd_ambiguities) do
    n = 3 + length(sd_ids)

    if length(information) == n * n do
      {:ok,
       %{
         baseline: baseline,
         ambiguities: sd_ids |> Enum.zip(sd_ambiguities) |> Map.new()
       }, unflatten_matrix(information, n), Map.new(header_refs), Map.new(fixed_cycles),
       Map.new(fixed_m), sd_ids}
    else
      {:error, {:invalid_rust_filter_state, :information_dimension}}
    end
  end

  defp rust_filter_state_to_elixir(_state), do: {:error, {:invalid_rust_filter_state, :shape}}

  defp unflatten_matrix(values, n), do: Enum.chunk_every(values, n)

  defp rust_sd_fixed_map_to_dd_ids(
         values,
         dd_ambiguity_pairs,
         dd_ambiguity_satellites,
         refs,
         header_refs
       ) do
    values
    |> Enum.reduce_while({:ok, %{}}, fn {sat_sd_id, value}, {:ok, acc} ->
      case rust_sd_id_to_dd_id(
             sat_sd_id,
             dd_ambiguity_pairs,
             dd_ambiguity_satellites,
             refs,
             header_refs
           ) do
        {:ok, dd_id} -> {:cont, {:ok, Map.put(acc, dd_id, value)}}
        {:error, _reason} = err -> {:halt, err}
      end
    end)
  end

  defp rust_sd_ids_to_dd_ids(ids, dd_ambiguity_pairs, dd_ambiguity_satellites, refs, header_refs) do
    ids
    |> Enum.reduce_while({:ok, []}, fn sat_sd_id, {:ok, acc} ->
      case rust_sd_id_to_dd_id(
             sat_sd_id,
             dd_ambiguity_pairs,
             dd_ambiguity_satellites,
             refs,
             header_refs
           ) do
        {:ok, dd_id} -> {:cont, {:ok, [dd_id | acc]}}
        {:error, _reason} = err -> {:halt, err}
      end
    end)
    |> case do
      {:ok, dd_ids} -> {:ok, Enum.sort(dd_ids)}
      {:error, _reason} = err -> err
    end
  end

  # Kernel SD ambiguity ids map back to the reference DD ids against their OWN
  # system's reference (the kernel header carries each system's reference SD
  # ambiguity arc, fixed for the run by the ReferenceChanged guard).
  defp rust_sd_id_to_dd_id(
         sat_sd_id,
         dd_ambiguity_pairs,
         dd_ambiguity_satellites,
         refs,
         header_refs
       ) do
    system = String.first(sat_sd_id)
    reference_sat = Map.fetch!(refs, system)
    ref_sd_id = Map.fetch!(header_refs, system)

    case Enum.find(dd_ambiguity_pairs, fn {_dd_id, pair} ->
           pair.sat_sd_id == sat_sd_id and pair.ref_sd_id == ref_sd_id
         end) do
      {dd_id, _pair} ->
        {:ok, dd_id}

      nil ->
        sat =
          Enum.find_value(dd_ambiguity_pairs, fn {dd_id, pair} ->
            if pair.sat_sd_id == sat_sd_id do
              Map.fetch!(dd_ambiguity_satellites, dd_id)
            end
          end)

        if sat do
          {:ok,
           double_difference_ambiguity_id(sat, sat_sd_id, %{
             satellite_id: reference_sat,
             ambiguity_id: ref_sd_id
           })}
        else
          {:error, {:unknown_rust_ambiguity_id, sat_sd_id}}
        end
    end
  end

  defp iterate_sequential_filter_epoch(
         base,
         epoch,
         refs,
         physical_sats,
         sd_ambiguity_ids,
         dd_ambiguity_pairs,
         prior_center,
         prior_state,
         prior_information,
         fixed_m,
         weights,
         solve_opts,
         filter_opts,
         row_mask,
         iter
       ) do
    with {:ok, rows} <-
           build_epoch_sequential_baseline_rows(
             base,
             epoch,
             refs,
             physical_sats,
             sd_ambiguity_ids,
             dd_ambiguity_pairs,
             prior_state,
             weights
           ),
         # Innovation-screen accept/reject decisions are made once, at the
         # epoch's predicted prior (kernel `prepare_innovation_screen`); the
         # index mask is applied to the FOLD INPUT on every relinearization
         # iteration. The full row set still flows to the search/residual paths.
         {measurement_information, measurement_rhs} <-
           baseline_normal_equations(
             apply_innovation_row_mask(rows, row_mask),
             baseline_unknown_count(sd_ambiguity_ids)
           ),
         {hold_information, hold_rhs} <-
           sequential_hold_normal_equations(
             sd_ambiguity_ids,
             dd_ambiguity_pairs,
             prior_state,
             fixed_m,
             filter_opts.hold_sigma_m
           ),
         prior_rhs =
           sequential_prior_rhs(sd_ambiguity_ids, prior_information, prior_center, prior_state),
         information =
           prior_information
           |> matrix_add(measurement_information)
           |> matrix_add(hold_information),
         rhs = measurement_rhs |> vector_add(hold_rhs) |> vector_add(prior_rhs),
         {information, rhs} =
           apply_reference_sd_gauge(
             information,
             rhs,
             epoch,
             refs,
             sd_ambiguity_ids,
             prior_center,
             prior_state,
             filter_opts.hold_sigma_m
           ),
         {:ok, dx} <- solve_linear(information, rhs),
         next_state = apply_baseline_delta(prior_state, sd_ambiguity_ids, dx),
         {baseline_step, ambiguity_step} <- baseline_step_norms(dx) do
      cond do
        baseline_step <= solve_opts.position_tolerance_m and
            ambiguity_step <= solve_opts.ambiguity_tolerance_m ->
          {:ok, next_state, information, rows}

        iter >= solve_opts.max_iterations ->
          {:ok, next_state, information, rows}

        true ->
          iterate_sequential_filter_epoch(
            base,
            epoch,
            refs,
            physical_sats,
            sd_ambiguity_ids,
            dd_ambiguity_pairs,
            prior_center,
            next_state,
            prior_information,
            fixed_m,
            weights,
            solve_opts,
            filter_opts,
            row_mask,
            iter + 1
          )
      end
    end
  end

  defp apply_innovation_row_mask(rows, nil), do: rows

  defp apply_innovation_row_mask(rows, mask) do
    rows
    |> Enum.zip(mask)
    |> Enum.filter(fn {_row, accepted?} -> accepted? end)
    |> Enum.map(fn {row, _accepted?} -> row end)
  end

  # Elixir reference of the kernel `prepare_innovation_screen`: build the DD
  # rows at the epoch's predicted prior, normalize each prefit residual by the
  # row's diagonal weight, and reject rows whose |normalized innovation|
  # exceeds the threshold. Returns `{nil, nil}` when the screen is disabled,
  # else `{meta, mask}` where `meta` matches the kernel's per-epoch
  # `innovation_screen` metadata shape exactly and `mask` is the per-row
  # accept list (row order is deterministic on both sides).
  defp sequential_innovation_screen(
         _base,
         _epoch,
         _refs,
         _physical_sats,
         _sd_ambiguity_ids,
         _dd_ambiguity_pairs,
         _state,
         _weights,
         %{enabled?: false}
       ), do: {:ok, nil, nil}

  defp sequential_innovation_screen(
         base,
         epoch,
         refs,
         physical_sats,
         sd_ambiguity_ids,
         dd_ambiguity_pairs,
         state,
         weights,
         %{enabled?: true, threshold_sigma: threshold, min_rows: min_rows}
       ) do
    with {:ok, rows} <-
           build_epoch_sequential_baseline_rows(
             base,
             epoch,
             refs,
             physical_sats,
             sd_ambiguity_ids,
             dd_ambiguity_pairs,
             state,
             weights
           ) do
      {mask_reversed, accepted, rejected, rejected_code, rejected_phase, max_norm, max_rej} =
        Enum.reduce(rows, {[], 0, 0, 0, 0, nil, nil}, fn row,
                                                         {mask, acc, rej, rej_c, rej_p, max_n,
                                                          max_r} ->
          # Kernel op order exactly: the screen normalizes by the DIAGONAL
          # information weight 1/(sd_var + ref_var) (the kernel DdRow.weight),
          # NOT this row's 1/sigma fold weight.
          weight = 1.0 / (row.sd_variance_m2 + row.ref_sd_variance_m2)
          normalized = abs(row.y * weight)
          max_n = if max_n, do: max(max_n, normalized), else: normalized

          if normalized > threshold do
            rej_c = if row.kind == :code, do: rej_c + 1, else: rej_c
            rej_p = if row.kind == :phase, do: rej_p + 1, else: rej_p
            max_r = if max_r, do: max(max_r, normalized), else: normalized
            {[false | mask], acc, rej + 1, rej_c, rej_p, max_n, max_r}
          else
            {[true | mask], acc + 1, rej, rej_c, rej_p, max_n, max_r}
          end
        end)

      meta = %{
        threshold_sigma: threshold,
        min_rows: min_rows,
        input_rows: length(rows),
        accepted_rows: accepted,
        rejected_rows: rejected,
        rejected_code_rows: rejected_code,
        rejected_phase_rows: rejected_phase,
        max_abs_normalized_innovation: max_norm,
        max_rejected_abs_normalized_innovation: max_rej,
        coasted?: accepted < min_rows
      }

      {:ok, meta, Enum.reverse(mask_reversed)}
    end
  end

  # Kinematic process-noise time update (predict step). Between epochs, inflate
  # the baseline (position) block of the prior covariance by Q = sigma^2 while
  # leaving the carried ambiguities tight. This is the information<->covariance
  # round-trip: Lambda -> P, P[baseline] += Q, P -> Lambda'. The prior mean is
  # unchanged for the default constant-position model; the velocity-propagated
  # model advances it immediately before this covariance update.
  #
  # Skipped on the first epoch (no prior motion) and when process noise is off
  # (sigma == 0, the static default). A singular round-trip falls back to the
  # un-inflated prior rather than failing the epoch.
  defp time_update_information(acc, epoch, filter_opts) do
    sigma = Map.get(filter_opts, :process_noise_baseline_sigma_m, 0.0)
    mean_acc = propagate_baseline_mean(acc, epoch, filter_opts)

    if mean_acc.epochs != [] and sigma > 0.0 do
      case inflate_baseline_information(mean_acc.information, sigma * sigma) do
        {:ok, information} -> {%{mean_acc | information: information}, mean_acc}
        :error -> {mean_acc, nil}
      end
    else
      {mean_acc, nil}
    end
  end

  defp propagate_baseline_mean(
         %{epochs: [_prev | _]} = acc,
         %{velocity_mps: {vx, vy, vz}, prediction_dt_s: dt_s},
         %{dynamics_model: :velocity_propagated}
       )
       when is_number(dt_s) and dt_s > 0.0 do
    if finite_number?(dt_s) and finite_number?(vx) and finite_number?(vy) and finite_number?(vz) do
      update_in(acc.state.baseline, &add3(&1, scale3({vx, vy, vz}, dt_s)))
    else
      acc
    end
  end

  defp propagate_baseline_mean(acc, _epoch, _filter_opts), do: acc

  # Rank-3 (Woodbury) baseline process-noise inflation in information space:
  #
  #   Λ' = Λ - Λ[:,0..2] · (Q⁻¹ + Λ₀₀)⁻¹ · Λ[0..2,:],   Q = q·I₃
  #
  # Mathematically equal to adding q to the baseline covariance diagonal and
  # re-inverting, but via a single 3×3 inverse instead of two full n×n inversions.
  # The double full-invert corrupts near-singular ambiguity directions when q does
  # not dominate (the filter then goes singular at small process noise); this
  # rank-3 form is stable for any q and preserves the ambiguity block exactly.
  # Returns `:error` on a singular 3×3 system (caller keeps the un-inflated prior).
  defp inflate_baseline_information(information, q) do
    inv_q = 1.0 / q

    m =
      for i <- 0..2 do
        for j <- 0..2 do
          at2(information, i, j) + if(i == j, do: inv_q, else: 0.0)
        end
      end

    case invert_3x3(m) do
      {:ok, m_inv} ->
        # cols = Λ[:,0..2] (n×3); Λ symmetric ⇒ Λ[0..2,:] = colsᵀ.
        cols = Enum.map(information, &Enum.take(&1, 3))
        # w = cols · M⁻¹ (n×3); Λ' = Λ - w · colsᵀ.
        w =
          Enum.map(cols, fn ci ->
            for a <- 0..2 do
              Enum.reduce(0..2, 0.0, fn b, acc -> acc + Enum.at(ci, b) * at2(m_inv, b, a) end)
            end
          end)

        updated =
          [information, w]
          |> Enum.zip()
          |> Enum.map(fn {row, wi} ->
            [row, cols]
            |> Enum.zip()
            |> Enum.map(fn {value, cj} ->
              value -
                Enum.reduce(0..2, 0.0, fn a, acc -> acc + Enum.at(wi, a) * Enum.at(cj, a) end)
            end)
          end)

        {:ok, updated}

      :error ->
        :error
    end
  end

  defp at2(matrix, i, j), do: matrix |> Enum.at(i) |> Enum.at(j)

  # Closed-form 3×3 inverse (adjugate / determinant). Deterministic arithmetic so
  # the Elixir reference and the Rust kernel match bit-for-bit. `<= 1e-12` matches
  # the singular-guard convention used elsewhere in the filter.
  defp invert_3x3([[a, b, c], [d, e, f], [g, h, i]]) do
    det = a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g)

    if abs(det) <= 1.0e-12 do
      :error
    else
      inv_det = 1.0 / det

      {:ok,
       [
         [(e * i - f * h) * inv_det, (c * h - b * i) * inv_det, (b * f - c * e) * inv_det],
         [(f * g - d * i) * inv_det, (a * i - c * g) * inv_det, (c * d - a * f) * inv_det],
         [(d * h - e * g) * inv_det, (b * g - a * h) * inv_det, (a * e - b * d) * inv_det]
       ]}
    end
  end

  # Ambiguity-conditioned ("fixed") baseline for the reported solution. The float
  # state carried between epochs is the Kalman/information posterior; when AR
  # newly succeeds this epoch, re-solve THIS epoch with the new integers held so
  # the reported baseline reflects the fix in the same epoch — matching RTKLIB's
  # fixed solution. Without this, the first fixed epoch (no prior hold yet)
  # reports its float baseline while claiming `:fixed` (the cold-start false
  # confidence: ~0.9 m at epoch 0, recovering only once the next epoch applies the
  # hold). The float `state`/`information` are still what gets carried forward.
  #
  # No newly-fixed ambiguities => the float state already incorporates every held
  # integer (the prior iterate used the same fixed_m), so it is already the
  # conditioned solution. A failed re-solve falls back to the float state.
  defp conditioned_fixed_state(
         float_state,
         [],
         _b,
         _e,
         _r,
         _p,
         _sd,
         _dd,
         _acc,
         _fm,
         _w,
         _so,
         _fo,
         _mask
       ), do: {:ok, float_state}

  defp conditioned_fixed_state(
         float_state,
         _newly_fixed,
         base,
         epoch,
         refs,
         physical_sats,
         sd_ambiguity_ids,
         dd_ambiguity_pairs,
         acc,
         fixed_m,
         weights,
         solve_opts,
         filter_opts,
         row_mask
       ) do
    case iterate_sequential_filter_epoch(
           base,
           epoch,
           refs,
           physical_sats,
           sd_ambiguity_ids,
           dd_ambiguity_pairs,
           acc.state,
           acc.state,
           acc.information,
           fixed_m,
           weights,
           solve_opts,
           filter_opts,
           row_mask,
           1
         ) do
      {:ok, conditioned, _information, _rows} -> {:ok, conditioned}
      _ -> {:ok, float_state}
    end
  end

  defp sequential_initial_information(n, baseline_sigma_m, ambiguity_sigma_m) do
    for i <- 0..(n - 1) do
      for j <- 0..(n - 1) do
        cond do
          i != j -> 0.0
          i < 3 -> 1.0 / (baseline_sigma_m * baseline_sigma_m)
          true -> 1.0 / (ambiguity_sigma_m * ambiguity_sigma_m)
        end
      end
    end
  end

  # Multi-system gauge fixing. Every measurement row, hold, and reported
  # quantity depends only on within-system single-difference DIFFERENCES, so
  # the common-mode level of each system's SD ambiguity block is unobservable.
  # In the single-system filter that gauge direction is carried by the initial
  # ambiguity prior alone; once per-epoch measurement information accumulates,
  # the prior's contribution falls below float precision, and with a
  # two-satellite system the gauge pivot cancels exactly to zero — the filter
  # dies with `:singular_geometry`. Multi-system runs therefore pin each
  # per-system reference SD ambiguity at its prior value with the hold weight
  # (a pure gauge constraint: double differences and the baseline are invariant
  # to it). The single-system path is left untouched — same objects, no
  # arithmetic — to stay bit-identical to the deployed reference and its Rust
  # kernel.
  defp apply_reference_sd_gauge(
         information,
         rhs,
         epoch,
         refs,
         sd_ambiguity_ids,
         prior_center,
         prior_state,
         hold_sigma_m
       )
       when map_size(refs) > 1 do
    weight = 1.0 / (hold_sigma_m * hold_sigma_m)
    available = epoch_sats(epoch)

    refs
    |> Enum.sort()
    |> Enum.reduce({information, rhs}, fn {_system, ref}, {information, rhs} ->
      if MapSet.member?(available, ref) do
        sd_id = single_difference(epoch, ref).ambiguity_id
        col = 3 + Enum.find_index(sd_ambiguity_ids, &(&1 == sd_id))

        residual =
          Map.fetch!(prior_center.ambiguities, sd_id) -
            Map.fetch!(prior_state.ambiguities, sd_id)

        information =
          List.update_at(information, col, fn row ->
            List.update_at(row, col, &(&1 + weight))
          end)

        {information, List.update_at(rhs, col, &(&1 + weight * residual))}
      else
        {information, rhs}
      end
    end)
  end

  defp apply_reference_sd_gauge(
         information,
         rhs,
         _epoch,
         _refs,
         _sd_ambiguity_ids,
         _prior_center,
         _prior_state,
         _hold_sigma_m
       ), do: {information, rhs}

  defp sequential_prior_rhs(ambiguity_ids, prior_information, prior_center, current_state) do
    prior_center
    |> state_vector(ambiguity_ids)
    |> vector_sub(state_vector(current_state, ambiguity_ids))
    |> then(&matvec(prior_information, &1))
  end

  defp sequential_hold_normal_equations(
         sd_ambiguity_ids,
         dd_ambiguity_pairs,
         state,
         fixed_m,
         hold_sigma_m
       ) do
    n = baseline_unknown_count(sd_ambiguity_ids)
    hold_weight = 1.0 / (hold_sigma_m * hold_sigma_m)
    zero = zero_matrix(n)
    rhs = zero_vector(n)

    fixed_m
    |> Enum.reduce({zero, rhs}, fn {dd_id, fixed}, {normal, rhs} ->
      %{sat_sd_id: sat_sd_id, ref_sd_id: ref_sd_id} = Map.fetch!(dd_ambiguity_pairs, dd_id)
      current = dd_state_ambiguity_m(state, sat_sd_id, ref_sd_id)

      h =
        design_single_difference_baseline_row(
          {0.0, 0.0, 0.0},
          sat_sd_id,
          ref_sd_id,
          sd_ambiguity_ids
        )

      weighted_h = Enum.map(h, &(&1 * hold_weight))

      normal =
        Enum.with_index(h)
        |> Enum.reduce(normal, fn {hi, i}, normal ->
          Enum.with_index(weighted_h)
          |> Enum.reduce(normal, fn {weighted_hj, j}, normal ->
            add_matrix_value(normal, i, j, hi * weighted_hj)
          end)
        end)

      rhs =
        Enum.with_index(weighted_h)
        |> Enum.reduce(rhs, fn {weighted_hi, i}, rhs ->
          add_vector_value(rhs, i, (fixed - current) * weighted_hi)
        end)

      {normal, rhs}
    end)
  end

  defp sequential_search_and_hold(
         state,
         covariance,
         rows,
         sd_ambiguity_ids,
         dd_ambiguity_ids,
         dd_ambiguity_pairs,
         float_only_dd_ids,
         fixed_cycles,
         wavelengths,
         offsets,
         integer_opts
       ) do
    observed_ids =
      rows
      |> Enum.map(& &1.ambiguity_id)
      |> Enum.uniq()
      |> Enum.sort()

    fixed_set = fixed_cycles |> Map.keys() |> MapSet.new()

    # Float-only systems' ambiguities contribute measurement rows but are never
    # entered into the LAMBDA search set; they stay float.
    search_ids =
      Enum.reject(
        observed_ids,
        &(MapSet.member?(fixed_set, &1) or MapSet.member?(float_only_dd_ids, &1))
      )

    if search_ids == [] do
      {:ok, fixed_cycles, fixed_ambiguities_m(fixed_cycles, wavelengths, offsets), %{}}
    else
      sd_covariance =
        submatrix(covariance, 3, length(sd_ambiguity_ids), 3, length(sd_ambiguity_ids))

      covariance_cycles =
        dd_covariance_m_to_cycles(
          sd_covariance,
          sd_ambiguity_ids,
          dd_ambiguity_ids,
          dd_ambiguity_pairs,
          wavelengths
        )

      subset_covariance = covariance_submatrix(dd_ambiguity_ids, covariance_cycles, search_ids)

      float_cycles =
        Map.new(search_ids, fn id ->
          %{sat_sd_id: sat_sd_id, ref_sd_id: ref_sd_id} = Map.fetch!(dd_ambiguity_pairs, id)
          ambiguity_m = dd_state_ambiguity_m(state, sat_sd_id, ref_sd_id)
          offset_m = Map.fetch!(offsets, id)
          {id, (ambiguity_m - offset_m) / Map.fetch!(wavelengths, id)}
        end)

      case IntegerLeastSquares.search(float_cycles, subset_covariance, integer_opts) do
        {:ok, new_fixed, %{integer_status: :fixed} = meta} ->
          meta = Map.put(meta, :ambiguity_offsets_m, Map.take(offsets, search_ids))

          fixed_cycles =
            Map.merge(fixed_cycles, new_fixed)

          {:ok, fixed_cycles, fixed_ambiguities_m(fixed_cycles, wavelengths, offsets), meta}

        {:ok, _new_fixed, meta} ->
          meta = Map.put(meta, :ambiguity_offsets_m, Map.take(offsets, search_ids))
          {:ok, fixed_cycles, fixed_ambiguities_m(fixed_cycles, wavelengths, offsets), meta}

        {:error, _reason} = err ->
          err
      end
    end
  end

  defp iterate_fixed_baseline(
         base,
         epochs,
         refs,
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
             refs,
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
            refs,
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
            refs,
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
            refs,
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
         refs,
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
             refs,
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
         refs,
         physical_sats,
         _ambiguity_ids,
         free_ambiguity_ids,
         fixed_m,
         state,
         weights
       ) do
    ref_data = epoch_reference_data(epoch, refs)

    epoch
    |> epoch_available_nonrefs(refs, physical_sats)
    |> Enum.reduce({:ok, []}, fn sat, {:ok, acc} ->
      %{sd: ref_dd, pos: ref_pos, base_pos: ref_base_pos, rover_pos: ref_rover_pos} =
        Map.fetch!(ref_data, satellite_system(sat))

      obs_dd = double_difference_measurement(epoch, sat, ref_dd)
      sat_pos = Map.fetch!(epoch.positions, sat)
      sat_base_pos = Map.fetch!(epoch.base_positions, sat)
      sat_rover_pos = Map.fetch!(epoch.rover_positions, sat)

      {geom_dd, deriv} =
        geometry_double_difference(
          base,
          state.baseline,
          sat_base_pos,
          sat_rover_pos,
          ref_base_pos,
          ref_rover_pos,
          weights
        )

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
        ref_sat: ref_dd.satellite_id,
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
        ref_sat: ref_dd.satellite_id,
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
         refs,
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
             refs,
             physical_sats,
             ambiguity_ids,
             free_ambiguity_ids,
             fixed_m,
             state,
             weights
           ),
         {:ok, residuals} <- baseline_residuals(rows) do
      code_residuals = Enum.map(residuals, & &1.code_m)
      phase_residuals = Enum.map(residuals, & &1.phase_m)
      rover = add3(base, state.baseline)

      {:ok,
       %FixedBaselineSolution{
         baseline_m: ecef_map(state.baseline),
         rover_position_m: ecef_map(rover),
         reference_satellite_id: reference_satellite_report(refs),
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
             reference_satellites: refs,
             ambiguity_satellites: Map.fetch!(float_sol.metadata, :ambiguity_satellites),
             dropped_cycle_slip_sats: Map.get(float_sol.metadata, :dropped_cycle_slip_sats, []),
             elevation_mask_deg: Map.get(float_sol.metadata, :elevation_mask_deg),
             elevation_masked_sats: Map.get(float_sol.metadata, :elevation_masked_sats, []),
             split_cycle_slip_arcs: Map.get(float_sol.metadata, :split_cycle_slip_arcs, []),
             measurement_covariance: %{
               model: :double_difference,
               code_sigma_m: weights.code_sigma_m,
               phase_sigma_m: weights.phase_sigma_m,
               stochastic_model: weights.stochastic_model,
               elevation_weighting: weights.elevation_weighting?,
               sagnac: weights.sagnac?,
               min_elevation_sin: @min_elevation_sin
             }
           })
       }}
    end
  end

  defp baseline_residuals(rows) do
    rows
    |> Enum.group_by(&{&1.epoch, &1.sat, &1.ambiguity_id})
    |> Enum.reduce_while({:ok, []}, fn {{epoch, sat, ambiguity_id}, grouped}, {:ok, acc} ->
      code = Enum.find(grouped, &(&1.kind == :code))
      phase = Enum.find(grouped, &(&1.kind == :phase))

      case {code, phase} do
        {%{y: code_y, weight: code_weight}, %{y: phase_y, weight: phase_weight}} ->
          residual = %{
            epoch: epoch,
            satellite_id: sat,
            reference_satellite_id: code.ref_sat,
            ambiguity_id: ambiguity_id,
            code_m: code_y,
            phase_m: phase_y,
            code_sigma_m: 1.0 / code_weight,
            phase_sigma_m: 1.0 / phase_weight,
            code_normalized: code_y * code_weight,
            phase_normalized: phase_y * phase_weight
          }

          {:cont, {:ok, [residual | acc]}}

        {nil, nil} ->
          {:halt,
           {:error, {:incomplete_residual_pair, epoch, sat, ambiguity_id, [:code, :phase]}}}

        {nil, _phase} ->
          {:halt, {:error, {:incomplete_residual_pair, epoch, sat, ambiguity_id, [:code]}}}

        {_code, nil} ->
          {:halt, {:error, {:incomplete_residual_pair, epoch, sat, ambiguity_id, [:phase]}}}
      end
    end)
    |> case do
      {:ok, residuals} ->
        {:ok, Enum.sort_by(residuals, &{inspect(&1.epoch), &1.satellite_id, &1.ambiguity_id})}

      {:error, _reason} = err ->
        err
    end
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
    # Rows correlate only through their shared reference single difference, so
    # the covariance is block-diagonal per {epoch, kind, reference}: each
    # system's double differences form their own block and there is no
    # cross-system correlation. With one system this reduces to the historical
    # {epoch, kind} blocking (the reference key is constant).
    rows
    |> Enum.group_by(&{&1.epoch_idx, &1.kind, &1.ref_sat})
    |> Enum.sort_by(fn {{epoch_idx, kind, ref_sat}, _rows} -> {epoch_idx, kind, ref_sat} end)
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

  defp geometric_range(sat_pos, receiver, %{sagnac?: true}) do
    {sx, sy, _sz} = sat_pos
    {rx, ry, _rz} = receiver

    range(sat_pos, receiver) +
      Constants.earth_rotation_rate_rad_s() * (sx * ry - sy * rx) /
        Constants.speed_of_light_m_s()
  end

  defp geometric_range(sat_pos, receiver, %{sagnac?: false}), do: range(sat_pos, receiver)

  # RTKLIB's first-order Sagnac model corrects the scalar range but keeps the
  # design row as the Euclidean line-of-sight unit vector.
  defp geometric_range_derivative(receiver, sat_pos, _weights),
    do: range_derivative(receiver, sat_pos)

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

  # Per-system reference selection for the two-receiver double-difference
  # primitive. The default keeps the historical deterministic choice — the
  # lexicographically first common satellite — applied per system.
  defp reference_satellites(common, opts) do
    systems = satellite_systems(common)
    common_by_system = Enum.group_by(common, &satellite_system/1)

    case Keyword.get(opts, :reference_satellite_id) do
      nil ->
        {:ok,
         Map.new(systems, fn system ->
           {system, common_by_system |> Map.fetch!(system) |> hd()}
         end)}

      sat when is_binary(sat) ->
        case systems do
          [system] ->
            if sat in Map.fetch!(common_by_system, system) do
              {:ok, %{system => sat}}
            else
              {:error, {:reference_satellite_missing, sat}}
            end

          _multiple ->
            {:error, {:reference_satellite_single_system, sat}}
        end

      refs when is_map(refs) ->
        Enum.reduce_while(systems, {:ok, %{}}, fn system, {:ok, acc} ->
          case Map.fetch(refs, system) do
            :error ->
              {:halt, {:error, {:reference_satellite_missing_system, system}}}

            {:ok, sat} when is_binary(sat) ->
              if sat in Map.fetch!(common_by_system, system) do
                {:cont, {:ok, Map.put(acc, system, sat)}}
              else
                {:halt, {:error, {:reference_satellite_missing, sat}}}
              end

            {:ok, _other} ->
              {:halt, {:error, {:invalid_option, :reference_satellite_id}}}
          end
        end)

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

  defp vector_add(a, b) do
    a
    |> Enum.zip(b)
    |> Enum.map(fn {x, y} -> x + y end)
  end

  defp vector_sub(a, b) do
    a
    |> Enum.zip(b)
    |> Enum.map(fn {x, y} -> x - y end)
  end

  defp state_vector(state, ambiguity_ids) do
    {x, y, z} = state.baseline
    [x, y, z | Enum.map(ambiguity_ids, &Map.fetch!(state.ambiguities, &1))]
  end

  defp add_matrix_value(matrix, row_idx, col_idx, value) do
    matrix
    |> Enum.with_index()
    |> Enum.map(fn {row, i} ->
      if i == row_idx do
        add_vector_value(row, col_idx, value)
      else
        row
      end
    end)
  end

  defp add_vector_value(vector, idx, value) do
    vector
    |> Enum.with_index()
    |> Enum.map(fn {x, col_idx} -> if col_idx == idx, do: x + value, else: x end)
  end

  defp dot3({ax, ay, az}, {bx, by, bz}), do: ax * bx + ay * by + az * bz
end
