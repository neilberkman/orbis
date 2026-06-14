defmodule Orbis.GNSS.Positioning do
  @moduledoc """
  GNSS single-point positioning (SPP): recover a receiver position, clock bias,
  and geometry diagnostics from one epoch of pseudorange observations against a
  precise SP3 ephemeris or a broadcast navigation product.

  This is the Elixir surface over the `astrodynamics-gnss` SPP solver. Given an
  ephemeris source — an `Orbis.GNSS.SP3` product or an `Orbis.GNSS.Broadcast`
  handle (GPS / Galileo / BeiDou / GLONASS) — a set of single-frequency
  pseudoranges, the receive epoch, and the broadcast/atmosphere parameters, it
  runs the transmit-time iteration and trust-region least-squares solve in the
  crate and returns an `Orbis.GNSS.Positioning.Solution`. A mixed-constellation
  set is solved together with one receiver clock per system. No positioning math
  lives on the Elixir side; this module marshals units and epoch arguments and
  decodes the result.

  ## Units at the boundary

    * pseudoranges and the initial guess position/clock are **meters**;
    * the recovered `position` is ITRF/IGS ECEF **meters**, matching the SP3 frame;
    * `geodetic` latitude/longitude are **radians** and height is **meters**;
    * `rx_clock_s` is **seconds**;
    * pressure is **hPa**, temperature is **kelvin**, relative humidity is a
      fraction in `[0, 1]`;
    * the Klobuchar `alpha`/`beta` coefficients are passed in their broadcast
      units.

  The epoch is interpreted in the SP3 product's own time scale (typically GPST);
  no leap-second shifting is applied. The seconds-since-J2000, second-of-day, and
  fractional day-of-year arguments the crate needs are derived from the supplied
  epoch via `Orbis.GNSS.Time`.

  ## Example

      {:ok, sp3} = Orbis.GNSS.SP3.load("igs.sp3")

      observations = [{"G01", 2.41e7}, {"G02", 2.49e7}, {"G05", 2.05e7}, {"G07", 2.30e7}]

      {:ok, solution} =
        Orbis.GNSS.Positioning.solve(sp3, observations, ~N[2020-06-24 12:00:00],
          ionosphere: true,
          troposphere: true,
          klobuchar_alpha: {1.0e-8, 2.2e-8, -6.0e-8, -1.2e-7},
          klobuchar_beta: {96_256.0, 131_072.0, -65_536.0, -589_824.0}
        )

      solution.position.x_m
      solution.rx_clock_s
  """

  alias Orbis.GNSS.Broadcast
  alias Orbis.GNSS.SP3
  alias Orbis.GNSS.Time
  alias Orbis.NIF

  defmodule Solution do
    @moduledoc """
    A single-point-positioning solution at one receive epoch.

    `position` is the converged ITRF/IGS ECEF position in meters. `geodetic` is
    the same point as `%{lat_rad, lon_rad, height_m}` when geodetic output was
    requested (the default), otherwise `nil`. `rx_clock_s` is the reference-system
    receiver clock bias in seconds; `system_clocks_s` is a map of GNSS letter
    (e.g. `"G"`, `"E"`) to that system's **absolute** receiver clock in seconds (a
    single entry for a one-system solve, one per constellation for a mixed solve).
    These are per-system clocks, not biases: the inter-system bias of a system is
    its clock minus the reference system's (`rx_clock_s`). `dop` carries the
    dilution-of-precision scalars for any full-rank geometry — a single-system
    solve uses the bit-exact four-state cofactor, a multi-system solve a general
    inverse with one clock column per constellation — and is `nil` only when the
    geometry is rank-deficient. `residuals_m` are the post-fit
    pseudorange residuals in meters, in `used_sats` order. `used_sats` are the
    contributing satellite id strings (e.g. `"G01"`); `rejected_sats` pairs each
    excluded satellite id with its reason atom (`:no_ephemeris` or
    `:low_elevation`). `metadata` reports solver iterations, convergence, the
    corrections applied, and the geometry redundancy: `used_count`, the distinct
    `systems`, the `redundancy` (degrees of freedom, `used_count - (3 + systems)`),
    and `raim_checkable?` (`redundancy >= 1`). An exactly determined fix has
    `redundancy < 1`, forcing the residuals near zero and leaving the fix
    unverifiable by RAIM.
    """
    @enforce_keys [
      :position,
      :geodetic,
      :rx_clock_s,
      :system_clocks_s,
      :dop,
      :residuals_m,
      :used_sats,
      :rejected_sats,
      :metadata
    ]
    defstruct [
      :position,
      :geodetic,
      :rx_clock_s,
      :system_clocks_s,
      :dop,
      :residuals_m,
      :used_sats,
      :rejected_sats,
      :metadata
    ]

    @type position :: %{x_m: float(), y_m: float(), z_m: float()}
    @type geodetic :: %{lat_rad: float(), lon_rad: float(), height_m: float()}
    @type dop :: %{
            gdop: float(),
            pdop: float(),
            hdop: float(),
            vdop: float(),
            tdop: float()
          }
    @type metadata :: %{
            :iterations => non_neg_integer(),
            :converged => boolean(),
            :status => atom(),
            :ionosphere_applied => boolean(),
            :troposphere_applied => boolean(),
            :used_count => non_neg_integer(),
            :systems => [String.t()],
            :redundancy => integer(),
            :raim_checkable? => boolean(),
            optional(:fde) => %{
              excluded: [{String.t(), :raim_excluded}],
              iterations: non_neg_integer()
            }
          }

    @type t :: %__MODULE__{
            position: position(),
            geodetic: geodetic() | nil,
            rx_clock_s: float(),
            system_clocks_s: %{String.t() => float()},
            dop: dop() | nil,
            residuals_m: [float()],
            used_sats: [String.t()],
            rejected_sats: [{String.t(), :no_ephemeris | :low_elevation}],
            metadata: metadata()
          }
  end

  @typedoc "A `{satellite_id, pseudorange_m}` pseudorange observation."
  @type observation :: {String.t(), number()}

  @typedoc "An epoch as a `NaiveDateTime` or `{{y, m, d}, {h, min, s}}` tuple."
  @type epoch :: NaiveDateTime.t() | tuple()

  @default_initial_guess {0.0, 0.0, 0.0, 0.0}
  # Mean Earth radius, used to place coarse cold-start seeds on a near-surface
  # shell. The crate freezes its elevation mask and weights at the seed geometry,
  # so a seed must be a near-surface point for the visible-satellite geometry to
  # be physical; the seed is only an input and is never gated.
  @mean_earth_radius_m 6_371_000.0
  # Default number of golden-spiral surface seeds when `:coarse_search` is on
  # without an explicit count. Pinned by the cold-start measurement: this many
  # seeds cover the sphere finely enough that at least one lands in the
  # convergence basin on the ESBC00DNK GPS-L1 arc, with margin under the bar.
  @default_coarse_seeds 24
  # A genuinely converged single-point fix has post-fit residuals at the
  # measurement-noise scale (metres to tens of metres). A converged-flagged
  # solution whose post-fit residual RMS exceeds this sanity bound did not
  # actually converge (for example a degenerate first step from the earth-centre
  # default seed that trips the step-tolerance test without moving), so it is
  # refused rather than returned as a plausible position.
  @max_converged_residual_rms_m 1.0e4
  # Position plausibility band, as geocentric radius. A real receiver fix sits
  # between roughly the polar surface (minus a margin for deep points) and a
  # generous low-Earth-orbit ceiling. A fix outside this band did not converge to
  # a physical receiver position (the earth-center default seed lands near radius
  # zero; a wrong-root least-squares fix lands far out), and is caught even when
  # the residuals are forced to zero by an exactly determined geometry.
  @min_plausible_radius_m 6_344_752.0
  @max_plausible_radius_m 8_378_137.0
  @default_alpha {0.0, 0.0, 0.0, 0.0}
  @default_beta {0.0, 0.0, 0.0, 0.0}
  # Standard-atmosphere surface meteorology, used when the troposphere term is
  # enabled and the caller does not override it.
  @default_pressure_hpa 1013.25
  @default_temperature_k 288.15
  @default_relative_humidity 0.5

  @doc """
  Solve single-point positioning for one receive epoch.

  `source` is a loaded ephemeris product — an `Orbis.GNSS.SP3` precise product or an
  `Orbis.GNSS.Broadcast` broadcast-navigation product. `observations` is a
  list of `{satellite_id, pseudorange_m}` pairs (ids like `"G01"`, pseudoranges
  in meters), and `epoch` is a `NaiveDateTime` or `{{y, m, d}, {h, min, s}}`
  tuple in the product's time scale.

  ## Options

    * `:ionosphere` - apply the broadcast Klobuchar ionosphere correction (default
      `false`); the L1 delay is scaled to each satellite's carrier by `(f_L1/f)^2`,
      covering GPS L1, Galileo E1, and BeiDou B1I (a system with no modeled
      single-frequency carrier, e.g. GLONASS, is rejected)
    * `:troposphere` - apply the Saastamoinen/Niell troposphere correction (default `false`)
    * `:klobuchar_alpha` - broadcast alpha coefficients, 4-tuple (default zeros)
    * `:klobuchar_beta` - broadcast beta coefficients, 4-tuple (default zeros)
    * `:pressure_hpa` - surface pressure, hPa (default `1013.25`)
    * `:temperature_k` - surface temperature, kelvin (default `288.15`)
    * `:relative_humidity` - relative humidity fraction `[0, 1]` (default `0.5`)
    * `:initial_guess` - `{x_m, y_m, z_m, b_m}` start point (default all zeros)
    * `:with_geodetic` - also return the geodetic position (default `true`)
    * `:max_pdop` - optional positive PDOP ceiling. When set, a fix whose
      geometry is rank-deficient or whose PDOP exceeds the ceiling is refused
      with `{:error, {:degenerate_geometry, pdop}}` (a non-positive value is
      `{:error, {:invalid_option, :max_pdop}}`); default unset.

  ### Robustness (opt-in, default off)

    * `:robust` - when `true` (default `false`), route the solve through the
      `Orbis.GNSS.QC` leave-one-out fault-detection-and-exclusion (FDE) loop: a
      chi-square RAIM test on the post-fit residuals flags an inconsistent set,
      the satellite with the largest normalized residual is excluded, and the
      position is re-solved, repeating until the set is self-consistent or the
      geometry runs out of redundancy. A clean set excludes nothing; the cleaned
      re-solve still runs through every integrity gate above, so `:robust`
      composes with the rank, `:max_pdop`, plausibility, and convergence
      refusals. The returned `Solution` carries the FDE summary at
      `solution.metadata.fde` as `%{excluded: [{sat, :raim_excluded}],
      iterations: n}`. If the loop hits its iteration cap (or runs out of
      satellites) while RAIM still flags a fault, the set could not be made
      self-consistent and the solve returns
      `{:error, {:fault_unresolved, statistic}}` rather than a still-faulted fix.
    * `:weights` - RAIM detection weights for the `:robust` test, a
      `%{sat => inverse_variance_weight}` map (`1 / sigma_i^2`) as built by
      `Orbis.GNSS.QC.weight_vector/2` from per-satellite elevation and/or C/N0.
      When given, it must carry a positive weight for every observed satellite;
      a partial map (missing an observed satellite) is rejected with
      `{:error, {:robust_requires_noise_model, :incomplete_weights}}`, and a
      non-positive, non-finite, or absurdly large observed weight with
      `{:error, {:invalid_option, :weights}}`, rather than silently falling back
      to unit weight or overflowing the detection statistic.
      A realistic noise model is REQUIRED for `:robust`: with unit weights the
      RAIM test reads ordinary code noise (several metres on a phone) as faults
      and over-excludes, degrading the fix. `:robust` therefore refuses to run
      without a noise basis: `:robust true` with no `:weights` and no
      `:unsafe_unit_weights` returns
      `{:error, {:robust_requires_noise_model, :no_weights}}` before any solve.
      Ignored when `:robust` is not set.
    * `:unsafe_unit_weights` - explicit escape hatch (default `false`). Set to
      `true` to run `:robust` FDE with unit weights (sigma = 1 m assumed), the
      only route to unit-weight FDE. Named to be self-documenting and grep-able
      because unit-weight FDE is the harmful mode on noisy real data; it is
      appropriate only when the measurements are genuinely sigma-1 clean (for
      example a synthetic oracle). Ignored when `:robust` is not set.
    * `:p_fa` - false-alarm probability for the `:robust` RAIM test (default
      `1.0e-3`); ignored when `:robust` is not set. See `Orbis.GNSS.QC.raim/2`.
    * `:coarse_search` - cold-start convergence-basin widening for a degraded or
      absent position prior (default `nil` = off, exact single solve from
      `:initial_guess`). The crate freezes its elevation mask and weights at the
      seed geometry, so a seed far from the true surface point (the `{0,0,0,0}`
      earth-center default, an antipodal last-known fix) either starves on the
      horizon or is refused by the integrity gates. When set, the solve runs
      once from each of a deterministic golden-spiral lattice of near-surface
      seeds (plus the caller's `:initial_guess`), routes every per-seed candidate
      through the same integrity gates, and selects the best redundant
      (`redundancy >= 1`) converged fix, so no hardcoded prior is needed.
      Accepts `true` (the default seed count), a positive integer seed count, or
      a keyword list `[seeds: n]`. Each seed is one extra crate solve, so the
      cost scales with the seed count; leave it off on the hot path where a good
      prior is known. Composes with the integrity gates (including `:max_pdop`)
      per candidate, but is mutually exclusive with `:robust`: setting both
      returns `{:error, {:incompatible_options, [:coarse_search, :robust]}}`. A
      non-positive or non-integer value returns
      `{:error, {:invalid_option, :coarse_search}}`.

  Regardless of options, a fix that did not converge to a physical receiver
  position is refused rather than returned: one whose geocentric radius is
  outside the plausible band (for example from the earth-center default seed, or
  a wrong-root least-squares fix) gives `{:error, {:implausible_position, radius_m}}`,
  and a converged-flagged fix whose post-fit residual RMS is physically
  implausible gives `{:error, {:no_convergence, rms_m}}`.

  A mixed GPS+Galileo+BeiDou+GLONASS observation set is solved together with one
  receiver clock per GNSS (an inter-system bias is the difference between a
  system's clock and the reference system's), and dilution of precision is
  reported for the combined geometry as well.

  Returns `{:ok, %Orbis.GNSS.Positioning.Solution{}}` or `{:error, reason}`,
  where `reason` is one of `{:too_few_satellites, used, required}` (`required` is
  `3 + n_systems`), `:singular_geometry`, `{:duplicate_observation, sat}`,
  `{:ephemeris_lost, sat}`, `{:ionosphere_unsupported, sat}` (the ionosphere
  correction was requested for a system with no modeled single-frequency
  carrier), `{:degenerate_geometry, reason}` (the geometry is rank-deficient, so
  `reason` is `:rank_deficient`, or exceeds the optional `:max_pdop` ceiling, so
  `reason` is the PDOP), `{:implausible_position, radius_m}` (the fix is outside
  the plausible geocentric-radius band), `{:no_convergence, rms_m}` (a
  converged-flagged fix with physically implausible post-fit residual RMS),
  `{:invalid_option, :max_pdop}`, `{:invalid_option, :coarse_search}`,
  `{:incompatible_options, [:coarse_search, :robust]}`, or, for `:robust`,
  `{:robust_requires_noise_model, :no_weights}` / `{:robust_requires_noise_model, :incomplete_weights}`
  (robust was requested with no realistic noise model, or a weights map that does
  not cover every observed satellite with a positive weight, and no explicit
  unsafe opt-in), `{:invalid_option, :p_fa}`, or `{:fault_unresolved, statistic}`
  (the FDE loop exhausted its iterations with a fault still flagged).
  """
  @spec solve(SP3.t() | Broadcast.t(), [observation()], epoch(), keyword()) ::
          {:ok, Solution.t()} | {:error, term()}
  def solve(source, observations, epoch, opts \\ [])

  def solve(%SP3{handle: handle} = source, observations, epoch, opts)
      when is_list(observations) do
    dispatch(source, :sp3, handle, observations, epoch, opts)
  end

  def solve(%Broadcast{handle: handle} = source, observations, epoch, opts)
      when is_list(observations) do
    dispatch(source, :broadcast, handle, observations, epoch, opts)
  end

  # Opt-in robustness layer over the bare crate solve. When neither :robust nor
  # :coarse_search is set, this is byte-identical to the pre-robustness path: a
  # single run_solve through the same post_process integrity gates. The robust
  # and coarse options are additive Elixir-side wrappers; they never change the
  # bare-path numerics or the gates.
  defp dispatch(source, source_tag, handle, observations, epoch, opts) do
    robust? = Keyword.get(opts, :robust, false)
    coarse = coarse_search_count(Keyword.get(opts, :coarse_search))

    cond do
      coarse == :invalid ->
        {:error, {:invalid_option, :coarse_search}}

      robust? and is_integer(coarse) ->
        # The robust FDE path and the per-seed coarse search are not composed;
        # rejecting the combination is honest rather than silently honoring one.
        {:error, {:incompatible_options, [:coarse_search, :robust]}}

      robust? ->
        robust_solve(source, observations, epoch, opts)

      is_integer(coarse) ->
        coarse_search_solve(source_tag, handle, observations, epoch, opts)

      true ->
        run_solve(source_tag, handle, observations, epoch, opts)
    end
  end

  # :robust routes through the QC leave-one-out FDE loop. A realistic noise model
  # is required: unit-weight FDE treats ordinary code noise as faults and
  # over-excludes, so it is refused unless the caller passes :weights or the
  # explicit :unsafe_unit_weights escape hatch. The cleaned re-solve still flows
  # through run_solve/post_process, so :robust composes with every integrity
  # gate. The FDE exclusion ledger is folded into metadata.fde.
  defp robust_solve(source, observations, epoch, opts) do
    weights = Keyword.get(opts, :weights)
    unsafe? = Keyword.get(opts, :unsafe_unit_weights, false)

    # Validate the robust contract BEFORE any solve, so solve/4 returns a tagged
    # error rather than raising from deep in the FDE/RAIM path. A weights map is
    # a valid noise model only when it carries a positive weight for EVERY
    # observed satellite: an empty or partial map would otherwise fall back to
    # unit weight per satellite inside RAIM, the exact silent-degrade this
    # contract exists to prevent.
    with :ok <- validate_p_fa(Keyword.get(opts, :p_fa)),
         :ok <- validate_unsafe_unit_weights(unsafe?),
         {:ok, raim_weights} <- robust_weights(weights, unsafe?, observations) do
      run_fde(source, observations, epoch, opts, raim_weights)
    end
  end

  defp validate_p_fa(nil), do: :ok
  # Reject a probability so small that `1.0 - p_fa` rounds to exactly 1.0 in
  # f64, which the chi-square inverse rejects (it would raise after a solve).
  defp validate_p_fa(p) when is_number(p) and p > 0.0 and p < 1.0 and 1.0 - p < 1.0, do: :ok
  defp validate_p_fa(_), do: {:error, {:invalid_option, :p_fa}}

  defp validate_unsafe_unit_weights(v) when is_boolean(v), do: :ok
  defp validate_unsafe_unit_weights(_), do: {:error, {:invalid_option, :unsafe_unit_weights}}

  # robust_weights is the single validation point for the FDE weight map: it
  # filters to the observed satellites before RAIM, so RAIM never sees an
  # unvalidated weight and cannot raise out of solve/4. A weight is 1/sigma^2;
  # the upper bound caps sigma at ~1 micrometre (far below any real GNSS sigma)
  # and keeps RAIM's r^2 * weight term finite rather than overflowing.
  @max_robust_weight 1.0e12

  defp robust_weights(weights, _unsafe?, observations) when is_map(weights) do
    sats = Enum.map(observations, fn {sat, _pr} -> sat end)
    present = Enum.filter(sats, &Map.has_key?(weights, &1))

    cond do
      length(present) < length(sats) ->
        {:error, {:robust_requires_noise_model, :incomplete_weights}}

      not Enum.all?(present, fn sat -> valid_robust_weight?(Map.get(weights, sat)) end) ->
        {:error, {:invalid_option, :weights}}

      true ->
        {:ok, Map.take(weights, sats)}
    end
  end

  defp robust_weights(_weights, true, _observations), do: {:ok, :unit}

  defp robust_weights(_weights, false, _observations),
    do: {:error, {:robust_requires_noise_model, :no_weights}}

  defp valid_robust_weight?(w), do: is_number(w) and w > 0.0 and w <= @max_robust_weight

  defp run_fde(source, observations, epoch, opts, raim_weights) do
    # Drop the robust-only options and RAIM-only options that are not part of the
    # public solve/4 surface (:n_systems would otherwise reach RAIM arithmetic
    # unvalidated and raise); RAIM derives :n_systems from the used satellites.
    core_opts =
      Keyword.drop(opts, [
        :robust,
        :weights,
        :unsafe_unit_weights,
        :p_fa,
        :coarse_search,
        :n_systems
      ])

    fde_opts =
      core_opts
      |> Keyword.put(:weights, raim_weights)
      |> maybe_put(:p_fa, Keyword.get(opts, :p_fa))

    case Orbis.GNSS.QC.fde(source, observations, epoch, fde_opts) do
      {:ok, %{solution: solution, excluded: excluded, iterations: iterations}} ->
        metadata = Map.put(solution.metadata, :fde, %{excluded: excluded, iterations: iterations})
        {:ok, %{solution | metadata: metadata}}

      {:error, _reason} = error ->
        error
    end
  end

  defp maybe_put(opts, _key, nil), do: opts
  defp maybe_put(opts, key, value), do: Keyword.put(opts, key, value)

  # Coarse cold-start: solve once per deterministic near-surface seed (plus the
  # caller's own :initial_guess), route EVERY per-seed candidate through the same
  # run_solve/post_process gates the single path uses, and select the best
  # redundant converged fix. The 0.24.0 plausibility/convergence gates drop the
  # never-iterated earth-center seed pass-through (it lands near radius zero), so
  # no residual-RMS workaround is needed; the scorer only ever sees real fixes.
  defp coarse_search_solve(source_tag, handle, observations, epoch, opts) do
    n = coarse_search_count(Keyword.get(opts, :coarse_search))
    core_opts = Keyword.delete(opts, :coarse_search)
    caller_guess = Keyword.get(opts, :initial_guess, @default_initial_guess)
    seeds = [caller_guess | coarse_seeds(n)]

    {candidates, last_error} =
      Enum.reduce(seeds, {[], {:error, :no_coarse_solution}}, fn seed, {acc, last} ->
        case run_solve(
               source_tag,
               handle,
               observations,
               epoch,
               Keyword.put(core_opts, :initial_guess, seed)
             ) do
          {:ok, %Solution{} = sol} -> {[sol | acc], last}
          {:error, _reason} = err -> {acc, err}
        end
      end)

    case select_coarse_candidate(candidates) do
      nil -> last_error
      %Solution{} = sol -> {:ok, sol}
    end
  end

  # Normalize the :coarse_search option to a seed count, or nil when off.
  defp coarse_search_count(nil), do: nil
  defp coarse_search_count(false), do: nil
  defp coarse_search_count(true), do: @default_coarse_seeds
  defp coarse_search_count(n) when is_integer(n) and n >= 1, do: n

  defp coarse_search_count(opts) when is_list(opts) do
    # Only a proper keyword list is the [seeds: n] form; any other list (a
    # charlist, a positional list) is an invalid request, not a default.
    if Keyword.keyword?(opts) do
      coarse_search_count(Keyword.get(opts, :seeds, true))
    else
      :invalid
    end
  end

  # Any other value (0, negative, non-integer) is an invalid request, surfaced
  # as a tagged error by dispatch rather than crashing solve/4.
  defp coarse_search_count(_other), do: :invalid

  # Eligibility uses the 0.24.0 redundancy metadata that post_process populates
  # (redundancy = used_count - (3 + distinct_systems)): keep only converged,
  # redundant (redundancy >= 1) candidates. This is correct for mixed-
  # constellation solves, where two systems estimate 5 states so 5 sats is zero
  # redundancy, unlike a hard-coded min-satellite constant.
  #
  # Among the eligible fits the scorer ranks on satellites-used-first, tie-broken
  # by post-fit residual RMS then GDOP. This is the ratified rule (see the spec
  # ratification note): a pure min-RMS rule prefers a smaller, more biased
  # satellite subset, because dropping a satellite lowers post-fit RMS without
  # lowering true position error; ranking on n_used selects the fix that explains
  # the most observations, closest to truth on the over-determined oracle.
  defp select_coarse_candidate(candidates) do
    candidates
    |> Enum.filter(fn sol -> sol.metadata.converged and sol.metadata.redundancy >= 1 end)
    |> case do
      [] ->
        nil

      eligible ->
        Enum.min_by(eligible, &{-length(&1.used_sats), candidate_rms(&1), candidate_gdop(&1)})
    end
  end

  defp candidate_rms(%Solution{residuals_m: residuals}), do: residual_rms(residuals)

  defp candidate_gdop(%Solution{dop: %{gdop: g}}), do: g
  defp candidate_gdop(%Solution{dop: nil}), do: :infinity

  # Deterministic near-surface seeds on a golden-spiral (Fibonacci) lattice at
  # mean Earth radius, clock bias zero. The lattice depends only on the seed
  # count, never on the answer, so the search carries no hardcoded prior.
  defp coarse_seeds(n) when n >= 1 do
    golden = :math.pi() * (3.0 - :math.sqrt(5.0))

    for i <- 0..(n - 1) do
      z = 1.0 - 2.0 * (i + 0.5) / n
      r = :math.sqrt(max(0.0, 1.0 - z * z))
      theta = golden * i
      x = @mean_earth_radius_m * r * :math.cos(theta)
      y = @mean_earth_radius_m * r * :math.sin(theta)
      zc = @mean_earth_radius_m * z
      {x, y, zc, 0.0}
    end
  end

  defp run_solve(source, handle, observations, epoch, opts) do
    with :ok <- validate_max_pdop(Keyword.get(opts, :max_pdop)),
         {:ok, t_rx_j2000_s} <- j2000_seconds(epoch) do
      sod = Time.second_of_day(epoch)
      doy = Time.day_of_year(epoch)

      apply_iono = Keyword.get(opts, :ionosphere, false)
      apply_tropo = Keyword.get(opts, :troposphere, false)
      alpha = Keyword.get(opts, :klobuchar_alpha, @default_alpha)
      beta = Keyword.get(opts, :klobuchar_beta, @default_beta)
      pressure = Keyword.get(opts, :pressure_hpa, @default_pressure_hpa)
      temperature = Keyword.get(opts, :temperature_k, @default_temperature_k)
      humidity = Keyword.get(opts, :relative_humidity, @default_relative_humidity)
      initial_guess = Keyword.get(opts, :initial_guess, @default_initial_guess)
      with_geodetic = Keyword.get(opts, :with_geodetic, true)

      obs = Enum.map(observations, fn {sat, pr} -> {sat, pr / 1.0} end)

      args = [
        handle,
        obs,
        t_rx_j2000_s,
        sod,
        doy,
        to_tuple4(initial_guess),
        apply_iono,
        apply_tropo,
        to_tuple4(alpha),
        to_tuple4(beta),
        pressure / 1.0,
        temperature / 1.0,
        humidity / 1.0,
        with_geodetic
      ]

      result =
        case source do
          :sp3 -> apply(NIF, :spp_solve, args)
          :broadcast -> apply(NIF, :spp_solve_broadcast, args)
        end

      result |> decode() |> post_process(opts)
    end
  rescue
    e in ErlangError -> {:error, e.original}
  end

  defp validate_max_pdop(nil), do: :ok
  defp validate_max_pdop(value) when is_number(value) and value > 0.0, do: :ok
  defp validate_max_pdop(_value), do: {:error, {:invalid_option, :max_pdop}}

  # Enrich every solution with redundancy metadata, then apply the convergence
  # sanity gate and the optional `:max_pdop` geometry gate. A solution that fails
  # either gate is returned as a tagged error rather than a plausible-looking fix.
  defp post_process({:error, _reason} = error, _opts), do: error

  defp post_process({:ok, solution}, opts) do
    solution = %{
      solution
      | metadata: Map.merge(solution.metadata, redundancy_meta(solution.used_sats))
    }

    # Rank-deficient geometry first (no trustworthy fix exists, and it is what
    # lets a wrong-root mirror land on the plausible shell with zero residuals),
    # then the optional :max_pdop ceiling, then position plausibility (earth-center
    # seed and far-out wrong roots), then the residual-RMS convergence sanity gate.
    with :ok <- check_rank(solution),
         :ok <- check_max_pdop(solution, Keyword.get(opts, :max_pdop)),
         :ok <- check_plausible_position(solution),
         :ok <- check_converged(solution) do
      {:ok, solution}
    end
  end

  # Degrees of freedom for the least-squares fix: one position triple plus one
  # receiver clock per distinct constellation. With `redundancy < 1` the geometry
  # is exactly determined (or under-determined), so the residuals are forced near
  # zero and RAIM cannot test the fix; `raim_checkable?` surfaces that.
  defp redundancy_meta(used_sats) do
    systems = used_sats |> Enum.map(&String.first/1) |> Enum.uniq() |> Enum.sort()
    used_count = length(used_sats)
    redundancy = used_count - (3 + length(systems))

    %{
      used_count: used_count,
      systems: systems,
      redundancy: redundancy,
      raim_checkable?: redundancy >= 1
    }
  end

  defp check_plausible_position(%Solution{position: %{x_m: x, y_m: y, z_m: z}}) do
    radius = :math.sqrt(x * x + y * y + z * z)

    if radius >= @min_plausible_radius_m and radius <= @max_plausible_radius_m do
      :ok
    else
      {:error, {:implausible_position, radius}}
    end
  end

  # The residual-RMS sanity gate only judges fixes the kernel flagged converged;
  # a fix reported as not converged is passed through unchanged so the caller can
  # inspect `metadata.converged` itself, exactly as before this gate existed.
  defp check_converged(%Solution{metadata: %{converged: true}, residuals_m: residuals}) do
    rms = residual_rms(residuals)

    if rms > @max_converged_residual_rms_m do
      {:error, {:no_convergence, rms}}
    else
      :ok
    end
  end

  defp check_converged(%Solution{}), do: :ok

  # A rank-deficient geometry (the crate could not form the DOP cofactor inverse)
  # has no unique trustworthy fix; refuse it by default, regardless of :max_pdop.
  defp check_rank(%Solution{dop: nil}), do: {:error, {:degenerate_geometry, :rank_deficient}}
  defp check_rank(%Solution{}), do: :ok

  defp check_max_pdop(_solution, nil), do: :ok

  defp check_max_pdop(%Solution{dop: %{pdop: pdop}}, max_pdop) when pdop > max_pdop,
    do: {:error, {:degenerate_geometry, pdop}}

  defp check_max_pdop(_solution, _max_pdop), do: :ok

  defp residual_rms([]), do: 0.0

  defp residual_rms(residuals) do
    # `residuals_m` is a list of post-fit residual floats (Solution typespec).
    # Do not swallow unexpected shapes: a non-number would crash here rather than
    # be coerced to 0.0 and silently lower the RMS past the sanity gate.
    {sum_sq, count} =
      Enum.reduce(residuals, {0.0, 0}, fn r, {acc, n} when is_number(r) ->
        {acc + r * r, n + 1}
      end)

    if count == 0, do: 0.0, else: :math.sqrt(sum_sq / count)
  end

  # --- decoding ------------------------------------------------------------

  defp decode(
         {:ok,
          {position, rx_clock_s, geodetic, dop, residuals, used, rejected, metadata,
           system_clocks}}
       ) do
    {:ok,
     %Solution{
       position: position_map(position),
       geodetic: geodetic_map(geodetic),
       rx_clock_s: rx_clock_s,
       system_clocks_s: Map.new(system_clocks),
       dop: dop_map(dop),
       residuals_m: residuals,
       used_sats: used,
       rejected_sats: Enum.map(rejected, fn {sat, reason} -> {sat, reason} end),
       metadata: metadata_map(metadata)
     }}
  end

  defp decode(error), do: map_solve_error(error)

  @doc false
  # Map a raw NIF error tuple onto the public `{:error, reason}` contract.
  # Exposed (undocumented) so the mapping for every advertised reason can be
  # tested directly — including the defensive `:singular_geometry` and
  # `:ephemeris_lost` paths, which the crate does not naturally reach from real
  # SP3 inputs and so cannot be driven through a `solve/4` fixture.
  def map_solve_error({:error, :too_few_satellites, used, required}),
    do: {:error, {:too_few_satellites, used, required}}

  def map_solve_error({:error, :singular_geometry}), do: {:error, :singular_geometry}

  def map_solve_error({:error, :duplicate_observation, sat}),
    do: {:error, {:duplicate_observation, sat}}

  def map_solve_error({:error, :ephemeris_lost, sat}), do: {:error, {:ephemeris_lost, sat}}

  def map_solve_error({:error, :ionosphere_unsupported, sat}),
    do: {:error, {:ionosphere_unsupported, sat}}

  def map_solve_error(other), do: {:error, other}

  defp position_map({x, y, z}), do: %{x_m: x, y_m: y, z_m: z}

  defp geodetic_map(nil), do: nil
  defp geodetic_map({lat, lon, h}), do: %{lat_rad: lat, lon_rad: lon, height_m: h}

  defp dop_map(nil), do: nil

  defp dop_map({gdop, pdop, hdop, vdop, tdop}),
    do: %{gdop: gdop, pdop: pdop, hdop: hdop, vdop: vdop, tdop: tdop}

  defp metadata_map({iterations, converged, status, iono, tropo}) do
    %{
      iterations: iterations,
      converged: converged,
      status: status,
      ionosphere_applied: iono,
      troposphere_applied: tropo
    }
  end

  # --- helpers -------------------------------------------------------------

  # SPP receive time is a continuous f64 second-of-J2000, so sub-second epochs
  # are supported: the whole-second part comes from the integer J2000-second
  # conversion and the microsecond fraction is added on.
  defp j2000_seconds(%NaiveDateTime{} = ndt) do
    {micro, _precision} = ndt.microsecond

    {:ok, seconds} = Time.epoch_to_j2000_seconds(%{ndt | microsecond: {0, 0}})
    {:ok, seconds + micro / 1_000_000.0}
  end

  # A `{{y, m, d}, {h, min, s}}` tuple with a fractional second: split the whole
  # second off for the integer J2000 conversion and carry the fraction, so a
  # sub-second tuple epoch is accepted on the same footing as a NaiveDateTime.
  defp j2000_seconds({{_y, _mo, _d} = date, {hour, minute, second}}) when is_float(second) do
    whole = trunc(second)
    frac = second - whole

    {:ok, seconds} = Time.epoch_to_j2000_seconds({date, {hour, minute, whole}})
    {:ok, seconds + frac}
  end

  defp j2000_seconds(epoch) do
    case Time.epoch_to_j2000_seconds(epoch) do
      {:ok, seconds} -> {:ok, seconds / 1.0}
      {:error, _} = err -> err
    end
  end

  defp to_tuple4({_a, _b, _c, _d} = t), do: t

  defp to_tuple4([a, b, c, d]), do: {a / 1.0, b / 1.0, c / 1.0, d / 1.0}
end
