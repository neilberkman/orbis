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
            iterations: non_neg_integer(),
            converged: boolean(),
            status: atom(),
            ionosphere_applied: boolean(),
            troposphere_applied: boolean(),
            used_count: non_neg_integer(),
            systems: [String.t()],
            redundancy: integer(),
            raim_checkable?: boolean()
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
  converged-flagged fix with physically implausible post-fit residual RMS), or
  `{:invalid_option, :max_pdop}`.
  """
  @spec solve(SP3.t() | Broadcast.t(), [observation()], epoch(), keyword()) ::
          {:ok, Solution.t()} | {:error, term()}
  def solve(source, observations, epoch, opts \\ [])

  def solve(%SP3{handle: handle}, observations, epoch, opts) when is_list(observations) do
    run_solve(:sp3, handle, observations, epoch, opts)
  end

  def solve(%Broadcast{handle: handle}, observations, epoch, opts) when is_list(observations) do
    run_solve(:broadcast, handle, observations, epoch, opts)
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
