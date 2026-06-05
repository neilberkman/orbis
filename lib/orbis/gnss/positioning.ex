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
    `:low_elevation`). `metadata` reports solver iterations, convergence, and
    which corrections were applied.
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
            troposphere_applied: boolean()
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

  A mixed GPS+Galileo+BeiDou+GLONASS observation set is solved together with one
  receiver clock per GNSS (an inter-system bias is the difference between a
  system's clock and the reference system's), and dilution of precision is
  reported for the combined geometry as well.

  Returns `{:ok, %Orbis.GNSS.Positioning.Solution{}}` or `{:error, reason}`,
  where `reason` is one of `{:too_few_satellites, used, required}` (`required` is
  `3 + n_systems`), `:singular_geometry`, `{:duplicate_observation, sat}`,
  `{:ephemeris_lost, sat}`, or `{:ionosphere_unsupported, sat}` (the ionosphere
  correction was requested for a system with no modeled single-frequency carrier).
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
    with {:ok, t_rx_j2000_s} <- j2000_seconds(epoch) do
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

      decode(result)
    end
  rescue
    e in ErlangError -> {:error, e.original}
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

    case Time.epoch_to_j2000_seconds(%{ndt | microsecond: {0, 0}}) do
      {:ok, seconds} -> {:ok, seconds + micro / 1_000_000.0}
      {:error, _} = err -> err
    end
  end

  # A `{{y, m, d}, {h, min, s}}` tuple with a fractional second: split the whole
  # second off for the integer J2000 conversion and carry the fraction, so a
  # sub-second tuple epoch is accepted on the same footing as a NaiveDateTime.
  defp j2000_seconds({{_y, _mo, _d} = date, {hour, minute, second}}) when is_float(second) do
    whole = trunc(second)
    frac = second - whole

    case Time.epoch_to_j2000_seconds({date, {hour, minute, whole}}) do
      {:ok, seconds} -> {:ok, seconds + frac}
      {:error, _} = err -> err
    end
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
