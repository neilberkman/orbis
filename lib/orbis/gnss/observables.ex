defmodule Orbis.GNSS.Observables do
  @moduledoc """
  Predict the GNSS observables a receiver at a known ECEF position would see for
  a satellite, from a precise (SP3) ephemeris source.

  This is the forward model behind the question "is this measurement physically
  plausible?": given a receiver position, a satellite, and a receive epoch, it
  computes the geometric range, the line-of-sight range rate, the L1 Doppler,
  the topocentric azimuth/elevation, the satellite clock offset, and the signal
  transmit time. It reads satellite states only from `Orbis.GNSS.SP3.position/3` and
  uses standard textbook GNSS geometry; it never solves the inverse (positioning)
  problem.

  ## Algorithm (standard GNSS geometry)

  * **Light-time / transmit-time correction.** The signal seen at the receive
    epoch `t_rx` left the satellite earlier, at
    `t_tx = t_rx - |r_sat(t_tx) - r_rx| / c`. This is solved by fixed-point
    iteration starting from `t_tx = t_rx`; a couple of iterations converge to
    sub-millimetre level for a coarse receiver position. The satellite state is
    evaluated at the fractional epoch `t_tx` (the SP3 spline is sampled at
    sub-second precision).

  * **Sagnac / Earth-rotation correction.** During the travel time `tau` the
    Earth-fixed (ECEF) frame rotates by `omega_e * tau`. The satellite position
    computed in the ECEF frame at `t_tx` is rotated about the Z axis by
    `Rz(omega_e * tau)` into the receive-epoch ECEF frame before differencing,
    with `omega_e = 7.2921151467e-5 rad/s`. This is the Sagnac (Earth-rotation)
    correction.

  * **Geometric range** is `|r_sat_rot - r_rx|` in metres, and the
    line-of-sight unit vector points from the receiver to the satellite.

  * **Range rate.** The satellite velocity at `t_tx` is obtained by central
    finite difference of `Orbis.GNSS.SP3.position/3` (+/- 0.5 s). For a static
    receiver (`v_rx = 0`) the range rate is the LOS projection
    `los . (v_sat - v_rx)`, which equals `d(range)/dt`.

  * **Doppler (IS-GPS-200 L1 carrier).** `doppler_hz = -range_rate * f / c`
    with the L1 carrier `f = 1575.42 MHz` and `c = 299792458 m/s`.

  ## Sign conventions

  `range_rate_m_s` is the time derivative of the geometric range: it is
  **negative when the satellite is approaching** (range decreasing) and positive
  when receding. The Doppler shift is the negative of the (scaled) range rate, so
  an **approaching satellite gives a positive Doppler** and a receding satellite
  a negative one.

  ## Result map

      %{
        geometric_range_m: float(),    # metres
        range_rate_m_s:    float(),    # d(range)/dt; negative = approaching
        doppler_hz:        float(),    # = -range_rate * carrier / c; + = approaching
        sat_clock_s:       float() | nil,  # SP3 clock offset at transmit time
        elevation_deg:     float(),    # topocentric elevation
        azimuth_deg:       float(),    # topocentric azimuth, [0, 360)
        transmit_time:     NaiveDateTime.t(),  # t_tx
        los_unit:          {float(), float(), float()}  # receiver -> satellite, ECEF unit
      }
  """

  alias Orbis.GNSS.Core.Constants
  alias Orbis.GNSS.Core.Types
  alias Orbis.GNSS.SP3

  # Central finite-difference half-step for the satellite velocity, seconds.
  @fd_half_s 0.5

  @type vec3 :: {float(), float(), float()}

  @type observables :: %{
          geometric_range_m: float(),
          range_rate_m_s: float(),
          doppler_hz: float(),
          sat_clock_s: float() | nil,
          elevation_deg: float(),
          azimuth_deg: float(),
          transmit_time: NaiveDateTime.t(),
          los_unit: vec3()
        }

  @doc """
  Predict the observables for `satellite_id` seen from `receiver_ecef` at `epoch`.

  `receiver_ecef` is the static receiver position in ITRF/ECEF metres, given as
  `{x_m, y_m, z_m}` or `%{x_m: _, y_m: _, z_m: _}`. `epoch` is the receive epoch,
  a `NaiveDateTime` (interpreted in the SP3 file's own time scale).

  ## Options

    * `:carrier_hz` - carrier frequency for the Doppler, default the L1 carrier
      `1575.42 MHz`.
    * `:light_time` - apply the light-time / transmit-time correction, default
      `true`. When `false`, the satellite is evaluated at `epoch`.
    * `:sagnac` - apply the Sagnac / Earth-rotation correction, default `true`.

  Returns `{:ok, observables}`, `{:error, :invalid_receiver}` for a malformed
  receiver position, or propagates any `Orbis.GNSS.SP3.position/3` error (e.g. an
  unknown satellite or a malformed satellite token) verbatim as
  `{:error, reason}`. Never raises.

  Note that `Orbis.GNSS.SP3.position/3` extrapolates its spline rather than reporting
  a coverage error, so an epoch outside the file's span does not yield a tagged
  error here; it produces an obviously non-physical geometry instead.
  """
  @spec predict(SP3.t(), String.t(), vec3() | map(), NaiveDateTime.t(), keyword()) ::
          {:ok, observables()} | {:error, term()}
  def predict(%SP3{} = sp3, satellite_id, receiver_ecef, %NaiveDateTime{} = epoch, opts \\ [])
      when is_binary(satellite_id) do
    carrier_hz = Keyword.get(opts, :carrier_hz, Constants.gps_l1_hz())
    light_time? = Keyword.get(opts, :light_time, true)
    sagnac? = Keyword.get(opts, :sagnac, true)

    with {:ok, {rx, ry, rz}} <- Types.normalize_ecef(receiver_ecef),
         {:ok, t_tx, tau, state, r_sat_rot} <-
           solve_transmit_time(sp3, satellite_id, {rx, ry, rz}, epoch, light_time?, sagnac?),
         {:ok, v_sat} <- sat_velocity(sp3, satellite_id, t_tx) do
      {sxr, syr, szr} = r_sat_rot
      dx = sxr - rx
      dy = syr - ry
      dz = szr - rz
      range = :math.sqrt(dx * dx + dy * dy + dz * dz)

      los = {dx / range, dy / range, dz / range}
      {lx, ly, lz} = los

      # The velocity must live in the same (receive-epoch) ECEF frame as the
      # Sagnac-rotated position, so apply the same Rz(omega_e * tau) rotation to
      # the finite-difference velocity before projecting onto the line of sight.
      {vx, vy, vz} = sagnac_rotate(v_sat, tau, sagnac?)

      # Range rate = d(range)/dt = LOS . (v_sat - v_rx), v_rx = 0 (static).
      # Negative when the satellite is approaching (range decreasing).
      range_rate = lx * vx + ly * vy + lz * vz

      # Doppler: approaching (range decreasing) -> positive shift.
      doppler_hz = -range_rate * carrier_hz / Constants.speed_of_light_m_s()

      {elevation_deg, azimuth_deg} = topocentric(rx, ry, rz, dx, dy, dz, range)

      {:ok,
       %{
         geometric_range_m: range,
         range_rate_m_s: range_rate,
         doppler_hz: doppler_hz,
         sat_clock_s: state.clock_s,
         elevation_deg: elevation_deg,
         azimuth_deg: azimuth_deg,
         transmit_time: t_tx,
         los_unit: los
       }}
    end
  end

  @doc """
  Predict observables for every satellite in the product, seen from `receiver_ecef`.

  Returns a map `satellite_id => {:ok, observables} | {:error, reason}`, so one
  satellite failing (e.g. no estimate at this epoch) does not sink the batch.
  Options are the same as `predict/5`.
  """
  @spec predict_all(SP3.t(), vec3() | map(), NaiveDateTime.t(), keyword()) ::
          %{optional(String.t()) => {:ok, observables()} | {:error, term()}}
  def predict_all(%SP3{} = sp3, receiver_ecef, %NaiveDateTime{} = epoch, opts \\ []) do
    sp3
    |> SP3.satellite_ids()
    |> Map.new(fn sat_id -> {sat_id, predict(sp3, sat_id, receiver_ecef, epoch, opts)} end)
  end

  # --- light-time fixed point ----------------------------------------------

  defp solve_transmit_time(sp3, sat_id, _receiver, epoch, false, sagnac?) do
    # No light-time correction: evaluate the satellite at the receive epoch.
    # tau = 0, so the Sagnac rotation (omega_e * tau) is identity even when
    # :sagnac is true.
    with {:ok, state} <- SP3.position(sp3, sat_id, epoch) do
      r_sat_rot = sagnac_rotate({state.x_m, state.y_m, state.z_m}, 0.0, sagnac?)
      {:ok, epoch, 0.0, state, r_sat_rot}
    end
  end

  defp solve_transmit_time(sp3, sat_id, receiver, epoch, true, sagnac?) do
    iterate_transmit_time(sp3, sat_id, receiver, epoch, sagnac?, 0.0, 3)
  end

  # Fixed-point iteration on tau = range / c. Starting from tau = 0 (t_tx =
  # t_rx), 2-3 iterations converge to sub-millimetre for a coarse receiver.
  defp iterate_transmit_time(sp3, sat_id, receiver, epoch, sagnac?, tau, iters_left) do
    t_tx = NaiveDateTime.add(epoch, -round(tau * 1_000_000), :microsecond)

    with {:ok, state} <- SP3.position(sp3, sat_id, t_tx) do
      r_sat_rot = sagnac_rotate({state.x_m, state.y_m, state.z_m}, tau, sagnac?)
      {sxr, syr, szr} = r_sat_rot
      {rx, ry, rz} = receiver
      dx = sxr - rx
      dy = syr - ry
      dz = szr - rz
      range = :math.sqrt(dx * dx + dy * dy + dz * dz)
      new_tau = range / Constants.speed_of_light_m_s()

      if iters_left <= 1 do
        # Recompute the satellite state and rotation at the converged tau so the
        # returned position, clock, and t_tx are all consistent with new_tau.
        finalize_transmit_time(sp3, sat_id, receiver, epoch, sagnac?, new_tau)
      else
        iterate_transmit_time(sp3, sat_id, receiver, epoch, sagnac?, new_tau, iters_left - 1)
      end
    end
  end

  defp finalize_transmit_time(sp3, sat_id, _receiver, epoch, sagnac?, tau) do
    t_tx = NaiveDateTime.add(epoch, -round(tau * 1_000_000), :microsecond)

    with {:ok, state} <- SP3.position(sp3, sat_id, t_tx) do
      r_sat_rot = sagnac_rotate({state.x_m, state.y_m, state.z_m}, tau, sagnac?)
      {:ok, t_tx, tau, state, r_sat_rot}
    end
  end

  # --- Sagnac / Earth-rotation rotation ------------------------------------

  defp sagnac_rotate(pos, _tau, false), do: pos

  defp sagnac_rotate({x, y, z}, tau, true) do
    # Rotate the ECEF-at-t_tx position into the receive-epoch ECEF frame by
    # Rz(omega_e * tau): the Earth rotates +omega_e over the travel time tau, so
    # the satellite expressed in the t_rx frame is the t_tx vector rotated by
    # +theta about Z.
    theta = Constants.earth_rotation_rate_rad_s() * tau
    c = :math.cos(theta)
    s = :math.sin(theta)
    {c * x + s * y, -s * x + c * y, z}
  end

  # --- satellite velocity (central finite difference) ----------------------

  defp sat_velocity(sp3, sat_id, t_tx) do
    half_us = round(@fd_half_s * 1_000_000)
    t_plus = NaiveDateTime.add(t_tx, half_us, :microsecond)
    t_minus = NaiveDateTime.add(t_tx, -half_us, :microsecond)

    with {:ok, sp} <- SP3.position(sp3, sat_id, t_plus),
         {:ok, sm} <- SP3.position(sp3, sat_id, t_minus) do
      denom = 2.0 * @fd_half_s

      {:ok, {(sp.x_m - sm.x_m) / denom, (sp.y_m - sm.y_m) / denom, (sp.z_m - sm.z_m) / denom}}
    end
  end

  # --- topocentric ENU az/el -----------------------------------------------

  defp topocentric(rx, ry, rz, dx, dy, dz, range) do
    # Receiver geodetic latitude/longitude give the ENU basis. Only lat/lon are
    # needed for the rotation (not height), so divide receiver metres to km for
    # the WGS84 inverse.
    geo = Orbis.Coordinates.to_geodetic({rx / 1000.0, ry / 1000.0, rz / 1000.0})
    lat = geo.latitude * :math.pi() / 180.0
    lon = geo.longitude * :math.pi() / 180.0

    sl = :math.sin(lat)
    cl = :math.cos(lat)
    so = :math.sin(lon)
    co = :math.cos(lon)

    # ECEF line-of-sight -> local ENU (standard topocentric rotation).
    e = -so * dx + co * dy
    n = -sl * co * dx - sl * so * dy + cl * dz
    u = cl * co * dx + cl * so * dy + sl * dz

    azimuth_deg = :math.atan2(e, n) * 180.0 / :math.pi()
    azimuth_deg = if azimuth_deg < 0.0, do: azimuth_deg + 360.0, else: azimuth_deg

    elevation_deg = :math.asin(u / range) * 180.0 / :math.pi()

    {elevation_deg, azimuth_deg}
  end
end
