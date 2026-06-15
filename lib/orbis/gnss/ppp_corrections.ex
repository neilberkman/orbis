defmodule Orbis.GNSS.PPPCorrections do
  @moduledoc """
  Precomputed, state-independent PPP per-range corrections for the static-arc
  float/fixed solve: solid-earth tide, carrier-phase wind-up, and satellite
  antenna PCO/PCV.

  All three depend only on epoch geometry (Sun/Moon direction, satellite
  position, and a fixed reference receiver position), not on the estimated
  receiver coordinates, so they are computed ONCE before the Gauss-Newton loop
  at the seed/approx position and looked up by the row builders and the shared
  `range_corrections_m` chokepoint. The sub-100 m seed-vs-truth offset moves each
  correction by far less than its own measurement noise (tide gradient
  ~1e-7 m/m; wind-up and PCO geometry vary at the milli-arcsec level over 100 m),
  so the state-independent precomputation is exact at solve precision.

  Crate-backed pieces (Sun/Moon ECEF and the IERS DEHANTTIDEINEL kernel) come
  through `Orbis.NIF`. The wind-up and satellite-PCO vector algebra is here in
  Elixir against the crate Sun direction and the SP3 satellite geometry, where it
  is short dot/cross algebra and directly checkable against RTKLIB ppp.c /
  preceph.c.
  """

  alias Orbis.GNSS.Antex
  alias Orbis.GNSS.Core.Constants
  alias Orbis.GNSS.Observables

  @two_pi 2.0 * :math.pi()
  # WGS84 Earth angular velocity (rad/s), for the inertial-velocity term in the
  # nominal yaw (RTKLIB OMGE).
  @omge 7.2921151467e-5

  @type vec3 :: {float(), float(), float()}

  @type t :: %{
          tide: %{optional(NaiveDateTime.t()) => vec3()},
          windup_m: %{optional({String.t(), NaiveDateTime.t()}) => float()},
          sat_pco_ecef: %{optional({String.t(), NaiveDateTime.t()}) => vec3()},
          sat_pcv_m: %{optional({String.t(), NaiveDateTime.t()}) => float()}
        }

  @doc """
  Build the precomputed correction tables for the arc.

  `config` keys (all optional; absent = that correction is off):

    * `:solid_earth_tide` - `true` to compute the per-epoch tide displacement.
    * `:phase_windup` - `true` to compute per-(sat,epoch) wind-up metres.
    * `:satellite_antenna` - `%{antex: %Antex{}, freq1: "G01", freq2: "G02"}` to
      compute per-(sat,epoch) satellite PCO (ECEF vector) and nadir PCV (metres).

  `ref_pos` is the reference receiver ECEF `{x,y,z}` (the solve seed/approx).
  """
  @spec build(Orbis.GNSS.SP3.t(), [map()], vec3(), keyword() | map()) :: t()
  def build(sp3, epochs, ref_pos, config) do
    config = Map.new(config)
    tide? = Map.get(config, :solid_earth_tide, false)
    windup? = Map.get(config, :phase_windup, false)
    sat_ant = Map.get(config, :satellite_antenna)

    sun_moon = if tide? or windup? or sat_ant, do: sun_moon_by_epoch(epochs), else: %{}

    tide = if tide?, do: tide_by_epoch(epochs, ref_pos, sun_moon), else: %{}

    {windup_m, sat_pco_ecef, sat_pcv_m} =
      per_sat_epoch(sp3, epochs, ref_pos, sun_moon, windup?, sat_ant)

    %{tide: tide, windup_m: windup_m, sat_pco_ecef: sat_pco_ecef, sat_pcv_m: sat_pcv_m}
  end

  @doc "Empty correction tables (no precomputed corrections)."
  @spec empty() :: t()
  def empty, do: %{tide: %{}, windup_m: %{}, sat_pco_ecef: %{}, sat_pcv_m: %{}}

  # --- Sun/Moon per epoch (crate NIF) --------------------------------------

  defp sun_moon_by_epoch(epochs) do
    epochs
    |> Enum.map(& &1.epoch)
    |> Enum.uniq()
    |> Map.new(fn epoch ->
      {sun, moon} = Orbis.NIF.sun_moon_ecef(naive_to_tuple(epoch))
      {epoch, %{sun: sun, moon: moon}}
    end)
  end

  # --- Solid-earth tide per epoch (crate NIF, IERS DEHANTTIDEINEL) ----------

  defp tide_by_epoch(epochs, {rx, ry, rz}, sun_moon) do
    epochs
    |> Enum.map(& &1.epoch)
    |> Enum.uniq()
    |> Map.new(fn epoch ->
      %{sun: sun, moon: moon} = Map.fetch!(sun_moon, epoch)
      {{y, mo, d}, {h, mi, s}} = naive_to_ymd_hms(epoch)
      fhr = h + mi / 60.0 + s / 3600.0
      d_tide = Orbis.NIF.solid_earth_tide(rx, ry, rz, y, mo, d, fhr, sun, moon)
      {epoch, d_tide}
    end)
  end

  # --- Per-(sat,epoch) wind-up + satellite PCO/PCV -------------------------

  defp per_sat_epoch(sp3, epochs, ref_pos, sun_moon, windup?, sat_ant)
       when windup? or not is_nil(sat_ant) do
    # Process arcs in chronological order so the wind-up floor-unwrap continuity
    # carries per (sat). `prev_phw` holds the last unwrapped phw cycle per sat.
    epochs
    |> Enum.reduce({%{}, %{}, %{}, %{}}, fn epoch_row, {wmap, pcomap, pcvmap, prev_phw} ->
      epoch = epoch_row.epoch
      %{sun: sun} = Map.fetch!(sun_moon, epoch)

      Enum.reduce(epoch_row.observations, {wmap, pcomap, pcvmap, prev_phw}, fn o,
                                                                               {wm, pco, pcv,
                                                                                pphw} ->
        sat = o.satellite_id

        case Observables.predict(sp3, sat, ref_pos, epoch) do
          {:ok, pred} ->
            {wm2, pphw2} =
              if windup? do
                prev = Map.get(pphw, sat)

                case windup_cycles(pred, ref_pos, sun, prev) do
                  {:ok, phw} ->
                    wm_m = windup_metres(phw, sat_ant_or_default(sat_ant, o))
                    {Map.put(wm, {sat, epoch}, wm_m), Map.put(pphw, sat, phw)}

                  :skip ->
                    {wm, pphw}
                end
              else
                {wm, pphw}
              end

            {pco2, pcv2} =
              if sat_ant do
                case satellite_antenna_correction(pred, sun, sat, epoch, sat_ant) do
                  {:ok, dant_ecef, pcv_m} ->
                    {Map.put(pco, {sat, epoch}, dant_ecef), Map.put(pcv, {sat, epoch}, pcv_m)}

                  :skip ->
                    {pco, pcv}
                end
              else
                {pco, pcv}
              end

            {wm2, pco2, pcv2, pphw2}

          {:error, _} ->
            {wm, pco, pcv, pphw}
        end
      end)
    end)
    |> then(fn {wm, pco, pcv, _pphw} -> {wm, pco, pcv} end)
  end

  defp per_sat_epoch(_sp3, _epochs, _ref_pos, _sun_moon, _windup?, _sat_ant), do: {%{}, %{}, %{}}

  # When wind-up is requested without satellite-antenna config the iono-free
  # wavelengths come from the carriers on the original (raw) obs row.
  defp sat_ant_or_default(nil, o) do
    raw = Map.get(o, :raw, o)
    %{freq1_hz: Map.fetch!(raw, :f1_hz), freq2_hz: Map.fetch!(raw, :f2_hz)}
  end

  defp sat_ant_or_default(sat_ant, _o), do: sat_ant

  # --- Carrier-phase wind-up (demo5 ppp.c model_phw / sat_yaw / yaw_nominal) -

  # Returns {:ok, phw_cycles} with floor()-based cross-epoch continuity unwrap,
  # or :skip on any zero-length normalize (eclipse/collinearity guard).
  defp windup_cycles(pred, {rx, ry, rz}, sun, prev_phw) do
    rs = pred.sat_pos_ecef_m
    vs = pred.sat_velocity_m_s

    # ek = unit(rr - rs): satellite -> receiver (RTKLIB model_phw line 309).
    with {:ok, {exs, eys}} <- sat_yaw(rs, vs, sun),
         {:ok, ek} <- unit(sub({rx, ry, rz}, rs)) do
      # Receiver-fixed unit vectors (ENU). exr = north, eyr = west.
      {n, e, _u} = enu_basis({rx, ry, rz})
      exr = n
      eyr = neg(e)

      eks = cross(ek, eys)
      ekr = cross(ek, eyr)

      ds = sub3(exs, add(scale(ek, dot(ek, exs)), eks))
      dr = sub3(exr, sub(scale(ek, dot(ek, exr)), ekr))

      nds = norm(ds)
      ndr = norm(dr)

      if nds == 0.0 or ndr == 0.0 do
        :skip
      else
        cosp = dot(ds, dr) / nds / ndr
        cosp = max(-1.0, min(1.0, cosp))
        ph = :math.acos(cosp) / @two_pi
        drs = cross(ds, dr)
        ph = if dot(ek, drs) < 0.0, do: -ph, else: ph

        phw =
          case prev_phw do
            nil -> ph
            p -> ph + Float.floor(p - ph + 0.5)
          end

        {:ok, phw}
      end
    else
      _ -> :skip
    end
  end

  # demo5 sat_yaw: builds the satellite-body x/y unit vectors (exs, eys) from the
  # nominal yaw with the noon/midnight singularity guard.
  defp sat_yaw(rs, vs, sun) do
    {sx, sy, _sz} = rs
    {vx, vy, vz} = vs
    # Inertial velocity (RTKLIB ri[3]-=OMGE*ri[1]; ri[4]+=OMGE*ri[0]).
    ri_v = {vx - @omge * sy, vy + @omge * sx, vz}
    n = cross(rs, ri_v)
    p = cross(sun, n)

    with {:ok, es} <- unit(rs),
         {:ok, esun} <- unit(sun),
         {:ok, en} <- unit(n),
         {:ok, ep} <- unit(p) do
      beta = :math.pi() / 2.0 - :math.acos(clamp(dot(esun, en)))
      ee = :math.acos(clamp(dot(es, ep)))
      mu = :math.pi() / 2.0 + if(dot(es, esun) <= 0.0, do: -ee, else: ee)

      mu =
        cond do
          mu < -:math.pi() / 2.0 -> mu + @two_pi
          mu >= :math.pi() / 2.0 -> mu - @two_pi
          true -> mu
        end

      yaw = yaw_nominal(beta, mu)
      ex = cross(en, es)
      cosy = :math.cos(yaw)
      siny = :math.sin(yaw)

      exs = add(scale(en, -siny), scale(ex, cosy))
      eys = add(scale(en, -cosy), scale(ex, -siny))
      {:ok, {exs, eys}}
    else
      _ -> :skip
    end
  end

  defp yaw_nominal(beta, mu) do
    if abs(beta) < 1.0e-12 and abs(mu) < 1.0e-12 do
      :math.pi()
    else
      :math.atan2(-:math.tan(beta), :math.sin(mu)) + :math.pi()
    end
  end

  # Iono-free wind-up in metres: meas_phase -= (c1*lam1 + c2*lam2)*phw, with the
  # IF combination c1=gamma, c2=-(gamma-1) of the per-carrier wavelengths.
  defp windup_metres(phw, ant) do
    f1 = ant_freq(ant, :freq1)
    f2 = ant_freq(ant, :freq2)
    c = Constants.speed_of_light_m_s()
    lam1 = c / f1
    lam2 = c / f2
    gamma = f1 * f1 / (f1 * f1 - f2 * f2)
    (gamma * lam1 - (gamma - 1.0) * lam2) * phw
  end

  defp ant_freq(%{freq1_hz: f1, freq2_hz: f2}, which), do: if(which == :freq1, do: f1, else: f2)

  defp ant_freq(%{freq1: f1, freq2: f2}, which) do
    label = if which == :freq1, do: f1, else: f2
    band_hz(label)
  end

  defp band_hz("G01"), do: Constants.gps_l1_hz()
  defp band_hz("G02"), do: Constants.gps_l2_hz()

  # --- Satellite antenna PCO/PCV (preceph.c satantoff + demo5 satantpcv) ----

  # Returns {:ok, dant_pco_ecef, pcv_if_m} or :skip on eclipse/collinearity.
  defp satellite_antenna_correction(pred, sun, sat, epoch, %{antex: antex, freq1: f1l, freq2: f2l}) do
    rs = pred.sat_pos_ecef_m
    prn = sat

    with %Antex.Antenna{} = ant <- Antex.satellite_antenna(antex, prn, epoch),
         {:ok, ez} <- unit(neg(rs)),
         {:ok, es} <- unit(sub(sun, rs)),
         {:ok, ey} <- unit(cross(ez, es)) do
      ex = cross(ey, ez)

      off1 = pco_or_nil(ant, f1l)
      off2 = pco_or_nil(ant, f2l)

      if off1 == nil or off2 == nil do
        :skip
      else
        f1 = band_hz(f1l)
        f2 = band_hz(f2l)
        gamma = f1 * f1 / (f1 * f1 - f2 * f2)

        # PCO body-frame NEU -> ECEF projection per freq, IF-combined.
        dant1 = body_to_ecef(off1, ex, ey, ez)
        dant2 = body_to_ecef(off2, ex, ey, ez)
        dant_ecef = sub(scale(dant1, gamma), scale(dant2, gamma - 1.0))

        # Nadir PCV (1-deg sat grid), IF-combined scalar.
        pcv_m = nadir_pcv_if(ant, pred, f1l, f2l, gamma)

        {:ok, dant_ecef, pcv_m}
      end
    else
      _ -> :skip
    end
  end

  # ANTEX satellite PCO is stored NEU in the SATELLITE BODY frame (x=NEU[0],
  # y=NEU[1], z=NEU[2]) consumed over (ex,ey,ez): dant = off0*ex + off1*ey +
  # off2*ez (preceph.c satantoff).
  defp body_to_ecef({o0, o1, o2}, ex, ey, ez) do
    add(add(scale(ex, o0), scale(ey, o1)), scale(ez, o2))
  end

  defp pco_or_nil(ant, label) do
    Antex.pco(ant, label)
  rescue
    _ -> nil
  end

  defp nadir_pcv_if(ant, pred, f1l, f2l, gamma) do
    rs = pred.sat_pos_ecef_m
    # eu = unit(rr - rs) but we only have los (rx->sat); nadir uses eu·ez with
    # eu = unit(rs->rx) = -los_unit, ez = unit(-rs). RTKLIB: eu=unit(rr-rs),
    # ez=unit(-rs), nadir=acos(eu·ez).
    eu =
      case unit(neg(pred.los_unit)) do
        {:ok, v} -> v
        :skip -> nil
      end

    ez =
      case unit(neg(rs)) do
        {:ok, v} -> v
        :skip -> nil
      end

    if eu == nil or ez == nil do
      0.0
    else
      cosa = clamp(dot(eu, ez))
      nadir_deg = :math.acos(cosa) * 180.0 / :math.pi()
      p1 = pcv_or_zero(ant, f1l, nadir_deg)
      p2 = pcv_or_zero(ant, f2l, nadir_deg)
      gamma * p1 - (gamma - 1.0) * p2
    end
  end

  defp pcv_or_zero(ant, label, nadir_deg) do
    Antex.pcv(ant, label, nadir_deg)
  rescue
    _ -> 0.0
  end

  # --- vector helpers ------------------------------------------------------

  defp enu_basis({x, y, z}) do
    {:ok, u} = up_unit({x, y, z})
    e = east_unit(u)
    {:ok, n} = north_unit(e, u)
    {n, e, u}
  end

  defp up_unit(v) do
    case norm(v) do
      n when n > 0.0 -> {:ok, scale(v, 1.0 / n)}
      _ -> {:ok, {0.0, 0.0, 1.0}}
    end
  end

  defp east_unit(up) do
    c = cross({0.0, 0.0, 1.0}, up)

    if c == {0.0, 0.0, 0.0} do
      {1.0, 0.0, 0.0}
    else
      case norm(c) do
        n when n > 0.0 -> scale(c, 1.0 / n)
        _ -> {1.0, 0.0, 0.0}
      end
    end
  end

  defp north_unit(east, up) do
    case unit(cross(east, up)) do
      {:ok, n} -> {:ok, n}
      :skip -> {:ok, {0.0, 0.0, 1.0}}
    end
  end

  defp unit(v) do
    case norm(v) do
      n when n > 0.0 -> {:ok, scale(v, 1.0 / n)}
      _ -> :skip
    end
  end

  defp clamp(x), do: max(-1.0, min(1.0, x))

  defp add({ax, ay, az}, {bx, by, bz}), do: {ax + bx, ay + by, az + bz}
  defp sub({ax, ay, az}, {bx, by, bz}), do: {ax - bx, ay - by, az - bz}
  # sub3 is sub with explicit naming to mirror the RTKLIB ds/dr assembly.
  defp sub3(a, b), do: sub(a, b)
  defp neg({x, y, z}), do: {-x, -y, -z}
  defp scale({x, y, z}, s), do: {x * s, y * s, z * s}
  defp dot({ax, ay, az}, {bx, by, bz}), do: ax * bx + ay * by + az * bz
  defp norm({x, y, z}), do: :math.sqrt(x * x + y * y + z * z)

  defp cross({ax, ay, az}, {bx, by, bz}) do
    {ay * bz - az * by, az * bx - ax * bz, ax * by - ay * bx}
  end

  defp naive_to_tuple(%NaiveDateTime{} = ndt) do
    {{ndt.year, ndt.month, ndt.day}, {ndt.hour, ndt.minute, ndt.second, elem(ndt.microsecond, 0)}}
  end

  defp naive_to_ymd_hms(%NaiveDateTime{} = ndt) do
    s = ndt.second + elem(ndt.microsecond, 0) / 1_000_000.0
    {{ndt.year, ndt.month, ndt.day}, {ndt.hour, ndt.minute, s}}
  end
end
