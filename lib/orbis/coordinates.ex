defmodule Orbis.Coordinates do
  @moduledoc """
  Coordinate frame transformations for satellite state vectors.

  Supports:
  - TEME → GCRS with bit-exact (0 ULP) Skyfield parity
  - GCRS → ITRS (Earth-fixed / ECEF)
  - ITRS → Geodetic (WGS84 lat/lon/alt)
  - Topocentric (azimuth/elevation/range) from a ground station
  """

  @doc """
  Convert a TEME state vector to GCRS.

  Accepts a map with `:position` and `:velocity` tuples (km and km/s),
  and a datetime (either `DateTime` or `{{y,m,d},{h,m,s}}` tuple).

  ## Options

    * `:skyfield_compat` - when `true`, uses Skyfield's exact computation
      path (AU scaling, Kahan triple product, FMA multiply) for 0 ULP
      parity. Default `false` uses direct km arithmetic which is more
      numerically precise.

  Returns a map with GCRS `:position` and `:velocity`.
  """
  def teme_to_gcrs(%{position: {x, y, z}, velocity: {vx, vy, vz}}, datetime, opts \\ []) do
    datetime_tuple = to_nif_datetime(datetime)
    skyfield_compat = Keyword.get(opts, :skyfield_compat, false)

    {{x_gcrs, y_gcrs, z_gcrs}, {vx_gcrs, vy_gcrs, vz_gcrs}} =
      Orbis.NIF.teme_to_gcrs(x, y, z, vx, vy, vz, datetime_tuple, skyfield_compat)

    %{
      position: {x_gcrs, y_gcrs, z_gcrs},
      velocity: {vx_gcrs, vy_gcrs, vz_gcrs}
    }
  end

  @doc """
  Convert a GCRS position to ITRS (Earth-fixed / ECEF).

  Accepts a map with a `:position` tuple (km) and a datetime.
  Pass `skyfield_compat: true` for 0 ULP Skyfield parity.

  Returns `{x, y, z}` in km.
  """
  def gcrs_to_itrs(%{position: {x, y, z}}, datetime, opts \\ []) do
    datetime_tuple = to_nif_datetime(datetime)
    skyfield_compat = Keyword.get(opts, :skyfield_compat, false)
    Orbis.NIF.gcrs_to_itrs(x, y, z, datetime_tuple, skyfield_compat)
  end

  @doc """
  Convert an ITRS/ECEF position to WGS84 geodetic coordinates.

  Accepts a position tuple `{x, y, z}` in km.

  Returns `%{latitude: degrees, longitude: degrees, altitude_km: km}`.
  """
  def to_geodetic({x, y, z}) do
    {lat, lon, alt} = Orbis.NIF.itrs_to_geodetic(x, y, z)
    %Orbis.Geodetic{latitude: lat, longitude: lon, altitude_km: alt}
  end

  @doc """
  Compute topocentric azimuth, elevation, and range from a ground station
  to a satellite given in GCRS.

  ## Parameters

    - `gcrs_state` - map with `:position` tuple (km) in GCRS
    - `datetime` - observation time
    - `station` - `%{latitude: deg, longitude: deg, altitude_m: meters}`

  Returns `%{azimuth: degrees, elevation: degrees, range_km: km}`.
  """
  def to_topocentric(%{position: {x, y, z}}, datetime, station, opts \\ []) do
    datetime_tuple = to_nif_datetime(datetime)
    alt_km = station.altitude_m / 1000.0
    skyfield_compat = Keyword.get(opts, :skyfield_compat, false)

    {az, el, range} =
      Orbis.NIF.gcrs_to_topocentric(
        x,
        y,
        z,
        station.latitude,
        station.longitude,
        alt_km,
        datetime_tuple,
        skyfield_compat
      )

    %Orbis.LookAngle{azimuth: az, elevation: el, range_km: range}
  end

  defp to_nif_datetime({{y, m, d}, {h, min, s}}) do
    {{y, m, d}, {h, min, s, 0}}
  end

  defp to_nif_datetime(%DateTime{} = dt) do
    {{dt.year, dt.month, dt.day}, {dt.hour, dt.minute, dt.second, elem(dt.microsecond, 0)}}
  end
end
