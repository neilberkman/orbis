defmodule Orbis do
  @moduledoc """
  Satellite toolkit for Elixir. SGP4 orbit propagation, coordinate
  transformations, and ground station pass prediction.
  """

  @doc """
  Parse a Two-Line Element set.

  Returns `{:ok, %Orbis.Elements{}}` or `{:error, reason}`.

  ## Examples

      iex> {:ok, el} = Orbis.parse_tle(
      ...>   "1 25544U 98067A   18184.80969102  .00001614  00000-0  31745-4 0  9993",
      ...>   "2 25544  51.6414 295.8524 0003435 262.6267 204.2868 15.54005638121106"
      ...> )
      iex> el.catalog_number
      "25544"

  """
  defdelegate parse_tle(line1, line2), to: Orbis.Format.TLE, as: :parse

  @doc """
  Propagate orbital elements to a specific datetime, returning TEME position and velocity.

  Returns `{:ok, %Orbis.TemeState{}}` with position in km and velocity in km/s,
  or `{:error, reason}`.

  ## Examples

      iex> {:ok, el} = Orbis.parse_tle(
      ...>   "1 25544U 98067A   18184.80969102  .00001614  00000-0  31745-4 0  9993",
      ...>   "2 25544  51.6414 295.8524 0003435 262.6267 204.2868 15.54005638121106"
      ...> )
      iex> {:ok, teme} = Orbis.propagate(el, ~U[2018-07-04 00:00:00Z])
      iex> {x, _y, _z} = teme.position
      iex> x > 3000 and x < 4000
      true

  """
  defdelegate propagate(tle, datetime), to: Orbis.SGP4

  @doc """
  Predict visible passes of a satellite over a ground station.

  See `Orbis.Passes.predict/5` for full documentation.
  """
  defdelegate predict_passes(tle, ground_station, start_time, end_time, opts \\ []),
    to: Orbis.Passes,
    as: :predict

  @doc """
  Compute the angle between satellite nadir and the Sun direction.

  See `Orbis.Angles.sun_angle/2` for details.
  """
  defdelegate sun_angle(satellite_gcrs_position, sun_position_from_earth),
    to: Orbis.Angles

  @doc """
  Compute the angle between satellite nadir and the Moon direction.

  See `Orbis.Angles.moon_angle/2` for details.
  """
  defdelegate moon_angle(satellite_gcrs_position, moon_position_from_earth),
    to: Orbis.Angles

  @doc """
  Convert a TEME state vector to GCRS (Geocentric Celestial Reference System).

  Pass `skyfield_compat: true` for IEEE 754 bit-identical output to
  Python Skyfield (0 ULP oracle test). Default uses direct km arithmetic
  which is more numerically precise.

  ## Example

      gcrs = Orbis.teme_to_gcrs(teme, datetime)
      gcrs = Orbis.teme_to_gcrs(teme, datetime, skyfield_compat: true)
  """
  def teme_to_gcrs(teme_state, datetime, opts \\ []) do
    Orbis.Coordinates.teme_to_gcrs(teme_state, datetime, opts)
  end

  @doc """
  Compute geodetic coordinates (lat/lon/alt) for a satellite at a given time.

  Propagates the TLE, transforms TEME -> GCRS -> ITRS, and converts to WGS84.

  Returns `{:ok, %{latitude: deg, longitude: deg, altitude_km: km}}`.

  ## Example

      {:ok, tle} = Orbis.parse_tle(line1, line2)
      {:ok, geo} = Orbis.geodetic(tle, datetime)
      geo.latitude  # => 51.23
  """
  def geodetic(%Orbis.Elements{} = tle, %DateTime{} = datetime) do
    with {:ok, teme} <- Orbis.SGP4.propagate(tle, datetime) do
      gcrs = Orbis.Coordinates.teme_to_gcrs(teme, datetime)
      itrs = Orbis.Coordinates.gcrs_to_itrs(gcrs, datetime)
      {:ok, Orbis.Coordinates.to_geodetic(itrs)}
    end
  end

  @doc """
  Check whether a satellite is in Earth's shadow (eclipse) at a given time.

  Propagates the TLE, transforms to GCRS, fetches the Sun position from the
  ephemeris, and returns the eclipse status.

  Returns `{:ok, :sunlit | :penumbra | :umbra}` or `{:error, reason}`.

  ## Example

      eph = Orbis.Ephemeris.load("de421.bsp")
      {:ok, status} = Orbis.eclipse(tle, datetime, eph)
  """
  defdelegate eclipse(tle, datetime, ephemeris), to: Orbis.Eclipse, as: :check

  @doc """
  Compute the look angle (azimuth/elevation/range) from a ground station
  to a satellite at a given time.

  The station is a map: `%{latitude: deg, longitude: deg, altitude_m: meters}`.

  Returns `{:ok, %{azimuth: deg, elevation: deg, range_km: km}}`.

  ## Example

      station = %{latitude: 40.0, longitude: -74.0, altitude_m: 0.0}
      {:ok, look} = Orbis.look_angle(tle, datetime, station)
      look.elevation  # => 25.7
  """
  def look_angle(%Orbis.Elements{} = tle, %DateTime{} = datetime, station) do
    with {:ok, teme} <- Orbis.SGP4.propagate(tle, datetime) do
      gcrs = Orbis.Coordinates.teme_to_gcrs(teme, datetime)
      {:ok, Orbis.Coordinates.to_topocentric(gcrs, datetime, station)}
    end
  end

  @doc """
  Compute Doppler shift for a satellite-ground link.

  Propagates the TLE, transforms to GCRS, and computes the range rate and
  Doppler shift at the given carrier frequency.

  The station is a map: `%{latitude: deg, longitude: deg, altitude_m: meters}`.

  Returns `{:ok, %{range_rate_km_s: float, doppler_hz: float, doppler_ratio: float}}`.

  ## Example

      station = %{latitude: 40.0, longitude: -74.0, altitude_m: 0.0}
      {:ok, d} = Orbis.doppler(tle, datetime, station, 437.0e6)
      d.doppler_hz  # => ~10_000.0
  """
  def doppler(%Orbis.Elements{} = tle, %DateTime{} = datetime, station, frequency_hz) do
    with {:ok, teme} <- Orbis.SGP4.propagate(tle, datetime) do
      gcrs = Orbis.Coordinates.teme_to_gcrs(teme, datetime)
      {:ok, Orbis.Doppler.shift(gcrs, datetime, station, frequency_hz)}
    end
  end
end
