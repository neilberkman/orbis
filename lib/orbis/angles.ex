defmodule Orbis.Angles do
  @moduledoc """
  Angular geometry calculations for satellites.

  Computes angular separations between a satellite and celestial bodies
  (Sun, Moon) using vector geometry in the GCRS frame. Useful for:

  - Solar panel pointing analysis
  - Lunar interference assessment
  - Optical brightness estimation (phase angle)
  - Eclipse geometry (Earth angular radius)

  All positions are expected in km in the GCRS (J2000/ICRF) frame.
  All returned angles are in degrees.

  ## Example

      eph = Orbis.Ephemeris.load("de421.bsp")
      {:ok, tle} = Orbis.parse_tle(line1, line2)
      result = Orbis.Angles.compute(tle, ~U[2024-06-21 12:00:00Z], eph)
      result.sun_angle      # degrees
      result.moon_angle     # degrees
      result.sun_elevation  # degrees (positive = sunlit side)
      result.earth_angle    # degrees (angular radius of Earth)
  """

  @earth_radius_km 6378.137

  @doc """
  Angle between satellite nadir (toward Earth) and the Sun direction.

  The nadir vector points from the satellite toward Earth's center,
  i.e., it is the negation of the satellite's GCRS position.

  ## Parameters

    - `satellite_gcrs_position` - `{x, y, z}` satellite position in GCRS (km)
    - `sun_position_from_earth` - `{x, y, z}` Sun position relative to Earth (km)

  Returns angle in degrees (0 = Sun is directly below satellite toward Earth,
  180 = Sun is directly above/away from Earth).
  """
  @spec sun_angle(
          {float(), float(), float()},
          {float(), float(), float()}
        ) :: float()
  def sun_angle(satellite_gcrs_position, sun_position_from_earth) do
    nadir = negate(satellite_gcrs_position)
    sun_from_sat = subtract(sun_position_from_earth, satellite_gcrs_position)
    angle_between(nadir, sun_from_sat)
  end

  @doc """
  Angle between satellite nadir (toward Earth) and the Moon direction.

  ## Parameters

    - `satellite_gcrs_position` - `{x, y, z}` satellite position in GCRS (km)
    - `moon_position_from_earth` - `{x, y, z}` Moon position relative to Earth (km)

  Returns angle in degrees.
  """
  @spec moon_angle(
          {float(), float(), float()},
          {float(), float(), float()}
        ) :: float()
  def moon_angle(satellite_gcrs_position, moon_position_from_earth) do
    nadir = negate(satellite_gcrs_position)
    moon_from_sat = subtract(moon_position_from_earth, satellite_gcrs_position)
    angle_between(nadir, moon_from_sat)
  end

  @doc """
  Sun elevation above or below the satellite's local horizontal plane.

  The local horizontal plane is perpendicular to the radial (nadir) vector.
  Positive elevation means the Sun is on the sunlit (away from Earth) side;
  negative means it is on the shadow (Earth) side.

  ## Parameters

    - `satellite_gcrs_position` - `{x, y, z}` satellite position in GCRS (km)
    - `sun_position_from_earth` - `{x, y, z}` Sun position relative to Earth (km)

  Returns elevation in degrees (-90 to +90).
  """
  @spec sun_elevation(
          {float(), float(), float()},
          {float(), float(), float()}
        ) :: float()
  def sun_elevation(satellite_gcrs_position, sun_position_from_earth) do
    sun_from_sat = subtract(sun_position_from_earth, satellite_gcrs_position)

    # The "up" direction from the satellite is along the radial (away from Earth),
    # which is simply the satellite position vector (since Earth is at origin).
    zenith = satellite_gcrs_position

    # Elevation = 90 - angle_between(zenith, sun_from_sat)
    # This gives positive when the Sun is above the local horizontal.
    zenith_angle = angle_between(zenith, sun_from_sat)
    90.0 - zenith_angle
  end

  @doc """
  Sun-satellite-observer phase angle.

  The phase angle is the angle at the satellite between the Sun and the
  observer. It determines the illumination geometry for optical brightness
  estimation:

  - 0 deg = full phase (Sun behind observer, satellite fully lit)
  - 180 deg = new phase (Sun behind satellite, satellite in shadow from observer)

  ## Parameters

    - `satellite_gcrs_position` - `{x, y, z}` satellite position in GCRS (km)
    - `sun_position_from_earth` - `{x, y, z}` Sun position relative to Earth (km)
    - `observer_position` - `{x, y, z}` observer position in GCRS (km)

  Returns phase angle in degrees (0 to 180).
  """
  @spec phase_angle(
          {float(), float(), float()},
          {float(), float(), float()},
          {float(), float(), float()}
        ) :: float()
  def phase_angle(satellite_gcrs_position, sun_position_from_earth, observer_position) do
    sun_from_sat = subtract(sun_position_from_earth, satellite_gcrs_position)
    observer_from_sat = subtract(observer_position, satellite_gcrs_position)
    angle_between(sun_from_sat, observer_from_sat)
  end

  @doc """
  Angular radius of the Earth as seen from the satellite.

  This is the half-angle of the cone that just encloses the Earth's disk
  as seen from the satellite's position: `asin(R_earth / |sat_position|)`.

  Useful for eclipse geometry — if the Sun is within this angular radius
  of the anti-nadir direction, the satellite may be in Earth's shadow.

  ## Parameters

    - `satellite_gcrs_position` - `{x, y, z}` satellite position in GCRS (km)

  Returns angular radius in degrees.
  """
  @spec earth_angular_radius({float(), float(), float()}) :: float()
  def earth_angular_radius(satellite_gcrs_position) do
    distance = magnitude(satellite_gcrs_position)
    ratio = @earth_radius_km / distance
    # Clamp to avoid domain errors if satellite is inside Earth (shouldn't happen)
    ratio = min(ratio, 1.0)
    :math.asin(ratio) * 180.0 / :math.pi()
  end

  @doc """
  Compute all standard angles for a satellite at a given time.

  Propagates the TLE, gets Sun and Moon positions from the ephemeris,
  and returns a map of angles.

  ## Parameters

    - `tle` - parsed `%Orbis.Elements{}` struct
    - `datetime` - `DateTime.t()` observation time
    - `ephemeris` - loaded `%Orbis.Ephemeris{}` handle

  ## Returns

      %{
        sun_angle: float(),       # nadir-to-Sun angle in degrees
        moon_angle: float(),      # nadir-to-Moon angle in degrees
        sun_elevation: float(),   # Sun elevation above local horizontal
        earth_angle: float()      # Earth angular radius from satellite
      }
  """
  @spec compute(Orbis.Elements.t(), DateTime.t(), Orbis.Ephemeris.t()) ::
          {:ok, map()} | {:error, String.t()}
  def compute(%Orbis.Elements{} = tle, %DateTime{} = datetime, %Orbis.Ephemeris{} = ephemeris) do
    with {:ok, teme} <- Orbis.SGP4.propagate(tle, datetime) do
      gcrs = Orbis.Coordinates.teme_to_gcrs(teme, datetime)
      sat_pos = gcrs.position

      sun_pos = Orbis.Ephemeris.position(ephemeris, :sun, :earth, datetime)
      moon_pos = Orbis.Ephemeris.position(ephemeris, :moon, :earth, datetime)

      {:ok,
       %{
         sun_angle: sun_angle(sat_pos, sun_pos),
         moon_angle: moon_angle(sat_pos, moon_pos),
         sun_elevation: sun_elevation(sat_pos, sun_pos),
         earth_angle: earth_angular_radius(sat_pos)
       }}
    end
  end

  # ------------------------------------------------------------------
  # Vector arithmetic helpers
  # ------------------------------------------------------------------

  defp subtract({ax, ay, az}, {bx, by, bz}) do
    {ax - bx, ay - by, az - bz}
  end

  defp negate({x, y, z}) do
    {-x, -y, -z}
  end

  defp dot({ax, ay, az}, {bx, by, bz}) do
    ax * bx + ay * by + az * bz
  end

  defp magnitude({x, y, z}) do
    :math.sqrt(x * x + y * y + z * z)
  end

  defp angle_between(a, b) do
    cos_theta = dot(a, b) / (magnitude(a) * magnitude(b))
    # Clamp to [-1, 1] for numerical safety
    cos_theta = max(-1.0, min(1.0, cos_theta))
    :math.acos(cos_theta) * 180.0 / :math.pi()
  end
end
