defmodule Orbis.Eclipse do
  @moduledoc """
  Earth shadow (eclipse) prediction for satellites.

  Determines whether a satellite is in full sunlight, penumbra, or umbra
  using a conical shadow model. This is critical for satellite power budgets,
  thermal analysis, and optical visibility.

  ## Shadow Model

  Uses the conical shadow model which computes the penumbra and umbra cones
  cast by Earth given the Sun's position. The satellite's position relative
  to these cones determines its illumination state.

  ## Example

      # Direct usage with pre-computed positions
      status = Orbis.Eclipse.status(satellite_gcrs_position, sun_position_from_earth)
      # => :sunlit | :penumbra | :umbra

      # Convenience: propagate + transform + check in one call
      status = Orbis.Eclipse.check(tle, datetime, ephemeris)
  """

  @r_earth 6_371.0
  @r_sun 696_340.0

  @doc """
  Determine the eclipse status of a satellite.

  ## Parameters

    * `sat_pos` - satellite GCRS position `{x, y, z}` in km
    * `sun_pos` - Sun position relative to Earth `{x, y, z}` in km
      (i.e., the vector from Earth center to the Sun)

  ## Returns

    * `:sunlit` - satellite is in full sunlight
    * `:penumbra` - satellite is partially shadowed
    * `:umbra` - satellite is in full shadow
  """
  @spec status(
          {float(), float(), float()},
          {float(), float(), float()}
        ) :: :sunlit | :penumbra | :umbra
  def status(sat_pos, sun_pos) do
    fraction = shadow_fraction(sat_pos, sun_pos)

    cond do
      fraction >= 1.0 -> :umbra
      fraction > 0.0 -> :penumbra
      true -> :sunlit
    end
  end

  @doc """
  Compute the shadow fraction for a satellite.

  Returns a value from `0.0` (full sunlight) to `1.0` (full umbra).
  Values between 0 and 1 indicate partial shadow (penumbra).

  ## Parameters

    * `sat_pos` - satellite GCRS position `{x, y, z}` in km
    * `sun_pos` - Sun position relative to Earth `{x, y, z}` in km
      (vector from Earth center to the Sun)
  """
  @spec shadow_fraction(
          {float(), float(), float()},
          {float(), float(), float()}
        ) :: float()
  def shadow_fraction({sx, sy, sz}, {sun_x, sun_y, sun_z}) do
    # Sun-Earth distance
    d_sun = :math.sqrt(sun_x * sun_x + sun_y * sun_y + sun_z * sun_z)

    # Unit vector from Earth to Sun
    sun_ux = sun_x / d_sun
    sun_uy = sun_y / d_sun
    sun_uz = sun_z / d_sun

    # Project satellite position onto the Earth-Sun line.
    # Positive = toward the Sun, negative = away from the Sun (shadow side).
    proj = sx * sun_ux + sy * sun_uy + sz * sun_uz

    # If the satellite is on the sunlit side (projection >= 0), it cannot
    # be in Earth's shadow.
    if proj >= 0 do
      0.0
    else
      # Perpendicular distance from the satellite to the Earth-Sun line
      # (the shadow axis).
      perp_x = sx - proj * sun_ux
      perp_y = sy - proj * sun_uy
      perp_z = sz - proj * sun_uz
      rho = :math.sqrt(perp_x * perp_x + perp_y * perp_y + perp_z * perp_z)

      # Distance behind Earth along the shadow axis (positive value).
      dist_behind = -proj

      # Half-angle of the umbra cone (converging cone).
      # The umbra cone vertex is beyond Earth; its half-angle is:
      #   sin(alpha_umbra) = (R_sun - R_earth) / d_sun
      alpha_umbra = :math.asin((@r_sun - @r_earth) / d_sun)

      # Half-angle of the penumbra cone (diverging cone).
      #   sin(alpha_penumbra) = (R_sun + R_earth) / d_sun
      alpha_penumbra = :math.asin((@r_sun + @r_earth) / d_sun)

      # Radius of the umbra cone at the satellite's distance behind Earth.
      # The umbra cone narrows as you go further from Earth. The apex is at
      # distance L_u = R_earth * d_sun / (R_sun - R_earth) from Earth center.
      # At distance `dist_behind` behind Earth along the axis:
      #   r_umbra = R_earth - dist_behind * tan(alpha_umbra)
      r_umbra = @r_earth - dist_behind * :math.tan(alpha_umbra)

      # Radius of the penumbra cone at the satellite's distance.
      # The penumbra cone expands:
      #   r_penumbra = R_earth + dist_behind * tan(alpha_penumbra)
      r_penumbra = @r_earth + dist_behind * :math.tan(alpha_penumbra)

      cond do
        # Beyond the penumbra cone: full sunlight
        rho >= r_penumbra ->
          0.0

        # Inside the umbra cone: full shadow
        r_umbra > 0 and rho <= r_umbra ->
          1.0

        # In the penumbra region: linearly interpolate
        # (This is a good approximation; the exact value depends on the
        # overlap area of the solar and Earth disks as seen from the satellite.)
        r_umbra > 0 ->
          # Satellite is between umbra and penumbra radii
          (r_penumbra - rho) / (r_penumbra - r_umbra)

        # Past the umbra tip (r_umbra <= 0): the umbra cone has fully
        # converged. We are in the antumbra region. Treat as penumbra
        # with shadow fraction decreasing from the axis.
        true ->
          # The penumbra still applies. Fraction based on distance from
          # axis relative to penumbra radius.
          max(0.0, (r_penumbra - rho) / r_penumbra)
      end
    end
  end

  @doc """
  Convenience function: propagate a TLE, transform to GCRS, fetch the Sun
  position from an ephemeris, and return the eclipse status.

  ## Parameters

    * `tle` - parsed `%Orbis.Elements{}` struct
    * `datetime` - `DateTime.t()` observation time
    * `ephemeris` - loaded `%Orbis.Ephemeris{}` handle

  ## Returns

    * `:sunlit`, `:penumbra`, or `:umbra`

  ## Example

      {:ok, tle} = Orbis.parse_tle(line1, line2)
      eph = Orbis.Ephemeris.load("de421.bsp")
      status = Orbis.Eclipse.check(tle, ~U[2024-06-21 12:00:00Z], eph)
  """
  @spec check(Orbis.Elements.t(), DateTime.t(), Orbis.Ephemeris.t()) ::
          {:ok, :sunlit | :penumbra | :umbra} | {:error, String.t()}
  def check(%Orbis.Elements{} = tle, %DateTime{} = datetime, %Orbis.Ephemeris{} = ephemeris) do
    with {:ok, teme} <- Orbis.SGP4.propagate(tle, datetime) do
      gcrs = Orbis.Coordinates.teme_to_gcrs(teme, datetime)
      sun_pos = Orbis.Ephemeris.position(ephemeris, :sun, :earth, datetime)

      {:ok, status(gcrs.position, sun_pos)}
    end
  end
end
