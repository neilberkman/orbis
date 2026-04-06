defmodule Orbis.Atmosphere do
  @moduledoc """
  Atmospheric density via the NRLMSISE-00 empirical model.

  NRLMSISE-00 computes total mass density and temperature from the surface
  to the lower exosphere (~1000 km). It is the standard model used for
  satellite drag prediction.

  ## Quick example

      position = %{latitude: 0.0, longitude: 0.0, altitude_km: 400.0}
      {:ok, result} = Orbis.Atmosphere.density(position, ~U[2024-06-20 12:00:00Z])
      result.density    # kg/m^3 (e.g., ~3e-12 at ISS altitude)
      result.temperature # K

  ## Solar/geomagnetic indices

  The model depends on space weather conditions. You can supply these via
  the `:f107`, `:f107a`, and `:ap` options. If omitted, moderate defaults
  are used (F10.7 = F10.7A = 150.0, Ap = 4.0).

  For operational use, fetch current indices from NOAA/SWPC.
  """

  @default_f107 150.0
  @default_f107a 150.0
  @default_ap 4.0

  @doc """
  Compute atmospheric density and temperature at a given position and time.

  ## Parameters

    - `position` - geodetic position as `%{latitude: deg, longitude: deg, altitude_km: km}`
    - `datetime` - observation time as `DateTime` or `{{y,m,d},{h,m,s}}` tuple
    - `opts` - keyword list of optional space weather indices:
      - `:f107` - daily 10.7 cm solar radio flux (default #{@default_f107})
      - `:f107a` - 81-day average of F10.7 (default #{@default_f107a})
      - `:ap` - daily geomagnetic Ap index (default #{@default_ap})

  ## Returns

  `{:ok, %{density: kg_per_m3, temperature: kelvin}}` or `{:error, reason}`.

  ## Examples

      # ISS altitude with default solar activity
      pos = %{latitude: 28.5, longitude: -80.6, altitude_km: 408.0}
      {:ok, result} = Orbis.Atmosphere.density(pos, ~U[2024-06-20 12:00:00Z])

      # With explicit solar indices
      {:ok, result} = Orbis.Atmosphere.density(pos, ~U[2024-06-20 12:00:00Z],
        f107: 180.0, f107a: 160.0, ap: 15.0)
  """
  def density(position, datetime, opts \\ [])

  def density(%{latitude: lat, longitude: lon, altitude_km: alt}, %DateTime{} = dt, opts) do
    f107 = Keyword.get(opts, :f107, @default_f107)
    f107a = Keyword.get(opts, :f107a, @default_f107a)
    ap = Keyword.get(opts, :ap, @default_ap)

    {year, doy} = datetime_to_year_doy(dt)
    sec = dt.hour * 3600 + dt.minute * 60 + dt.second + elem(dt.microsecond, 0) / 1_000_000

    call_nif(lat, lon, alt, year, doy, sec, f107, f107a, ap)
  end

  def density(
        %{latitude: lat, longitude: lon, altitude_km: alt},
        {{year, month, day}, {hour, minute, second}},
        opts
      ) do
    f107 = Keyword.get(opts, :f107, @default_f107)
    f107a = Keyword.get(opts, :f107a, @default_f107a)
    ap = Keyword.get(opts, :ap, @default_ap)

    doy = day_of_year(year, month, day)
    sec = hour * 3600 + minute * 60 + second

    call_nif(lat, lon, alt, year, doy, sec * 1.0, f107, f107a, ap)
  end

  defp call_nif(lat, lon, alt, year, doy, sec, f107, f107a, ap) do
    {density, temperature} =
      Orbis.NIF.atmosphere_density(
        lat * 1.0,
        lon * 1.0,
        alt * 1.0,
        year,
        doy,
        sec * 1.0,
        f107 * 1.0,
        f107a * 1.0,
        ap * 1.0
      )

    {:ok, %{density: density, temperature: temperature}}
  rescue
    e -> {:error, Exception.message(e)}
  end

  defp datetime_to_year_doy(%DateTime{} = dt) do
    doy = day_of_year(dt.year, dt.month, dt.day)
    {dt.year, doy}
  end

  defp day_of_year(year, month, day) do
    {:ok, date} = Date.new(year, month, day)
    Date.day_of_year(date)
  end
end
