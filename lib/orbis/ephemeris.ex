defmodule Orbis.Ephemeris do
  @moduledoc """
  JPL SPK/BSP ephemeris file reader.

  Computes positions of solar system bodies (Sun, Moon, planets) using
  JPL Development Ephemeris files (DE421, DE440, etc.) in SPK/BSP format.

  These files contain Chebyshev polynomial coefficients that are evaluated
  to produce high-precision positions in km.

  ## Example

      eph = Orbis.Ephemeris.load("de421.bsp")
      {x, y, z} = Orbis.Ephemeris.position(eph, :sun, :earth, ~U[2024-01-01 12:00:00Z])

  ## Body Names

  Supported body atoms: `:ssb` (Solar System Barycenter), `:mercury`,
  `:venus`, `:earth_moon_barycenter` (or `:emb`), `:mars`, `:jupiter`,
  `:saturn`, `:uranus`, `:neptune`, `:pluto`, `:sun`, `:moon`, `:earth`.
  """

  defstruct [:path]

  @type t :: %__MODULE__{path: String.t()}

  @body_atoms ~w(
    ssb solar_system_barycenter
    mercury mercury_barycenter
    venus venus_barycenter
    earth_moon_barycenter emb
    mars mars_barycenter
    jupiter jupiter_barycenter
    saturn saturn_barycenter
    uranus uranus_barycenter
    neptune neptune_barycenter
    pluto pluto_barycenter
    sun moon earth
  )a

  @doc """
  Load an SPK/BSP ephemeris file.

  Returns an ephemeris handle that can be passed to `position/4`.
  The file is not held open; it is read on each `position/4` call.

  ## Example

      eph = Orbis.Ephemeris.load("/path/to/de421.bsp")
  """
  @spec load(String.t()) :: t()
  def load(path) when is_binary(path) do
    if !File.exists?(path) do
      raise ArgumentError, "ephemeris file not found: #{path}"
    end

    %__MODULE__{path: Path.expand(path)}
  end

  @doc """
  Compute the position of `target` relative to `observer` at the given time.

  Returns `{x, y, z}` in km in the J2000/ICRF reference frame.

  The `target` and `observer` are body atoms (see module docs) or
  NAIF integer codes.

  The `datetime` can be a `DateTime`, a `NaiveDateTime`, or a Julian Date
  (TDB) as a float.

  ## Options

    * `:skyfield_compat` - when `true`, replicates Skyfield's exact
      VectorSum accumulation order for bit-identical (0 ULP) parity.
      Default `false` uses direct common-ancestor subtraction which
      is more numerically precise.

  ## Examples

      # Default (precise)
      {x, y, z} = Orbis.Ephemeris.position(eph, :moon, :earth, datetime)

      # Skyfield-compatible (0 ULP oracle test)
      {x, y, z} = Orbis.Ephemeris.position(eph, :moon, :earth, datetime, skyfield_compat: true)
  """
  @spec position(
          t(),
          atom() | integer(),
          atom() | integer(),
          DateTime.t() | NaiveDateTime.t() | float(),
          keyword()
        ) ::
          {float(), float(), float()}
  def position(%__MODULE__{path: path}, target, observer, datetime, opts \\ []) do
    target_name = resolve_body_name(target)
    observer_name = resolve_body_name(observer)
    {whole, fraction} = to_jd_tdb_split(datetime)
    skyfield_compat = Keyword.get(opts, :skyfield_compat, false)

    Orbis.NIF.get_body_position(
      path,
      target_name,
      observer_name,
      whole,
      fraction,
      skyfield_compat
    )
  end

  # ------------------------------------------------------------------
  # Private helpers
  # ------------------------------------------------------------------

  defp resolve_body_name(atom) when atom in @body_atoms do
    Atom.to_string(atom)
  end

  defp resolve_body_name(code) when is_integer(code) do
    code_to_name(code)
  end

  defp resolve_body_name(other) do
    raise ArgumentError, "invalid body: #{inspect(other)}"
  end

  defp code_to_name(0), do: "ssb"
  defp code_to_name(1), do: "mercury_barycenter"
  defp code_to_name(2), do: "venus_barycenter"
  defp code_to_name(3), do: "earth_moon_barycenter"
  defp code_to_name(4), do: "mars_barycenter"
  defp code_to_name(5), do: "jupiter_barycenter"
  defp code_to_name(6), do: "saturn_barycenter"
  defp code_to_name(7), do: "uranus_barycenter"
  defp code_to_name(8), do: "neptune_barycenter"
  defp code_to_name(9), do: "pluto_barycenter"
  defp code_to_name(10), do: "sun"
  defp code_to_name(301), do: "moon"
  defp code_to_name(399), do: "earth"
  defp code_to_name(code), do: raise(ArgumentError, "unknown NAIF body code: #{code}")

  # Convert UTC datetime to split TDB Julian Date {whole, fraction} using
  # the precise time scale NIF. The split form preserves full precision
  # for the Chebyshev argument computation inside the SPK reader.
  defp to_jd_tdb_split(%DateTime{} = dt) do
    second_with_micro = dt.second + elem(dt.microsecond, 0) / 1_000_000

    Orbis.NIF.utc_to_tdb_jd_split(
      dt.year,
      dt.month,
      dt.day,
      dt.hour,
      dt.minute,
      second_with_micro
    )
  end

  defp to_jd_tdb_split(%NaiveDateTime{} = dt) do
    second_with_micro = dt.second + elem(dt.microsecond, 0) / 1_000_000

    Orbis.NIF.utc_to_tdb_jd_split(
      dt.year,
      dt.month,
      dt.day,
      dt.hour,
      dt.minute,
      second_with_micro
    )
  end

  # If a float is passed, assume it's already TDB Julian Date.
  # Split into integer day + fraction for precision.
  defp to_jd_tdb_split(jd) when is_float(jd), do: {Float.floor(jd), jd - Float.floor(jd)}
  defp to_jd_tdb_split(jd) when is_integer(jd), do: {jd * 1.0, 0.0}
end
