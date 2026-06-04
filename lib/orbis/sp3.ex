defmodule Orbis.SP3 do
  @moduledoc """
  SP3-c / SP3-d precise-ephemeris products (IGS precise orbits + clocks).

  This is the Elixir surface over the `astrodynamics-gnss` SP3 parser and
  `scipy.interpolate`-matched position/clock interpolation. It is **not** the
  JPL-SPK reader (`Orbis.Ephemeris`): SP3 carries GNSS satellite states in the
  ITRF/IGS ECEF frame, in meters, tagged by a GNSS satellite id like `"G01"`.

  A file is parsed **once** into a resource handle held by the BEAM; evaluation
  operates on that handle and never re-reads the file.

  ## Example

      {:ok, sp3} = Orbis.SP3.load("/path/to/igs.sp3")
      {:ok, state} =
        Orbis.SP3.position(sp3, "G01", ~N[2020-06-24 00:00:00])

      state.x_m       # ITRF/IGS ECEF X, meters
      state.clock_s   # satellite clock offset, seconds (or nil if no estimate)

  ## Epochs

  The query epoch is interpreted in the file's **own** time scale (read from the
  SP3 header — typically GPST). Pass a `NaiveDateTime` or a
  `{{year, month, day}, {hour, minute, second}}` tuple; it is converted to the
  split Julian date with the same midnight-boundary convention the parser uses
  (no leap-second shifting — the epoch stays in the file's scale).
  """

  alias Orbis.NIF

  @enforce_keys [:handle, :time_scale]
  defstruct [:handle, :time_scale]

  @type t :: %__MODULE__{handle: reference(), time_scale: String.t()}

  defmodule State do
    @moduledoc """
    An interpolated SP3 satellite state at one epoch.

    Position is ITRF/IGS-realization ECEF, in meters (frame and unit are fixed
    in the field names per the spec's frames-in-the-type-system rule). `clock_s`
    is the satellite clock offset in seconds, or `nil` when the product carries
    no clock estimate for that satellite/epoch.
    """
    @enforce_keys [:x_m, :y_m, :z_m, :clock_s]
    defstruct [:x_m, :y_m, :z_m, :clock_s]

    @type t :: %__MODULE__{
            x_m: float(),
            y_m: float(),
            z_m: float(),
            clock_s: float() | nil
          }
  end

  @doc """
  Load and parse an SP3-c / SP3-d file into a product handle.

  Returns `{:ok, %Orbis.SP3{}}` or `{:error, reason}`. The file is read and
  parsed exactly once; the parsed product is held as a resource handle.
  """
  @spec load(String.t()) :: {:ok, t()} | {:error, term()}
  def load(path) when is_binary(path) do
    with {:ok, bytes} <- File.read(path) do
      parse_bytes(bytes)
    end
  end

  @doc """
  Like `load/1` but raises on failure.
  """
  @spec load!(String.t()) :: t()
  def load!(path) when is_binary(path) do
    case load(path) do
      {:ok, sp3} -> sp3
      {:error, reason} -> raise ArgumentError, "could not load SP3 #{path}: #{inspect(reason)}"
    end
  end

  @doc """
  Parse an in-memory SP3 byte buffer (already decompressed) into a handle.
  """
  @spec parse(binary()) :: {:ok, t()} | {:error, term()}
  def parse(bytes) when is_binary(bytes), do: parse_bytes(bytes)

  defp parse_bytes(bytes) do
    case NIF.sp3_parse(bytes) do
      handle when is_reference(handle) ->
        {:ok, %__MODULE__{handle: handle, time_scale: NIF.sp3_time_scale(handle)}}

      {:error, _} = err ->
        err

      other ->
        {:error, other}
    end
  rescue
    e in ErlangError -> {:error, e.original}
  end

  @doc """
  Return the SP3/RINEX satellite identifiers declared by the product header.

  These are canonical three-character tokens such as `"G01"`, `"E12"`, or
  `"C30"`. The list is read from the already-loaded SP3 handle; no file I/O or
  interpolation is performed.

  ## Examples

      {:ok, sp3} = Orbis.SP3.parse(sp3_bytes)
      ids = Orbis.SP3.satellite_ids(sp3)
      "G01" in ids
  """
  @spec satellite_ids(t()) :: [String.t()]
  def satellite_ids(%__MODULE__{handle: handle}) do
    NIF.sp3_satellite_ids(handle)
  rescue
    e in ErlangError ->
      raise ArgumentError, "could not read SP3 satellite ids: #{inspect(e.original)}"
  end

  @doc """
  Interpolate the state of satellite `sat_id` at `epoch`.

  `sat_id` is the canonical SP3/RINEX token, e.g. `"G01"` (GPS PRN 1), `"E12"`,
  `"C30"`. `epoch` is a `NaiveDateTime` or a
  `{{year, month, day}, {hour, minute, second}}` tuple, interpreted in the
  file's own time scale.

  Returns `{:ok, %Orbis.SP3.State{}}` or `{:error, reason}`.
  """
  @spec position(t(), String.t(), NaiveDateTime.t() | tuple()) ::
          {:ok, State.t()} | {:error, term()}
  def position(%__MODULE__{handle: handle, time_scale: scale}, sat_id, epoch)
      when is_binary(sat_id) do
    with {:ok, system_letter, prn} <- parse_sat_id(sat_id),
         {jd_whole, jd_fraction} <- epoch_to_split_jd(epoch) do
      case NIF.sp3_position(handle, system_letter, prn, scale, jd_whole, jd_fraction) do
        {x_m, y_m, z_m, clock} ->
          # `clock` is already `nil` (no estimate) or a float (seconds).
          {:ok, %State{x_m: x_m, y_m: y_m, z_m: z_m, clock_s: clock}}

        {:error, _} = err ->
          err

        other ->
          {:error, other}
      end
    end
  rescue
    e in ErlangError -> {:error, e.original}
  end

  # --- helpers -------------------------------------------------------------

  defp parse_sat_id(<<letter::binary-size(1), rest::binary>>) do
    case Integer.parse(rest) do
      {prn, ""} when prn >= 0 and prn <= 255 -> {:ok, String.upcase(letter), prn}
      _ -> {:error, {:bad_sat_id, letter <> rest}}
    end
  end

  defp parse_sat_id(other), do: {:error, {:bad_sat_id, other}}

  defp epoch_to_split_jd(%NaiveDateTime{} = ndt) do
    {micro, _precision} = ndt.microsecond
    seconds = ndt.second + micro / 1_000_000.0
    epoch_to_split_jd({{ndt.year, ndt.month, ndt.day}, {ndt.hour, ndt.minute, seconds}})
  end

  defp epoch_to_split_jd({{year, month, day}, {hour, minute, second}}) do
    # Fliegel-Van Flandern Gregorian -> Julian Day Number, in integer
    # arithmetic so the whole-day boundary is exact. Mirrors the crate's
    # `civil_to_julian_split` (sp3.rs): jd_whole is the *.5 midnight boundary of
    # the civil day, the time-of-day is the fraction. No leap-second shifting —
    # the epoch stays in the file's own scale.
    a = div(14 - month, 12)
    y = year + 4800 - a
    m = month + 12 * a - 3

    jdn =
      day + div(153 * m + 2, 5) + 365 * y + div(y, 4) - div(y, 100) + div(y, 400) - 32_045

    jd_whole = jdn - 0.5
    day_seconds = hour * 3600.0 + minute * 60.0 + second / 1.0
    fraction = day_seconds / 86_400.0
    {jd_whole, fraction}
  end
end
