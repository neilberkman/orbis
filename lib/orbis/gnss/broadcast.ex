defmodule Orbis.GNSS.Broadcast do
  @moduledoc """
  A parsed RINEX broadcast-navigation product (GPS LNAV, Galileo I/NAV+F/NAV,
  BeiDou D1/D2, GLONASS).

  Holds the broadcast Keplerian elements and clock terms as a resource handle,
  the broadcast-ephemeris counterpart to `Orbis.GNSS.SP3`. Pass a handle to
  `Orbis.GNSS.Positioning.solve/4` to position from broadcast ephemeris instead
  of a precise SP3 product. The navigation file is parsed exactly once; the
  parsed product is held as a reference, not re-parsed per call.

  Parsing covers RINEX 3.x and 4.xx files: GPS, Galileo, and BeiDou records
  (including BeiDou geostationary satellites), and GLONASS (a PZ-90.11
  state-vector model propagated by Runge-Kutta integration rather than Keplerian
  elements). Other constellations in a mixed file are skipped, as are version-4
  CNAV-family messages.

  The orbit and clock models follow IS-GPS-200 (GPS LNAV), the Galileo OS-SIS-ICD
  (I/NAV + F/NAV), and the BeiDou BDS-SIS-ICD (D1/D2), parsed from RINEX 3.x/4.xx
  navigation records.

  ## Epochs

  `position/3` interprets the query epoch in GPS time (GPST). A `NaiveDateTime` or
  `{{year, month, day}, {hour, minute, second}}` is converted to a continuous
  second-of-J2000 via `Orbis.GNSS.Time`; the crate maps that onto each system's
  own time scale (BDT for BeiDou, UTC-referenced for GLONASS) before selecting the
  governing record. No leap-second shifting is applied to the supplied epoch.
  """

  alias Orbis.GNSS.Time
  alias Orbis.NIF

  @enforce_keys [:handle]
  defstruct [:handle]

  @type t :: %__MODULE__{handle: reference()}

  defmodule State do
    @moduledoc """
    A broadcast-evaluated satellite state at one epoch.

    Position is ITRF/IGS-realization ECEF, in meters (frame and unit fixed in the
    field names). `clock_s` is the satellite clock offset in seconds, the
    broadcast clock-polynomial total including the relativistic eccentricity term
    and the broadcast group delay. The sign convention matches `Orbis.GNSS.SP3`:
    a positive `clock_s` means the satellite clock is **ahead** of system time, so
    the geometric range correction is `range + c * clock_s`.
    """
    @enforce_keys [:x_m, :y_m, :z_m, :clock_s]
    defstruct [:x_m, :y_m, :z_m, :clock_s]

    @type t :: %__MODULE__{
            x_m: float(),
            y_m: float(),
            z_m: float(),
            clock_s: float()
          }
  end

  @doc """
  Parse a RINEX 3.x or 4.xx navigation file from disk.

  Returns `{:ok, %Orbis.GNSS.Broadcast{}}` or `{:error, reason}`. The file
  is read and parsed once; the parsed product is held as a resource handle.
  """
  @spec load(String.t()) :: {:ok, t()} | {:error, term()}
  def load(path) when is_binary(path) do
    with {:ok, text} <- File.read(path) do
      parse(text)
    end
  end

  @doc """
  Like `load/1` but raises on failure.
  """
  @spec load!(String.t()) :: t()
  def load!(path) when is_binary(path) do
    case load(path) do
      {:ok, eph} ->
        eph

      {:error, reason} ->
        raise ArgumentError, "could not load RINEX NAV #{path}: #{inspect(reason)}"
    end
  end

  @doc """
  Parse an in-memory RINEX 3.x or 4.xx navigation text buffer into a handle.
  """
  @spec parse(String.t()) :: {:ok, t()} | {:error, term()}
  def parse(text) when is_binary(text) do
    case NIF.broadcast_parse(text) do
      handle when is_reference(handle) -> {:ok, %__MODULE__{handle: handle}}
      {:error, _} = err -> err
    end
  rescue
    e in ErlangError -> {:error, e.original}
  end

  @doc """
  Evaluate the broadcast state of satellite `sat_id` at `epoch`.

  `sat_id` is the canonical RINEX token, e.g. `"G01"` (GPS PRN 1), `"E12"`,
  `"C30"`, `"R07"`. `epoch` is a `NaiveDateTime` or a
  `{{year, month, day}, {hour, minute, second}}` tuple, interpreted in GPS time.

  Returns `{:ok, %Orbis.GNSS.Broadcast.State{}}` with the ECEF position (meters)
  and satellite clock offset (seconds), `{:error, :no_ephemeris}` when no
  broadcast record covers that satellite at that epoch (the validity window has no
  match — this is **not** extrapolated), or `{:error, reason}` for a malformed
  satellite token or a non-integer-second tuple epoch.

  Evaluating the same satellite across a window reuses the parsed handle; the
  navigation file is never re-read.
  """
  @spec position(t(), String.t(), NaiveDateTime.t() | tuple()) ::
          {:ok, State.t()} | {:error, term()}
  def position(%__MODULE__{handle: handle}, sat_id, epoch) when is_binary(sat_id) do
    with {:ok, system_letter, prn} <- parse_sat_id(sat_id),
         {:ok, t_j2000_s} <- j2000_seconds(epoch) do
      case NIF.broadcast_position(handle, system_letter, prn, t_j2000_s) do
        {x_m, y_m, z_m, clock_s} ->
          {:ok, %State{x_m: x_m, y_m: y_m, z_m: z_m, clock_s: clock_s}}

        nil ->
          {:error, :no_ephemeris}

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

  # The broadcast eval takes a continuous second-of-J2000 (GPST-aligned). Carry a
  # sub-second NaiveDateTime/tuple fraction on the integer conversion so a
  # fractional epoch is accepted, mirroring `Orbis.GNSS.Positioning`.
  defp j2000_seconds(%NaiveDateTime{} = ndt) do
    {micro, _precision} = ndt.microsecond

    case Time.epoch_to_j2000_seconds(%{ndt | microsecond: {0, 0}}) do
      {:ok, seconds} -> {:ok, seconds + micro / 1_000_000.0}
      {:error, _} = err -> err
    end
  end

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
end
