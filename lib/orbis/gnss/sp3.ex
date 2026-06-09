defmodule Orbis.GNSS.SP3 do
  @moduledoc """
  SP3-c / SP3-d precise-ephemeris products (IGS precise orbits + clocks).

  This is the Elixir surface over the `astrodynamics-gnss` SP3 parser and
  `scipy.interpolate`-matched position/clock interpolation. It is **not** the
  JPL-SPK reader (`Orbis.Ephemeris`): SP3 carries GNSS satellite states in the
  ITRF/IGS ECEF frame, in meters, tagged by a GNSS satellite id like `"G01"`.

  A file is parsed **once** into a resource handle held by the BEAM; evaluation
  operates on that handle and never re-reads the file.

  ## Example

      {:ok, sp3} = Orbis.GNSS.SP3.load("/path/to/igs.sp3")
      {:ok, state} =
        Orbis.GNSS.SP3.position(sp3, "G01", ~N[2020-06-24 00:00:00])

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

  Returns `{:ok, %Orbis.GNSS.SP3{}}` or `{:error, reason}`. The file is read and
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

      {:ok, sp3} = Orbis.GNSS.SP3.parse(sp3_bytes)
      ids = Orbis.GNSS.SP3.satellite_ids(sp3)
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
  Serialize the product to standard SP3-c / SP3-d text as iodata. Pure — no I/O.

  This is the inverse of `load/1` / `parse/1`: a read → (`merge/2`) → write
  pipeline round-trips to a single standard SP3 file any reader consumes. The
  output is deterministic (same product → identical bytes). Header fields
  (version, epoch count, satellite list, time system, week / seconds-of-week /
  MJD / interval) are derived from the product. A satellite absent at an epoch is
  written as the SP3 missing-orbit sentinel — so a quarantined `merge/2` cell
  re-reads as missing, never a fabricated position.

  To write to disk (optionally gzipped, with an atomic commit), use
  `Orbis.GNSS.Data.write_sp3/3`.

  ## Examples

      {:ok, sp3} = Orbis.GNSS.SP3.load("igs.sp3")
      iodata = Orbis.GNSS.SP3.to_iodata(sp3)
      {:ok, ^sp3_equal} = Orbis.GNSS.SP3.parse(IO.iodata_to_binary(iodata))
  """
  @spec to_iodata(t(), keyword()) :: iodata()
  def to_iodata(%__MODULE__{handle: handle}, _opts \\ []) do
    NIF.sp3_to_iodata(handle)
  rescue
    e in ErlangError ->
      raise ArgumentError, "could not serialize SP3 product: #{inspect(e.original)}"
  end

  @doc """
  Interpolate the state of satellite `sat_id` at `epoch`.

  `sat_id` is the canonical SP3/RINEX token, e.g. `"G01"` (GPS PRN 1), `"E12"`,
  `"C30"`. `epoch` is a `NaiveDateTime` or a
  `{{year, month, day}, {hour, minute, second}}` tuple, interpreted in the
  file's own time scale.

  Returns `{:ok, %Orbis.GNSS.SP3.State{}}` or `{:error, reason}`.
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

  @doc """
  Merge several SP3 products from different analysis centers into one consistent
  precise-ephemeris dataset.

  `sources` is a list of loaded products **in precedence order** (earlier wins
  ties). This is orthogonal to time-stitching: it combines providers at the same
  epochs. For every `(epoch, satellite)` cell in the union of the inputs:

    * **Union coverage** — a satellite present in any input is present in the
      merged product for that epoch (filling a single center's dropouts).
    * **Consensus** — the largest subset of sources agreeing within tolerance is
      combined; sources outside it are recorded as outliers. A cell with no
      agreeing subset of `:min_agree` is *quarantined* (omitted), never averaged
      across disagreeing centers. A lone source is carried through.

  Returns `{:ok, %Orbis.GNSS.SP3{}, report}` or `{:error, reason}`, where
  `report` is a map with `:quarantined`, `:single_source`, and
  `:position_outliers` lists. Each entry is a map
  `%{satellite: "G03", jd_whole: float, jd_fraction: float, sources: [0, 2]}`
  (`sources` are zero-based indices into `sources`).

  ## Options

    * `:position_tolerance_m` — position agreement tolerance, meters (default `0.5`)
    * `:clock_tolerance_s` — clock agreement tolerance, seconds (default `5.0e-9`)
    * `:min_agree` — agreeing sources required to accept a contested cell (default `2`)
    * `:clock_min_common` — common clocked satellites for the clock-datum estimate (default `5`)
    * `:combine` — `:mean` (default), `:median`, or `:precedence`
  """
  @spec merge([t()], keyword()) :: {:ok, t(), map()} | {:error, term()}
  def merge(sources, opts \\ []) when is_list(sources) do
    handles = Enum.map(sources, fn %__MODULE__{handle: handle} -> handle end)
    position_tolerance_m = Keyword.get(opts, :position_tolerance_m, 0.5)
    clock_tolerance_s = Keyword.get(opts, :clock_tolerance_s, 5.0e-9)
    min_agree = Keyword.get(opts, :min_agree, 2)
    clock_min_common = Keyword.get(opts, :clock_min_common, 5)
    combine = opts |> Keyword.get(:combine, :mean) |> to_string()

    case NIF.sp3_merge(
           handles,
           position_tolerance_m,
           clock_tolerance_s,
           min_agree,
           clock_min_common,
           combine
         ) do
      {handle, {quarantined, single_source, position_outliers}} when is_reference(handle) ->
        report = %{
          quarantined: Enum.map(quarantined, &to_flag/1),
          single_source: Enum.map(single_source, &to_flag/1),
          position_outliers: Enum.map(position_outliers, &to_flag/1)
        }

        {:ok, %__MODULE__{handle: handle, time_scale: NIF.sp3_time_scale(handle)}, report}

      {:error, _} = err ->
        err

      other ->
        {:error, other}
    end
  rescue
    e in ErlangError -> {:error, e.original}
  end

  @doc """
  Estimate the per-epoch reference-clock offset of `other` relative to
  `reference` (the clock-datum primitive).

  Precise clock products from different centers are referenced to different
  station/ensemble clocks, so their raw clocks differ by a per-epoch common
  offset that drifts over the day. This returns that datum: a list of maps
  `%{jd_whole: float, jd_fraction: float, offset_s: float, satellites: integer}`,
  one per epoch where at least `:min_common` common clocked satellites let the
  (robust median) offset be estimated. Subtract `offset_s` from `other`'s clocks
  to put both products on `reference`'s datum. Orbit positions need no such
  treatment — every center reports ITRF center-of-mass coordinates.

  ## Options

    * `:min_common` — minimum common clocked satellites per epoch (default `5`)
  """
  @spec clock_reference_offset(t(), t(), keyword()) :: [map()]
  def clock_reference_offset(
        %__MODULE__{handle: reference},
        %__MODULE__{handle: other},
        opts \\ []
      ) do
    min_common = Keyword.get(opts, :min_common, 5)

    reference
    |> NIF.sp3_clock_reference_offset(other, min_common)
    |> Enum.map(fn {jd_whole, jd_fraction, offset_s, satellites} ->
      %{jd_whole: jd_whole, jd_fraction: jd_fraction, offset_s: offset_s, satellites: satellites}
    end)
  rescue
    e in ErlangError ->
      raise ArgumentError, "could not estimate clock reference offset: #{inspect(e.original)}"
  end

  @doc """
  Return a copy of `other` with its clocks shifted onto `reference`'s clock datum
  (the clock-datum primitive, applied).

  At every epoch the offset could be estimated, each clocked satellite's offset
  has the datum subtracted, so the result's clocks are directly comparable to
  `reference`'s. Positions are untouched. Epochs without an estimate are left
  unchanged. The returned product interpolates like any other SP3.

  Returns `{:ok, %Orbis.GNSS.SP3{}}` or `{:error, reason}`.

  ## Options

    * `:min_common` — minimum common clocked satellites per epoch (default `5`)
  """
  @spec align_clock_reference(t(), t(), keyword()) :: {:ok, t()} | {:error, term()}
  def align_clock_reference(
        %__MODULE__{handle: reference},
        %__MODULE__{handle: other},
        opts \\ []
      ) do
    min_common = Keyword.get(opts, :min_common, 5)

    case NIF.sp3_align_clock_reference(reference, other, min_common) do
      handle when is_reference(handle) ->
        {:ok, %__MODULE__{handle: handle, time_scale: NIF.sp3_time_scale(handle)}}

      {:error, _} = err ->
        err

      other_result ->
        {:error, other_result}
    end
  rescue
    e in ErlangError -> {:error, e.original}
  end

  # --- helpers -------------------------------------------------------------

  defp to_flag({satellite, jd_whole, jd_fraction, sources}) do
    %{satellite: satellite, jd_whole: jd_whole, jd_fraction: jd_fraction, sources: sources}
  end

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
