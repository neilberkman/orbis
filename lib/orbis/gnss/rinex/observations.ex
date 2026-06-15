defmodule Orbis.GNSS.RINEX.Observations do
  @moduledoc """
  RINEX 3 observation products: parse a station's observation file, expose its
  header (including the surveyed `APPROX POSITION XYZ`, optional antenna
  `DELTA H/E/N` offset, and carrier phase-shift records), and extract the
  single-frequency pseudoranges that `Orbis.GNSS.Positioning.solve/4` consumes.

  This is the Elixir surface over the `astrodynamics-gnss` RINEX observation
  parser and its Hatanaka (CRINEX) decoder. A file is parsed **once** into a
  resource handle held by the BEAM; accessors operate on that handle and never
  re-read the file.

  Both plain RINEX (`.rnx`) and Hatanaka-compressed CRINEX (`.crx`) text are
  accepted: `load/1` and `parse_auto/1` sniff the first line for the
  `CRINEX VERS / TYPE` marker and decode CRINEX before parsing. Gzip is handled
  upstream by `Orbis.GNSS.Data.fetch/2`, so this module only ever sees plain text
  or CRINEX text.

  ## Example

      {:ok, obs} = Orbis.GNSS.RINEX.Observations.load("ESBC00DNK_..._MO.crx")

      Orbis.GNSS.RINEX.Observations.approx_position(obs)
      # => {3_582_105.291, 532_589.7313, 5_232_754.8054}

      [%{index: i, epoch: epoch} | _] = Orbis.GNSS.RINEX.Observations.epochs(obs)
      {:ok, prs} = Orbis.GNSS.RINEX.Observations.pseudoranges(obs, i, codes: %{"G" => ["C1C"]})
      # prs :: [{"G01", range_m}, ...] — feeds solve/4 verbatim

  ## Default pseudorange codes

  The per-system defaults are version-aware: GPS `C1C`, Galileo `C1C` then
  `C1X`, BeiDou `C1I` for RINEX 3.02 / `C2I` for 3.01 and 3.03+ (the B1I label
  changed between minor versions), GLONASS `C1C`. Override per system with the
  `:codes` option, e.g. `codes: %{"G" => ["C1C"], "C" => ["C2I"]}`.
  """

  alias Orbis.NIF

  @enforce_keys [:handle]
  defstruct [:handle]

  @type t :: %__MODULE__{handle: reference()}

  @typedoc "A pseudorange observation `{satellite_id, range_m}`."
  @type observation :: {String.t(), float()}

  @typedoc "An epoch descriptor as returned by `epochs/1`."
  @type epoch_entry :: %{
          index: non_neg_integer(),
          epoch: {{integer(), integer(), integer()}, {integer(), integer(), float()}},
          flag: 0..255,
          sat_count: non_neg_integer()
        }

  @typedoc "GLONASS satellite id to FDMA frequency-channel number."
  @type glonass_slot_map :: %{String.t() => integer()}

  @doc """
  Load and parse a RINEX observation file from disk.

  The file may be plain RINEX (`.rnx`) or Hatanaka CRINEX (`.crx`); the first
  line is sniffed for the CRINEX marker and decoded if present. Returns
  `{:ok, %Orbis.GNSS.RINEX.Observations{}}` or `{:error, reason}`.
  """
  @spec load(String.t()) :: {:ok, t()} | {:error, term()}
  def load(path) when is_binary(path) do
    with {:ok, text} <- File.read(path) do
      parse_auto(text)
    end
  end

  @doc """
  Like `load/1` but raises on failure.
  """
  @spec load!(String.t()) :: t()
  def load!(path) when is_binary(path) do
    case load(path) do
      {:ok, obs} ->
        obs

      {:error, reason} ->
        raise ArgumentError, "could not load RINEX OBS #{path}: #{inspect(reason)}"
    end
  end

  @doc """
  Parse text, auto-detecting plain RINEX vs CRINEX from the first line.
  """
  @spec parse_auto(binary()) :: {:ok, t()} | {:error, term()}
  def parse_auto(text) when is_binary(text) do
    if crinex?(text), do: parse_crinex(text), else: parse(text)
  end

  @doc """
  Parse plain RINEX 3 observation text into a handle.
  """
  @spec parse(binary()) :: {:ok, t()} | {:error, term()}
  def parse(text) when is_binary(text) do
    wrap(NIF.rinex_obs_parse(text))
  rescue
    e in ErlangError -> {:error, e.original}
  end

  @doc """
  Decode Hatanaka CRINEX text and parse the result into a handle.
  """
  @spec parse_crinex(binary()) :: {:ok, t()} | {:error, term()}
  def parse_crinex(text) when is_binary(text) do
    wrap(NIF.crinex_obs_parse(text))
  rescue
    e in ErlangError -> {:error, e.original}
  end

  @doc """
  Decode Hatanaka CRINEX text into the plain RINEX observation text it expands
  to. Returns `{:ok, rinex_text}` or `{:error, reason}`.
  """
  @spec decode_crinex(binary()) :: {:ok, String.t()} | {:error, term()}
  def decode_crinex(text) when is_binary(text) do
    case NIF.crinex_decode(text) do
      rnx when is_binary(rnx) -> {:ok, rnx}
      {:error, _} = err -> err
      other -> {:error, other}
    end
  rescue
    e in ErlangError -> {:error, e.original}
  end

  @doc """
  The surveyed a-priori receiver position `{x_m, y_m, z_m}` (ECEF meters), or
  `nil` if the file carries no `APPROX POSITION XYZ`.
  """
  @spec approx_position(t()) :: {float(), float(), float()} | nil
  def approx_position(%__MODULE__{handle: handle}) do
    NIF.rinex_obs_approx_position(handle)
  rescue
    e in ErlangError ->
      raise ArgumentError, "could not read approx position: #{inspect(e.original)}"
  end

  @doc """
  The antenna reference-point offset from the marker `{height_m, east_m, north_m}`,
  or `nil` if the file carries no `ANTENNA: DELTA H/E/N` header record.

  RINEX stores this field in local height/east/north coordinates. For a station
  whose `APPROX POSITION XYZ` is the marker, add this local offset before
  comparing an observation-derived baseline to antenna-reference-point truth.
  """
  @spec antenna_delta_hen(t()) :: {float(), float(), float()} | nil
  def antenna_delta_hen(%__MODULE__{handle: handle}) do
    NIF.rinex_obs_antenna_delta_hen(handle)
  rescue
    e in ErlangError ->
      raise ArgumentError, "could not read antenna delta H/E/N: #{inspect(e.original)}"
  end

  @doc """
  Carrier phase-shift header records from `SYS / PHASE SHIFT`, in file order.

  Each row is a map with `:system`, `:code`, `:correction_cycles`, and
  `:satellites`. An empty satellite list means the correction applies to every
  satellite for that system/code.
  """
  @spec phase_shifts(t()) :: [
          %{
            system: String.t(),
            code: String.t(),
            correction_cycles: float(),
            satellites: [String.t()]
          }
        ]
  def phase_shifts(%__MODULE__{handle: handle}) do
    handle
    |> NIF.rinex_obs_phase_shifts()
    |> Enum.map(fn {system, code, correction_cycles, satellites} ->
      %{
        system: system,
        code: code,
        correction_cycles: correction_cycles,
        satellites: satellites
      }
    end)
  rescue
    e in ErlangError ->
      raise ArgumentError, "could not read phase shifts: #{inspect(e.original)}"
  end

  @doc """
  The per-constellation observation-code table as a map of system letter to the
  ordered code list, e.g. `%{"G" => ["C1C", ...], "E" => [...]}`.
  """
  @spec observation_codes(t()) :: %{String.t() => [String.t()]}
  def observation_codes(%__MODULE__{handle: handle}) do
    handle
    |> NIF.rinex_obs_codes()
    |> Map.new()
  rescue
    e in ErlangError ->
      raise ArgumentError, "could not read observation codes: #{inspect(e.original)}"
  end

  @doc """
  The GLONASS satellite slot/frequency-channel map from the optional
  `GLONASS SLOT / FRQ #` header records.

  The keys are RINEX satellite ids such as `"R01"` and values are the FDMA
  frequency-channel numbers used to derive GLONASS G1/G2 carrier frequencies.
  Returns `%{}` when the observation file does not carry the header records.
  """
  @spec glonass_slots(t()) :: glonass_slot_map()
  def glonass_slots(%__MODULE__{handle: handle}) do
    handle
    |> NIF.rinex_obs_glonass_slots()
    |> Map.new()
  rescue
    e in ErlangError ->
      raise ArgumentError, "could not read GLONASS slots: #{inspect(e.original)}"
  end

  @doc """
  The epoch list as `[%{index:, epoch:, flag:, sat_count:}]`. The `:epoch` is a
  `{{y, mo, d}, {h, mi, second_float}}` tuple in the file's time scale — exactly
  the form `Orbis.GNSS.Positioning.solve/4` accepts.
  """
  @spec epochs(t()) :: [epoch_entry()]
  def epochs(%__MODULE__{handle: handle}) do
    handle
    |> NIF.rinex_obs_epochs()
    |> Enum.with_index()
    |> Enum.map(fn {{epoch, flag, sat_count}, index} ->
      %{index: index, epoch: epoch, flag: flag, sat_count: sat_count}
    end)
  rescue
    e in ErlangError ->
      raise ArgumentError, "could not read epochs: #{inspect(e.original)}"
  end

  @doc """
  Extract single-frequency pseudoranges for one epoch.

  `epoch` is either the integer epoch index (from `epochs/1`) or an epoch tuple
  `{{y, mo, d}, {h, mi, s}}`, which is resolved to its index.

  Without `:codes`, the version-aware defaults are applied across every system.
  When `:codes` is given it **defines the whole policy**: only the listed systems
  are extracted, each with its given code preference, e.g. `codes: %{"G" =>
  ["C1C"]}` yields GPS-only pseudoranges and `codes: %{"G" => ["C1C"], "C" =>
  ["C2I"]}` yields GPS + BeiDou.

  Returns `{:ok, [{"G01", range_m}, ...]}` (ascending satellite id) or
  `{:error, :epoch_out_of_range}` / `{:error, :unknown_epoch}`.
  """
  @spec pseudoranges(t(), non_neg_integer() | tuple(), keyword()) ::
          {:ok, [observation()]} | {:error, term()}
  def pseudoranges(obs, epoch, opts \\ [])

  def pseudoranges(%__MODULE__{handle: handle}, index, opts) when is_integer(index) do
    overrides = codes_overrides(Keyword.get(opts, :codes, %{}))

    case NIF.rinex_obs_pseudoranges(handle, index, overrides) do
      {:ok, prs} -> {:ok, prs}
      {:error, reason} -> {:error, reason}
      other -> {:error, other}
    end
  rescue
    e in ErlangError -> {:error, e.original}
  end

  def pseudoranges(%__MODULE__{} = obs, {{_y, _mo, _d}, {_h, _mi, _s}} = epoch, opts) do
    case Enum.find(epochs(obs), fn e -> e.epoch == epoch end) do
      %{index: index} -> pseudoranges(obs, index, opts)
      nil -> {:error, :unknown_epoch}
    end
  end

  # Speed of light in vacuum, m/s (for the carrier-phase wavelength).
  @speed_of_light_m_s 299_792_458.0

  @doc """
  Every observation value for one epoch, keyed by satellite.

  `epoch` is the integer index (from `epochs/1`) or an epoch tuple
  `{{y, mo, d}, {h, mi, s}}`. Unlike `pseudoranges/3` this returns the raw RINEX
  observations across code types — pseudorange, carrier phase, Doppler, and
  signal strength — so callers can build carrier-phase combinations.

  Returns `{:ok, %{satellite_id => [obs]}}` where each `obs` is

      %{code: "L1C", kind: :carrier_phase, value: 1.23e8, units: :cycles,
        lli: 0 | nil, ssi: 7 | nil}

  `kind`/`units` follow the RINEX code's leading letter (`C` → `:pseudorange`/
  `:meters`, `L` → `:carrier_phase`/`:cycles`, `D` → `:doppler`/`:hz`, `S` →
  `:signal_strength`/`:db_hz`). A blank observation has a `nil` value. Returns
  `{:error, :epoch_out_of_range}` / `{:error, :unknown_epoch}`.

  ## Options

    * `:codes` — a per-system code filter, e.g. `%{"G" => ["L1C", "L2W"]}`. By
      default every code for every satellite is returned; a non-empty filter
      restricts the result (and the data crossing the NIF boundary) to the listed
      systems, and within each to the listed codes. A system mapped to `[]` keeps
      all of that system's codes — e.g. `%{"G" => []}` is GPS-only, all codes.
  """
  @spec values(t(), non_neg_integer() | tuple(), keyword()) ::
          {:ok, %{String.t() => [map()]}} | {:error, term()}
  def values(obs, epoch, opts \\ [])

  def values(%__MODULE__{handle: handle}, index, opts) when is_integer(index) do
    overrides = codes_overrides(Keyword.get(opts, :codes, %{}))

    case NIF.rinex_obs_values(handle, index, overrides) do
      {:ok, rows} ->
        {:ok,
         Map.new(rows, fn {sat, code_values} ->
           {sat,
            Enum.map(code_values, fn {code, value, lli, ssi} ->
              decode_obs(code, value, lli, ssi)
            end)}
         end)}

      {:error, reason} ->
        {:error, reason}

      other ->
        {:error, other}
    end
  rescue
    e in ErlangError -> {:error, e.original}
  end

  def values(%__MODULE__{} = obs, {{_y, _mo, _d}, {_h, _mi, _s}} = epoch, opts) do
    case Enum.find(epochs(obs), fn e -> e.epoch == epoch end) do
      %{index: index} -> values(obs, index, opts)
      nil -> {:error, :unknown_epoch}
    end
  end

  @doc """
  Carrier-phase observations for one epoch (the `L*` codes), keyed by satellite.

  Convenience over `values/3` that keeps only carrier phase and adds the
  wavelength and the phase expressed in metres when the carrier frequency is
  known for the satellite's system and band (GPS, Galileo, BeiDou, and GLONASS
  when the file carries `GLONASS SLOT / FRQ #` records). GLONASS G1/G2
  wavelengths are derived from the satellite's parsed FDMA frequency-channel
  number; a GLONASS satellite without a channel map entry keeps
  `wavelength_m`/`value_m` as `nil`.

  ## `SYS / PHASE SHIFT` correction

  When the file carries `SYS / PHASE SHIFT` header records with a non-zero
  `correction_cycles`, that fractional-cycle bias is added to `value_cycles`
  (and folded into `value_m`) so the carrier phase is aligned to a common
  reference, which is what an integer-ambiguity resolver requires. A record's
  satellite list scopes the correction: an empty list applies to every
  satellite of that system/code, a non-empty list applies only to the listed
  satellites. The applied offset is reported on each phase row as
  `:phase_shift_cycles` (0.0 when no record matches), and `:value_cycles`
  already includes it. With the common all-zero `SYS / PHASE SHIFT` records the
  output is bit-identical to the uncorrected values.

  Returns `{:ok, %{satellite_id => [phase]}}` where each `phase` is

      %{code: "L1C", value_cycles: 1.23e8, lli: 0 | nil, ssi: 7 | nil,
        frequency_hz: 1.57542e9 | nil, wavelength_m: 0.1903 | nil,
        value_m: 2.34e7 | nil, phase_shift_cycles: 0.0}
  """
  @spec phases(t(), non_neg_integer() | tuple(), keyword()) ::
          {:ok, %{String.t() => [map()]}} | {:error, term()}
  def phases(%__MODULE__{} = obs, epoch, opts \\ []) do
    with {:ok, by_sat} <- values(obs, epoch, opts) do
      glonass_slots = glonass_slots(obs)
      shift_table = phase_shift_table(obs)

      phases =
        Map.new(by_sat, fn {sat, observations} ->
          rows =
            observations
            |> Enum.filter(&(&1.kind == :carrier_phase))
            |> Enum.map(&phase_row(sat, &1, glonass_slots, shift_table))

          {sat, rows}
        end)

      {:ok, phases}
    end
  end

  @doc """
  Carrier frequency in hertz for a system letter and RINEX band digit.

  The two-argument form covers fixed-frequency systems (`"G"`, `"E"`, `"C"`)
  and returns `nil` for GLONASS because its G1/G2 carriers are FDMA
  channel-dependent. Use the three-argument form with the parsed GLONASS
  frequency-channel number:

      Orbis.GNSS.RINEX.Observations.band_frequency_hz("R", "1", 1)
      # => 1602562500.0

  For GLONASS, band `"1"` is G1 (`1602 MHz + k * 562.5 kHz`) and band `"2"` is
  G2 (`1246 MHz + k * 437.5 kHz`), where `k` is the frequency-channel number.
  Unknown bands return `nil`.
  """
  @spec band_frequency_hz(String.t(), String.t()) :: float() | nil
  def band_frequency_hz(system, band), do: band_frequency_hz(system, band, nil)

  @spec band_frequency_hz(String.t(), String.t(), integer() | nil) :: float() | nil
  def band_frequency_hz("G", band, _channel) do
    case band do
      "1" -> 1_575_420_000.0
      "2" -> 1_227_600_000.0
      "5" -> 1_176_450_000.0
      _ -> nil
    end
  end

  def band_frequency_hz("E", band, _channel) do
    case band do
      "1" -> 1_575_420_000.0
      "5" -> 1_176_450_000.0
      "6" -> 1_278_750_000.0
      "7" -> 1_207_140_000.0
      "8" -> 1_191_795_000.0
      _ -> nil
    end
  end

  def band_frequency_hz("C", band, _channel) do
    case band do
      "1" -> 1_575_420_000.0
      "2" -> 1_561_098_000.0
      "5" -> 1_176_450_000.0
      "6" -> 1_268_520_000.0
      "7" -> 1_207_140_000.0
      "8" -> 1_191_795_000.0
      _ -> nil
    end
  end

  def band_frequency_hz("R", band, channel) when is_integer(channel) do
    case band do
      "1" -> 1_602_000_000.0 + channel * 562_500.0
      "2" -> 1_246_000_000.0 + channel * 437_500.0
      _ -> nil
    end
  end

  def band_frequency_hz(_system, _band, _channel), do: nil

  # --- helpers -------------------------------------------------------------

  defp decode_obs(code, value, lli, ssi) do
    kind = obs_kind(code)
    %{code: code, kind: kind, value: value, units: obs_units(kind), lli: lli, ssi: ssi}
  end

  defp obs_kind(<<?C, _::binary>>), do: :pseudorange
  defp obs_kind(<<?L, _::binary>>), do: :carrier_phase
  defp obs_kind(<<?D, _::binary>>), do: :doppler
  defp obs_kind(<<?S, _::binary>>), do: :signal_strength
  defp obs_kind(_), do: :unknown

  defp obs_units(:pseudorange), do: :meters
  defp obs_units(:carrier_phase), do: :cycles
  defp obs_units(:doppler), do: :hz
  defp obs_units(:signal_strength), do: :db_hz
  defp obs_units(:unknown), do: :unknown

  defp phase_row(sat, obs, glonass_slots, shift_table) do
    freq =
      band_frequency_hz(
        String.first(sat),
        String.at(obs.code, 1),
        Map.get(glonass_slots, sat)
      )

    shift_cycles = phase_shift_cycles(shift_table, sat, obs.code)

    value_cycles =
      case obs.value do
        v when is_number(v) -> v + shift_cycles
        _ -> obs.value
      end

    {wavelength_m, value_m} =
      cond do
        is_number(freq) and is_number(value_cycles) ->
          lambda = @speed_of_light_m_s / freq
          {lambda, value_cycles * lambda}

        is_number(freq) ->
          {@speed_of_light_m_s / freq, nil}

        true ->
          {nil, nil}
      end

    %{
      code: obs.code,
      value_cycles: value_cycles,
      lli: obs.lli,
      ssi: obs.ssi,
      frequency_hz: freq,
      wavelength_m: wavelength_m,
      value_m: value_m,
      phase_shift_cycles: shift_cycles
    }
  end

  # Build the `SYS / PHASE SHIFT` lookup keyed by {system, code}. Each value is
  # the list of {correction_cycles, satellites_set} records for that key; a
  # record with an empty satellite list applies to every satellite of the
  # system/code, otherwise only to the listed satellites. Empty when the file
  # carries no records, in which case every lookup returns a 0.0 correction.
  defp phase_shift_table(%__MODULE__{} = obs) do
    obs
    |> phase_shifts()
    |> Enum.reduce(%{}, fn %{system: system, code: code} = row, acc ->
      entry = {row.correction_cycles, MapSet.new(row.satellites)}
      Map.update(acc, {system, code}, [entry], &[entry | &1])
    end)
  end

  # Resolve the phase-shift correction (cycles) for a satellite/code. Prefers a
  # per-satellite record (non-empty satellite list containing the sat) over a
  # system-wide record (empty satellite list); returns 0.0 when none matches.
  defp phase_shift_cycles(table, sat, code) do
    system = String.first(sat)

    case Map.get(table, {system, code}) do
      nil ->
        0.0

      records ->
        per_sat =
          Enum.find_value(records, fn {cycles, sats} ->
            if MapSet.size(sats) > 0 and MapSet.member?(sats, sat), do: cycles
          end)

        if is_number(per_sat) do
          per_sat
        else
          Enum.find_value(records, 0.0, fn {cycles, sats} ->
            if MapSet.size(sats) == 0, do: cycles
          end)
        end
    end
  end

  defp wrap(handle) when is_reference(handle), do: {:ok, %__MODULE__{handle: handle}}
  defp wrap({:error, _} = err), do: err
  defp wrap(other), do: {:error, other}

  defp crinex?(text) do
    case String.split(text, "\n", parts: 2) do
      [first | _] -> String.contains?(first, "CRINEX VERS")
      _ -> false
    end
  end

  # Normalize a `%{"G" => ["C1C"]}` override map into the NIF's
  # `[{"G", ["C1C"]}]` form. An empty map means "use the crate defaults".
  defp codes_overrides(map) when is_map(map) do
    Enum.map(map, fn {sys, codes} -> {to_string(sys), Enum.map(codes, &to_string/1)} end)
  end

  defp codes_overrides(list) when is_list(list) do
    Enum.map(list, fn {sys, codes} -> {to_string(sys), Enum.map(codes, &to_string/1)} end)
  end
end
