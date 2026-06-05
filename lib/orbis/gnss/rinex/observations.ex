defmodule Orbis.GNSS.RINEX.Observations do
  @moduledoc """
  RINEX 3 observation products: parse a station's observation file, expose its
  header (including the surveyed `APPROX POSITION XYZ`), and extract the
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

  # --- helpers -------------------------------------------------------------

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
