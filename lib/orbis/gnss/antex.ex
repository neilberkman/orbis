defmodule Orbis.GNSS.Antex do
  @moduledoc """
  Parser and lookup helpers for ANTEX 1.4 receiver and satellite antenna blocks.

  The module parses the relevant parts of ANTEX text into Elixir structs and
  exposes minimal helpers for Phase 3a:

    * `load/1`, `load!/1` — read and parse an ANTEX file.
    * `antenna/2` — resolve an antenna by the `TYPE / SERIAL` key.
    * `pco/2` — retrieve frequency PCO values in meters.
    * `pcv/4` — interpolate phase-center variation in meters.

  Interpolation is linear in zenith and, when azimuth data exists, linear in
  azimuth with wrap handling at 0/360°.
  """

  alias Orbis.GNSS.Antex.Antenna
  alias Orbis.GNSS.Antex.Frequency

  @enforce_keys [:antennas]
  defstruct [:antennas]

  defmodule Antenna do
    @moduledoc false
    @enforce_keys [
      :id,
      :kind,
      :type,
      :serial,
      :dazi_deg,
      :zenith_start_deg,
      :zenith_end_deg,
      :zenith_step_deg,
      :sinex_code,
      :frequencies
    ]

    defstruct [
      :id,
      :kind,
      :type,
      :serial,
      :dazi_deg,
      :zenith_start_deg,
      :zenith_end_deg,
      :zenith_step_deg,
      :sinex_code,
      :valid_from,
      :valid_until,
      :frequencies
    ]
  end

  defmodule Frequency do
    @moduledoc false
    @enforce_keys [:frequency, :pco_m, :pcv_samples]
    defstruct [:frequency, :pco_m, :pcv_samples]
  end

  @type t :: %__MODULE__{antennas: %{optional(String.t()) => Antenna.t()}}

  @type parse_error :: {:error, term()}

  defmodule State do
    defstruct antennas: %{}, current_antenna: nil, current_frequency: nil
  end

  @doc """
  Load and parse an ANTEX file from `path`.
  """
  @spec load(String.t()) :: {:ok, t()} | parse_error
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
      {:ok, antex} ->
        antex

      {:error, reason} ->
        raise ArgumentError, "could not load ANTEX #{path}: #{inspect(reason)}"
    end
  end

  @doc """
  Parse ANTEX text already in memory.
  """
  @spec parse(binary()) :: {:ok, t()} | parse_error
  def parse(text) when is_binary(text) do
    state =
      text
      |> String.split(["\r\n", "\n", "\r"], trim: false)
      |> Enum.reduce(%State{}, &step/2)
      |> finalize_antenna()

    {:ok, %__MODULE__{antennas: state.antennas}}
  rescue
    e ->
      {:error, e}
  end

  @doc """
  Return an antenna by its `TYPE / SERIAL` key.
  """
  @spec antenna(t(), String.t()) :: Antenna.t() | nil
  def antenna(%__MODULE__{antennas: antennas}, id) when is_binary(id) do
    Map.get(antennas, String.trim(id))
  end

  @doc """
  Return the satellite antenna block for PRN `prn` (e.g. `"G05"`) valid at the
  given epoch, or `nil` if none.

  ANTEX carries one satellite-antenna block per spacecraft per validity window
  (the same PRN is reused across spacecraft over time), so the block must be
  selected by `VALID FROM`/`VALID UNTIL`. An unbounded `valid_until` (the active
  spacecraft) is treated as open-ended.
  """
  @spec satellite_antenna(t(), String.t(), NaiveDateTime.t()) :: Antenna.t() | nil
  def satellite_antenna(%__MODULE__{antennas: antennas}, prn, %NaiveDateTime{} = epoch)
      when is_binary(prn) do
    antennas
    |> Map.values()
    |> Enum.find(fn a ->
      a.kind == :satellite and String.trim(a.serial) == String.trim(prn) and
        valid_at?(a, epoch)
    end)
  end

  defp valid_at?(%Antenna{valid_from: from, valid_until: until}, epoch) do
    after_from? = is_nil(from) or NaiveDateTime.compare(epoch, from) != :lt
    before_until? = is_nil(until) or NaiveDateTime.compare(epoch, until) != :gt
    after_from? and before_until?
  end

  @doc """
  Frequency-dependent PCO (north/east/up in meters).
  """
  @spec pco(Antenna.t(), String.t()) :: {float(), float(), float()}
  def pco(%Antenna{} = antenna, frequency) when is_binary(frequency) do
    case Map.get(antenna.frequencies, String.trim(frequency)) do
      %Frequency{pco_m: pco_m} ->
        pco_m

      nil ->
        raise ArgumentError, "unknown frequency #{inspect(frequency)} for #{inspect(antenna.id)}"
    end
  end

  @doc """
  Frequency-dependent phase-center variation in meters.

  Interpolation is linear in zenith and azimuth. Azimuth is optional: when not
  given (or when the antenna has no azimuth-dependent rows), the NOAZI row is
  used.
  """
  @spec pcv(Antenna.t(), String.t(), float(), float() | nil) :: float()
  def pcv(%Antenna{} = antenna, frequency, zenith_deg, azimuth_deg \\ nil)
      when is_number(zenith_deg) do
    frequency_block = Map.get(antenna.frequencies, String.trim(frequency))

    if frequency_block == nil do
      raise ArgumentError, "unknown frequency #{inspect(frequency)} for #{inspect(antenna.id)}"
    end

    samples = frequency_block.pcv_samples

    noazi_samples =
      Enum.filter(samples, fn sample -> sample.grid == :noazi end)
      |> Enum.map(fn sample -> {sample.zenith_deg, sample.value_m} end)

    if azimuth_deg == nil or Enum.all?(samples, fn sample -> sample.grid == :noazi end) do
      interpolate(noazi_samples, zenith_deg)
    else
      azimuth_samples =
        samples
        |> Enum.filter(fn sample -> sample.grid == :azi end)
        |> Enum.group_by(& &1.azimuth_deg, fn sample -> {sample.zenith_deg, sample.value_m} end)

      if azimuth_samples == %{} do
        interpolate(noazi_samples, zenith_deg)
      else
        interpolate_azimuth(azimuth_samples, azimuth_deg, zenith_deg)
      end
    end
  end

  # --- parser --------------------------------------------------------------

  defp step(line, state) do
    case tag(line) do
      "START OF ANTENNA" ->
        state
        |> finalize_antenna()
        |> begin_antenna()

      "END OF ANTENNA" ->
        finalize_antenna(state)

      "TYPE / SERIAL NO" ->
        parse_type_serial(line, state)

      "DAZI" ->
        parse_dazi(line, state)

      "ZEN1 / ZEN2 / DZEN" ->
        parse_zenith_grid(line, state)

      "SINEX CODE" ->
        parse_sinex_code(line, state)

      "VALID FROM" ->
        parse_valid(line, state, :valid_from)

      "VALID UNTIL" ->
        parse_valid(line, state, :valid_until)

      "START OF FREQUENCY" ->
        begin_frequency(line, state)

      "END OF FREQUENCY" ->
        finalize_frequency(state)

      "NORTH / EAST / UP" ->
        parse_pco(line, state)

      _ ->
        parse_pcv_row(line, state)
    end
  end

  defp parse_type_serial(line, state) do
    %{state | current_antenna: decode_antenna_header(line), current_frequency: nil}
  end

  defp parse_dazi(_line, %{current_antenna: nil} = state) do
    state
  end

  defp parse_dazi(line, %{current_antenna: current} = state) do
    case parse_floats_from_prefix(line) do
      [dazi | _] ->
        %{state | current_antenna: %{current | dazi_deg: dazi}}

      _ ->
        state
    end
  end

  defp parse_zenith_grid(_line, %{current_antenna: nil} = state), do: state

  defp parse_zenith_grid(line, %{current_antenna: current} = state) do
    case parse_floats_from_prefix(line) do
      [zen1, zen2, dzen] ->
        %{
          state
          | current_antenna: %{
              current
              | zenith_start_deg: zen1,
                zenith_end_deg: zen2,
                zenith_step_deg: dzen
            }
        }

      _ ->
        state
    end
  end

  defp parse_sinex_code(_line, %{current_antenna: nil} = state), do: state

  defp parse_sinex_code(line, %{current_antenna: current} = state) do
    code = String.trim(String.slice(line, 0, 60))

    if code == "" do
      state
    else
      %{state | current_antenna: %{current | sinex_code: code}}
    end
  end

  defp parse_valid(_line, %{current_antenna: nil} = state, _field), do: state

  defp parse_valid(line, %{current_antenna: current} = state, field) do
    case parse_floats_from_prefix(line) do
      [y, mo, d, h, mi, s | _] ->
        ndt =
          NaiveDateTime.new!(
            trunc(y),
            trunc(mo),
            trunc(d),
            trunc(h),
            trunc(mi),
            trunc(s)
          )

        %{state | current_antenna: Map.put(current, field, ndt)}

      _ ->
        state
    end
  end

  defp begin_antenna(state), do: %{state | current_antenna: nil, current_frequency: nil}

  defp decode_antenna_header(line) do
    type_serial = String.trim(String.slice(line, 0, 60))
    type_field = String.trim(String.slice(line, 0, 20))
    serial_field = String.trim(String.slice(line, 20, 20))

    kind =
      if Regex.match?(~r/^[A-Z][0-9]{2}$/, serial_field) do
        :satellite
      else
        :receiver
      end

    %Antenna{
      id: type_serial,
      kind: kind,
      type: type_field,
      serial: serial_field,
      dazi_deg: 0.0,
      zenith_start_deg: 0.0,
      zenith_end_deg: 0.0,
      zenith_step_deg: 0.0,
      sinex_code: nil,
      valid_from: nil,
      valid_until: nil,
      frequencies: %{}
    }
  end

  defp begin_frequency(_line, %{current_antenna: nil} = state), do: state

  defp begin_frequency(line, %{current_antenna: _current} = state) do
    frequency = String.trim(String.slice(line, 0, 20))

    %{state | current_frequency: %{frequency: frequency, phase: :pco, pco_m: nil, samples: []}}
  end

  defp parse_pco(_line, %{current_frequency: nil} = state), do: state

  defp parse_pco(line, %{current_frequency: %{phase: :pco}} = state) do
    case parse_floats_from_prefix(line) do
      [north, east, up] ->
        %{
          state
          | current_frequency: %{
              state.current_frequency
              | pco_m: {north / 1000.0, east / 1000.0, up / 1000.0},
                phase: :pcv
            }
        }

      _ ->
        state
    end
  end

  defp parse_pco(_line, state), do: state

  defp parse_pcv_row(line, %{current_frequency: %{phase: :pcv}} = state) do
    tokens = parse_tokens(line)

    case tokens do
      ["NOAZI" | values] ->
        add_pcv_values(nil, values, state)

      [first | values] ->
        case parse_float(first) do
          {:ok, azimuth} ->
            add_pcv_values(azimuth, values, state)

          :error ->
            state
        end

      _ ->
        state
    end
  end

  defp parse_pcv_row(_line, state), do: state

  defp add_pcv_values(
         azimuth_deg,
         values,
         %{current_antenna: current, current_frequency: current_freq} = state
       )
       when is_list(values) do
    grid_start = current.zenith_start_deg
    grid_step = current.zenith_step_deg

    samples =
      values
      |> Enum.with_index()
      |> Enum.reduce([], fn {value_text, index}, acc ->
        case parse_float(value_text) do
          {:ok, value} ->
            zenith = if grid_step == 0.0, do: grid_start, else: grid_start + grid_step * index

            sample = %{
              grid: if(azimuth_deg == nil, do: :noazi, else: :azi),
              azimuth_deg: azimuth_deg,
              zenith_deg: zenith,
              value_m: value / 1000.0
            }

            [sample | acc]

          :error ->
            acc
        end
      end)
      |> Enum.reverse()

    %{state | current_frequency: %{current_freq | samples: current_freq.samples ++ samples}}
  end

  defp finalize_frequency(%{current_frequency: nil} = state), do: state

  defp finalize_frequency(%{current_frequency: current_freq, current_antenna: current} = state) do
    frequency = %Frequency{
      frequency: current_freq.frequency,
      pco_m: current_freq.pco_m || {0.0, 0.0, 0.0},
      pcv_samples: current_freq.samples
    }

    %{
      state
      | current_antenna: %{
          current
          | frequencies: Map.put(current.frequencies, frequency.frequency, frequency)
        },
        current_frequency: nil
    }
  end

  defp finalize_antenna(%State{current_antenna: nil} = state), do: state

  defp finalize_antenna(%State{antennas: antennas, current_antenna: current} = state) do
    %{
      state
      | antennas: Map.put(antennas, current.id, current),
        current_antenna: nil,
        current_frequency: nil
    }
  end

  # --- interpolation -------------------------------------------------------

  defp interpolate_azimuth(azimuth_samples, azimuth_deg, zenith_deg) do
    azimuths = azimuth_samples |> Map.keys() |> Enum.sort()
    azimuth = normalize_azimuth(azimuth_deg)

    {low_deg, high_deg} = azimuth_bracket(azimuths, azimuth)

    low_samples = Map.fetch!(azimuth_samples, low_deg)
    high_samples = Map.fetch!(azimuth_samples, high_deg)

    low = low_deg
    high = if high_deg < low, do: high_deg + 360.0, else: high_deg
    target = if azimuth < low, do: azimuth + 360.0, else: azimuth

    low_value = interpolate(low_samples, zenith_deg)
    high_value = interpolate(high_samples, zenith_deg)

    if high == low do
      low_value
    else
      t = (target - low) / (high - low)
      low_value + (high_value - low_value) * t
    end
  end

  defp azimuth_bracket(azimuths, azimuth) do
    first = hd(azimuths)
    last = List.last(azimuths)

    low = Enum.find(Enum.reverse(azimuths), fn a -> a <= azimuth end)
    high = Enum.find(azimuths, fn a -> a >= azimuth end)

    case {low, high} do
      {nil, _} ->
        {last, first}

      {_, nil} ->
        {last, first}

      {low, high} ->
        {low, high}
    end
  end

  defp interpolate(samples, zenith_deg) do
    samples = Enum.sort_by(samples, fn {z, _} -> z end)
    first = hd(samples)
    last = List.last(samples)

    cond do
      zenith_deg <= elem(first, 0) ->
        elem(first, 1)

      zenith_deg >= elem(last, 0) ->
        elem(last, 1)

      true ->
        {low_zen, low_value} = lower_sample(samples, zenith_deg)
        {high_zen, high_value} = upper_sample(samples, zenith_deg)

        if high_zen == low_zen do
          low_value
        else
          t = (zenith_deg - low_zen) / (high_zen - low_zen)
          low_value + (high_value - low_value) * t
        end
    end
  end

  defp lower_sample(samples, zenith_deg) do
    samples
    |> Enum.reduce_while(nil, fn {z, v}, acc ->
      if z <= zenith_deg do
        {:cont, {z, v}}
      else
        {:halt, acc}
      end
    end)
  end

  defp upper_sample(samples, zenith_deg) do
    samples
    |> Enum.reduce_while(nil, fn {z, v}, acc ->
      if z >= zenith_deg do
        {:halt, {z, v}}
      else
        {:cont, acc}
      end
    end)
  end

  # --- formatting / helpers -----------------------------------------------

  defp tag(line), do: String.trim(String.slice(line, 60, 20))

  defp parse_tokens(line) do
    line
    |> String.trim()
    |> String.split(~r/\s+/, trim: true)
  end

  defp parse_floats_from_prefix(line) do
    line
    |> parse_tokens()
    |> Enum.map(&parse_float/1)
    |> Enum.take_while(&match?({:ok, _}, &1))
    |> Enum.map(fn {:ok, value} -> value end)
  end

  defp parse_float(token) do
    case Float.parse(token) do
      {value, _rest} -> {:ok, value}
      :error -> :error
    end
  end

  defp normalize_azimuth(azimuth_deg) when is_number(azimuth_deg) do
    wrapped = :math.fmod(azimuth_deg, 360.0)

    if wrapped < 0.0 do
      wrapped + 360.0
    else
      wrapped
    end
  end
end
