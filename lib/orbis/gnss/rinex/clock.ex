defmodule Orbis.GNSS.RINEX.Clock do
  @moduledoc """
  RINEX clock (`.CLK`) reader for satellite clock-bias records.

  Precise clock products are distributed as RINEX clock files alongside the SP3
  orbit. The SP3 orbit carries satellite clocks too, but only at the SP3 epoch
  spacing (15 minutes for IGS final), whereas the companion `.CLK` file carries
  the same clocks at a much finer cadence (30 seconds for IGS final). Linearly
  interpolating a 15-minute clock across the gap is a metre-level error on the
  faster satellite oscillators; the 30s clock removes almost all of it.

  This reader parses the `AS` (satellite) records of a RINEX clock file
  (versions 2 and 3) into a per-satellite, time-ordered series of clock biases in
  seconds, and interpolates linearly between the two bracketing records at a
  requested epoch:

      AS G05  2026 05 13 00 00  0.000000  2   -2.329120317895e-04  4.4959e-11

  The fields are the record type (`AS`), the satellite id, the epoch
  (year month day hour minute second), the value count, then the clock bias in
  seconds and an optional bias sigma. `AR` (receiver) records are ignored.

  Use `clock_s/3` to read a satellite's interpolated clock bias (seconds) at an
  epoch, matching the convention of `Orbis.GNSS.SP3.State.clock_s`.
  """

  @enforce_keys [:series]
  defstruct [:series]

  @type t :: %__MODULE__{series: %{String.t() => [{float(), float()}]}}

  @doc """
  Load a RINEX clock file, raising on error.
  """
  @spec load!(String.t()) :: t()
  def load!(path) when is_binary(path) do
    case load(path) do
      {:ok, clock} ->
        clock

      {:error, reason} ->
        raise ArgumentError, "could not load RINEX clock #{path}: #{inspect(reason)}"
    end
  end

  @doc """
  Load a RINEX clock file.

  Returns `{:ok, %Orbis.GNSS.RINEX.Clock{}}` or `{:error, reason}`. The series is
  per-satellite, sorted ascending by GPS-seconds time tag, with each entry
  `{gps_seconds, clock_bias_s}`.
  """
  @spec load(String.t()) :: {:ok, t()} | {:error, term()}
  def load(path) when is_binary(path) do
    case File.read(path) do
      {:ok, contents} -> {:ok, parse(contents)}
      {:error, reason} -> {:error, reason}
    end
  end

  @doc """
  Interpolated satellite clock bias in seconds at `epoch`.

  Returns `{:ok, bias_s}` when the satellite has records bracketing the epoch (or
  an exact-match record), `{:error, :no_clock}` when the satellite is unknown or
  the epoch lies outside its record span. Linear interpolation between the two
  nearest records; no extrapolation past the first/last record.
  """
  @spec clock_s(t(), String.t(), NaiveDateTime.t()) :: {:ok, float()} | {:error, :no_clock}
  def clock_s(%__MODULE__{series: series}, satellite_id, %NaiveDateTime{} = epoch)
      when is_binary(satellite_id) do
    case Map.get(series, satellite_id) do
      nil -> {:error, :no_clock}
      records -> interpolate(records, naive_to_gps_seconds(epoch))
    end
  end

  # --- parsing -------------------------------------------------------------

  defp parse(contents) do
    body =
      contents
      |> String.split(["\r\n", "\n"], trim: false)
      |> drop_header()

    series =
      body
      |> Enum.flat_map(&parse_record/1)
      |> Enum.group_by(fn {sat, _t, _bias} -> sat end, fn {_sat, t, bias} -> {t, bias} end)
      |> Map.new(fn {sat, points} ->
        {sat, points |> Enum.sort_by(&elem(&1, 0)) |> dedup_by_time()}
      end)

    %__MODULE__{series: series}
  end

  defp drop_header(lines) do
    Enum.drop_while(lines, fn line -> not String.contains?(line, "END OF HEADER") end)
    |> case do
      [] -> lines
      [_end_of_header | rest] -> rest
    end
  end

  defp parse_record(line) do
    case String.split(line, ~r/\s+/, trim: true) do
      ["AS", sat, y, mo, d, h, mi, s, _count, bias | _rest] ->
        with {:ok, t} <- record_time(y, mo, d, h, mi, s),
             {bias_s, ""} <- Float.parse(bias) do
          [{sat, t, bias_s}]
        else
          _ -> []
        end

      _ ->
        []
    end
  end

  defp record_time(y, mo, d, h, mi, s) do
    with {yi, ""} <- Integer.parse(y),
         {moi, ""} <- Integer.parse(mo),
         {di, ""} <- Integer.parse(d),
         {hi, ""} <- Integer.parse(h),
         {mii, ""} <- Integer.parse(mi),
         {sf, ""} <- Float.parse(s),
         ws = trunc(sf),
         us = round((sf - ws) * 1_000_000),
         {:ok, date} <- Date.new(yi, moi, di),
         {:ok, time} <- Time.new(hi, mii, ws, {us, 6}),
         {:ok, ndt} <- NaiveDateTime.new(date, time) do
      {:ok, naive_to_gps_seconds(ndt)}
    else
      _ -> :error
    end
  end

  # Records at the same time tag (should not happen in a well-formed file) keep
  # the last; the series stays strictly monotone in time for interpolation.
  defp dedup_by_time(points) do
    points
    |> Enum.reduce([], fn
      {t, _bias} = point, [{t, _prev} | rest] -> [point | rest]
      point, acc -> [point | acc]
    end)
    |> Enum.reverse()
  end

  # --- interpolation -------------------------------------------------------

  defp interpolate(records, t) do
    case bracket(records, t) do
      {:exact, bias} -> {:ok, bias}
      {:between, {t0, b0}, {t1, b1}} -> {:ok, b0 + (b1 - b0) * (t - t0) / (t1 - t0)}
      :out_of_range -> {:error, :no_clock}
    end
  end

  defp bracket(records, t) do
    Enum.reduce_while(records, {:start, nil}, fn {ti, bi} = point, {_state, prev} ->
      cond do
        ti == t -> {:halt, {:exact, bi}}
        ti > t and is_nil(prev) -> {:halt, :out_of_range}
        ti > t -> {:halt, {:between, prev, point}}
        true -> {:cont, {:scan, point}}
      end
    end)
    |> case do
      {:scan, _last} -> :out_of_range
      {:start, _} -> :out_of_range
      other -> other
    end
  end

  defp naive_to_gps_seconds(%NaiveDateTime{} = ndt) do
    # GPS clock files are in GPS time; the absolute epoch zero is arbitrary here
    # since interpolation only uses differences. Use NaiveDateTime gregorian
    # seconds (microsecond precision) as a monotone time tag.
    NaiveDateTime.diff(ndt, ~N[1980-01-06 00:00:00], :microsecond) / 1_000_000.0
  end
end
