defmodule Orbis.Format.TLE do
  @moduledoc """
  Parse and encode Two-Line Element sets.

  TLE is the legacy fixed-width format for satellite orbital elements,
  designed for 80-column punch cards in the 1960s. Despite its age,
  it remains the most widely used format for distributing orbital data.

  ## Parsing

  The parser is liberal in what it accepts:
  - Trailing whitespace and extra characters are trimmed
  - Leading dots in floats (`.123` → `0.123`)
  - Spaces in numeric fields

  Checksum validation is performed and reported but does not prevent parsing.

  ## Examples

      {:ok, elements} = Orbis.Format.TLE.parse(line1, line2)
      {line1, line2} = Orbis.Format.TLE.encode(elements)
  """

  alias Orbis.Elements

  require Logger

  @microseconds_per_day 86_400 * 1_000_000

  @doc """
  Parse a two-line element set into an `%Orbis.Elements{}` struct.

  Returns `{:ok, elements}` or `{:error, reason}`.
  Logs a warning if checksums are invalid but still parses.

  ## Examples

      iex> {:ok, el} = Orbis.Format.TLE.parse(
      ...>   "1 25544U 98067A   18184.80969102  .00001614  00000-0  31745-4 0  9993",
      ...>   "2 25544  51.6414 295.8524 0003435 262.6267 204.2868 15.54005638121106"
      ...> )
      iex> el.catalog_number
      "25544"
      iex> el.inclination_deg
      51.6414

  """
  @spec parse(String.t(), String.t()) :: {:ok, Elements.t()} | {:error, String.t()}
  def parse(longstr1, longstr2) do
    with :ok <- validate_ascii(longstr1, longstr2),
         line1 = clean_line(longstr1, 69),
         line2 = clean_line(longstr2, 69),
         :ok <- validate_format(line1, line2),
         {:ok, fields} <- extract_fields(line1, line2) do
      validate_checksums(line1, line2)
      {:ok, build_elements(fields)}
    end
  end

  @doc """
  Encode an `%Orbis.Elements{}` struct as TLE-format strings.

  Returns `{line1, line2}` — two 69-character strings with valid checksums.
  Round-trips are character-exact for standard TLEs.

  ## Examples

      iex> l1 = "1 25544U 98067A   18184.80969102  .00001614  00000-0  31745-4 0  9993"
      iex> l2 = "2 25544  51.6414 295.8524 0003435 262.6267 204.2868 15.54005638121106"
      iex> {:ok, el} = Orbis.Format.TLE.parse(l1, l2)
      iex> {gen_l1, gen_l2} = Orbis.Format.TLE.encode(el)
      iex> gen_l1 == l1
      true
      iex> gen_l2 == l2
      true

  """
  @spec encode(Elements.t()) :: {String.t(), String.t()}
  def encode(%Elements{} = el) do
    validate_encode_fields!(el)
    {epochyr, epochdays} = epoch_to_tle(el.epoch)

    cat = String.pad_leading(String.trim(el.catalog_number || "0"), 5)
    cls = el.classification || "U"
    intl = String.pad_trailing(el.international_designator || "", 8)

    l1_body =
      "1 #{cat}#{cls} #{intl} #{fmt_epoch(epochyr, epochdays)} #{fmt_ndot(el.mean_motion_dot)} #{fmt_exp(el.mean_motion_double_dot)} #{fmt_exp(el.bstar)} #{el.ephemeris_type || 0} #{String.pad_leading(to_string(el.elset_number || 999), 4)}"

    line1 = pad_and_checksum(l1_body)

    l2_body =
      "2 #{cat} #{fmt_deg(el.inclination_deg, 8)} #{fmt_deg(el.raan_deg, 8)} #{fmt_ecc(el.eccentricity)} #{fmt_deg(el.arg_perigee_deg, 8)} #{fmt_deg(el.mean_anomaly_deg, 8)} #{fmt_mm(el.mean_motion)}#{String.pad_leading(to_string(el.rev_number || 0), 5)}"

    line2 = pad_and_checksum(l2_body)

    {line1, line2}
  end

  # -- Parsing internals --

  defp validate_ascii(line1, line2) do
    if String.to_charlist(line1) |> Enum.all?(&(&1 <= 127)) and
         String.to_charlist(line2) |> Enum.all?(&(&1 <= 127)) do
      :ok
    else
      {:error, "TLE lines contain non-ASCII characters"}
    end
  end

  defp clean_line(line, max_len) do
    line
    |> String.trim_trailing()
    |> then(fn l -> if byte_size(l) > max_len, do: String.slice(l, 0, max_len), else: l end)
  end

  defp validate_format(line1, line2) do
    with :ok <- validate_line(line1, "1", 64, line1_positions()),
         :ok <- validate_line(line2, "2", 68, line2_positions()) do
      if String.slice(line1, 2..6) == String.slice(line2, 2..6) do
        :ok
      else
        {:error, "Satellite numbers in lines 1 and 2 do not match"}
      end
    end
  end

  defp validate_line(line, prefix, min_len, positions) do
    cond do
      String.length(line) < min_len ->
        {:error, tle_format_error()}

      not String.starts_with?(line, prefix <> " ") ->
        {:error, tle_format_error()}

      not Enum.all?(positions, fn {pos, ch} -> String.at(line, pos) == ch end) ->
        {:error, tle_format_error()}

      true ->
        :ok
    end
  end

  defp line1_positions,
    do: [{8, " "}, {23, "."}, {32, " "}, {34, "."}, {43, " "}, {52, " "}, {61, " "}, {63, " "}]

  defp line2_positions,
    do: [
      {7, " "},
      {11, "."},
      {16, " "},
      {20, "."},
      {25, " "},
      {33, " "},
      {37, "."},
      {42, " "},
      {46, "."},
      {51, " "}
    ]

  defp validate_checksums(line1, line2) do
    check_one_checksum(line1, "line 1")
    check_one_checksum(line2, "line 2")
  end

  defp check_one_checksum(line, label) do
    if String.length(line) >= 69 do
      expected = String.to_integer(String.at(line, 68))
      computed = compute_checksum(line)

      if expected != computed do
        Logger.warning(
          "TLE #{label} checksum mismatch: expected #{expected}, computed #{computed}"
        )
      end
    end
  end

  defp compute_checksum(line) do
    line
    |> String.to_charlist()
    |> Enum.take(68)
    |> Enum.reduce(0, fn
      c, acc when c >= ?0 and c <= ?9 -> acc + (c - ?0)
      ?-, acc -> acc + 1
      _, acc -> acc
    end)
    |> rem(10)
  end

  defp extract_fields(line1, line2) do
    fields = %{
      catalog_number: String.slice(line1, 2..6) |> String.trim(),
      classification: String.at(line1, 7) || "U",
      intldesg: String.trim_trailing(String.slice(line1, 9..16)),
      two_digit_year: String.slice(line1, 18..19) |> String.trim() |> String.to_integer(),
      epochdays: String.slice(line1, 20..31) |> parse_float(),
      ndot: String.slice(line1, 33..42) |> parse_float(),
      nddot: parse_exp_field(line1, 44, 45..49, 50..51),
      bstar: parse_exp_field(line1, 53, 54..58, 59..60),
      ephtype: String.at(line1, 62) |> String.to_integer(),
      elnum: String.slice(line1, 64..67) |> String.trim() |> String.to_integer(),
      inclo: String.slice(line2, 8..15) |> parse_float(),
      nodeo: String.slice(line2, 17..24) |> parse_float(),
      ecco: ("0." <> String.replace(String.slice(line2, 26..32), " ", "0")) |> String.to_float(),
      argpo: String.slice(line2, 34..41) |> parse_float(),
      mo: String.slice(line2, 43..50) |> parse_float(),
      no_kozai: String.slice(line2, 52..62) |> parse_float(),
      revnum: String.slice(line2, 63..67) |> String.trim() |> String.to_integer()
    }

    {:ok, fields}
  rescue
    e -> {:error, "TLE parse error: #{Exception.message(e)}"}
  end

  defp parse_exp_field(line, sign_pos, mant_range, exp_range) do
    sign = if String.at(line, sign_pos) == "-", do: -1, else: 1
    mantissa = ("0." <> String.slice(line, mant_range)) |> String.trim() |> String.to_float()
    exp = String.slice(line, exp_range) |> String.trim() |> String.to_integer()
    sign * mantissa * :math.pow(10.0, exp)
  end

  defp build_elements(fields) do
    epoch_year =
      if fields.two_digit_year < 57,
        do: 2000 + fields.two_digit_year,
        else: 1900 + fields.two_digit_year

    %Elements{
      catalog_number: fields.catalog_number,
      classification: fields.classification,
      international_designator: fields.intldesg,
      epoch: calculate_epoch(epoch_year, fields.epochdays),
      mean_motion_dot: fields.ndot,
      mean_motion_double_dot: fields.nddot,
      bstar: fields.bstar,
      ephemeris_type: fields.ephtype,
      elset_number: fields.elnum,
      inclination_deg: fields.inclo,
      raan_deg: fields.nodeo,
      eccentricity: fields.ecco,
      arg_perigee_deg: fields.argpo,
      mean_anomaly_deg: fields.mo,
      mean_motion: fields.no_kozai,
      rev_number: fields.revnum
    }
  end

  defp calculate_epoch(year, epochdays) do
    days_from_jan1 = epochdays - 1
    whole_days = trunc(days_from_jan1)
    fractional_day = days_from_jan1 - whole_days

    start = DateTime.new!(Date.new!(year, 1, 1), Time.new!(0, 0, 0, 0), "Etc/UTC")
    with_days = DateTime.add(start, whole_days, :day)
    microseconds = round(fractional_day * @microseconds_per_day)
    DateTime.add(with_days, microseconds, :microsecond)
  end

  defp parse_float(str) do
    trimmed = String.trim(str)

    normalized =
      case trimmed do
        "+" <> rest -> rest
        _ -> trimmed
      end

    normalized =
      case normalized do
        "." <> _ -> "0" <> normalized
        "-." <> rest -> "-0." <> rest
        _ -> normalized
      end

    String.to_float(normalized)
  end

  defp validate_encode_fields!(el) do
    cat = String.trim(el.catalog_number || "0")

    if String.length(cat) > 5 do
      raise ArgumentError, "catalog_number #{inspect(cat)} exceeds 5 characters"
    end

    if (el.rev_number || 0) > 99_999 do
      raise ArgumentError, "rev_number #{el.rev_number} exceeds 5 digits"
    end

    if (el.elset_number || 0) > 9999 do
      raise ArgumentError, "elset_number #{el.elset_number} exceeds 4 digits"
    end

    :ok
  end

  defp tle_format_error do
    "TLE format error: line does not match the Two-Line Element fixed-width format"
  end

  # -- Encoding internals --

  defp epoch_to_tle(epoch) do
    year = epoch.year
    jan1 = DateTime.new!(Date.new!(year, 1, 1), Time.new!(0, 0, 0, 0), "Etc/UTC")
    diff_us = DateTime.diff(epoch, jan1, :microsecond)
    epochdays = 1.0 + diff_us / @microseconds_per_day
    {rem(year, 100), epochdays}
  end

  defp fmt_epoch(yr, days) do
    yr_str = String.pad_leading(to_string(yr), 2, "0")
    days_str = :erlang.float_to_binary(days, decimals: 8)
    yr_str <> String.pad_leading(days_str, 12, "0")
  end

  defp fmt_ndot(val) do
    # ndot: 10 chars (cols 33-42), format " .XXXXXXXX" or "-.XXXXXXXX"
    sign = if val < 0, do: "-", else: " "
    str = :erlang.float_to_binary(abs(val), decimals: 8)
    str = String.replace_prefix(str, "0", "")
    sign <> String.pad_leading(str, 9)
  end

  # TLE "assumed decimal" format for nddot and bstar:
  # [sign][5 mantissa digits][exp sign][1 exp digit]
  # Value = 0.XXXXX * 10^N, e.g., " 31745-4" means 0.31745e-4 = 3.1745e-5
  defp fmt_exp(val) when val == 0.0, do: " 00000-0"

  defp fmt_exp(val) do
    sign = if val < 0, do: "-", else: " "
    av = abs(val)
    raw_exp = Float.floor(:math.log10(av)) |> trunc()
    # Assumed decimal: express as 0.XXXXX * 10^(raw_exp + 1)
    exp = raw_exp + 1
    mantissa = av / :math.pow(10.0, exp)
    mant_str = :erlang.float_to_binary(mantissa, decimals: 5) |> String.slice(2, 5)
    exp_sign = if exp >= 0, do: "+", else: "-"
    sign <> mant_str <> exp_sign <> to_string(abs(exp))
  end

  defp fmt_ecc(ecc) do
    :erlang.float_to_binary(ecc, decimals: 7)
    |> String.replace_prefix("0.", "")
    |> String.pad_leading(7, "0")
  end

  defp fmt_deg(val, width) do
    :erlang.float_to_binary(val, decimals: 4)
    |> String.pad_leading(width)
  end

  defp fmt_mm(val) do
    :erlang.float_to_binary(val, decimals: 8)
    |> String.pad_leading(11)
  end

  defp pad_and_checksum(body) do
    padded = String.pad_trailing(String.slice(body, 0, 68), 68)
    padded <> to_string(compute_checksum(padded))
  end
end
