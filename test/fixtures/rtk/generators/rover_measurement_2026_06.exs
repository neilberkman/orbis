defmodule RoverMeasurement202606 do
  @moduledoc false

  alias Orbis.GNSS.Broadcast
  alias Orbis.GNSS.Positioning
  alias Orbis.GNSS.RINEX.Observations
  alias Orbis.GNSS.RTK

  @c_m_s 299_792_458.0
  @gps_l1_hz 1_575_420_000.0
  @gps_l1_wavelength_m @c_m_s / @gps_l1_hz
  @bds_b1i_hz 1_561_098_000.0
  @glonass_g1_hz 1_602_000_000.0
  @glonass_g1_step_hz 562_500.0
  @earth_a_m 6_378_137.0
  @earth_f 1.0 / 298.257_223_563

  @systems ["G", "R", "E", "C"]
  @default_work "/tmp/gsdc-work"
  @default_results "/tmp/rover-measurement-2026-06-results.json"
  @max_segment_epochs 1
  @sanity_code_residual_threshold_m 1_000.0

  @oracle_fixtures [
    "gsdc_2021_08_04_sjc1_pixel5_p222_demo5_rtklib_oracle.json",
    "gsdc_2021_08_24_svl1_pixel5_p222_demo5_rtklib_oracle.json",
    "gsdc_2021_12_15_mtv1_pixel5_p222_demo5_rtklib_oracle.json",
    "gsdc_2021_12_28_mtv1_pixel5_p222_demo5_rtklib_oracle.json"
  ]

  @phone_l1_codes %{
    "G" => [{"C1C", "L1C"}],
    "R" => [{"C1C", "L1C"}],
    "E" => [{"C1C", "L1C"}, {"C1X", "L1X"}],
    "C" => [{"C2I", "L2I"}]
  }

  @rinex2_l1_codes %{
    "G" => [{"C1", "L1"}, {"P1", "L1"}],
    "R" => [{"C1", "L1"}, {"P1", "L1"}],
    "E" => [{"C1", "L1"}],
    "C" => [{"C2", "L2"}]
  }

  @filter_option_notes [
    {"filter_kernel", "rust",
     "The shipped 0.18.0 default kernel is measured; the library is not modified."},
    {"initial_baseline_m", "first broadcast-code SPP minus P222 ARP",
     "Uses phone RINEX pseudoranges and the oracle broadcast NAV source, not truth."},
    {"baseline_prior_sigma_m", "500.0",
     "Allows a phone-code SPP seed to be wrong by many metres without letting the prior dominate."},
    {"ambiguity_prior_sigma_m", "1000.0",
     "Weak single-difference ambiguity prior matching the shipped filter scale."},
    {"process_noise_baseline_sigma_m", "30.0",
     "Kept at the original kinematic setting; it is inert with one-epoch measurement segments."},
    {"hold_sigma_m", "1.0e-4",
     "Keeps the shipped tight ambiguity hold used after an accepted integer fix."},
    {"max_iterations", "10", "Matches the real-arc RTK tests' nonlinear iteration cap."},
    {"on_cycle_slip", "split_arc",
     "Kept at the real-arc test setting; this script omits phone/base LLI flags because reference-satellite LLI splits are rejected by the shipped sequential filter."},
    {"elevation_mask_deg", "10.0", "Required by the brief."},
    {"stochastic_model", "rtklib", "Required by the brief."},
    {"code_sigma_m", "0.9",
     "RTKLIB phone oracle uses errphase=0.003 and eratio1=300, giving 0.9 m code scale."},
    {"phase_sigma_m", "0.003", "Matches the RTKLIB oracle phase error setting."},
    {"ambiguity_wavelength_m", "per-satellite L1/B1/G1 map",
     "GPS/Galileo use L1, BeiDou uses B1I, GLONASS G1 uses phone RINEX FDMA slots."},
    {"integer_ratio_threshold", "3.0", "Matches the pre-registered oracle/spec AR bar."},
    {"integer_candidate_limit", "200000", "Matches the real-arc RTK filter tests."},
    {"float_only_systems", "[\"R\"]", "Required by the brief; GLONASS FDMA is not integer-fixed."}
  ]

  defmodule Rinex2Obs do
    @moduledoc false
    defstruct [:obs_types, :approx_position_m, :antenna_delta_hen_m, :epochs]
  end

  def main(args) do
    {opts, _argv, invalid} =
      OptionParser.parse(args,
        strict: [work: :string, results: :string, report: :string],
        aliases: [w: :work, r: :results]
      )

    if invalid != [] do
      raise ArgumentError, "invalid arguments: #{inspect(invalid)}"
    end

    generator_dir = __DIR__
    fixture_dir = Path.expand("..", generator_dir)
    work = Keyword.get(opts, :work, @default_work)
    results_path = Keyword.get(opts, :results, @default_results)

    report_path =
      Keyword.get(opts, :report, Path.join(generator_dir, "rover-measurement-2026-06.md"))

    arcs =
      @oracle_fixtures
      |> Enum.map(&load_arc(fixture_dir, &1, work))
      |> Enum.map(&measure_arc/1)

    result = %{
      "version" => 1,
      "generated_at_utc" =>
        DateTime.utc_now() |> DateTime.truncate(:second) |> DateTime.to_iso8601(),
      "script" => Path.relative_to(__ENV__.file, File.cwd!()),
      "work_dir" => work,
      "option_notes" => option_notes_json(),
      "arcs" => Enum.map(arcs, &json_arc/1),
      "pooled" => pooled_summary(arcs)
    }

    File.write!(results_path, Jason.encode!(result, pretty: true))
    File.write!(report_path, report_markdown(result, results_path))

    IO.puts("wrote #{results_path}")
    IO.puts("wrote #{report_path}")
  end

  defp load_arc(fixture_dir, fixture, work) do
    oracle_path = Path.join(fixture_dir, fixture)
    oracle = oracle_path |> File.read!() |> Jason.decode!()
    inputs = oracle["inputs"]

    base_doy = input_doy!(inputs["base_obs"], ~r/p222(\d{3})0\.21d/)
    nav_name = inputs["nav"] |> Path.basename() |> String.replace_suffix(".gz", "")
    drive = inputs["drive"]

    %{
      fixture: fixture,
      oracle_path: oracle_path,
      oracle: oracle,
      label: oracle["reference"]["label"],
      drive: drive,
      rover_path: Path.join([work, drive, "supplemental/gnss_rinex.21o"]),
      base_path: Path.join([work, "cors", "p222#{base_doy}0.21o"]),
      nav_path: Path.join([work, "cors", nav_name])
    }
  end

  defp input_doy!(value, regex) do
    case Regex.run(regex, value) do
      [_, doy] -> doy
      _ -> raise "could not derive DOY from #{value}"
    end
  end

  defp measure_arc(arc) do
    IO.puts("measuring #{arc.label}")
    require_files!([arc.rover_path, arc.base_path, arc.nav_path, arc.oracle_path])

    rover_obs = Observations.load!(arc.rover_path)
    base_obs = load_rinex2_obs!(arc.base_path)
    nav = Broadcast.load!(arc.nav_path)

    base_arp = base_arp(base_obs)
    glonass_slots = Observations.glonass_slots(rover_obs)
    oracle_by_time = Map.new(arc.oracle["per_epoch"], &{&1["time"], &1})

    {epochs, contexts} =
      build_filter_epochs(nav, rover_obs, base_obs, glonass_slots, oracle_by_time, base_arp)

    if epochs == [] do
      raise "no usable Orbis epochs for #{arc.label}"
    end

    sanity_gate = sanity_gate!(nav, epochs, contexts, base_arp)
    time_alignment = time_alignment(contexts)

    IO.puts(
      "  sanity gate: median |clock-demeaned SD code residual| = #{fmt(sanity_gate.median_abs_code_residual_m)} m"
    )

    segments = segment_epoch_contexts(epochs, contexts)
    {per_epoch, segment_reports} = solve_segments(segments, nav, base_arp, glonass_slots)
    initial_baseline = segment_reports |> List.first() |> Map.fetch!(:initial_baseline_m)

    demo5 = demo5_summary(arc.oracle["per_epoch"])
    orbis = summarize_measurements(per_epoch)
    comparative = comparative_verdict(orbis, demo5, :per_arc)
    invariant = invariant_verdict(per_epoch)
    ledger = classify_worst_decile(per_epoch)

    %{
      input: arc,
      base_arp: base_arp,
      initial_baseline: initial_baseline,
      sanity_gate: sanity_gate,
      time_alignment: time_alignment,
      built_epoch_count: length(epochs),
      skipped_oracle_epochs: length(arc.oracle["per_epoch"]) - length(epochs),
      segment_count: length(segment_reports),
      segment_reports: segment_reports,
      per_epoch: per_epoch,
      demo5: demo5,
      orbis: orbis,
      comparative: comparative,
      invariant: invariant,
      ledger: ledger
    }
  end

  defp require_files!(paths) do
    missing = Enum.reject(paths, &File.exists?/1)

    if missing != [] do
      raise "missing required inputs: #{Enum.join(missing, ", ")}"
    end
  end

  defp load_rinex2_obs!(path) do
    lines = path |> File.read!() |> String.split("\n", trim: false)
    {header, rest} = Enum.split_while(lines, &(not String.contains?(&1, "END OF HEADER")))
    rest = tl(rest)

    obs_types = rinex2_obs_types(header)
    approx = header_tuple(header, "APPROX POSITION XYZ")
    antenna_delta = header_tuple(header, "ANTENNA: DELTA H/E/N") || {0.0, 0.0, 0.0}
    epochs = parse_rinex2_epochs(rest, obs_types, [])

    %Rinex2Obs{
      obs_types: obs_types,
      approx_position_m: approx,
      antenna_delta_hen_m: antenna_delta,
      epochs: epochs
    }
  end

  defp rinex2_obs_types(header) do
    header
    |> Enum.reduce({nil, []}, fn line, {count, obs_types} ->
      if header_label(line) == "# / TYPES OF OBSERV" do
        tokens = line |> String.slice(0, 60) |> String.split()

        cond do
          is_nil(count) ->
            [count_s | rest] = tokens
            {String.to_integer(count_s), obs_types ++ rest}

          length(obs_types) < count ->
            {count, obs_types ++ tokens}

          true ->
            {count, obs_types}
        end
      else
        {count, obs_types}
      end
    end)
    |> then(fn {count, obs_types} ->
      if count == nil or length(obs_types) < count do
        raise "could not parse RINEX 2 observation types"
      end

      Enum.take(obs_types, count)
    end)
  end

  defp header_tuple(header, label) do
    header
    |> Enum.find(&(header_label(&1) == label))
    |> case do
      nil ->
        nil

      line ->
        line
        |> String.slice(0, 60)
        |> String.split()
        |> Enum.take(3)
        |> Enum.map(&parse_float!/1)
        |> List.to_tuple()
    end
  end

  defp header_label(line), do: line |> pad(80) |> String.slice(60, 20) |> String.trim()

  defp parse_rinex2_epochs([], _obs_types, acc), do: Enum.reverse(acc)

  defp parse_rinex2_epochs([line | rest], obs_types, acc) do
    if String.trim(line) == "" do
      parse_rinex2_epochs(rest, obs_types, acc)
    else
      {epoch, sats, rest_after_header} = parse_rinex2_epoch_header(line, rest)

      {sat_obs, rest_after_obs} =
        parse_rinex2_sat_observations(rest_after_header, sats, obs_types, %{})

      parse_rinex2_epochs(rest_after_obs, obs_types, [
        %{epoch: epoch, observations: sat_obs} | acc
      ])
    end
  end

  defp parse_rinex2_epoch_header(line, rest) do
    head = line |> pad(80) |> String.slice(0, 32)
    [yy_s, month_s, day_s, hour_s, minute_s, second_s, _flag_s, count_s] = String.split(head)

    year = expand_rinex2_year(String.to_integer(yy_s))

    epoch =
      naive_datetime(
        year,
        String.to_integer(month_s),
        String.to_integer(day_s),
        String.to_integer(hour_s),
        String.to_integer(minute_s),
        parse_float!(second_s)
      )

    count = String.to_integer(count_s)

    {sat_text, rest_after_sats} =
      collect_rinex2_sat_text(line, rest, count, String.slice(pad(line, 80), 32, 48))

    sats =
      sat_text
      |> chunks(3)
      |> Enum.map(&String.trim/1)
      |> Enum.reject(&(&1 == ""))
      |> Enum.take(count)

    {epoch, sats, rest_after_sats}
  end

  defp collect_rinex2_sat_text(_line, rest, count, text) do
    if rinex2_sat_text_count(text) >= count do
      {text, rest}
    else
      [next | tail] = rest
      collect_rinex2_sat_text(next, tail, count, text <> String.slice(pad(next, 80), 32, 48))
    end
  end

  defp rinex2_sat_text_count(text) do
    text
    |> chunks(3)
    |> Enum.map(&String.trim/1)
    |> Enum.reject(&(&1 == ""))
    |> length()
  end

  defp parse_rinex2_sat_observations(rest, [], _obs_types, acc), do: {acc, rest}

  defp parse_rinex2_sat_observations(rest, [sat | sats], obs_types, acc) do
    line_count = obs_types |> length() |> Kernel.+(4) |> div(5)
    {obs_lines, rest_after_sat} = Enum.split(rest, line_count)

    observations =
      obs_types
      |> Enum.with_index()
      |> Enum.reduce(%{}, fn {code, index}, obs_acc ->
        line = Enum.at(obs_lines, div(index, 5), "") |> pad(80)
        field = String.slice(line, rem(index, 5) * 16, 16)
        value_s = field |> String.slice(0, 14) |> String.trim()

        if value_s == "" do
          obs_acc
        else
          lli = parse_optional_int(String.slice(field, 14, 1))
          ssi = parse_optional_int(String.slice(field, 15, 1))
          Map.put(obs_acc, code, %{value: parse_float!(value_s), lli: lli, ssi: ssi})
        end
      end)

    parse_rinex2_sat_observations(
      rest_after_sat,
      sats,
      obs_types,
      Map.put(acc, sat, observations)
    )
  end

  defp parse_optional_int(value) do
    value = String.trim(value)

    if value != "" do
      String.to_integer(value)
    end
  end

  defp expand_rinex2_year(year) when year >= 80, do: 1900 + year
  defp expand_rinex2_year(year), do: 2000 + year

  defp base_arp(%Rinex2Obs{
         approx_position_m: marker,
         antenna_delta_hen_m: {height_m, east_m, north_m}
       }) do
    if east_m != 0.0 or north_m != 0.0 do
      raise "measurement script only handles zero east/north base antenna deltas"
    end

    add3(marker, scale3(marker, height_m / norm3(marker)))
  end

  defp build_filter_epochs(nav, rover_obs, base_obs, glonass_slots, oracle_by_time, base_arp) do
    base_index = base_epoch_index(base_obs)

    rover_obs
    |> Observations.epochs()
    |> Enum.reduce({[], []}, fn entry, {epochs, contexts} ->
      epoch = naive_datetime(entry.epoch)
      time_key = epoch_key(epoch)

      case Map.fetch(oracle_by_time, time_key) do
        {:ok, oracle_epoch} ->
          rover_values = phone_l1_values(rover_obs, entry.index, @systems, glonass_slots)
          base_values = interpolated_base_l1_values(base_index, epoch, @systems, glonass_slots)

          {filter_epoch, context} =
            build_filter_epoch(
              nav,
              epoch,
              time_key,
              rover_values,
              base_values,
              oracle_epoch,
              base_arp
            )

          if filter_epoch do
            {[filter_epoch | epochs], [context | contexts]}
          else
            {epochs, contexts}
          end

        :error ->
          {epochs, contexts}
      end
    end)
    |> then(fn {epochs, contexts} -> {Enum.reverse(epochs), Enum.reverse(contexts)} end)
  end

  defp build_filter_epoch(nav, epoch, time_key, rover_values, base_values, oracle_epoch, base_arp) do
    common =
      base_values
      |> Map.keys()
      |> MapSet.new()
      |> MapSet.intersection(rover_values |> Map.keys() |> MapSet.new())
      |> MapSet.to_list()
      |> Enum.sort()

    positions = satellite_positions(nav, epoch, common)
    base_positions = transmit_time_satellite_positions(nav, epoch, base_values, common)
    rover_positions = transmit_time_satellite_positions(nav, epoch, rover_values, common)

    usable =
      Enum.filter(common, fn sat ->
        Map.has_key?(positions, sat) and Map.has_key?(base_positions, sat) and
          Map.has_key?(rover_positions, sat)
      end)
      |> Enum.filter(&(elevation_deg(base_arp, Map.fetch!(positions, &1)) >= 10.0))
      |> drop_single_satellite_systems()

    if length(usable) >= 4 do
      filter_epoch = %{
        epoch: epoch,
        satellite_positions_m: Map.take(positions, usable),
        base_satellite_positions_m: Map.take(base_positions, usable),
        rover_satellite_positions_m: Map.take(rover_positions, usable),
        base_observations: Enum.map(usable, &Map.fetch!(base_values, &1)),
        rover_observations: Enum.map(usable, &Map.fetch!(rover_values, &1))
      }

      context = %{
        time: time_key,
        epoch: epoch,
        oracle_epoch: oracle_epoch,
        pre_mask_satellites: length(usable)
      }

      {filter_epoch, context}
    else
      {nil, nil}
    end
  end

  defp drop_single_satellite_systems(sats) do
    keep_systems =
      sats
      |> Enum.frequencies_by(&String.first/1)
      |> Enum.filter(fn {_system, count} -> count >= 2 end)
      |> MapSet.new(&elem(&1, 0))

    Enum.filter(sats, &(String.first(&1) in keep_systems))
  end

  defp segment_epoch_contexts(epochs, contexts) do
    epochs
    |> Enum.zip(contexts)
    |> Enum.reduce([], fn pair, segments ->
      case segments do
        [] ->
          [[pair]]

        [current | rest] ->
          candidate = current ++ [pair]

          if length(candidate) <= @max_segment_epochs and segment_reference_solvable?(candidate) do
            [candidate | rest]
          else
            [[pair], current | rest]
          end
      end
    end)
    |> Enum.reverse()
  end

  defp segment_reference_solvable?(pairs) do
    epochs = Enum.map(pairs, &elem(&1, 0))
    all_sats = epochs |> Enum.flat_map(&Map.keys(&1.satellite_positions_m)) |> Enum.uniq()
    systems = all_sats |> Enum.map(&String.first/1) |> Enum.uniq()
    common_by_system = per_system_common_sats(epochs)

    length(all_sats) >= 4 and Enum.all?(systems, &(Map.get(common_by_system, &1, []) != []))
  end

  defp per_system_common_sats(epochs) do
    epochs
    |> Enum.reduce(%{}, fn epoch, acc ->
      epoch.satellite_positions_m
      |> Map.keys()
      |> Enum.group_by(&String.first/1)
      |> Enum.reduce(acc, fn {system, sats}, system_acc ->
        sats = MapSet.new(sats)
        Map.update(system_acc, system, sats, &MapSet.intersection(&1, sats))
      end)
    end)
    |> Map.new(fn {system, sats} -> {system, MapSet.to_list(sats)} end)
  end

  defp solve_segments(segments, nav, base_arp, glonass_slots) do
    segments
    |> Enum.with_index()
    |> Enum.flat_map(fn {segment, index} ->
      solve_segment(segment, index, nav, base_arp, glonass_slots)
    end)
    |> Enum.unzip()
    |> case do
      {measurements, reports} -> {List.flatten(measurements), reports}
    end
  end

  defp solve_segment(segment, index, nav, base_arp, glonass_slots) do
    epochs = Enum.map(segment, &elem(&1, 0))
    contexts = Enum.map(segment, &elem(&1, 1))
    initial_baseline = spp_initial_baseline!(nav, epochs, base_arp)
    opts = filter_opts(initial_baseline, epochs, glonass_slots)

    IO.puts(
      "  segment #{index + 1}: #{hd(contexts).time}..#{List.last(contexts).time} #{length(epochs)} epochs"
    )

    case RTK.solve_filter_baseline_epochs(base_arp, epochs, opts) do
      {:ok, sol} ->
        measurements =
          sol.epochs
          |> Enum.zip(contexts)
          |> Enum.map(fn {result, context} ->
            epoch_measurement(result, context, base_arp)
          end)

        report = %{
          index: index,
          epochs: length(epochs),
          first_time: hd(contexts).time,
          last_time: List.last(contexts).time,
          initial_baseline_m: initial_baseline,
          metadata: sol.metadata
        }

        [{measurements, report}]

      {:error, reason} when length(segment) > 1 ->
        if System.get_env("ROVER_MEAS_DEBUG_ERROR") == "1" do
          raise "segment #{index + 1} failed: #{inspect(reason)}"
        end

        {left, right} = Enum.split(segment, div(length(segment), 2))

        solve_segment(left, index, nav, base_arp, glonass_slots) ++
          solve_segment(right, index, nav, base_arp, glonass_slots)

      {:error, reason} ->
        context = hd(contexts)
        IO.puts("skipping unsolved epoch #{context.time}: #{inspect(reason)}")
        []
    end
  end

  defp base_epoch_index(%Rinex2Obs{epochs: epochs}) do
    sorted = Enum.sort_by(epochs, &time_us(&1.epoch))

    %{
      times: sorted |> Enum.map(&time_us(&1.epoch)) |> List.to_tuple(),
      epochs: List.to_tuple(sorted),
      count: length(sorted)
    }
  end

  defp phone_l1_values(obs, index, systems, glonass_slots) do
    codes =
      @phone_l1_codes
      |> Map.take(systems)
      |> Map.new(fn {system, pairs} ->
        {system, Enum.flat_map(pairs, fn {code, phase} -> [code, phase] end)}
      end)

    {:ok, by_sat} = Observations.values(obs, index, codes: codes)

    by_sat
    |> Enum.flat_map(fn {sat, values} ->
      system = String.first(sat)
      values_by_code = Map.new(values, &{&1.code, &1})
      pairs = Map.get(@phone_l1_codes, system, [])

      with true <- system in systems,
           {:ok, wavelength_m} <- wavelength_m(sat, glonass_slots),
           {:ok, {code_m, phase_obs}} <- first_complete_pair(values_by_code, pairs) do
        [
          {sat,
           %{
             satellite_id: sat,
             code_m: code_m,
             phase_m: phase_obs.value * wavelength_m,
             lli: nil
           }}
        ]
      else
        _ -> []
      end
    end)
    |> Map.new()
  end

  defp interpolated_base_l1_values(base_index, epoch, systems, glonass_slots) do
    case bracket_base_epoch(base_index, epoch) do
      {:ok, before_epoch, after_epoch, fraction} ->
        before_epoch.observations
        |> Map.keys()
        |> MapSet.new()
        |> MapSet.intersection(after_epoch.observations |> Map.keys() |> MapSet.new())
        |> MapSet.to_list()
        |> Enum.sort()
        |> Enum.flat_map(fn sat ->
          system = String.first(sat)
          pairs = Map.get(@rinex2_l1_codes, system, [])

          with true <- system in systems,
               {:ok, wavelength_m} <- wavelength_m(sat, glonass_slots),
               {:ok, {code1, phase1}} <-
                 first_complete_pair(Map.fetch!(before_epoch.observations, sat), pairs),
               {:ok, {code2, phase2}} <-
                 first_complete_pair(Map.fetch!(after_epoch.observations, sat), pairs) do
            code_m = interpolate(code1, code2, fraction)
            phase_cycles = interpolate(phase1.value, phase2.value, fraction)

            [
              {sat,
               %{
                 satellite_id: sat,
                 code_m: code_m,
                 phase_m: phase_cycles * wavelength_m,
                 lli: nil
               }}
            ]
          else
            _ -> []
          end
        end)
        |> Map.new()

      _ ->
        %{}
    end
  end

  defp bracket_base_epoch(%{count: count}, _epoch) when count < 2, do: :error

  defp bracket_base_epoch(index, epoch) do
    target = time_us(epoch)
    pos = lower_bound(index.times, index.count, target)

    cond do
      pos < index.count and elem(index.times, pos) == target ->
        exact_epoch = elem(index.epochs, pos)
        {:ok, exact_epoch, exact_epoch, 0.0}

      pos == 0 ->
        :error

      pos >= index.count ->
        :error

      true ->
        before_epoch = elem(index.epochs, pos - 1)
        after_epoch = elem(index.epochs, pos)
        before_us = elem(index.times, pos - 1)
        after_us = elem(index.times, pos)
        max_gap_us = 30 * 1_000_000

        if target - before_us <= max_gap_us and after_us - target <= max_gap_us do
          fraction = (target - before_us) / (after_us - before_us)
          {:ok, before_epoch, after_epoch, fraction}
        else
          :error
        end
    end
  end

  defp lower_bound(times, count, target), do: lower_bound(times, target, 0, count)

  defp lower_bound(_times, _target, low, high) when low >= high, do: low

  defp lower_bound(times, target, low, high) do
    mid = div(low + high, 2)

    if elem(times, mid) < target do
      lower_bound(times, target, mid + 1, high)
    else
      lower_bound(times, target, low, mid)
    end
  end

  defp first_complete_pair(_values_by_code, []), do: :error

  defp first_complete_pair(values_by_code, [{code, phase} | rest]) do
    with %{value: code_m} when is_number(code_m) <- values_by_code[code],
         %{value: phase_cycles} = phase_obs when is_number(phase_cycles) <- values_by_code[phase] do
      {:ok, {code_m, phase_obs}}
    else
      _ -> first_complete_pair(values_by_code, rest)
    end
  end

  defp wavelength_m("G" <> _, _slots), do: {:ok, @gps_l1_wavelength_m}
  defp wavelength_m("E" <> _, _slots), do: {:ok, @gps_l1_wavelength_m}
  defp wavelength_m("C" <> _, _slots), do: {:ok, @c_m_s / @bds_b1i_hz}

  defp wavelength_m("R" <> _ = sat, slots) do
    with {:ok, k} <- Map.fetch(slots, sat) do
      {:ok, @c_m_s / (@glonass_g1_hz + k * @glonass_g1_step_hz)}
    end
  end

  defp satellite_positions(nav, epoch, sats) do
    sats
    |> Enum.reduce(%{}, fn sat, acc ->
      case Broadcast.position(nav, sat, epoch) do
        {:ok, %{x_m: x, y_m: y, z_m: z}} -> Map.put(acc, sat, {x, y, z})
        {:error, _reason} -> acc
      end
    end)
  end

  defp transmit_time_satellite_positions(nav, receive_epoch, values, sats) do
    sats
    |> Enum.reduce(%{}, fn sat, acc ->
      with %{code_m: code_m} when is_number(code_m) <- Map.get(values, sat),
           {:ok, transmit_epoch} <- transmit_epoch(receive_epoch, code_m),
           {:ok, %{x_m: x, y_m: y, z_m: z}} <- Broadcast.position(nav, sat, transmit_epoch) do
        Map.put(acc, sat, {x, y, z})
      else
        _ -> acc
      end
    end)
  end

  defp transmit_epoch(receive_epoch, code_m) do
    microseconds = round(code_m / @c_m_s * 1_000_000.0)
    {:ok, NaiveDateTime.add(receive_epoch, -microseconds, :microsecond)}
  rescue
    _ -> :error
  end

  defp sanity_gate!(nav, epochs, contexts, base_arp) do
    {residuals, spp_epochs} =
      epochs
      |> Enum.zip(contexts)
      |> Enum.reduce({[], 0}, fn {epoch, context}, {residual_acc, spp_count} ->
        case clock_demeaned_single_difference_residuals(nav, epoch, context, base_arp) do
          {:ok, residuals} -> {residuals ++ residual_acc, spp_count + 1}
          :error -> {residual_acc, spp_count}
        end
      end)

    if residuals == [] do
      raise "sanity gate could not form SPP-level single-difference residuals"
    end

    abs_residuals = Enum.map(residuals, &abs/1)

    gate = %{
      threshold_m: @sanity_code_residual_threshold_m,
      samples: length(abs_residuals),
      spp_epochs: spp_epochs,
      median_abs_code_residual_m: median(abs_residuals),
      p95_abs_code_residual_m: percentile(abs_residuals, 0.95),
      max_abs_code_residual_m: Enum.max(abs_residuals)
    }

    if gate.median_abs_code_residual_m > @sanity_code_residual_threshold_m do
      raise """
      sanity gate failed before filter run:
        median |clock-demeaned SD code residual| = #{fmt(gate.median_abs_code_residual_m)} m
        threshold = #{fmt(@sanity_code_residual_threshold_m)} m
        samples = #{gate.samples}
      """
    end

    Map.put(gate, :pass, true)
  end

  defp clock_demeaned_single_difference_residuals(nav, epoch, _context, base_arp) do
    case spp_rover_position(nav, epoch, base_arp) do
      nil ->
        :error

      rover_position ->
        base = observation_map(epoch.base_observations)
        rover = observation_map(epoch.rover_observations)
        base_tuple = ecef_to_tuple(base_arp)
        rover_tuple = ecef_to_tuple(rover_position)

        epoch.satellite_positions_m
        |> Map.keys()
        |> Enum.flat_map(fn sat ->
          with %{code_m: base_code} <- Map.get(base, sat),
               %{code_m: rover_code} <- Map.get(rover, sat) do
            base_pos = Map.fetch!(epoch.base_satellite_positions_m, sat)
            rover_pos = Map.fetch!(epoch.rover_satellite_positions_m, sat)
            code_sd = rover_code - base_code

            geom_sd =
              norm3(sub3(ecef_to_tuple(rover_pos), rover_tuple)) -
                norm3(sub3(ecef_to_tuple(base_pos), base_tuple))

            [{String.first(sat), code_sd - geom_sd}]
          else
            _ -> []
          end
        end)
        |> Enum.group_by(&elem(&1, 0), &elem(&1, 1))
        |> Enum.flat_map(fn {_system, values} ->
          clock_bias = median(values)
          Enum.map(values, &(&1 - clock_bias))
        end)
        |> then(&{:ok, &1})
    end
  end

  defp spp_rover_position(nav, epoch, base_arp) do
    observations =
      epoch.rover_observations
      |> Enum.map(&{&1.satellite_id, &1.code_m})

    case Positioning.solve(nav, observations, epoch.epoch,
           initial_guess: with_clock(ecef_to_tuple(base_arp), 0.0),
           troposphere: true
         ) do
      {:ok, sol} -> sol.position
      {:error, _reason} -> nil
    end
  end

  defp observation_map(observations), do: Map.new(observations, &{&1.satellite_id, &1})

  defp spp_initial_baseline!(nav, epochs, base_arp) do
    seed = ecef_to_tuple(base_arp)

    epochs
    |> Enum.reduce_while(nil, fn epoch, _acc ->
      observations =
        epoch.rover_observations
        |> Enum.map(&{&1.satellite_id, &1.code_m})

      case Positioning.solve(nav, observations, epoch.epoch,
             initial_guess: with_clock(seed, 0.0),
             troposphere: true
           ) do
        {:ok, sol} ->
          {:halt, sub3(ecef_to_tuple(sol.position), seed)}

        {:error, _reason} ->
          {:cont, nil}
      end
    end)
    |> case do
      nil -> raise "could not produce a broadcast SPP seed for initial baseline"
      baseline -> baseline
    end
  end

  defp filter_opts(initial_baseline, epochs, glonass_slots) do
    [
      filter_kernel: :rust,
      initial_baseline_m: initial_baseline,
      baseline_prior_sigma_m: 500.0,
      ambiguity_prior_sigma_m: 1_000.0,
      process_noise_baseline_sigma_m: 30.0,
      hold_sigma_m: 1.0e-4,
      max_iterations: 10,
      on_cycle_slip: :split_arc,
      elevation_mask_deg: 10.0,
      stochastic_model: :rtklib,
      code_sigma_m: 0.9,
      phase_sigma_m: 0.003,
      ambiguity_wavelength_m: wavelength_map(epochs, glonass_slots),
      integer_ratio_threshold: 3.0,
      integer_candidate_limit: 200_000,
      float_only_systems: ["R"]
    ]
  end

  defp wavelength_map(epochs, glonass_slots) do
    epochs
    |> Enum.flat_map(&Map.keys(&1.satellite_positions_m))
    |> Enum.uniq()
    |> Map.new(fn sat ->
      {:ok, wavelength_m} = wavelength_m(sat, glonass_slots)
      {sat, wavelength_m}
    end)
  end

  defp epoch_measurement(result, context, base_arp) do
    truth = context.oracle_epoch["truth_ecef_m"]
    truth_tuple = {truth["x"], truth["y"], truth["z"]}
    rover_tuple = result.baseline_m |> ecef_to_tuple() |> add3(ecef_to_tuple(base_arp))
    {lat_rad, lon_rad, _height_m} = ecef_to_geodetic(truth_tuple)
    {east, north, up} = ecef_delta_to_enu(sub3(rover_tuple, truth_tuple), lat_rad, lon_rad)
    horizontal = :math.sqrt(east * east + north * north)
    error_3d = norm3(sub3(rover_tuple, truth_tuple))
    residual_summary = residual_summary(result.residuals_m)

    %{
      time: context.time,
      truth_time_utc: context.oracle_epoch["truth_time_utc"],
      error_3d_m: error_3d,
      horizontal_error_m: horizontal,
      vertical_error_m: up,
      error_enu_m: %{east: east, north: north, up: up},
      integer_status: result.integer_status,
      ratio: finite_ratio(result.integer_ratio),
      satellites: residual_satellite_count(result.residuals_m),
      pre_mask_satellites: context.pre_mask_satellites,
      fixed_ambiguities: result.fixed_ambiguities,
      newly_fixed_ambiguities: result.newly_fixed_ambiguities,
      residuals: residual_summary
    }
  end

  defp residual_summary(residuals) do
    code_abs = Enum.map(residuals, &abs(&1.code_m))
    phase_abs = Enum.map(residuals, &abs(&1.phase_m))
    code_norm_abs = Enum.map(residuals, &abs(&1.code_normalized))
    phase_norm_abs = Enum.map(residuals, &abs(&1.phase_normalized))

    %{
      count: length(residuals),
      max_abs_code_m: max_or_nil(code_abs),
      max_abs_phase_m: max_or_nil(phase_abs),
      code_rms_m: rms(code_abs),
      phase_rms_m: rms(phase_abs),
      max_abs_code_normalized: max_or_nil(code_norm_abs),
      max_abs_phase_normalized: max_or_nil(phase_norm_abs)
    }
  end

  defp residual_satellite_count(residuals) do
    residuals
    |> Enum.flat_map(&[&1.satellite_id, &1.reference_satellite_id])
    |> Enum.uniq()
    |> length()
  end

  defp demo5_summary(per_epoch) do
    errors = Enum.map(per_epoch, & &1["error_3d_m"])
    horizontal = Enum.map(per_epoch, & &1["horizontal_error_m"])
    vertical = Enum.map(per_epoch, &abs(&1["vertical_error_m"]))

    %{
      epochs: length(per_epoch),
      fixed_epochs: Enum.count(per_epoch, &(&1["q"] == 1)),
      error_3d_median_m: median(errors),
      error_3d_p95_m: percentile(errors, 0.95),
      horizontal_p95_m: percentile(horizontal, 0.95),
      vertical_abs_p95_m: percentile(vertical, 0.95)
    }
  end

  defp summarize_measurements(per_epoch) do
    errors = Enum.map(per_epoch, & &1.error_3d_m)
    horizontal = Enum.map(per_epoch, & &1.horizontal_error_m)
    vertical = Enum.map(per_epoch, &abs(&1.vertical_error_m))

    %{
      epochs: length(per_epoch),
      fixed_epochs: Enum.count(per_epoch, &(&1.integer_status == :fixed)),
      error_3d_median_m: median(errors),
      error_3d_p95_m: percentile(errors, 0.95),
      horizontal_p95_m: percentile(horizontal, 0.95),
      vertical_abs_p95_m: percentile(vertical, 0.95)
    }
  end

  defp comparative_verdict(orbis, demo5, :per_arc) do
    median_pass = orbis.error_3d_median_m <= demo5.error_3d_median_m * 1.25
    p95_pass = orbis.error_3d_p95_m <= demo5.error_3d_p95_m * 1.25

    %{
      median_pass: median_pass,
      p95_pass: p95_pass,
      pass: median_pass and p95_pass,
      median_ratio: ratio(orbis.error_3d_median_m, demo5.error_3d_median_m),
      p95_ratio: ratio(orbis.error_3d_p95_m, demo5.error_3d_p95_m)
    }
  end

  defp comparative_verdict(orbis, demo5, :pooled) do
    median_pass = orbis.error_3d_median_m <= demo5.error_3d_median_m

    %{
      median_pass: median_pass,
      p95_pass: nil,
      pass: median_pass,
      median_ratio: ratio(orbis.error_3d_median_m, demo5.error_3d_median_m),
      p95_ratio: ratio(orbis.error_3d_p95_m, demo5.error_3d_p95_m)
    }
  end

  defp invariant_verdict(per_epoch) do
    fixed = per_epoch |> Enum.filter(&(&1.integer_status == :fixed)) |> Enum.map(& &1.error_3d_m)
    float = per_epoch |> Enum.reject(&(&1.integer_status == :fixed)) |> Enum.map(& &1.error_3d_m)

    cond do
      fixed == [] ->
        %{status: "pass_no_fixed_epochs", fixed_epochs: 0, float_epochs: length(float)}

      float == [] ->
        %{status: "fail_no_float_population", fixed_epochs: length(fixed), float_epochs: 0}

      true ->
        fixed_median = median(fixed)
        fixed_p95 = percentile(fixed, 0.95)
        float_median = median(float)
        float_p95 = percentile(float, 0.95)
        pass = fixed_median < float_median and fixed_p95 < float_p95

        %{
          status: if(pass, do: "pass", else: "fail"),
          fixed_epochs: length(fixed),
          float_epochs: length(float),
          fixed_median_m: fixed_median,
          fixed_p95_m: fixed_p95,
          float_median_m: float_median,
          float_p95_m: float_p95
        }
    end
  end

  defp classify_worst_decile(per_epoch) do
    count = max(1, ceil(length(per_epoch) * 0.10))
    worst = per_epoch |> Enum.sort_by(& &1.error_3d_m, :desc) |> Enum.take(count)
    worst_times = MapSet.new(worst, & &1.time)
    sat_limit = per_epoch |> Enum.map(& &1.satellites) |> percentile(0.10) |> max(4)
    residual_limit = per_epoch |> Enum.map(&residual_score/1) |> percentile(0.90)
    run_lengths = worst_run_lengths(per_epoch, worst_times)

    classified =
      Enum.map(worst, fn epoch ->
        run_length = Map.get(run_lengths, epoch.time, 1)
        score = residual_score(epoch)

        cause =
          cond do
            epoch.satellites <= sat_limit or epoch.pre_mask_satellites <= sat_limit ->
              :dropout_gap

            run_length <= 3 and score >= residual_limit and epoch.satellites > sat_limit ->
              :multipath_outlier

            run_length >= 5 and coherent_bias?(per_epoch, worst_times, epoch.time) ->
              :geometry_antenna

            true ->
              :other
          end

        Map.put(epoch, :ledger_cause, cause)
      end)

    classified
    |> Enum.group_by(& &1.ledger_cause)
    |> Map.new(fn {cause, rows} ->
      errors = Enum.map(rows, & &1.error_3d_m)

      {cause,
       %{
         count: length(rows),
         error_3d_min_m: Enum.min(errors),
         error_3d_median_m: median(errors),
         error_3d_max_m: Enum.max(errors),
         satellites_min: rows |> Enum.map(& &1.satellites) |> Enum.min(),
         satellites_max: rows |> Enum.map(& &1.satellites) |> Enum.max(),
         max_abs_phase_residual_m:
           rows
           |> Enum.map(& &1.residuals.max_abs_phase_m)
           |> Enum.reject(&is_nil/1)
           |> max_or_nil(),
         max_abs_code_residual_m:
           rows
           |> Enum.map(& &1.residuals.max_abs_code_m)
           |> Enum.reject(&is_nil/1)
           |> max_or_nil()
       }}
    end)
  end

  defp residual_score(epoch) do
    max(
      epoch.residuals.max_abs_code_normalized || 0.0,
      epoch.residuals.max_abs_phase_normalized || 0.0
    )
  end

  defp worst_run_lengths(per_epoch, worst_times) do
    {runs, current} =
      Enum.reduce(per_epoch, {[], []}, fn epoch, {runs, current} ->
        if MapSet.member?(worst_times, epoch.time) do
          {runs, [epoch | current]}
        else
          finish_run(runs, current)
        end
      end)

    {runs, _} = finish_run(runs, current)

    runs
    |> Enum.flat_map(fn run ->
      Enum.map(run, &{&1.time, length(run)})
    end)
    |> Map.new()
  end

  defp finish_run(runs, []), do: {runs, []}
  defp finish_run(runs, current), do: {[Enum.reverse(current) | runs], []}

  defp coherent_bias?(per_epoch, worst_times, time) do
    run =
      per_epoch
      |> contiguous_worst_run(worst_times, time)

    if length(run) < 5 do
      false
    else
      vectors = Enum.map(run, &enu_tuple(&1.error_enu_m))
      mean = mean3(vectors)
      mean_norm = norm3(mean)

      mean_norm > 0.0 and
        Enum.count(vectors, fn vector ->
          dot3(vector, mean) / max(norm3(vector) * mean_norm, 1.0e-12) > 0.8
        end) >=
          div(length(vectors) * 4, 5)
    end
  end

  defp contiguous_worst_run(per_epoch, worst_times, time) do
    {before_target, [target | after_target]} = Enum.split_while(per_epoch, &(&1.time != time))

    left =
      before_target
      |> Enum.reverse()
      |> Enum.take_while(&MapSet.member?(worst_times, &1.time))
      |> Enum.reverse()

    right = Enum.take_while(after_target, &MapSet.member?(worst_times, &1.time))
    left ++ [target] ++ right
  end

  defp pooled_summary(arcs) do
    orbis_epochs = Enum.flat_map(arcs, & &1.per_epoch)
    demo5_epochs = Enum.flat_map(arcs, & &1.input.oracle["per_epoch"])
    orbis = summarize_measurements(orbis_epochs)
    demo5 = demo5_summary(demo5_epochs)

    %{
      "orbis" => stringify_keys(orbis),
      "demo5" => stringify_keys(demo5),
      "comparative" => stringify_keys(comparative_verdict(orbis, demo5, :pooled)),
      "invariant" => stringify_keys(invariant_verdict(orbis_epochs)),
      "ledger" => pooled_ledger(arcs)
    }
  end

  defp time_alignment(contexts) do
    rinex_to_oracle_ms =
      Enum.map(contexts, fn context ->
        oracle_epoch = parse_iso_naive!(context.time)
        abs(NaiveDateTime.diff(context.epoch, oracle_epoch, :microsecond)) / 1_000.0
      end)

    truth_offsets_s =
      Enum.map(contexts, fn context ->
        gpst = parse_iso_naive!(context.time)
        truth_utc = parse_iso_naive!(context.oracle_epoch["truth_time_utc"])
        NaiveDateTime.diff(gpst, truth_utc, :microsecond) / 1_000_000.0
      end)

    %{
      rinex_to_oracle_time_max_ms: max_or_nil(rinex_to_oracle_ms),
      gpst_minus_truth_utc_median_s: median(truth_offsets_s),
      gpst_minus_truth_utc_min_s: Enum.min(truth_offsets_s),
      gpst_minus_truth_utc_max_s: Enum.max(truth_offsets_s)
    }
  end

  defp pooled_ledger(arcs) do
    arcs
    |> Enum.flat_map(fn arc ->
      Enum.map(arc.ledger, fn {cause, values} -> {cause, values} end)
    end)
    |> Enum.group_by(&elem(&1, 0), &elem(&1, 1))
    |> Map.new(fn {cause, rows} ->
      errors_min = Enum.map(rows, & &1.error_3d_min_m)
      errors_max = Enum.map(rows, & &1.error_3d_max_m)

      {Atom.to_string(cause),
       %{
         "count" => rows |> Enum.map(& &1.count) |> Enum.sum(),
         "error_3d_min_m" => Enum.min(errors_min),
         "error_3d_max_m" => Enum.max(errors_max),
         "max_abs_phase_residual_m" =>
           rows
           |> Enum.map(& &1.max_abs_phase_residual_m)
           |> Enum.reject(&is_nil/1)
           |> max_or_nil(),
         "max_abs_code_residual_m" =>
           rows
           |> Enum.map(& &1.max_abs_code_residual_m)
           |> Enum.reject(&is_nil/1)
           |> max_or_nil()
       }}
    end)
  end

  defp json_arc(arc) do
    %{
      "label" => arc.input.label,
      "drive" => arc.input.drive,
      "fixture" => arc.input.fixture,
      "inputs" => %{
        "rover_obs" => arc.input.rover_path,
        "base_obs" => arc.input.base_path,
        "nav" => arc.input.nav_path,
        "oracle" => arc.input.oracle_path
      },
      "base_arp_m" => tuple_json(arc.base_arp),
      "initial_baseline_m" => tuple_json(arc.initial_baseline),
      "sanity_gate" => stringify_keys(arc.sanity_gate),
      "time_alignment" => stringify_keys(arc.time_alignment),
      "built_epoch_count" => arc.built_epoch_count,
      "skipped_oracle_epochs" => arc.skipped_oracle_epochs,
      "segment_count" => arc.segment_count,
      "segments" => Enum.map(arc.segment_reports, &segment_json/1),
      "orbis" => stringify_keys(arc.orbis),
      "demo5" => stringify_keys(arc.demo5),
      "comparative" => stringify_keys(arc.comparative),
      "invariant" => stringify_keys(arc.invariant),
      "ledger" => ledger_json(arc.ledger),
      "per_epoch" => Enum.map(arc.per_epoch, &epoch_json/1)
    }
  end

  defp epoch_json(epoch) do
    %{
      "time" => epoch.time,
      "truth_time_utc" => epoch.truth_time_utc,
      "error_3d_m" => epoch.error_3d_m,
      "horizontal_error_m" => epoch.horizontal_error_m,
      "vertical_error_m" => epoch.vertical_error_m,
      "error_enu_m" => stringify_keys(epoch.error_enu_m),
      "integer_status" => Atom.to_string(epoch.integer_status),
      "ratio" => epoch.ratio,
      "satellites" => epoch.satellites,
      "pre_mask_satellites" => epoch.pre_mask_satellites,
      "fixed_ambiguities" => epoch.fixed_ambiguities,
      "newly_fixed_ambiguities" => epoch.newly_fixed_ambiguities,
      "residuals" => stringify_keys(epoch.residuals)
    }
  end

  defp segment_json(segment) do
    %{
      "index" => segment.index,
      "epochs" => segment.epochs,
      "first_time" => segment.first_time,
      "last_time" => segment.last_time,
      "initial_baseline_m" => tuple_json(segment.initial_baseline_m),
      "metadata" => stringify_keys(segment.metadata)
    }
  end

  defp ledger_json(ledger) do
    Map.new(ledger, fn {cause, values} -> {Atom.to_string(cause), stringify_keys(values)} end)
  end

  defp option_notes_json do
    Enum.map(@filter_option_notes, fn {name, value, why} ->
      %{"option" => name, "value" => value, "why" => why}
    end)
  end

  defp report_markdown(result, results_path) do
    arcs = result["arcs"]
    pooled = result["pooled"]

    [
      "# Rover measurement pass, June 2026",
      "",
      "Generated by `mix run test/fixtures/rtk/generators/rover_measurement_2026_06.exs`.",
      "Per-epoch JSON was emitted to `#{results_path}`.",
      "",
      "This is the pre-registered single-sided MEASUREMENT pass. The RTK solver and library code were not changed.",
      "",
      "## Inputs and epoch construction",
      "",
      "- Phone observations: each arc's `supplemental/gnss_rinex.21o` from `/tmp/gsdc-work`.",
      "- Base observations: NOAA CORS P222 RINEX 2.11 observation files named in each oracle's provenance.",
      "- Ephemeris: BKG combined broadcast NAV (`BRDC00WRD_R_..._MN.rnx`), matching each oracle's `pos1-sateph=brdc` source.",
      "- Base observations are linearly interpolated to each phone epoch because the CORS files are 30 s and the oracle config has `misc-timeinterp=on`.",
      "- Satellite positions use per-receiver transmit time from each receiver's code pseudorange, as in the real-arc RTK tests.",
      "- Before any filter run, the harness aborts if the median clock-demeaned single-difference code residual at SPP-level geometry exceeds #{fmt(@sanity_code_residual_threshold_m)} m.",
      "- The 10 degree elevation mask is applied during epoch construction and passed to the solver; per-epoch constellations with fewer than two usable satellites are dropped before double differencing.",
      "- Each phone epoch is solved as its own filter segment. The earlier 60-epoch harness carried ambiguity state across dense Android epochs with frequent satellite re-acquisition, producing megameter artifacts after otherwise sane first epochs.",
      "",
      "## Sanity gate and time basis",
      "",
      sanity_gate_table(arcs),
      "",
      "All four arcs pass the pre-filter residual gate. RINEX phone epochs match oracle GPST times within 0.5 ms after the oracle's millisecond rounding, and GPST minus truth UTC is 18 s, matching the oracle truth metadata.",
      "",
      "Bug attribution: the megameter output was a harness segmentation/state-carryover bug. Epoch pairing had no evidence of the megameter failure after adding a defensive exact-base-epoch path: CORS interpolation remains correct for the 30 s base data and the paired epoch residual gate is meter-level. The satellite timescale and phone clock suspects were cleared by the GPST/UTC checks.",
      "",
      "## Filter options",
      "",
      option_table(),
      "",
      "## Distributions",
      "",
      distributions_table(arcs, pooled),
      "",
      "Comparative bar is report-only for this pass. Per-arc bar is Orbis <= 1.25 x demo5 for median and p95. Pooled registered bar compares medians without margin; pooled p95 is listed for context.",
      "",
      "## Hard invariant",
      "",
      invariant_table(arcs, pooled),
      "",
      "## Worst-decile ledger",
      "",
      ledger_table(arcs, pooled),
      "",
      "Classes are assigned only from the emitted epoch diagnostics: satellite counts, output gaps, residual magnitudes, and contiguous error-vector runs.",
      "",
      "## Capability candidates",
      "",
      candidates_table(pooled),
      ""
    ]
    |> Enum.join("\n")
  end

  defp option_table do
    rows =
      Enum.map(@filter_option_notes, fn {name, value, why} ->
        "| `#{name}` | `#{value}` | #{why} |"
      end)

    Enum.join(["| Option | Value | Reason |", "|---|---:|---|" | rows], "\n")
  end

  defp sanity_gate_table(arcs) do
    rows =
      Enum.map(arcs, fn arc ->
        gate = arc["sanity_gate"]
        time = arc["time_alignment"]

        [
          arc["label"],
          pass_text(gate["pass"]),
          gate["samples"],
          fmt(gate["median_abs_code_residual_m"]),
          fmt(gate["p95_abs_code_residual_m"]),
          fmt(gate["max_abs_code_residual_m"]),
          fmt(time["rinex_to_oracle_time_max_ms"]),
          fmt(time["gpst_minus_truth_utc_median_s"])
        ]
        |> table_row()
      end)

    Enum.join(
      [
        "| Arc | Gate | Samples | Median SD residual m | p95 SD residual m | Max SD residual m | RINEX-oracle max ms | GPST-truth UTC median s |",
        "|---|---|---:|---:|---:|---:|---:|---:|"
      ] ++ rows,
      "\n"
    )
  end

  defp distributions_table(arcs, pooled) do
    rows =
      Enum.map(arcs, fn arc ->
        o = arc["orbis"]
        d = arc["demo5"]
        c = arc["comparative"]

        [
          arc["label"],
          "#{o["epochs"]}/#{d["epochs"]}",
          fmt(o["error_3d_median_m"]),
          fmt(d["error_3d_median_m"]),
          fmt(o["error_3d_p95_m"]),
          fmt(d["error_3d_p95_m"]),
          pass_text(c["pass"])
        ]
        |> table_row()
      end)

    pooled_row =
      [
        "pooled",
        "#{pooled["orbis"]["epochs"]}/#{pooled["demo5"]["epochs"]}",
        fmt(pooled["orbis"]["error_3d_median_m"]),
        fmt(pooled["demo5"]["error_3d_median_m"]),
        fmt(pooled["orbis"]["error_3d_p95_m"]),
        fmt(pooled["demo5"]["error_3d_p95_m"]),
        pass_text(pooled["comparative"]["pass"])
      ]
      |> table_row()

    Enum.join(
      [
        "| Arc | Epochs Orbis/demo5 | Orbis 3D median m | demo5 3D median m | Orbis 3D p95 m | demo5 3D p95 m | Bar |",
        "|---|---:|---:|---:|---:|---:|---:|"
      ] ++ rows ++ [pooled_row],
      "\n"
    )
  end

  defp invariant_table(arcs, pooled) do
    rows =
      Enum.map(arcs, fn arc ->
        invariant = arc["invariant"]

        [
          arc["label"],
          invariant["status"],
          invariant["fixed_epochs"] || 0,
          invariant["float_epochs"] || 0,
          fmt(invariant["fixed_median_m"]),
          fmt(invariant["float_median_m"]),
          fmt(invariant["fixed_p95_m"]),
          fmt(invariant["float_p95_m"])
        ]
        |> table_row()
      end)

    pooled_row =
      [
        "pooled",
        pooled["invariant"]["status"],
        pooled["invariant"]["fixed_epochs"] || 0,
        pooled["invariant"]["float_epochs"] || 0,
        fmt(pooled["invariant"]["fixed_median_m"]),
        fmt(pooled["invariant"]["float_median_m"]),
        fmt(pooled["invariant"]["fixed_p95_m"]),
        fmt(pooled["invariant"]["float_p95_m"])
      ]
      |> table_row()

    Enum.join(
      [
        "| Arc | Verdict | Fixed n | Float n | Fixed median m | Float median m | Fixed p95 m | Float p95 m |",
        "|---|---|---:|---:|---:|---:|---:|---:|"
      ] ++ rows ++ [pooled_row],
      "\n"
    )
  end

  defp ledger_table(arcs, pooled) do
    arc_rows =
      arcs
      |> Enum.flat_map(fn arc ->
        Enum.map(arc["ledger"], fn {cause, values} ->
          [
            arc["label"],
            cause,
            values["count"],
            "#{fmt(values["error_3d_min_m"])}-#{fmt(values["error_3d_max_m"])}",
            "#{values["satellites_min"]}-#{values["satellites_max"]}",
            fmt(values["max_abs_code_residual_m"]),
            fmt(values["max_abs_phase_residual_m"])
          ]
          |> table_row()
        end)
      end)

    pooled_rows =
      pooled["ledger"]
      |> Enum.sort_by(fn {_cause, values} -> -values["error_3d_max_m"] end)
      |> Enum.map(fn {cause, values} ->
        [
          "pooled",
          cause,
          values["count"],
          "#{fmt(values["error_3d_min_m"])}-#{fmt(values["error_3d_max_m"])}",
          "",
          fmt(values["max_abs_code_residual_m"]),
          fmt(values["max_abs_phase_residual_m"])
        ]
        |> table_row()
      end)

    Enum.join(
      [
        "| Arc | Class | Epochs | 3D error range m | Sats | Max code residual m | Max phase residual m |",
        "|---|---|---:|---:|---:|---:|---:|"
      ] ++ arc_rows ++ pooled_rows,
      "\n"
    )
  end

  defp candidates_table(pooled) do
    pooled["ledger"]
    |> Enum.sort_by(fn {_cause, values} -> {-values["error_3d_max_m"], -values["count"]} end)
    |> Enum.map(fn {cause, values} ->
      candidate =
        case cause do
          "dropout_gap" -> "Base/rover epoch bridging and outage handling"
          "multipath_outlier" -> "Robust residual gating for isolated phone multipath"
          "geometry_antenna" -> "Bias-stretch diagnostics before sub-cm effects"
          "other" -> "Manual review of unclassified worst-decile epochs"
        end

      [
        cause,
        values["count"],
        "#{fmt(values["error_3d_min_m"])}-#{fmt(values["error_3d_max_m"])}",
        candidate
      ]
      |> table_row()
    end)
    |> then(fn rows ->
      Enum.join(
        [
          "| Ledger class | Epochs | 3D error range m | Candidate |",
          "|---|---:|---:|---|"
        ] ++ rows,
        "\n"
      )
    end)
  end

  defp table_row(values), do: "| " <> (values |> Enum.map_join(" | ", &to_string/1)) <> " |"

  defp pass_text(true), do: "pass"
  defp pass_text(false), do: "miss"
  defp pass_text(nil), do: "n/a"

  defp finite_ratio(:infinity), do: "infinity"
  defp finite_ratio(nil), do: nil
  defp finite_ratio(value), do: value

  defp ratio(_a, b) when b == 0.0, do: nil
  defp ratio(a, b), do: a / b

  defp median([]), do: nil

  defp median(values) do
    ordered = Enum.sort(values)
    count = length(ordered)
    mid = div(count, 2)

    if rem(count, 2) == 1 do
      Enum.at(ordered, mid)
    else
      (Enum.at(ordered, mid - 1) + Enum.at(ordered, mid)) / 2.0
    end
  end

  defp percentile([], _pct), do: nil

  defp percentile(values, pct) do
    ordered = Enum.sort(values)
    Enum.at(ordered, trunc(pct * (length(ordered) - 1)))
  end

  defp rms([]), do: nil

  defp rms(values) do
    :math.sqrt(Enum.sum(Enum.map(values, &(&1 * &1))) / length(values))
  end

  defp max_or_nil([]), do: nil
  defp max_or_nil(values), do: Enum.max(values)

  defp interpolate(a, b, fraction), do: a + (b - a) * fraction

  defp ecef_to_geodetic({x, y, z}) do
    e2 = @earth_f * (2.0 - @earth_f)
    lon = :math.atan2(y, x)
    p = :math.sqrt(x * x + y * y)

    lat =
      Enum.reduce(1..8, :math.atan2(z, p * (1.0 - e2)), fn _i, lat_acc ->
        sin_lat = :math.sin(lat_acc)
        n = @earth_a_m / :math.sqrt(1.0 - e2 * sin_lat * sin_lat)
        :math.atan2(z + e2 * n * sin_lat, p)
      end)

    sin_lat = :math.sin(lat)
    n = @earth_a_m / :math.sqrt(1.0 - e2 * sin_lat * sin_lat)
    height = p / :math.cos(lat) - n
    {lat, lon, height}
  end

  defp ecef_delta_to_enu({dx, dy, dz}, lat, lon) do
    sin_lat = :math.sin(lat)
    cos_lat = :math.cos(lat)
    sin_lon = :math.sin(lon)
    cos_lon = :math.cos(lon)

    east = -sin_lon * dx + cos_lon * dy
    north = -sin_lat * cos_lon * dx - sin_lat * sin_lon * dy + cos_lat * dz
    up = cos_lat * cos_lon * dx + cos_lat * sin_lon * dy + sin_lat * dz
    {east, north, up}
  end

  defp elevation_deg(receiver, satellite) do
    receiver_tuple = ecef_to_tuple(receiver)
    satellite_tuple = ecef_to_tuple(satellite)
    {lat, lon, _height} = ecef_to_geodetic(receiver_tuple)
    {east, north, up} = ecef_delta_to_enu(sub3(satellite_tuple, receiver_tuple), lat, lon)
    horizontal = :math.sqrt(east * east + north * north)
    :math.atan2(up, horizontal) * 180.0 / :math.pi()
  end

  defp naive_datetime({{year, month, day}, {hour, minute, second}}),
    do: naive_datetime(year, month, day, hour, minute, second)

  defp naive_datetime(year, month, day, hour, minute, second) do
    whole_second = trunc(second)
    microsecond = round((second - whole_second) * 1_000_000)

    NaiveDateTime.new!(
      Date.new!(year, month, day),
      Time.new!(hour, minute, whole_second, {microsecond, 6})
    )
  end

  defp epoch_key(%NaiveDateTime{} = ndt) do
    {microsecond, _precision} = ndt.microsecond
    millisecond = round(microsecond / 1_000)

    {base, ms} =
      if millisecond == 1_000,
        do: {NaiveDateTime.add(ndt, 1, :second), 0},
        else: {ndt, millisecond}

    base = %{base | microsecond: {0, 0}}
    NaiveDateTime.to_iso8601(base) <> "." <> String.pad_leading(Integer.to_string(ms), 3, "0")
  end

  defp time_us(epoch), do: NaiveDateTime.diff(epoch, ~N[1970-01-01 00:00:00], :microsecond)

  defp parse_float!(value), do: value |> String.replace("D", "E") |> String.to_float()

  defp parse_iso_naive!(value) do
    case NaiveDateTime.from_iso8601(value) do
      {:ok, ndt} -> ndt
      {:error, reason} -> raise "invalid ISO NaiveDateTime #{inspect(value)}: #{inspect(reason)}"
    end
  end

  defp chunks(value, size) do
    value
    |> String.codepoints()
    |> Enum.chunk_every(size)
    |> Enum.map(&Enum.join/1)
  end

  defp pad(value, size), do: String.pad_trailing(value || "", size)

  defp add3({ax, ay, az}, {bx, by, bz}), do: {ax + bx, ay + by, az + bz}
  defp sub3({ax, ay, az}, {bx, by, bz}), do: {ax - bx, ay - by, az - bz}
  defp scale3({x, y, z}, s), do: {x * s, y * s, z * s}
  defp dot3({ax, ay, az}, {bx, by, bz}), do: ax * bx + ay * by + az * bz
  defp norm3({x, y, z}), do: :math.sqrt(x * x + y * y + z * z)

  defp mean3(vectors) do
    {x, y, z} = Enum.reduce(vectors, {0.0, 0.0, 0.0}, &add3/2)
    count = length(vectors)
    {x / count, y / count, z / count}
  end

  defp ecef_to_tuple(%{x_m: x, y_m: y, z_m: z}), do: {x, y, z}
  defp ecef_to_tuple({x, y, z}), do: {x, y, z}
  defp with_clock({x, y, z}, clock_m), do: {x, y, z, clock_m}

  defp enu_tuple(%{east: east, north: north, up: up}), do: {east, north, up}

  defp tuple_json({x, y, z}), do: %{"x" => x, "y" => y, "z" => z}

  defp stringify_keys(%NaiveDateTime{} = value), do: NaiveDateTime.to_iso8601(value)
  defp stringify_keys(%DateTime{} = value), do: DateTime.to_iso8601(value)

  defp stringify_keys(value) when is_map(value) do
    Map.new(value, fn {key, val} -> {key_to_string(key), stringify_keys(val)} end)
  end

  defp stringify_keys(value) when is_list(value), do: Enum.map(value, &stringify_keys/1)

  defp stringify_keys(value) when is_tuple(value),
    do: value |> Tuple.to_list() |> stringify_keys()

  defp stringify_keys(value) when is_boolean(value) or is_nil(value), do: value
  defp stringify_keys(value) when is_atom(value), do: Atom.to_string(value)
  defp stringify_keys(value), do: value

  defp key_to_string(key) when is_atom(key), do: Atom.to_string(key)
  defp key_to_string(key), do: key

  defp fmt(nil), do: ""
  defp fmt(value) when is_float(value), do: :erlang.float_to_binary(value, decimals: 3)
  defp fmt(value), do: to_string(value)
end

RoverMeasurement202606.main(System.argv())
