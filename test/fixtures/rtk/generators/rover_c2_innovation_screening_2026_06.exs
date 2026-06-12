defmodule RoverC2InnovationScreening202606 do
  @moduledoc false
  alias Orbis.GNSS.Broadcast
  alias Orbis.GNSS.CarrierPhase
  alias Orbis.GNSS.Positioning
  alias Orbis.GNSS.RINEX.Observations
  alias Orbis.GNSS.RTK
  alias Orbis.NIF

  @compile :nowarn_unused_function

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
  @default_results "/tmp/c2-innovation-screening-2026-06-results.json"
  @max_segment_epochs 60
  @comparison_segment_epochs 1
  @sanity_code_residual_threshold_m 1_000.0
  @divergence_threshold_m 1_000.0
  @rtk_filter_state_version 3
  @stable_reference_median_m 9.533
  @demo5_reference_median_m 4.007
  @innovation_screen_sigmas [5.0, 8.0, 10.0, 15.0, 25.0]
  @innovation_screen_min_rows 8
  @process_noise_sweep_sigmas [1.0, 3.0, 5.5, 10.0, 30.0]
  @continuity_gap_s 2.5
  @code_phase_min_threshold_m 15.0
  @code_phase_sigma_multiplier 8.0
  @doppler_min_threshold_m 1.0
  @doppler_sigma_multiplier 8.0
  @gf_threshold_m 0.05
  @mw_threshold_cycles 4.0
  @detector_segment_max_sd_ambiguity_ids 80
  @detector_max_injected_events_per_satellite 5

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

  @phone_detector_codes %{
    "G" => ["C1C", "L1C", "D1C", "S1C", "C5X", "L5X", "D5X", "S5X"],
    "R" => ["C1C", "L1C", "D1C", "S1C"],
    "E" => ["C1C", "L1C", "D1C", "S1C", "C5X", "L5X", "D5X", "S5X"],
    "C" => ["C2I", "L2I", "D2I", "S2I"]
  }

  @primary_detector_codes %{
    "G" => %{code: "C1C", phase: "L1C", doppler: "D1C"},
    "R" => %{code: "C1C", phase: "L1C", doppler: "D1C"},
    "E" => %{code: "C1C", phase: "L1C", doppler: "D1C"},
    "C" => %{code: "C2I", phase: "L2I", doppler: "D2I"}
  }

  @dual_detector_codes %{
    "G" => %{code: "C5X", phase: "L5X"},
    "E" => %{code: "C5X", phase: "L5X"}
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
     "Kept at the original kinematic setting; it applies between carried-state epochs."},
    {"hold_sigma_m", "1.0e-4",
     "Keeps the shipped tight ambiguity hold used after an accepted integer fix."},
    {"max_iterations", "10", "Matches the real-arc RTK tests' nonlinear iteration cap."},
    {"on_cycle_slip", "split_arc",
     "Kept at the real-arc test setting; detector segmentation is injected with pre-segmented rover ambiguity ids rather than library changes."},
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
    if System.get_env("ORBIS_C2_KEEP_C1_HELPER_REFS") == "1" do
      retained_c1_helper_references()
    end

    {opts, _argv, invalid} =
      OptionParser.parse(args,
        strict: [
          work: :string,
          results: :string,
          report: :string,
          limit: :integer,
          epoch_limit: :integer
        ],
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
      Keyword.get(
        opts,
        :report,
        Path.join(generator_dir, "c2-innovation-screening-2026-06.md")
      )

    run_c2_exploration(opts, fixture_dir, work, results_path, report_path)
  end

  # This runner is forked from the larger C1 exploration harness. Keep captures
  # for retained C1 report helpers so script compilation stays quiet.
  defp retained_c1_helper_references do
    [
      &arc_diagnosis/4,
      &availability_systems/1,
      &best_completed_process_noise_row/1,
      &build_filter_epochs/7,
      &candidates_table/1,
      &classify_worst_decile/1,
      &clock_demeaned_single_difference_residuals/4,
      &coherent_bias?/3,
      &comparative_verdict/3,
      &comparison_json/2,
      &condition_range_text/1,
      &contiguous_worst_run/3,
      &costs_table/2,
      &detector_agreement_table/1,
      &detector_availability_table/1,
      &detector_json/1,
      &detector_stats_table/1,
      &detector_threshold_table/1,
      &diagnosis_summary/1,
      &diagnosis_table/1,
      &distributions_table/2,
      &enu_tuple/1,
      &epoch2_condition/1,
      &epoch2_condition_text/1,
      &epoch_count_text/1,
      &epoch_json/1,
      &filter_case_json/2,
      &find_pooled_process_noise_cell/3,
      &find_process_noise_cell/3,
      &finish_run/2,
      &first_bad_changes/2,
      &first_bad_mechanism/2,
      &first_bad_verdict/2,
      &fmt_ratio/1,
      &fmt_sci/1,
      &gauge_case_json/1,
      &gauge_epoch2_condition_text/1,
      &gauge_report_markdown/1,
      &gauge_result_table/2,
      &gauge_verdict/2,
      &gauge_verdict_text/1,
      &input_consistency_summary/4,
      &invariant_table/2,
      &json_arc/1,
      &key_experiment_table/2,
      &ledger_table/2,
      &mean3/1,
      &measure_arc/1,
      &measure_gauge_case/3,
      &nil_safe_sub/2,
      &option_notes_json/0,
      &option_table/0,
      &parse_iso_naive!/1,
      &pass_text/1,
      &pooled_completed_metric/1,
      &pooled_cost/1,
      &pooled_process_noise_sweep/1,
      &pooled_summary/1,
      &process_noise_epoch2/1,
      &process_noise_label/1,
      &process_noise_sweep_cell/4,
      &process_noise_sweep_json/1,
      &process_noise_sweep_table/2,
      &process_noise_sweep_verdict_text/1,
      &promotion_verdict/3,
      &promotion_verdict_text/1,
      &report_markdown/2,
      &residual_score/1,
      &run_full_exploration/5,
      &run_gauge_test/4,
      &run_process_noise_sweep/7,
      &same_sigma?/2,
      &sanity_gate!/4,
      &sanity_gate_table/1,
      &segment_column_summary/1,
      &segment_for_epoch/2,
      &segment_json/1,
      &segmentation_costs/5,
      &segmented_arc_ids_present?/1,
      &short_arc_label/1,
      &solve_segments/4,
      &spp_rover_position/3,
      &time_alignment/1,
      &tuple_json/1,
      &worst_run_lengths/2
    ]
  end

  defp run_c2_exploration(opts, fixture_dir, work, results_path, report_path) do
    arcs =
      @oracle_fixtures
      |> maybe_limit(Keyword.get(opts, :limit))
      |> Enum.map(&load_arc(fixture_dir, &1, work))
      |> Enum.map(&measure_c2_arc(&1, Keyword.get(opts, :epoch_limit)))

    pooled = c2_pooled_summary(arcs)
    best_k = pooled.verdict.best_k

    result = %{
      "version" => 1,
      "generated_at_utc" =>
        DateTime.utc_now() |> DateTime.truncate(:second) |> DateTime.to_iso8601(),
      "script" => Path.relative_to(__ENV__.file, File.cwd!()),
      "work_dir" => work,
      "brief" => "/tmp/c2-kernel-brief.md",
      "epoch_limit" => Keyword.get(opts, :epoch_limit),
      "filter_kernel" => "rust",
      "innovation_screen_sigmas" => @innovation_screen_sigmas,
      "innovation_screen_min_rows" => @innovation_screen_min_rows,
      "detector_config" => detector_config_json(),
      "arcs" => Enum.map(arcs, &c2_arc_json/1),
      "pooled" => stringify_keys(pooled),
      "best_k_error_series" => best_k_error_series(arcs, best_k)
    }

    File.write!(results_path, Jason.encode!(result, pretty: true))
    File.write!(report_path, c2_report_markdown(result, results_path))

    IO.puts("wrote #{results_path}")
    IO.puts("wrote #{report_path}")
  end

  defp measure_c2_arc(arc, epoch_limit) do
    IO.puts("C2 innovation screening: #{arc.label}")
    require_files!([arc.rover_path, arc.base_path, arc.nav_path, arc.oracle_path])

    rover_obs = Observations.load!(arc.rover_path)
    base_obs = load_rinex2_obs!(arc.base_path)
    nav = Broadcast.load!(arc.nav_path)

    base_arp = base_arp(base_obs)
    glonass_slots = Observations.glonass_slots(rover_obs)
    oracle_by_time = Map.new(arc.oracle["per_epoch"], &{&1["time"], &1})
    detector = detect_rover_slips(rover_obs, glonass_slots)

    {segmented_epochs, segmented_contexts} =
      build_filter_epochs(
        nav,
        rover_obs,
        base_obs,
        glonass_slots,
        oracle_by_time,
        base_arp,
        detector.ambiguity_ids
      )
      |> maybe_limit_epoch_contexts(epoch_limit)

    if segmented_epochs == [] do
      raise "no usable Orbis epochs for #{arc.label}"
    end

    comparison_segments =
      segment_epoch_contexts(segmented_epochs, segmented_contexts, @comparison_segment_epochs)

    {comparison_per_epoch, comparison_segment_reports} =
      solve_segments(comparison_segments, nav, base_arp, glonass_slots,
        diagnostics?: false,
        log?: false
      )

    comparison_orbis = summarize_measurements(comparison_per_epoch)
    demo5 = demo5_summary(arc.oracle["per_epoch"])

    sweep =
      Enum.map(@innovation_screen_sigmas, fn sigma ->
        IO.puts("  k=#{fmt(sigma)}")

        filter_case =
          run_filter_case(
            "detector_segmented_innovation_k_#{innovation_screen_label(sigma)}",
            segmented_epochs,
            segmented_contexts,
            nav,
            base_arp,
            glonass_slots,
            @max_segment_epochs,
            diagnostics?: false,
            log?: false,
            max_sd_ambiguity_ids: @detector_segment_max_sd_ambiguity_ids,
            filter_kernel: :rust,
            innovation_screen_sigma: sigma,
            innovation_screen_min_rows: @innovation_screen_min_rows
          )

        innovation_screen_cell(sigma, filter_case, length(segmented_epochs))
      end)

    %{
      input: arc,
      base_arp: base_arp,
      built_epoch_count: length(segmented_epochs),
      skipped_oracle_epochs: length(arc.oracle["per_epoch"]) - length(segmented_epochs),
      detector: detector,
      comparison: %{
        mode: "per_epoch_segments_comparison_only",
        max_segment_epochs: @comparison_segment_epochs,
        segment_count: length(comparison_segment_reports),
        segment_reports: comparison_segment_reports,
        per_epoch: comparison_per_epoch,
        orbis: comparison_orbis
      },
      demo5: demo5,
      innovation_screen_sweep: sweep
    }
  end

  defp maybe_limit_epoch_contexts({epochs, contexts}, limit)
       when is_integer(limit) and limit > 0 do
    {Enum.take(epochs, limit), Enum.take(contexts, limit)}
  end

  defp maybe_limit_epoch_contexts(epoch_contexts, _limit), do: epoch_contexts

  defp innovation_screen_label(sigma_m) do
    sigma_m
    |> :erlang.float_to_binary(decimals: 1)
    |> String.replace(".", "_")
  end

  defp run_full_exploration(opts, fixture_dir, work, results_path, report_path) do
    arcs =
      @oracle_fixtures
      |> maybe_limit(Keyword.get(opts, :limit))
      |> Enum.map(&load_arc(fixture_dir, &1, work))
      |> Enum.map(&measure_arc/1)

    result = %{
      "version" => 1,
      "generated_at_utc" =>
        DateTime.utc_now() |> DateTime.truncate(:second) |> DateTime.to_iso8601(),
      "script" => Path.relative_to(__ENV__.file, File.cwd!()),
      "work_dir" => work,
      "detector_config" => detector_config_json(),
      "option_notes" => option_notes_json(),
      "arcs" => Enum.map(arcs, &json_arc/1),
      "pooled" => pooled_summary(arcs)
    }

    File.write!(results_path, Jason.encode!(result, pretty: true))
    File.write!(report_path, report_markdown(result, results_path))

    IO.puts("wrote #{results_path}")
    IO.puts("wrote #{report_path}")
  end

  defp run_gauge_test(fixture_dir, work, results_path, report_path) do
    gps_cases =
      @oracle_fixtures
      |> Enum.map(&load_arc(fixture_dir, &1, work))
      |> Enum.map(&measure_gauge_case(&1, "gps_only", ["G"]))

    control =
      @oracle_fixtures
      |> List.first()
      |> then(&load_arc(fixture_dir, &1, work))
      |> measure_gauge_case("multi_system_control", @systems)

    result = %{
      "version" => 1,
      "generated_at_utc" =>
        DateTime.utc_now() |> DateTime.truncate(:second) |> DateTime.to_iso8601(),
      "script" => Path.relative_to(__ENV__.file, File.cwd!()),
      "work_dir" => work,
      "results_path" => results_path,
      "brief" => "/tmp/gauge-test-brief.md",
      "gps_only" => Enum.map(gps_cases, &gauge_case_json/1),
      "control" => gauge_case_json(control),
      "verdict" => stringify_keys(gauge_verdict(gps_cases, control))
    }

    File.write!(results_path, Jason.encode!(result, pretty: true))
    File.write!(report_path, "\n" <> gauge_report_markdown(result), [:append])

    IO.puts("wrote #{results_path}")
    IO.puts("appended #{report_path}")
  end

  defp measure_gauge_case(arc, mode, systems) do
    IO.puts("gauge test #{mode}: #{arc.label} systems=#{Enum.join(systems, "")}")
    require_files!([arc.rover_path, arc.base_path, arc.nav_path, arc.oracle_path])

    rover_obs = Observations.load!(arc.rover_path)
    base_obs = load_rinex2_obs!(arc.base_path)
    nav = Broadcast.load!(arc.nav_path)

    base_arp = base_arp(base_obs)
    glonass_slots = Observations.glonass_slots(rover_obs)
    oracle_by_time = Map.new(arc.oracle["per_epoch"], &{&1["time"], &1})

    {epochs, contexts} =
      build_filter_epochs(
        nav,
        rover_obs,
        base_obs,
        glonass_slots,
        oracle_by_time,
        base_arp,
        %{},
        systems
      )

    if epochs == [] do
      raise "no usable Orbis epochs for #{arc.label} systems=#{inspect(systems)}"
    end

    filter_case =
      run_filter_case(
        "#{mode}_sigma_30",
        epochs,
        contexts,
        nav,
        base_arp,
        glonass_slots,
        @max_segment_epochs,
        diagnostics?: true,
        diagnostic_epoch_limit: 2,
        log?: false,
        process_noise_baseline_sigma_m: 30.0
      )

    completed = filter_case.solved_epoch_count == length(epochs)

    %{
      input: arc,
      mode: mode,
      systems: systems,
      built_epoch_count: length(epochs),
      filter_case: filter_case,
      completed_arc: completed,
      epoch2: process_noise_epoch2(filter_case.per_epoch),
      completed_arc_median_m: if(completed, do: filter_case.orbis.error_3d_median_m),
      completed_arc_p95_m: if(completed, do: filter_case.orbis.error_3d_p95_m)
    }
  end

  defp maybe_limit(fixtures, nil), do: fixtures

  defp maybe_limit(fixtures, limit) when is_integer(limit) and limit > 0,
    do: Enum.take(fixtures, limit)

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
    IO.puts("exploring #{arc.label}")
    require_files!([arc.rover_path, arc.base_path, arc.nav_path, arc.oracle_path])

    rover_obs = Observations.load!(arc.rover_path)
    base_obs = load_rinex2_obs!(arc.base_path)
    nav = Broadcast.load!(arc.nav_path)

    base_arp = base_arp(base_obs)
    glonass_slots = Observations.glonass_slots(rover_obs)
    oracle_by_time = Map.new(arc.oracle["per_epoch"], &{&1["time"], &1})
    detector = detect_rover_slips(rover_obs, glonass_slots)

    IO.puts(
      "  detector: #{detector.agreement.event_keys} event keys, #{detector.segmentation.eligible_event_keys} eligible keys, #{detector.segmentation.injected_event_keys} injected keys, #{detector.segmentation.split_ambiguity_id_count} split ids"
    )

    {epochs, contexts} =
      build_filter_epochs(nav, rover_obs, base_obs, glonass_slots, oracle_by_time, base_arp)

    if epochs == [] do
      raise "no usable Orbis epochs for #{arc.label}"
    end

    {segmented_epochs, segmented_contexts} =
      build_filter_epochs(
        nav,
        rover_obs,
        base_obs,
        glonass_slots,
        oracle_by_time,
        base_arp,
        detector.ambiguity_ids
      )

    sanity_gate = sanity_gate!(nav, epochs, contexts, base_arp)
    time_alignment = time_alignment(contexts)

    IO.puts(
      "  sanity gate: median |clock-demeaned SD code residual| = #{fmt(sanity_gate.median_abs_code_residual_m)} m"
    )

    stock =
      run_filter_case(
        "stock_carried_state",
        epochs,
        contexts,
        nav,
        base_arp,
        glonass_slots,
        @max_segment_epochs,
        diagnostics?: false
      )

    comparison_segments = segment_epoch_contexts(epochs, contexts, @comparison_segment_epochs)

    {comparison_per_epoch, comparison_segment_reports} =
      solve_segments(comparison_segments, nav, base_arp, glonass_slots,
        diagnostics?: false,
        log?: false
      )

    segmented =
      run_filter_case(
        "detector_segmented_carried_state",
        segmented_epochs,
        segmented_contexts,
        nav,
        base_arp,
        glonass_slots,
        @max_segment_epochs,
        diagnostics?: false,
        log?: false,
        max_sd_ambiguity_ids: @detector_segment_max_sd_ambiguity_ids
      )

    initial_baseline = stock.segment_reports |> List.first() |> Map.fetch!(:initial_baseline_m)

    demo5 = demo5_summary(arc.oracle["per_epoch"])
    comparison_orbis = summarize_measurements(comparison_per_epoch)

    comparison = %{
      mode: "per_epoch_segments_comparison_only",
      max_segment_epochs: @comparison_segment_epochs,
      segment_count: length(comparison_segment_reports),
      segment_reports: comparison_segment_reports,
      per_epoch: comparison_per_epoch,
      orbis: comparison_orbis,
      comparative: comparative_verdict(comparison_orbis, demo5, :per_arc),
      invariant: invariant_verdict(comparison_per_epoch)
    }

    process_noise_sweep =
      run_process_noise_sweep(
        epochs,
        contexts,
        segmented_epochs,
        segmented_contexts,
        nav,
        base_arp,
        glonass_slots
      )

    cost =
      segmentation_costs(
        detector,
        length(segmented_epochs),
        segmented,
        stock,
        comparison
      )

    %{
      input: arc,
      base_arp: base_arp,
      initial_baseline: initial_baseline,
      sanity_gate: sanity_gate,
      time_alignment: time_alignment,
      built_epoch_count: length(epochs),
      skipped_oracle_epochs: length(arc.oracle["per_epoch"]) - length(epochs),
      detector: detector,
      stock: stock,
      segmented: segmented,
      comparison: comparison,
      comparison_per_epoch: comparison_per_epoch,
      demo5: demo5,
      process_noise_sweep: process_noise_sweep,
      cost: cost
    }
  end

  defp run_process_noise_sweep(
         epochs,
         contexts,
         segmented_epochs,
         segmented_contexts,
         nav,
         base_arp,
         glonass_slots
       ) do
    [
      {"stock", epochs, contexts, []},
      {"detector", segmented_epochs, segmented_contexts,
       [max_sd_ambiguity_ids: @detector_segment_max_sd_ambiguity_ids]}
    ]
    |> Enum.flat_map(fn {mode, case_epochs, case_contexts, extra_opts} ->
      Enum.map(@process_noise_sweep_sigmas, fn sigma_m ->
        IO.puts("  process-noise sweep #{mode}: sigma=#{fmt(sigma_m)} m")

        opts =
          [
            diagnostics?: true,
            log?: false,
            process_noise_baseline_sigma_m: sigma_m
          ] ++ extra_opts

        filter_case =
          run_filter_case(
            "#{mode}_process_noise_#{process_noise_label(sigma_m)}",
            case_epochs,
            case_contexts,
            nav,
            base_arp,
            glonass_slots,
            @max_segment_epochs,
            opts
          )

        process_noise_sweep_cell(mode, sigma_m, filter_case, length(case_epochs))
      end)
    end)
  end

  defp process_noise_label(sigma_m) do
    sigma_m
    |> :erlang.float_to_binary(decimals: 1)
    |> String.replace(".", "_")
  end

  defp process_noise_sweep_cell(mode, sigma_m, filter_case, built_epoch_count) do
    completed_arc = filter_case.solved_epoch_count == built_epoch_count

    %{
      mode: mode,
      process_noise_baseline_sigma_m: sigma_m,
      built_epoch_count: built_epoch_count,
      solved_epoch_count: filter_case.solved_epoch_count,
      completed_arc: completed_arc,
      segment_count: filter_case.segment_count,
      longest_stable_prefix_epochs: filter_case.longest_stable_prefix_epochs,
      first_unstable_time: filter_case.first_unstable_time,
      epoch2: process_noise_epoch2(filter_case.per_epoch),
      orbis: filter_case.orbis,
      completed_arc_median_m: if(completed_arc, do: filter_case.orbis.error_3d_median_m),
      completed_arc_p95_m: if(completed_arc, do: filter_case.orbis.error_3d_p95_m),
      errors_3d_m: Enum.map(filter_case.per_epoch, & &1.error_3d_m)
    }
  end

  defp innovation_screen_cell(sigma, filter_case, built_epoch_count) do
    completed_arc = filter_case.solved_epoch_count == built_epoch_count

    screen_rows =
      Enum.map(filter_case.per_epoch, & &1.innovation_screen) |> Enum.reject(&is_nil/1)

    rejected_rows = Enum.map(screen_rows, & &1.rejected_rows)
    accepted_rows = Enum.map(screen_rows, & &1.accepted_rows)
    input_rows = Enum.map(screen_rows, & &1.input_rows)
    coasted_epochs = Enum.count(screen_rows, & &1.coasted?)
    total_rejected_rows = Enum.sum(rejected_rows)
    total_input_rows = Enum.sum(input_rows)

    %{
      innovation_screen_sigma: sigma,
      innovation_screen_min_rows: @innovation_screen_min_rows,
      built_epoch_count: built_epoch_count,
      solved_epoch_count: filter_case.solved_epoch_count,
      completed_arc: completed_arc,
      segment_count: filter_case.segment_count,
      longest_stable_prefix_epochs: filter_case.longest_stable_prefix_epochs,
      first_unstable_time: filter_case.first_unstable_time,
      orbis: filter_case.orbis,
      completed_arc_median_m: if(completed_arc, do: filter_case.orbis.error_3d_median_m),
      completed_arc_p95_m: if(completed_arc, do: filter_case.orbis.error_3d_p95_m),
      coasted_epochs: coasted_epochs,
      coasting_fraction: ratio(coasted_epochs, max(length(screen_rows), 1)),
      rejected_row_count: total_rejected_rows,
      input_row_count: total_input_rows,
      row_rejection_fraction: ratio(total_rejected_rows, max(total_input_rows, 1)),
      rows_rejected_per_epoch: distribution(rejected_rows),
      rows_accepted_per_epoch: distribution(accepted_rows),
      rows_rejected_counts: rejected_rows,
      rows_accepted_counts: accepted_rows,
      errors_3d_m: Enum.map(filter_case.per_epoch, & &1.error_3d_m),
      error_series: Enum.map(filter_case.per_epoch, &error_series_epoch/1)
    }
  end

  defp error_series_epoch(epoch) do
    %{
      time: epoch.time,
      error_3d_m: epoch.error_3d_m,
      horizontal_error_m: epoch.horizontal_error_m,
      vertical_error_m: epoch.vertical_error_m,
      integer_status: epoch.integer_status
    }
  end

  defp distribution(values) do
    %{
      min: min_or_nil(values),
      median: median(values),
      p95: percentile(values, 0.95),
      max: max_or_nil(values)
    }
  end

  defp process_noise_epoch2(per_epoch) do
    case Enum.at(per_epoch, 1) do
      nil ->
        nil

      epoch ->
        state = epoch.state_diagnostics || %{}

        %{
          time: epoch.time,
          stable: stable_epoch?(epoch),
          error_3d_m: epoch.error_3d_m,
          carried_baseline_error_3d_m: Map.get(state, :carried_baseline_error_3d_m),
          information_condition_estimate: Map.get(state, :information_condition_estimate),
          sd_ambiguity_columns: Map.get(state, :sd_ambiguity_columns),
          hold_count: Map.get(state, :hold_count)
        }
    end
  end

  defp run_filter_case(
         name,
         epochs,
         contexts,
         nav,
         base_arp,
         glonass_slots,
         max_segment_epochs,
         opts
       ) do
    segments =
      segment_epoch_contexts(
        epochs,
        contexts,
        max_segment_epochs,
        Keyword.get(opts, :max_sd_ambiguity_ids)
      )

    {per_epoch, segment_reports} = solve_segments(segments, nav, base_arp, glonass_slots, opts)
    summary = summarize_measurements(per_epoch)

    %{
      mode: name,
      max_segment_epochs: max_segment_epochs,
      segment_count: length(segment_reports),
      segment_reports: segment_reports,
      per_epoch: per_epoch,
      orbis: summary,
      comparative: nil,
      invariant: invariant_verdict(per_epoch),
      longest_stable_prefix_epochs: longest_stable_prefix(per_epoch),
      first_unstable_time: first_unstable_time(per_epoch),
      solved_epoch_count: length(per_epoch)
    }
  end

  defp longest_stable_prefix(per_epoch) do
    per_epoch
    |> Enum.take_while(&stable_epoch?/1)
    |> length()
  end

  defp first_unstable_time(per_epoch) do
    per_epoch
    |> Enum.find(&(not stable_epoch?(&1)))
    |> case do
      nil -> nil
      epoch -> epoch.time
    end
  end

  defp stable_epoch?(epoch), do: epoch.error_3d_m < @divergence_threshold_m

  defp segmentation_costs(detector, built_epoch_count, segmented, stock, comparison) do
    stock_columns = segment_column_summary(stock.segment_reports)
    segmented_columns = segment_column_summary(segmented.segment_reports)

    new_failure_mode =
      cond do
        segmented.solved_epoch_count < stock.solved_epoch_count ->
          "detector segmentation lost #{stock.solved_epoch_count - segmented.solved_epoch_count} epochs versus stock"

        segmented.solved_epoch_count < comparison.orbis.epochs ->
          "detector segmentation solved fewer epochs than one-epoch comparison"

        true ->
          "none observed in this run"
      end

    %{
      detector_event_keys: detector.agreement.event_keys,
      detector_eligible_event_keys: detector.segmentation.eligible_event_keys,
      detector_injected_event_keys: detector.segmentation.injected_event_keys,
      detector_events: length(detector.events),
      split_ambiguity_ids: detector.segmentation.split_ambiguity_id_count,
      stock_segments: stock.segment_count,
      detector_segments: segmented.segment_count,
      stock_solved_epochs: stock.solved_epoch_count,
      detector_solved_epochs: segmented.solved_epoch_count,
      one_epoch_solved_epochs: comparison.orbis.epochs,
      detector_epochs_lost: built_epoch_count - segmented.solved_epoch_count,
      stock_epochs_lost: built_epoch_count - stock.solved_epoch_count,
      stock_sd_columns: stock_columns,
      detector_sd_columns: segmented_columns,
      extra_max_sd_columns:
        nil_safe_sub(
          segmented_columns.max_single_difference_columns,
          stock_columns.max_single_difference_columns
        ),
      new_failure_mode: new_failure_mode
    }
  end

  defp segment_column_summary(segment_reports) do
    columns =
      segment_reports
      |> Enum.map(fn report ->
        get_in(report, [:metadata, :single_difference_ambiguity_count]) ||
          get_in(report, [:diagnostics, :max_sd_ambiguity_columns])
      end)
      |> Enum.reject(&is_nil/1)

    %{
      max_single_difference_columns: max_or_nil(columns),
      median_single_difference_columns: median(columns),
      p95_single_difference_columns: percentile(columns, 0.95)
    }
  end

  defp nil_safe_sub(a, b) when is_number(a) and is_number(b), do: a - b
  defp nil_safe_sub(_a, _b), do: nil

  defp detect_rover_slips(rover_obs, glonass_slots) do
    rows = detector_rows(rover_obs, glonass_slots)
    thresholds = detector_thresholds(rows)
    events = detector_events(rows, thresholds)
    ambiguity_ids = detector_ambiguity_ids(rows, events)
    agreement = detector_agreement(events)

    %{
      raw_epoch_count: length(Observations.epochs(rover_obs)),
      observed_satellite_count: rows |> Enum.map(& &1.satellite_id) |> Enum.uniq() |> length(),
      dual_frequency: detector_availability(rows, :dual),
      doppler: detector_availability(rows, :doppler),
      thresholds: thresholds,
      events: events,
      stats: detector_stats(rows, events),
      agreement: agreement,
      ambiguity_ids: ambiguity_ids,
      segmentation: detector_segmentation_summary(ambiguity_ids, events)
    }
  end

  defp detector_rows(rover_obs, glonass_slots) do
    codes = Map.take(@phone_detector_codes, @systems)

    rover_obs
    |> Observations.epochs()
    |> Enum.flat_map(fn entry ->
      epoch = naive_datetime(entry.epoch)
      time_key = epoch_key(epoch)
      {:ok, by_sat} = Observations.values(rover_obs, entry.index, codes: codes)

      by_sat
      |> Enum.flat_map(fn {sat, values} ->
        detector_row(sat, values, glonass_slots, epoch, time_key, entry.index)
      end)
    end)
  end

  defp detector_row(sat, values, glonass_slots, epoch, time_key, index) do
    system = String.first(sat)
    primary_codes = Map.get(@primary_detector_codes, system)
    values_by_code = Map.new(values, &{&1.code, &1})

    with %{code: code_code, phase: phase_code, doppler: doppler_code} <- primary_codes,
         %{value: code_m} when is_number(code_m) <- values_by_code[code_code],
         %{value: phase_cycles} = phase_obs when is_number(phase_cycles) <-
           values_by_code[phase_code],
         primary_frequency_hz when is_number(primary_frequency_hz) <-
           detector_frequency(sat, phase_code, glonass_slots) do
      primary_wavelength_m = @c_m_s / primary_frequency_hz

      row = %{
        epoch: epoch,
        time: time_key,
        epoch_index: index,
        satellite_id: sat,
        system: system,
        primary_code: code_code,
        primary_phase_code: phase_code,
        primary_frequency_hz: primary_frequency_hz,
        primary_wavelength_m: primary_wavelength_m,
        primary_code_m: code_m,
        primary_phase_cycles: phase_cycles,
        primary_phase_m: phase_cycles * primary_wavelength_m,
        phase_minus_code_m: phase_cycles * primary_wavelength_m - code_m,
        primary_doppler_hz: detector_value(values_by_code, doppler_code),
        primary_lli: phase_obs.lli,
        dual: detector_dual(values_by_code, system, sat, glonass_slots, code_m, phase_cycles)
      }

      [row]
    else
      _ -> []
    end
  end

  defp detector_dual(
         values_by_code,
         system,
         sat,
         glonass_slots,
         primary_code_m,
         primary_phase_cycles
       ) do
    with %{code: code_code, phase: phase_code} <- Map.get(@dual_detector_codes, system),
         %{value: code_m} when is_number(code_m) <- values_by_code[code_code],
         %{value: phase_cycles} = phase_obs when is_number(phase_cycles) <-
           values_by_code[phase_code],
         frequency_hz when is_number(frequency_hz) <-
           detector_frequency(sat, phase_code, glonass_slots),
         primary_frequency_hz when is_number(primary_frequency_hz) <-
           detector_frequency(
             sat,
             Map.fetch!(@primary_detector_codes, system).phase,
             glonass_slots
           ) do
      wavelength_m = @c_m_s / frequency_hz
      primary_wavelength_m = @c_m_s / primary_frequency_hz
      primary_phase_m = primary_phase_cycles * primary_wavelength_m
      phase_m = phase_cycles * wavelength_m

      {:ok, mw_m} =
        CarrierPhase.melbourne_wubbena(
          primary_phase_cycles,
          phase_cycles,
          primary_code_m,
          code_m,
          primary_frequency_hz,
          frequency_hz
        )

      {:ok, lambda_wl_m} = CarrierPhase.wide_lane_wavelength(primary_frequency_hz, frequency_hz)

      %{
        code: code_code,
        phase_code: phase_code,
        frequency_hz: frequency_hz,
        wavelength_m: wavelength_m,
        code_m: code_m,
        phase_cycles: phase_cycles,
        phase_m: phase_m,
        lli: phase_obs.lli,
        geometry_free_m: primary_phase_m - phase_m,
        melbourne_wubbena_m: mw_m,
        wide_lane_wavelength_m: lambda_wl_m
      }
    else
      _ -> nil
    end
  end

  defp detector_frequency(sat, code, glonass_slots) do
    Observations.band_frequency_hz(
      String.first(sat),
      String.at(code, 1),
      Map.get(glonass_slots, sat)
    )
  end

  defp detector_value(values_by_code, code) do
    case Map.get(values_by_code, code) do
      %{value: value} when is_number(value) -> value
      _ -> nil
    end
  end

  defp detector_thresholds(rows) do
    steps = detector_steps(rows)
    code_phase_values = Enum.map(steps, &abs(&1.code_phase_step_m))

    doppler_values =
      steps
      |> Enum.map(& &1.doppler_residual_m)
      |> Enum.reject(&is_nil/1)
      |> Enum.map(&abs/1)

    %{
      code_phase:
        robust_threshold(
          code_phase_values,
          @code_phase_min_threshold_m,
          @code_phase_sigma_multiplier
        ),
      doppler_phase:
        robust_threshold(doppler_values, @doppler_min_threshold_m, @doppler_sigma_multiplier),
      geometry_free_m: @gf_threshold_m,
      melbourne_wubbena_cycles: @mw_threshold_cycles,
      calibration: %{
        code_phase: robust_calibration(code_phase_values),
        doppler_phase: robust_calibration(doppler_values),
        continuity_gap_s: @continuity_gap_s
      }
    }
  end

  defp robust_threshold(values, minimum, multiplier) do
    case robust_calibration(values) do
      %{median: nil} ->
        minimum

      %{median: median_value, mad_sigma: mad_sigma} ->
        max(minimum, median_value + multiplier * mad_sigma)
    end
  end

  defp robust_calibration([]), do: %{samples: 0, median: nil, mad: nil, mad_sigma: nil}

  defp robust_calibration(values) do
    median_value = median(values)
    mad = values |> Enum.map(&abs(&1 - median_value)) |> median()

    %{
      samples: length(values),
      median: median_value,
      mad: mad,
      mad_sigma: 1.4826 * mad
    }
  end

  defp detector_events(rows, thresholds) do
    rows
    |> detector_steps()
    |> Enum.flat_map(&events_for_step(&1, thresholds))
    |> Enum.sort_by(&{&1.time, &1.satellite_id, Atom.to_string(&1.detector)})
  end

  defp detector_steps(rows) do
    rows
    |> Enum.group_by(& &1.satellite_id)
    |> Enum.flat_map(fn {_sat, sat_rows} ->
      sat_rows
      |> Enum.sort_by(&time_us(&1.epoch))
      |> Enum.chunk_every(2, 1, :discard)
      |> Enum.flat_map(fn [previous, current] ->
        dt_s = NaiveDateTime.diff(current.epoch, previous.epoch, :microsecond) / 1_000_000.0

        if dt_s > 0.0 and dt_s <= @continuity_gap_s do
          [detector_step(previous, current, dt_s)]
        else
          []
        end
      end)
    end)
  end

  defp detector_step(previous, current, dt_s) do
    code_phase_step_m = current.phase_minus_code_m - previous.phase_minus_code_m

    doppler_residual_m =
      if is_number(previous.primary_doppler_hz) and is_number(current.primary_doppler_hz) do
        current.primary_wavelength_m *
          (current.primary_phase_cycles - previous.primary_phase_cycles +
             0.5 * (previous.primary_doppler_hz + current.primary_doppler_hz) * dt_s)
      end

    dual_step =
      if previous.dual && current.dual do
        gf_step_m = current.dual.geometry_free_m - previous.dual.geometry_free_m
        mw_step_m = current.dual.melbourne_wubbena_m - previous.dual.melbourne_wubbena_m

        %{
          geometry_free_step_m: gf_step_m,
          melbourne_wubbena_step_m: mw_step_m,
          melbourne_wubbena_step_cycles: abs(mw_step_m) / abs(current.dual.wide_lane_wavelength_m)
        }
      end

    %{
      previous: previous,
      current: current,
      dt_s: dt_s,
      code_phase_step_m: code_phase_step_m,
      doppler_residual_m: doppler_residual_m,
      dual_step: dual_step
    }
  end

  defp events_for_step(step, thresholds) do
    current = step.current
    lambda = current.primary_wavelength_m

    []
    |> maybe_event(
      abs(step.code_phase_step_m) > thresholds.code_phase,
      current,
      :code_phase,
      step.code_phase_step_m,
      abs(step.code_phase_step_m) / lambda,
      step.dt_s
    )
    |> maybe_event(
      is_number(step.doppler_residual_m) and
        abs(step.doppler_residual_m) > thresholds.doppler_phase,
      current,
      :doppler_phase,
      step.doppler_residual_m,
      if(is_number(step.doppler_residual_m), do: abs(step.doppler_residual_m) / lambda),
      step.dt_s
    )
    |> maybe_dual_events(step, thresholds, current, lambda)
    |> maybe_event(
      lli_set?(current.primary_lli),
      current,
      :rinex_lli,
      0.0,
      0.0,
      step.dt_s
    )
  end

  defp maybe_dual_events(events, %{dual_step: nil}, _thresholds, _current, _lambda), do: events

  defp maybe_dual_events(events, step, thresholds, current, lambda) do
    dual = step.dual_step

    events
    |> maybe_event(
      abs(dual.geometry_free_step_m) > thresholds.geometry_free_m,
      current,
      :geometry_free,
      dual.geometry_free_step_m,
      abs(dual.geometry_free_step_m) / lambda,
      step.dt_s
    )
    |> maybe_event(
      dual.melbourne_wubbena_step_cycles > thresholds.melbourne_wubbena_cycles,
      current,
      :melbourne_wubbena,
      dual.melbourne_wubbena_step_m,
      dual.melbourne_wubbena_step_cycles,
      step.dt_s
    )
  end

  defp maybe_event(events, false, _row, _detector, _magnitude_m, _cycles, _dt_s), do: events

  defp maybe_event(events, true, row, detector, magnitude_m, cycles, dt_s) do
    [
      %{
        time: row.time,
        epoch: row.epoch,
        satellite_id: row.satellite_id,
        system: row.system,
        detector: detector,
        magnitude_m: abs(magnitude_m),
        signed_magnitude_m: magnitude_m,
        magnitude_cycles: cycles,
        dt_s: dt_s
      }
      | events
    ]
  end

  defp detector_ambiguity_ids(rows, events) do
    slip_keys = segmentation_event_keys(events)

    rows
    |> Enum.group_by(& &1.satellite_id)
    |> Enum.reduce(%{}, fn {sat, sat_rows}, acc ->
      {entries, _segment} =
        sat_rows
        |> Enum.sort_by(&time_us(&1.epoch))
        |> Enum.map_reduce(0, fn row, segment ->
          segment =
            if MapSet.member?(slip_keys, {row.time, sat}) do
              segment + 1
            else
              segment
            end

          id = if segment == 0, do: sat, else: "#{sat}@det##{segment}"
          {{row.time, sat, id}, segment}
        end)

      Enum.reduce(entries, acc, fn {time, sat_id, ambiguity_id}, id_acc ->
        Map.put(id_acc, {time, sat_id}, ambiguity_id)
      end)
    end)
  end

  defp detector_stats(rows, events) do
    exposure = detector_exposure_minutes(rows)
    by_detector = Enum.group_by(events, & &1.detector)

    [:code_phase, :geometry_free, :melbourne_wubbena, :doppler_phase, :rinex_lli]
    |> Map.new(fn detector ->
      detector_events = Map.get(by_detector, detector, [])
      exposure_minutes = Map.get(exposure, detector, 0.0)

      {detector,
       %{
         events: length(detector_events),
         exposure_satellite_minutes: exposure_minutes,
         events_per_satellite_minute: ratio(length(detector_events), exposure_minutes),
         magnitude_m_median: detector_events |> Enum.map(& &1.magnitude_m) |> median(),
         magnitude_m_p95: detector_events |> Enum.map(& &1.magnitude_m) |> percentile(0.95),
         magnitude_cycles_median:
           detector_events
           |> Enum.map(& &1.magnitude_cycles)
           |> Enum.reject(&is_nil/1)
           |> median(),
         magnitude_cycles_p95:
           detector_events
           |> Enum.map(& &1.magnitude_cycles)
           |> Enum.reject(&is_nil/1)
           |> percentile(0.95)
       }}
    end)
  end

  defp detector_exposure_minutes(rows) do
    steps = detector_steps(rows)
    primary_minutes = steps |> Enum.map(& &1.dt_s) |> Enum.sum() |> Kernel./(60.0)

    doppler_minutes =
      steps
      |> Enum.filter(&is_number(&1.doppler_residual_m))
      |> Enum.map(& &1.dt_s)
      |> Enum.sum()
      |> Kernel./(60.0)

    dual_minutes =
      steps
      |> Enum.filter(& &1.dual_step)
      |> Enum.map(& &1.dt_s)
      |> Enum.sum()
      |> Kernel./(60.0)

    %{
      code_phase: primary_minutes,
      rinex_lli: primary_minutes,
      doppler_phase: doppler_minutes,
      geometry_free: dual_minutes,
      melbourne_wubbena: dual_minutes
    }
  end

  defp detector_agreement(events) do
    by_key = Enum.group_by(events, &{&1.time, &1.satellite_id})
    detectors = [:code_phase, :geometry_free, :melbourne_wubbena, :doppler_phase, :rinex_lli]

    detector_sets =
      Map.new(by_key, fn {key, rows} ->
        {key, rows |> Enum.map(& &1.detector) |> Enum.uniq() |> Enum.sort()}
      end)

    pairwise =
      for a <- detectors, b <- detectors, Atom.to_string(a) < Atom.to_string(b), into: %{} do
        count =
          Enum.count(detector_sets, fn {_key, set} ->
            a in set and b in set
          end)

        {"#{a}+#{b}", count}
      end

    %{
      event_keys: map_size(detector_sets),
      multi_detector_keys: Enum.count(detector_sets, fn {_key, set} -> length(set) > 1 end),
      pairwise_overlap: pairwise
    }
  end

  defp detector_availability(rows, :dual) do
    detector_availability(rows, &(&1.dual != nil))
  end

  defp detector_availability(rows, :doppler) do
    detector_availability(rows, &is_number(&1.primary_doppler_hz))
  end

  defp detector_availability(rows, predicate) do
    rows
    |> Enum.group_by(& &1.system)
    |> Map.new(fn {system, system_rows} ->
      available = Enum.count(system_rows, predicate)

      {system,
       %{
         observations: length(system_rows),
         available_observations: available,
         available_fraction: ratio(available, length(system_rows))
       }}
    end)
  end

  defp detector_segmentation_summary(ambiguity_ids, events) do
    ids = ambiguity_ids |> Map.values() |> Enum.uniq()
    eligible_keys = eligible_segmentation_event_keys(events)
    injected_keys = segmentation_event_keys(events)

    %{
      event_keys: detector_agreement(events).event_keys,
      eligible_event_keys: MapSet.size(eligible_keys),
      injected_event_keys: MapSet.size(injected_keys),
      injection_policy:
        "raw LLI or independent agreement, capped at #{@detector_max_injected_events_per_satellite} injected keys per satellite for the runnable key experiment",
      split_ambiguity_id_count: Enum.count(ids, &String.contains?(&1, "@det#")),
      total_ambiguity_id_count: length(ids)
    }
  end

  defp segmentation_event_keys(events) do
    eligible = eligible_segmentation_event_keys(events)

    events
    |> Enum.map(&{&1.time, &1.satellite_id})
    |> Enum.uniq()
    |> Enum.filter(&MapSet.member?(eligible, &1))
    |> Enum.group_by(fn {_time, sat} -> sat end)
    |> Enum.flat_map(fn {_sat, keys} ->
      keys
      |> Enum.sort_by(fn {time, _sat} -> time end)
      |> Enum.take(@detector_max_injected_events_per_satellite)
    end)
    |> MapSet.new()
  end

  defp eligible_segmentation_event_keys(events) do
    events
    |> Enum.group_by(&{&1.time, &1.satellite_id})
    |> Enum.filter(fn {_key, rows} ->
      detectors = rows |> MapSet.new(& &1.detector)

      MapSet.member?(detectors, :rinex_lli) or
        (MapSet.member?(detectors, :code_phase) and
           MapSet.member?(detectors, :doppler_phase)) or
        (MapSet.member?(detectors, :doppler_phase) and
           MapSet.member?(detectors, :geometry_free)) or
        (MapSet.member?(detectors, :doppler_phase) and
           MapSet.member?(detectors, :melbourne_wubbena))
    end)
    |> MapSet.new(fn {key, _rows} -> key end)
  end

  defp lli_set?(lli) when is_integer(lli), do: Bitwise.band(lli, 1) == 1
  defp lli_set?(_lli), do: false

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

  defp build_filter_epochs(
         nav,
         rover_obs,
         base_obs,
         glonass_slots,
         oracle_by_time,
         base_arp,
         ambiguity_ids \\ %{},
         systems \\ @systems
       ) do
    base_index = base_epoch_index(base_obs)

    rover_obs
    |> Observations.epochs()
    |> Enum.reduce({[], []}, fn entry, {epochs, contexts} ->
      epoch = naive_datetime(entry.epoch)
      time_key = epoch_key(epoch)

      case Map.fetch(oracle_by_time, time_key) do
        {:ok, oracle_epoch} ->
          rover_values =
            phone_l1_values(
              rover_obs,
              entry.index,
              systems,
              glonass_slots,
              time_key,
              ambiguity_ids
            )

          base_values = interpolated_base_l1_values(base_index, epoch, systems, glonass_slots)

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

  defp segment_epoch_contexts(epochs, contexts, max_segment_epochs, max_sd_ambiguity_ids \\ nil) do
    epochs
    |> Enum.zip(contexts)
    |> Enum.reduce([], fn pair, segments ->
      case segments do
        [] ->
          [[pair]]

        [current | rest] ->
          candidate = current ++ [pair]

          if length(candidate) <= max_segment_epochs and segment_reference_solvable?(candidate) and
               segment_state_size_ok?(candidate, max_sd_ambiguity_ids) do
            [candidate | rest]
          else
            [[pair], current | rest]
          end
      end
    end)
    |> Enum.reverse()
  end

  defp segment_state_size_ok?(_pairs, nil), do: true

  defp segment_state_size_ok?(pairs, max_sd_ambiguity_ids) do
    pairs
    |> segment_single_difference_ids()
    |> length()
    |> Kernel.<=(max_sd_ambiguity_ids)
  end

  defp segment_single_difference_ids(pairs) do
    pairs
    |> Enum.flat_map(fn {epoch, _context} ->
      base = observation_map(epoch.base_observations)
      rover = observation_map(epoch.rover_observations)

      epoch.satellite_positions_m
      |> Map.keys()
      |> Enum.flat_map(fn sat ->
        with %{ambiguity_id: base_id} <- Map.get(base, sat),
             %{ambiguity_id: rover_id} <- Map.get(rover, sat) do
          [single_difference_ambiguity_id(sat, base_id, rover_id)]
        else
          _ -> []
        end
      end)
    end)
    |> Enum.uniq()
  end

  defp single_difference_ambiguity_id(sat, base_id, rover_id) do
    case {base_id, rover_id} do
      {^sat, ^sat} -> sat
      {^sat, rover_id} -> rover_id
      {base_id, ^sat} -> base_id
      {base_id, rover_id} when base_id == rover_id -> base_id
      {base_id, rover_id} -> "#{sat}:base=#{base_id},rover=#{rover_id}"
    end
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

  defp solve_segments(segments, nav, base_arp, glonass_slots, opts \\ []) do
    segments
    |> Enum.with_index()
    |> Enum.flat_map(fn {segment, index} ->
      solve_segment(segment, index, nav, base_arp, glonass_slots, opts)
    end)
    |> Enum.unzip()
    |> case do
      {measurements, reports} -> {List.flatten(measurements), reports}
    end
  end

  defp solve_segment(segment, index, nav, base_arp, glonass_slots, solve_opts) do
    epochs = Enum.map(segment, &elem(&1, 0))
    contexts = Enum.map(segment, &elem(&1, 1))
    initial_baseline = spp_initial_baseline!(nav, epochs, base_arp)

    opts =
      initial_baseline
      |> filter_opts(epochs, glonass_slots)
      |> maybe_override_filter_options(solve_opts)

    diagnostics? = Keyword.get(solve_opts, :diagnostics?, true)
    log? = Keyword.get(solve_opts, :log?, true)

    if log? do
      IO.puts(
        "  segment #{index + 1}: #{hd(contexts).time}..#{List.last(contexts).time} #{length(epochs)} epochs"
      )
    end

    case RTK.solve_filter_baseline_epochs(base_arp, epochs, opts) do
      {:ok, sol} ->
        state_diagnostics =
          if diagnostics?,
            do:
              filter_state_diagnostics(
                segment,
                sol,
                initial_baseline,
                opts,
                base_arp,
                Keyword.get(solve_opts, :diagnostic_epoch_limit)
              ),
            else: []

        state_by_time = Map.new(state_diagnostics, &{&1.time, &1})

        measurements =
          sol.epochs
          |> Enum.zip(contexts)
          |> Enum.map(fn {result, context} ->
            result
            |> epoch_measurement(context, base_arp)
            |> Map.put(:state_diagnostics, Map.get(state_by_time, context.time))
          end)

        report = %{
          index: index,
          epochs: length(epochs),
          first_time: hd(contexts).time,
          last_time: List.last(contexts).time,
          initial_baseline_m: initial_baseline,
          diagnostics: segment_diagnostics_summary(measurements, state_diagnostics),
          metadata: sol.metadata
        }

        [{measurements, report}]

      {:error, reason} when length(segment) > 1 ->
        if System.get_env("ROVER_MEAS_DEBUG_ERROR") == "1" do
          raise "segment #{index + 1} failed: #{inspect(reason)}"
        end

        {left, right} = Enum.split(segment, div(length(segment), 2))

        solve_segment(left, index, nav, base_arp, glonass_slots, solve_opts) ++
          solve_segment(right, index, nav, base_arp, glonass_slots, solve_opts)

      {:error, reason} ->
        context = hd(contexts)
        IO.puts("skipping unsolved epoch #{context.time}: #{inspect(reason)}")
        []
    end
  end

  defp maybe_override_filter_options(opts, solve_opts) do
    opts
    |> maybe_put_filter_option(solve_opts, :process_noise_baseline_sigma_m)
    |> maybe_put_filter_option(solve_opts, :filter_kernel)
    |> maybe_put_filter_option(solve_opts, :innovation_screen_sigma)
    |> maybe_put_filter_option(solve_opts, :innovation_screen_min_rows)
  end

  defp maybe_put_filter_option(opts, solve_opts, key) do
    case Keyword.fetch(solve_opts, key) do
      {:ok, value} -> Keyword.put(opts, key, value)
      :error -> opts
    end
  end

  defp filter_state_diagnostics(segment, sol, initial_baseline, opts, base_arp, epoch_limit) do
    epochs = Enum.map(segment, &elem(&1, 0))
    contexts = Enum.map(segment, &elem(&1, 1))

    diagnostic_pairs =
      case epoch_limit do
        limit when is_integer(limit) and limit > 0 -> Enum.take(Enum.zip(epochs, contexts), limit)
        _ -> Enum.zip(epochs, contexts)
      end

    diagnostic_epochs = Enum.map(diagnostic_pairs, &elem(&1, 0))
    diagnostic_contexts = Enum.map(diagnostic_pairs, &elem(&1, 1))

    all_sats =
      epochs |> Enum.flat_map(&Map.keys(&1.satellite_positions_m)) |> Enum.uniq() |> Enum.sort()

    refs = diagnostic_reference_satellites(base_arp, epochs, all_sats)
    sd_ids = all_sats
    initial_ambiguities = diagnostic_initial_ambiguities(epochs, sd_ids)

    state =
      diagnostic_initial_state(sd_ids, refs, epochs, initial_baseline, initial_ambiguities, opts)

    rust_epochs = Enum.map(diagnostic_epochs, &diagnostic_rust_epoch(&1, refs, all_sats))

    wavelengths =
      opts
      |> Keyword.fetch!(:ambiguity_wavelength_m)
      |> Map.drop(Map.values(refs))
      |> Enum.sort()

    offsets = Enum.map(wavelengths, fn {sat, _wavelength} -> {sat, 0.0} end)

    case NIF.rtk_filter_update_epochs(
           state,
           rust_epochs,
           ecef_to_tuple(base_arp),
           {Keyword.fetch!(opts, :code_sigma_m), Keyword.fetch!(opts, :phase_sigma_m), "rtklib",
            false, true},
           wavelengths,
           offsets,
           {
             Keyword.fetch!(opts, :hold_sigma_m),
             1.0e-4,
             1.0e-4,
             Keyword.fetch!(opts, :max_iterations),
             Keyword.fetch!(opts, :process_noise_baseline_sigma_m),
             Keyword.fetch!(opts, :integer_ratio_threshold),
             {Keyword.fetch!(opts, :float_only_systems), 0.0, 8}
           }
         ) do
      {:ok, updates} ->
        diagnostic_contexts
        |> Enum.zip(diagnostic_epochs)
        |> Enum.zip(updates)
        |> Enum.zip(sol.epochs)
        |> Enum.reduce({[], nil}, fn {{{context, epoch}, update}, result}, {acc, previous} ->
          diagnostic =
            diagnostic_epoch(
              context,
              epoch,
              update,
              result,
              refs,
              base_arp,
              previous
            )

          {[diagnostic | acc], diagnostic}
        end)
        |> elem(0)
        |> Enum.reverse()

      {:error, epoch_index, reason} ->
        raise "diagnostic NIF failed at epoch #{epoch_index}: #{inspect(reason)}"

      {:error, reason} ->
        raise "diagnostic NIF failed: #{inspect(reason)}"
    end
  end

  defp diagnostic_reference_satellites(base_arp, epochs, all_sats) do
    systems = all_sats |> Enum.map(&String.first/1) |> Enum.uniq() |> Enum.sort()
    common_by_system = per_system_common_sats(epochs)

    Map.new(systems, fn system ->
      sats = Map.fetch!(common_by_system, system) |> Enum.sort()

      reference =
        Enum.min_by(sats, fn sat ->
          {-average_reference_score(base_arp, sat, epochs), sat}
        end)

      {system, reference}
    end)
  end

  defp average_reference_score(base_arp, sat, epochs) do
    up = unit3(base_arp) || {0.0, 0.0, 1.0}

    values =
      epochs
      |> Enum.filter(&Map.has_key?(&1.satellite_positions_m, sat))
      |> Enum.map(fn epoch ->
        sat_pos = Map.fetch!(epoch.satellite_positions_m, sat)

        case unit3(sub3(sat_pos, base_arp)) do
          nil -> -1.0
          los -> dot3(los, up)
        end
      end)

    Enum.sum(values) / length(values)
  end

  defp unit3(v) do
    case norm3(v) do
      n when n > 0.0 -> scale3(v, 1.0 / n)
      _zero -> nil
    end
  end

  defp diagnostic_initial_ambiguities(epochs, sd_ids) do
    zero = Map.new(sd_ids, &{&1, 0.0})
    wanted = MapSet.new(sd_ids)

    seeded =
      Enum.reduce_while(epochs, %{}, fn epoch, acc ->
        if map_size(acc) == length(sd_ids) do
          {:halt, acc}
        else
          epoch_map = observation_map(epoch.base_observations)
          rover_map = observation_map(epoch.rover_observations)

          seeded_epoch =
            epoch.satellite_positions_m
            |> Map.keys()
            |> Enum.reduce(acc, fn sat, sat_acc ->
              with true <- MapSet.member?(wanted, sat),
                   false <- Map.has_key?(sat_acc, sat),
                   %{code_m: base_code, phase_m: base_phase} <- Map.get(epoch_map, sat),
                   %{code_m: rover_code, phase_m: rover_phase} <- Map.get(rover_map, sat) do
                code_sd = rover_code - base_code
                phase_sd = rover_phase - base_phase
                Map.put(sat_acc, sat, phase_sd - code_sd)
              else
                _ -> sat_acc
              end
            end)

          {:cont, seeded_epoch}
        end
      end)

    Map.merge(zero, seeded)
  end

  defp diagnostic_initial_state(sd_ids, refs, epochs, initial_baseline, initial_ambiguities, opts) do
    n = 3 + length(sd_ids)
    information = diagnostic_initial_information(n, opts)

    header_refs =
      refs
      |> Enum.sort()
      |> Enum.map(fn {system, reference_sat} ->
        epoch = Enum.find(epochs, &Map.has_key?(&1.satellite_positions_m, reference_sat))
        {system, reference_satellite_id(epoch, reference_sat)}
      end)

    {
      {@rtk_filter_state_version, header_refs, sd_ids,
       Keyword.fetch!(opts, :ambiguity_prior_sigma_m), 0},
      initial_baseline,
      Enum.map(sd_ids, &Map.fetch!(initial_ambiguities, &1)),
      List.flatten(information),
      [],
      []
    }
  end

  defp diagnostic_initial_information(n, opts) do
    baseline_sigma_m = Keyword.fetch!(opts, :baseline_prior_sigma_m)
    ambiguity_sigma_m = Keyword.fetch!(opts, :ambiguity_prior_sigma_m)

    for i <- 0..(n - 1) do
      for j <- 0..(n - 1) do
        cond do
          i != j -> 0.0
          i < 3 -> 1.0 / (baseline_sigma_m * baseline_sigma_m)
          true -> 1.0 / (ambiguity_sigma_m * ambiguity_sigma_m)
        end
      end
    end
  end

  defp diagnostic_rust_epoch(epoch, refs, all_sats) do
    available = epoch.satellite_positions_m |> Map.keys() |> MapSet.new()
    reference_set = refs |> Map.values() |> MapSet.new()

    references =
      refs
      |> Enum.sort()
      |> Enum.filter(fn {_system, sat} -> MapSet.member?(available, sat) end)
      |> Enum.map(fn {_system, sat} -> diagnostic_rust_sat(epoch, sat) end)

    nonrefs =
      all_sats
      |> Enum.reject(&MapSet.member?(reference_set, &1))
      |> Enum.filter(&MapSet.member?(available, &1))
      |> Enum.map(&diagnostic_rust_sat(epoch, &1))

    {references, nonrefs, nil, 0.0}
  end

  defp diagnostic_rust_sat(epoch, sat) do
    base = observation_map(epoch.base_observations) |> Map.fetch!(sat)
    rover = observation_map(epoch.rover_observations) |> Map.fetch!(sat)

    {
      {sat, reference_satellite_id(epoch, sat)},
      {base.code_m, base.phase_m, rover.code_m, rover.phase_m},
      {
        Map.fetch!(epoch.base_satellite_positions_m, sat),
        Map.fetch!(epoch.rover_satellite_positions_m, sat),
        Map.fetch!(epoch.satellite_positions_m, sat)
      }
    }
  end

  defp reference_satellite_id(_epoch, sat), do: sat

  defp diagnostic_epoch(context, epoch, update, result, refs, base_arp, previous) do
    {state_term, reported_baseline, ratio, fixed?, newly_fixed, fixed_ids, _screen} = update

    {{@rtk_filter_state_version, header_refs, sd_ids, _ambiguity_sigma_m, state_epoch_count},
     carried_baseline, _sd_ambiguities, information_flat, fixed_cycles, _fixed_m} = state_term

    n = 3 + length(sd_ids)
    information = Enum.chunk_every(information_flat, n)
    truth_tuple = truth_tuple(context)
    carried_rover = carried_baseline |> add3(ecef_to_tuple(base_arp))
    reported_rover = reported_baseline |> add3(ecef_to_tuple(base_arp))
    sats = epoch.satellite_positions_m |> Map.keys() |> Enum.sort()
    previous_sats = (previous && previous.satellite_ids) || []
    previous_set = MapSet.new(previous_sats)
    current_set = MapSet.new(sats)
    current_systems = sats |> Enum.map(&String.first/1) |> Enum.uniq()

    %{
      time: context.time,
      segment_epoch_index: result.index,
      state_epoch_count: state_epoch_count,
      reference_satellites: Map.new(header_refs),
      expected_reference_satellites: refs,
      reference_satellites_present?:
        Enum.all?(current_systems, fn system ->
          refs |> Map.fetch!(system) |> then(&MapSet.member?(current_set, &1))
        end),
      satellite_ids: sats,
      satellites_added:
        current_set |> MapSet.difference(previous_set) |> MapSet.to_list() |> Enum.sort(),
      satellites_removed:
        previous_set |> MapSet.difference(current_set) |> MapSet.to_list() |> Enum.sort(),
      gap_s:
        previous && NaiveDateTime.diff(context.epoch, previous.epoch, :microsecond) / 1_000_000.0,
      epoch: context.epoch,
      sd_ambiguity_columns: length(sd_ids),
      hold_count: length(fixed_cycles),
      fixed_sd_ids: Enum.sort(fixed_ids),
      newly_fixed_sd_ids: Enum.sort(newly_fixed),
      information_condition_estimate: information_condition_estimate(information),
      carried_baseline_error_3d_m: norm3(sub3(carried_rover, truth_tuple)),
      reported_baseline_error_3d_m: norm3(sub3(reported_rover, truth_tuple)),
      nif_integer_fixed?: fixed?,
      nif_integer_ratio: finite_ratio(ratio)
    }
  end

  defp truth_tuple(context) do
    truth = context.oracle_epoch["truth_ecef_m"]
    {truth["x"], truth["y"], truth["z"]}
  end

  defp information_condition_estimate(matrix) do
    row_sums =
      matrix
      |> Enum.map(fn row -> row |> Enum.map(&abs/1) |> Enum.sum() end)
      |> Enum.reject(&(&1 == 0.0))

    case row_sums do
      [] -> nil
      values -> Enum.max(values) / Enum.min(values)
    end
  end

  defp segment_diagnostics_summary(measurements, state_diagnostics) do
    first_bad = first_bad_epoch(measurements)

    %{
      first_bad_epoch: first_bad && first_bad_excerpt(first_bad),
      max_reported_error_3d_m: measurements |> Enum.map(& &1.error_3d_m) |> max_or_nil(),
      max_carried_error_3d_m:
        state_diagnostics
        |> Enum.map(& &1.carried_baseline_error_3d_m)
        |> max_or_nil(),
      max_information_condition_estimate:
        state_diagnostics
        |> Enum.map(& &1.information_condition_estimate)
        |> reject_infinity()
        |> max_or_nil(),
      max_sd_ambiguity_columns:
        state_diagnostics
        |> Enum.map(& &1.sd_ambiguity_columns)
        |> max_or_nil(),
      max_hold_count: state_diagnostics |> Enum.map(& &1.hold_count) |> max_or_nil()
    }
  end

  defp reject_infinity(values), do: Enum.reject(values, &(&1 == :infinity))

  defp base_epoch_index(%Rinex2Obs{epochs: epochs}) do
    sorted = Enum.sort_by(epochs, &time_us(&1.epoch))

    %{
      times: sorted |> Enum.map(&time_us(&1.epoch)) |> List.to_tuple(),
      epochs: List.to_tuple(sorted),
      count: length(sorted)
    }
  end

  defp phone_l1_values(obs, index, systems, glonass_slots, time_key, ambiguity_ids) do
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
             ambiguity_id: Map.get(ambiguity_ids, {time_key, sat}, sat),
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
      innovation_screen: Map.get(result, :innovation_screen),
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

  defp arc_diagnosis(per_epoch, segment_reports, sanity_gate, time_alignment) do
    first_bad = first_bad_epoch(per_epoch)

    if first_bad do
      index = Enum.find_index(per_epoch, &(&1.time == first_bad.time))
      previous = if index && index > 0, do: Enum.at(per_epoch, index - 1)
      changes = first_bad_changes(first_bad, previous)
      verdict = first_bad_verdict(first_bad, changes)

      %{
        threshold_m: @divergence_threshold_m,
        verdict: verdict,
        mechanism: first_bad_mechanism(first_bad, changes),
        input_consistency:
          input_consistency_summary(sanity_gate, time_alignment, first_bad, changes),
        first_bad_epoch: first_bad_excerpt(first_bad),
        previous_epoch: previous && first_bad_excerpt(previous),
        changes_at_first_bad: changes,
        segment: segment_for_epoch(segment_reports, first_bad.time)
      }
    else
      %{
        threshold_m: @divergence_threshold_m,
        verdict: "no_megameter_divergence",
        mechanism: "no epoch crossed the divergence threshold in the completed multi-epoch run",
        input_consistency: input_consistency_summary(sanity_gate, time_alignment, nil, %{}),
        first_bad_epoch: nil,
        previous_epoch: nil,
        changes_at_first_bad: %{},
        segment: nil
      }
    end
  end

  defp first_bad_epoch(per_epoch), do: Enum.find(per_epoch, &bad_epoch?/1)

  defp bad_epoch?(epoch) do
    state = epoch.state_diagnostics || %{}

    epoch.error_3d_m >= @divergence_threshold_m or
      Map.get(state, :carried_baseline_error_3d_m, 0.0) >= @divergence_threshold_m
  end

  defp first_bad_excerpt(epoch) do
    state = epoch.state_diagnostics || %{}

    %{
      time: epoch.time,
      error_3d_m: epoch.error_3d_m,
      carried_baseline_error_3d_m: Map.get(state, :carried_baseline_error_3d_m),
      information_condition_estimate: Map.get(state, :information_condition_estimate),
      segment_epoch_index: Map.get(state, :segment_epoch_index),
      sd_ambiguity_columns: Map.get(state, :sd_ambiguity_columns),
      hold_count: Map.get(state, :hold_count),
      fixed_ambiguities: epoch.fixed_ambiguities,
      newly_fixed_ambiguities: epoch.newly_fixed_ambiguities,
      satellites: epoch.satellites,
      pre_mask_satellites: epoch.pre_mask_satellites,
      max_abs_code_residual_m: epoch.residuals.max_abs_code_m,
      max_abs_phase_residual_m: epoch.residuals.max_abs_phase_m
    }
  end

  defp first_bad_changes(first_bad, previous) do
    state = first_bad.state_diagnostics || %{}
    previous_state = (previous && previous.state_diagnostics) || %{}

    %{
      gap_s: Map.get(state, :gap_s),
      satellites_added: Map.get(state, :satellites_added, []),
      satellites_removed: Map.get(state, :satellites_removed, []),
      reference_satellites_present?: Map.get(state, :reference_satellites_present?),
      references: Map.get(state, :reference_satellites, %{}),
      newly_fixed_sd_ids: Map.get(state, :newly_fixed_sd_ids, []),
      previous_newly_fixed_sd_ids: Map.get(previous_state, :newly_fixed_sd_ids, []),
      hold_count: Map.get(state, :hold_count),
      previous_hold_count: Map.get(previous_state, :hold_count),
      sd_ambiguity_columns: Map.get(state, :sd_ambiguity_columns),
      previous_sd_ambiguity_columns: Map.get(previous_state, :sd_ambiguity_columns),
      segmented_arc_ids_present?: segmented_arc_ids_present?(state)
    }
  end

  defp first_bad_verdict(_first_bad, %{reference_satellites_present?: false}), do: "harness_bug"

  defp first_bad_verdict(_first_bad, _changes), do: "filter_behavior"

  defp first_bad_mechanism(_first_bad, %{reference_satellites_present?: false}) do
    "segment admitted an epoch without its selected reference satellite"
  end

  defp first_bad_mechanism(_first_bad, %{segmented_arc_ids_present?: true}) do
    "segmented ambiguity ids grew inside the carried filter state"
  end

  defp first_bad_mechanism(_first_bad, changes) do
    added = Map.get(changes, :satellites_added, [])
    removed = Map.get(changes, :satellites_removed, [])
    hold_count = Map.get(changes, :hold_count) || 0
    previous_hold_count = Map.get(changes, :previous_hold_count) || 0
    newly_fixed = Map.get(changes, :newly_fixed_sd_ids, [])
    previous_newly_fixed = Map.get(changes, :previous_newly_fixed_sd_ids, [])

    cond do
      hold_count > 0 and (added != [] or removed != []) ->
        "tight ambiguity holds remain active across satellite-set churn"

      newly_fixed != [] ->
        "integer hold accepted at the first divergent epoch"

      previous_newly_fixed != [] ->
        "integer hold accepted immediately before the first divergent epoch"

      hold_count > previous_hold_count ->
        "hold set expanded immediately before divergence"

      added != [] or removed != [] ->
        "carried float state diverged at a satellite-set change before any harness inconsistency"

      true ->
        "carried float state diverged without an input inconsistency marker"
    end
  end

  defp segmented_arc_ids_present?(state) do
    state
    |> Map.get(:satellite_ids, [])
    |> Enum.any?(&String.contains?(&1, "~ra"))
  end

  defp input_consistency_summary(sanity_gate, time_alignment, _first_bad, changes) do
    cond do
      Map.get(changes, :reference_satellites_present?) == false ->
        "failed: selected reference absent at the first bad epoch"

      sanity_gate.pass and time_alignment.rinex_to_oracle_time_max_ms <= 0.5 ->
        "passed: meter-level SPP residual gate and sub-ms RINEX/oracle alignment"

      sanity_gate.pass ->
        "passed residual gate; time alignment should be reviewed"

      true ->
        "failed residual gate"
    end
  end

  defp segment_for_epoch(segment_reports, time) do
    Enum.find_value(segment_reports, fn segment ->
      if segment.first_time <= time and time <= segment.last_time do
        %{
          index: segment.index,
          epochs: segment.epochs,
          first_time: segment.first_time,
          last_time: segment.last_time,
          diagnostics: segment.diagnostics
        }
      end
    end)
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
    demo5_epochs = Enum.flat_map(arcs, & &1.input.oracle["per_epoch"])
    stock_epochs = Enum.flat_map(arcs, & &1.stock.per_epoch)
    segmented_epochs = Enum.flat_map(arcs, & &1.segmented.per_epoch)
    comparison_epochs = Enum.flat_map(arcs, & &1.comparison.per_epoch)

    demo5 = demo5_summary(demo5_epochs)
    stock = summarize_measurements(stock_epochs)
    segmented = summarize_measurements(segmented_epochs)
    comparison = summarize_measurements(comparison_epochs)

    promotion = promotion_verdict(arcs, segmented, comparison)

    %{
      "demo5" => stringify_keys(demo5),
      "stock" =>
        stringify_keys(%{
          orbis: stock,
          comparative: comparative_verdict(stock, demo5, :pooled),
          invariant: invariant_verdict(stock_epochs)
        }),
      "segmented" =>
        stringify_keys(%{
          orbis: segmented,
          comparative: comparative_verdict(segmented, demo5, :pooled),
          invariant: invariant_verdict(segmented_epochs)
        }),
      "comparison" =>
        stringify_keys(%{
          mode: "per_epoch_segments_comparison_only",
          orbis: comparison,
          comparative: comparative_verdict(comparison, demo5, :pooled),
          invariant: invariant_verdict(comparison_epochs)
        }),
      "process_noise_sweep" => stringify_keys(pooled_process_noise_sweep(arcs)),
      "cost" => stringify_keys(pooled_cost(arcs)),
      "promotion" => stringify_keys(promotion)
    }
  end

  defp pooled_process_noise_sweep(arcs) do
    for mode <- ["stock", "detector"], sigma_m <- @process_noise_sweep_sigmas do
      cells =
        Enum.map(arcs, fn arc ->
          Enum.find(arc.process_noise_sweep, fn cell ->
            cell.mode == mode and cell.process_noise_baseline_sigma_m == sigma_m
          end)
        end)

      completed_cells = Enum.filter(cells, & &1.completed_arc)
      completed_errors = Enum.flat_map(completed_cells, & &1.errors_3d_m)

      %{
        mode: mode,
        process_noise_baseline_sigma_m: sigma_m,
        completed_arcs: length(completed_cells),
        total_arcs: length(cells),
        all_arcs_completed: length(completed_cells) == length(cells),
        min_stable_prefix_epochs:
          cells |> Enum.map(& &1.longest_stable_prefix_epochs) |> Enum.min(),
        max_stable_prefix_epochs:
          cells |> Enum.map(& &1.longest_stable_prefix_epochs) |> Enum.max(),
        completed_pooled_median_m: median(completed_errors),
        completed_pooled_p95_m: percentile(completed_errors, 0.95),
        epoch2_condition_min:
          cells |> Enum.map(&epoch2_condition/1) |> Enum.reject(&is_nil/1) |> min_or_nil(),
        epoch2_condition_max:
          cells |> Enum.map(&epoch2_condition/1) |> Enum.reject(&is_nil/1) |> max_or_nil()
      }
    end
  end

  defp epoch2_condition(%{epoch2: nil}), do: nil

  defp epoch2_condition(cell) do
    cell.epoch2.information_condition_estimate
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

  defp pooled_cost(arcs) do
    %{
      detector_event_keys: arcs |> Enum.map(& &1.cost.detector_event_keys) |> Enum.sum(),
      detector_eligible_event_keys:
        arcs |> Enum.map(& &1.cost.detector_eligible_event_keys) |> Enum.sum(),
      detector_injected_event_keys:
        arcs |> Enum.map(& &1.cost.detector_injected_event_keys) |> Enum.sum(),
      detector_events: arcs |> Enum.map(& &1.cost.detector_events) |> Enum.sum(),
      split_ambiguity_ids: arcs |> Enum.map(& &1.cost.split_ambiguity_ids) |> Enum.sum(),
      detector_epochs_lost: arcs |> Enum.map(& &1.cost.detector_epochs_lost) |> Enum.sum(),
      stock_epochs_lost: arcs |> Enum.map(& &1.cost.stock_epochs_lost) |> Enum.sum(),
      max_detector_sd_columns:
        arcs
        |> Enum.map(& &1.cost.detector_sd_columns.max_single_difference_columns)
        |> Enum.reject(&is_nil/1)
        |> max_or_nil(),
      max_stock_sd_columns:
        arcs
        |> Enum.map(& &1.cost.stock_sd_columns.max_single_difference_columns)
        |> Enum.reject(&is_nil/1)
        |> max_or_nil()
    }
  end

  defp promotion_verdict(arcs, segmented, comparison) do
    min_prefix = arcs |> Enum.map(& &1.segmented.longest_stable_prefix_epochs) |> Enum.min()
    stock_min_prefix = arcs |> Enum.map(& &1.stock.longest_stable_prefix_epochs) |> Enum.min()
    beats_reference = segmented.error_3d_median_m < @stable_reference_median_m
    beats_one_epoch = segmented.error_3d_median_m < comparison.error_3d_median_m
    survives = min_prefix > stock_min_prefix and min_prefix > 1
    promote = survives and beats_reference

    %{
      verdict: if(promote, do: "promote_c1", else: "do_not_promote_c1"),
      survives_beyond_stock_prefix: survives,
      stock_min_stable_prefix_epochs: stock_min_prefix,
      segmented_min_stable_prefix_epochs: min_prefix,
      segmented_pooled_median_m: segmented.error_3d_median_m,
      one_epoch_pooled_median_m: comparison.error_3d_median_m,
      reference_single_epoch_median_m: @stable_reference_median_m,
      beats_reference_single_epoch_median: beats_reference,
      beats_current_one_epoch_median: beats_one_epoch
    }
  end

  defp c2_pooled_summary(arcs) do
    demo5_epochs = Enum.flat_map(arcs, & &1.input.oracle["per_epoch"])
    comparison_epochs = Enum.flat_map(arcs, & &1.comparison.per_epoch)
    demo5 = demo5_summary(demo5_epochs)
    comparison = summarize_measurements(comparison_epochs)
    sweep = c2_pooled_sweep(arcs, demo5, comparison)

    %{
      demo5: demo5,
      comparison: comparison,
      innovation_screen_sweep: sweep,
      verdict: c2_verdict(sweep, demo5, comparison)
    }
  end

  defp c2_pooled_sweep(arcs, demo5, comparison) do
    for sigma <- @innovation_screen_sigmas do
      cells =
        Enum.map(arcs, fn arc ->
          Enum.find(arc.innovation_screen_sweep, &(&1.innovation_screen_sigma == sigma))
        end)

      completed_cells = Enum.filter(cells, & &1.completed_arc)
      completed_errors = Enum.flat_map(completed_cells, & &1.errors_3d_m)
      rejected_counts = Enum.flat_map(cells, & &1.rows_rejected_counts)
      coasted_epochs = cells |> Enum.map(& &1.coasted_epochs) |> Enum.sum()
      solved_epochs = cells |> Enum.map(& &1.solved_epoch_count) |> Enum.sum()
      rejected_row_count = cells |> Enum.map(& &1.rejected_row_count) |> Enum.sum()
      input_row_count = cells |> Enum.map(& &1.input_row_count) |> Enum.sum()
      median_m = median(completed_errors)

      %{
        innovation_screen_sigma: sigma,
        completed_arcs: length(completed_cells),
        total_arcs: length(cells),
        all_arcs_completed: length(completed_cells) == length(cells),
        min_stable_prefix_epochs:
          cells |> Enum.map(& &1.longest_stable_prefix_epochs) |> Enum.min(),
        max_stable_prefix_epochs:
          cells |> Enum.map(& &1.longest_stable_prefix_epochs) |> Enum.max(),
        completed_pooled_median_m: median_m,
        completed_pooled_p95_m: percentile(completed_errors, 0.95),
        coasted_epochs: coasted_epochs,
        solved_epoch_count: solved_epochs,
        coasting_fraction: ratio(coasted_epochs, max(solved_epochs, 1)),
        rejected_row_count: rejected_row_count,
        input_row_count: input_row_count,
        row_rejection_fraction: ratio(rejected_row_count, max(input_row_count, 1)),
        rows_rejected_per_epoch: distribution(rejected_counts),
        beats_one_epoch_bar: is_number(median_m) and median_m <= @stable_reference_median_m,
        beats_current_one_epoch: is_number(median_m) and median_m <= comparison.error_3d_median_m,
        beats_demo5_x_1_25: is_number(median_m) and median_m <= demo5.error_3d_median_m * 1.25
      }
    end
  end

  defp c2_verdict(sweep, demo5, comparison) do
    best =
      sweep
      |> Enum.filter(&is_number(&1.completed_pooled_median_m))
      |> Enum.min_by(& &1.completed_pooled_median_m, fn -> nil end)

    %{
      best_k: best && best.innovation_screen_sigma,
      best_completed_pooled_median_m: best && best.completed_pooled_median_m,
      one_epoch_reference_median_m: @stable_reference_median_m,
      measured_one_epoch_pooled_median_m: comparison.error_3d_median_m,
      demo5_reference_median_m: @demo5_reference_median_m,
      measured_demo5_pooled_median_m: demo5.error_3d_median_m,
      beats_report_bar: best && best.completed_pooled_median_m <= @stable_reference_median_m,
      promotion_grade: best && best.completed_pooled_median_m <= demo5.error_3d_median_m * 1.25
    }
  end

  defp c2_arc_json(arc) do
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
      "built_epoch_count" => arc.built_epoch_count,
      "skipped_oracle_epochs" => arc.skipped_oracle_epochs,
      "detector_segmentation" => stringify_keys(arc.detector.segmentation),
      "comparison" => %{
        "mode" => arc.comparison.mode,
        "segment_count" => arc.comparison.segment_count,
        "orbis" => stringify_keys(arc.comparison.orbis)
      },
      "demo5" => stringify_keys(arc.demo5),
      "innovation_screen_sweep" =>
        Enum.map(arc.innovation_screen_sweep, &innovation_screen_sweep_json/1)
    }
  end

  defp innovation_screen_sweep_json(cell) do
    %{
      "innovation_screen_sigma" => cell.innovation_screen_sigma,
      "innovation_screen_min_rows" => cell.innovation_screen_min_rows,
      "built_epoch_count" => cell.built_epoch_count,
      "solved_epoch_count" => cell.solved_epoch_count,
      "completed_arc" => cell.completed_arc,
      "segment_count" => cell.segment_count,
      "longest_stable_prefix_epochs" => cell.longest_stable_prefix_epochs,
      "first_unstable_time" => cell.first_unstable_time,
      "orbis" => stringify_keys(cell.orbis),
      "completed_arc_median_m" => cell.completed_arc_median_m,
      "completed_arc_p95_m" => cell.completed_arc_p95_m,
      "coasted_epochs" => cell.coasted_epochs,
      "coasting_fraction" => cell.coasting_fraction,
      "rejected_row_count" => cell.rejected_row_count,
      "input_row_count" => cell.input_row_count,
      "row_rejection_fraction" => cell.row_rejection_fraction,
      "rows_rejected_per_epoch" => stringify_keys(cell.rows_rejected_per_epoch),
      "rows_accepted_per_epoch" => stringify_keys(cell.rows_accepted_per_epoch),
      "rows_rejected_counts" => cell.rows_rejected_counts,
      "rows_accepted_counts" => cell.rows_accepted_counts,
      "errors_3d_m" => cell.errors_3d_m,
      "error_series" => Enum.map(cell.error_series, &stringify_keys/1)
    }
  end

  defp best_k_error_series(_arcs, nil), do: []

  defp best_k_error_series(arcs, best_k) do
    Enum.map(arcs, fn arc ->
      cell =
        Enum.find(arc.innovation_screen_sweep, &(&1.innovation_screen_sigma == best_k))

      one_epoch_median_m = arc.comparison.orbis.error_3d_median_m
      series = cell.error_series
      below = Enum.filter(series, &(&1.error_3d_m <= one_epoch_median_m))
      first_below = Enum.find(series, &(&1.error_3d_m <= one_epoch_median_m))
      tail = Enum.take(series, -min(length(series), 100))
      tail_errors = Enum.map(tail, & &1.error_3d_m)
      all_errors = Enum.map(series, & &1.error_3d_m)

      %{
        "label" => arc.input.label,
        "innovation_screen_sigma" => best_k,
        "one_epoch_median_m" => one_epoch_median_m,
        "demo5_median_m" => arc.demo5.error_3d_median_m,
        "epoch_count" => length(series),
        "below_one_epoch_count" => length(below),
        "below_one_epoch_fraction" => ratio(length(below), max(length(series), 1)),
        "first_below_one_epoch_index" =>
          first_below && Enum.find_index(series, &(&1 == first_below)),
        "first_below_one_epoch_time" => first_below && first_below.time,
        "tail_100_median_m" => median(tail_errors),
        "tail_100_p95_m" => percentile(tail_errors, 0.95),
        "min_error_m" => min_or_nil(all_errors),
        "last_error_m" => List.last(all_errors),
        "series" => Enum.map(series, &stringify_keys/1)
      }
    end)
  end

  defp c2_report_markdown(result, results_path) do
    arcs = result["arcs"]
    pooled = result["pooled"]

    [
      "# C2 kernel innovation screening exploration, June 2026",
      "",
      "Generated by `mix run test/fixtures/rtk/generators/rover_c2_innovation_screening_2026_06.exs`.",
      "Per-run JSON was emitted to `#{results_path}`.",
      "",
      "This is a single-sided exploration artifact. It ports per-epoch innovation screening into the Rust RTK filter kernel on branch `explore/c2-kernel-screen`.",
      "",
      "## Setup",
      "",
      "- Inputs match C1: phone RINEX from `/tmp/gsdc-work`, P222 CORS base observations, and BKG combined broadcast NAV.",
      "- Detector segmentation is ON via the C1 detector-injected rover ambiguity ids.",
      epoch_limit_note(result),
      "- The measured path uses the Rust sequential filter (`filter_kernel: :rust`) with `innovation_screen_sigma` in #{Enum.map_join(@innovation_screen_sigmas, ", ", &fmt/1)} and `innovation_screen_min_rows: #{@innovation_screen_min_rows}`.",
      "- The screen rejects individual code/phase DD rows whose absolute normalized predicted residual is greater than k before the measurement update; if fewer than #{@innovation_screen_min_rows} rows survive, the epoch coasts on the prediction.",
      "",
      "## K Sweep",
      "",
      c2_k_sweep_table(arcs, pooled),
      "",
      "Median and p95 columns are reported only for arcs that solved every built epoch. Pooled values combine the completed arcs for that k.",
      "",
      "## Verdict",
      "",
      c2_verdict_text(pooled),
      "",
      "## Convergence",
      "",
      c2_convergence_text(result),
      "",
      c2_convergence_table(result),
      "",
      "## Costs",
      "",
      c2_cost_table(pooled),
      "",
      "Rows rejected per epoch is pooled across all solved epochs for each k. Coasting fraction is coasted epochs divided by solved epochs. Rejected row % is rejected rows divided by input screen rows.",
      ""
    ]
    |> Enum.join("\n")
  end

  defp epoch_limit_note(%{"epoch_limit" => limit}) when is_integer(limit) do
    "- Runtime boundary: this report uses the first #{limit} built epochs per arc. Full-arc Elixir sweeps were not tractable in this exploration worktree."
  end

  defp epoch_limit_note(_result) do
    "- Runtime boundary: this report uses every built epoch in each arc."
  end

  defp c2_k_sweep_table(arcs, pooled) do
    pooled_by_k = Map.new(pooled["innovation_screen_sweep"], &{&1["innovation_screen_sigma"], &1})

    rows =
      Enum.map(@innovation_screen_sigmas, fn sigma ->
        cells =
          Enum.map(arcs, fn arc ->
            Enum.find(arc["innovation_screen_sweep"], &(&1["innovation_screen_sigma"] == sigma))
          end)

        pooled_cell = Map.fetch!(pooled_by_k, sigma)

        [
          fmt(sigma),
          "#{pooled_cell["completed_arcs"]}/#{pooled_cell["total_arcs"]}",
          Enum.map_join(cells, "/", & &1["longest_stable_prefix_epochs"]),
          Enum.map_join(cells, "/", &completed_metric(&1, "completed_arc_median_m")),
          Enum.map_join(cells, "/", &completed_metric(&1, "completed_arc_p95_m")),
          "#{fmt(pooled_cell["completed_pooled_median_m"])}/#{fmt(pooled_cell["completed_pooled_p95_m"])}",
          "#{pooled_cell["coasted_epochs"]}/#{fmt(100.0 * pooled_cell["coasting_fraction"])}%",
          "#{fmt(100.0 * pooled_cell["row_rejection_fraction"])}%",
          rejected_distribution_text(pooled_cell["rows_rejected_per_epoch"]),
          yes_no(pooled_cell["beats_one_epoch_bar"]),
          yes_no(pooled_cell["beats_demo5_x_1_25"])
        ]
        |> table_row()
      end)

    Enum.join(
      [
        "| k | Complete arcs | Stable prefix SJC1/SVL1/MTV15/MTV28 | Completed median m SJC1/SVL1/MTV15/MTV28 | Completed p95 m SJC1/SVL1/MTV15/MTV28 | Pooled med/p95 m | Coasted epochs/% | Rejected row % | Rejected rows/epoch med/p95/max | <=9.5m | <=demo5*1.25 |",
        "|---:|---:|---|---|---|---:|---:|---:|---:|---|---|"
      ] ++ rows,
      "\n"
    )
  end

  defp c2_cost_table(pooled) do
    rows =
      Enum.map(pooled["innovation_screen_sweep"], fn cell ->
        [
          fmt(cell["innovation_screen_sigma"]),
          cell["coasted_epochs"],
          fmt(100.0 * cell["coasting_fraction"]),
          fmt(100.0 * cell["row_rejection_fraction"]),
          rejected_distribution_text(cell["rows_rejected_per_epoch"])
        ]
        |> table_row()
      end)

    Enum.join(
      [
        "| k | Coasted epochs | Coasting % | Rejected row % | Rejected rows/epoch med/p95/max |",
        "|---:|---:|---:|---:|---:|"
      ] ++ rows,
      "\n"
    )
  end

  defp rejected_distribution_text(%{"median" => median, "p95" => p95, "max" => max}) do
    "#{fmt(median)}/#{fmt(p95)}/#{fmt(max)}"
  end

  defp c2_convergence_text(result) do
    series = result["best_k_error_series"]
    verdict = result["pooled"]["verdict"]
    best_k = verdict["best_k"]

    if series == [] do
      "No best-k time series was emitted because no k cell completed."
    else
      arc_passes =
        Enum.count(series, fn arc ->
          arc["tail_100_median_m"] <= arc["one_epoch_median_m"]
        end)

      tail_errors =
        Enum.flat_map(series, fn arc ->
          arc["series"]
          |> Enum.take(-min(length(arc["series"]), 100))
          |> Enum.map(& &1["error_3d_m"])
        end)

      pooled_one_epoch = verdict["measured_one_epoch_pooled_median_m"]
      pooled_tail_median = median(tail_errors)
      pooled_tail_p95 = percentile(tail_errors, 0.95)

      answer =
        if arc_passes == length(series) and pooled_tail_median <= pooled_one_epoch do
          "yes"
        else
          "no"
        end

      "Best k=#{fmt(best_k)} time series are emitted in `best_k_error_series` in the JSON. Convergence below the one-epoch baseline after the transient: #{answer}. Final-100 pooled median/p95 is #{fmt(pooled_tail_median)}/#{fmt(pooled_tail_p95)} m versus measured one-epoch pooled median #{fmt(pooled_one_epoch)} m; #{arc_passes}/#{length(series)} arcs have final-100 median at or below their own one-epoch median."
    end
  end

  defp c2_convergence_table(result) do
    rows =
      Enum.map(result["best_k_error_series"], fn arc ->
        [
          short_arc_label(arc),
          arc["epoch_count"],
          fmt(arc["one_epoch_median_m"]),
          first_below_one_epoch_text(arc),
          "#{fmt(100.0 * arc["below_one_epoch_fraction"])}%",
          "#{fmt(arc["tail_100_median_m"])}/#{fmt(arc["tail_100_p95_m"])}",
          "#{fmt(arc["min_error_m"])}/#{fmt(arc["last_error_m"])}"
        ]
        |> table_row()
      end)

    Enum.join(
      [
        "| Arc | Epochs | One-epoch median m | First <= one-epoch | Below one-epoch % | Final-100 med/p95 m | Min/last m |",
        "|---|---:|---:|---|---:|---:|---:|"
      ] ++ rows,
      "\n"
    )
  end

  defp first_below_one_epoch_text(%{"first_below_one_epoch_index" => nil}), do: "never"

  defp first_below_one_epoch_text(%{
         "first_below_one_epoch_index" => index,
         "first_below_one_epoch_time" => time
       }) do
    "#{index} @ #{time}"
  end

  defp c2_verdict_text(pooled) do
    verdict = pooled["verdict"]
    best_k = verdict["best_k"]
    best_median = verdict["best_completed_pooled_median_m"]

    cond do
      is_nil(best_k) ->
        "Verdict: no completed k cell. The prototype did not produce a usable four-arc carried-state trajectory."

      verdict["promotion_grade"] ->
        "Verdict: promotion-grade. k=#{fmt(best_k)} produced pooled median #{fmt(best_median)} m, beating the 9.5 m one-epoch bar and the demo5 x 1.25 bar (demo5 pooled median #{fmt(verdict["measured_demo5_pooled_median_m"])} m)."

      verdict["beats_report_bar"] ->
        "Verdict: passes the report bar but not promotion-grade. k=#{fmt(best_k)} produced pooled median #{fmt(best_median)} m, beating 9.5 m but not demo5 x 1.25 (demo5 pooled median #{fmt(verdict["measured_demo5_pooled_median_m"])} m)."

      true ->
        "Verdict: no. Best k=#{fmt(best_k)} produced pooled median #{fmt(best_median)} m versus the 9.5 m one-epoch bar and demo5 pooled median #{fmt(verdict["measured_demo5_pooled_median_m"])} m."
    end
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
      "detector" => detector_json(arc.detector),
      "stock" => filter_case_json(arc.stock, arc.demo5),
      "segmented" => filter_case_json(arc.segmented, arc.demo5),
      "comparison" => comparison_json(arc.comparison, arc.demo5),
      "demo5" => stringify_keys(arc.demo5),
      "process_noise_sweep" => Enum.map(arc.process_noise_sweep, &process_noise_sweep_json/1),
      "cost" => stringify_keys(arc.cost)
    }
  end

  defp process_noise_sweep_json(cell) do
    %{
      "mode" => cell.mode,
      "process_noise_baseline_sigma_m" => cell.process_noise_baseline_sigma_m,
      "built_epoch_count" => cell.built_epoch_count,
      "solved_epoch_count" => cell.solved_epoch_count,
      "completed_arc" => cell.completed_arc,
      "segment_count" => cell.segment_count,
      "longest_stable_prefix_epochs" => cell.longest_stable_prefix_epochs,
      "first_unstable_time" => cell.first_unstable_time,
      "epoch2" => stringify_keys(cell.epoch2),
      "orbis" => stringify_keys(cell.orbis),
      "completed_arc_median_m" => cell.completed_arc_median_m,
      "completed_arc_p95_m" => cell.completed_arc_p95_m
    }
  end

  defp gauge_case_json(gauge_case) do
    %{
      "mode" => gauge_case.mode,
      "label" => gauge_case.input.label,
      "drive" => gauge_case.input.drive,
      "fixture" => gauge_case.input.fixture,
      "systems" => gauge_case.systems,
      "built_epoch_count" => gauge_case.built_epoch_count,
      "solved_epoch_count" => gauge_case.filter_case.solved_epoch_count,
      "completed_arc" => gauge_case.completed_arc,
      "longest_stable_prefix_epochs" => gauge_case.filter_case.longest_stable_prefix_epochs,
      "first_unstable_time" => gauge_case.filter_case.first_unstable_time,
      "epoch2" => stringify_keys(gauge_case.epoch2),
      "orbis" => stringify_keys(gauge_case.filter_case.orbis),
      "completed_arc_median_m" => gauge_case.completed_arc_median_m,
      "completed_arc_p95_m" => gauge_case.completed_arc_p95_m
    }
  end

  defp gauge_verdict(gps_cases, control) do
    gps_conditions =
      gps_cases
      |> Enum.map(&epoch2_condition/1)
      |> Enum.reject(&is_nil/1)

    control_condition = epoch2_condition(control)
    gps_max_condition = max_or_nil(gps_conditions)

    drop_ratio =
      if is_number(control_condition) and is_number(gps_max_condition) and
           gps_max_condition != 0.0 do
        control_condition / gps_max_condition
      end

    min_prefix =
      gps_cases |> Enum.map(& &1.filter_case.longest_stable_prefix_epochs) |> Enum.min()

    convicted =
      is_number(drop_ratio) and drop_ratio >= 1.0e9 and min_prefix > 2

    status = if(convicted, do: "convicted", else: "acquitted")

    rationale =
      if convicted do
        "GPS-only removes the multi-system reference-SD gauge rows, drops the epoch-2 condition estimate by about ten orders, and survives beyond epoch 2."
      else
        "GPS-only either retained the high epoch-2 condition estimate or failed to survive beyond epoch 2."
      end

    %{
      status: status,
      gps_min_stable_prefix_epochs: min_prefix,
      gps_epoch2_condition_min: min_or_nil(gps_conditions),
      gps_epoch2_condition_max: gps_max_condition,
      control_epoch2_condition: control_condition,
      condition_drop_ratio: drop_ratio,
      rationale: rationale
    }
  end

  defp detector_json(detector) do
    %{
      "raw_epoch_count" => detector.raw_epoch_count,
      "observed_satellite_count" => detector.observed_satellite_count,
      "dual_frequency" => stringify_keys(detector.dual_frequency),
      "doppler" => stringify_keys(detector.doppler),
      "thresholds" => stringify_keys(detector.thresholds),
      "stats" => stringify_keys(detector.stats),
      "agreement" => stringify_keys(detector.agreement),
      "segmentation" => stringify_keys(detector.segmentation),
      "events" => Enum.map(detector.events, &stringify_keys/1)
    }
  end

  defp filter_case_json(filter_case, demo5) do
    %{
      "mode" => filter_case.mode,
      "max_segment_epochs" => filter_case.max_segment_epochs,
      "segment_count" => filter_case.segment_count,
      "solved_epoch_count" => filter_case.solved_epoch_count,
      "longest_stable_prefix_epochs" => filter_case.longest_stable_prefix_epochs,
      "first_unstable_time" => filter_case.first_unstable_time,
      "orbis" => stringify_keys(filter_case.orbis),
      "comparative" => stringify_keys(comparative_verdict(filter_case.orbis, demo5, :per_arc)),
      "invariant" => stringify_keys(filter_case.invariant),
      "segments" => Enum.map(filter_case.segment_reports, &segment_json/1),
      "per_epoch" => Enum.map(filter_case.per_epoch, &epoch_json/1)
    }
  end

  defp comparison_json(comparison, demo5) do
    %{
      "mode" => comparison.mode,
      "max_segment_epochs" => comparison.max_segment_epochs,
      "segment_count" => comparison.segment_count,
      "orbis" => stringify_keys(comparison.orbis),
      "comparative" => stringify_keys(comparative_verdict(comparison.orbis, demo5, :per_arc)),
      "invariant" => stringify_keys(comparison.invariant),
      "segments" => Enum.map(comparison.segment_reports, &segment_json/1),
      "per_epoch" => Enum.map(comparison.per_epoch, &epoch_json/1)
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
      "state_diagnostics" => stringify_keys(epoch.state_diagnostics),
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
      "diagnostics" => stringify_keys(segment.diagnostics),
      "metadata" => stringify_keys(segment.metadata)
    }
  end

  defp option_notes_json do
    Enum.map(@filter_option_notes, fn {name, value, why} ->
      %{"option" => name, "value" => value, "why" => why}
    end)
  end

  defp detector_config_json do
    %{
      "systems" => @systems,
      "primary_codes" => stringify_keys(@primary_detector_codes),
      "dual_codes" => stringify_keys(@dual_detector_codes),
      "continuity_gap_s" => @continuity_gap_s,
      "code_phase_min_threshold_m" => @code_phase_min_threshold_m,
      "code_phase_sigma_multiplier" => @code_phase_sigma_multiplier,
      "doppler_min_threshold_m" => @doppler_min_threshold_m,
      "doppler_sigma_multiplier" => @doppler_sigma_multiplier,
      "geometry_free_threshold_m" => @gf_threshold_m,
      "melbourne_wubbena_threshold_cycles" => @mw_threshold_cycles,
      "detector_segment_max_sd_ambiguity_ids" => @detector_segment_max_sd_ambiguity_ids,
      "detector_max_injected_events_per_satellite" => @detector_max_injected_events_per_satellite,
      "injection_policy" =>
        "raw LLI or independent agreement, capped per satellite for the runnable key experiment"
    }
  end

  defp gauge_report_markdown(result) do
    verdict = result["verdict"]

    [
      "## Multi-system gauge micro-experiment",
      "",
      "Generated by `mix run test/fixtures/rtk/generators/rover_c1_slip_exploration_2026_06.exs --gauge-test`.",
      "Per-run JSON was emitted to `#{result["results_path"]}`.",
      "",
      "Setup: stock carried-state filter, sigma 30, no detector segmentation. GPS-only rows filter the epoch builder to `G` satellites; the control row uses the same code path with `G/R/E/C`.",
      "",
      gauge_result_table(result["gps_only"], result["control"]),
      "",
      gauge_verdict_text(verdict),
      ""
    ]
    |> Enum.join("\n")
  end

  defp gauge_result_table(gps_rows, control) do
    rows =
      (gps_rows ++ [control])
      |> Enum.map(fn row ->
        [
          row["mode"],
          row["label"],
          Enum.join(row["systems"], ""),
          row["built_epoch_count"],
          row["solved_epoch_count"],
          yes_no(row["completed_arc"]),
          row["longest_stable_prefix_epochs"],
          gauge_epoch2_condition_text(row),
          completed_metric(row, "completed_arc_median_m"),
          completed_metric(row, "completed_arc_p95_m")
        ]
        |> table_row()
      end)

    Enum.join(
      [
        "| Case | Arc | Systems | Built | Solved | Complete | Longest stable prefix | Epoch-2 cond est | Completed median m | Completed p95 m |",
        "|---|---|---:|---:|---:|---|---:|---:|---:|---:|"
      ] ++ rows,
      "\n"
    )
  end

  defp gauge_epoch2_condition_text(%{"epoch2" => nil}), do: "n/a"

  defp gauge_epoch2_condition_text(row) do
    row
    |> get_in(["epoch2", "information_condition_estimate"])
    |> fmt_sci()
  end

  defp gauge_verdict_text(%{"status" => "convicted"} = verdict) do
    "Verdict: convicted. GPS-only epoch-2 condition estimates span #{fmt_sci(verdict["gps_epoch2_condition_min"])}-#{fmt_sci(verdict["gps_epoch2_condition_max"])} versus the control at #{fmt_sci(verdict["control_epoch2_condition"])}, a #{fmt_ratio(verdict["condition_drop_ratio"])}x drop, and the GPS-only minimum stable prefix is #{epoch_count_text(verdict["gps_min_stable_prefix_epochs"])}. The multi-system reference-SD gauge weight is the structural driver of the 2e14 conditioning spike. Capability candidate: well-conditioned gauge, with the REFERENCE design changing the gauge weight to dominate the nullspace without crushing the ambiguity prior."
  end

  defp gauge_verdict_text(verdict) do
    "Verdict: acquitted. GPS-only epoch-2 condition estimates span #{fmt_sci(verdict["gps_epoch2_condition_min"])}-#{fmt_sci(verdict["gps_epoch2_condition_max"])} versus the control at #{fmt_sci(verdict["control_epoch2_condition"])}, a #{fmt_ratio(verdict["condition_drop_ratio"])}x drop, and the GPS-only minimum stable prefix is #{epoch_count_text(verdict["gps_min_stable_prefix_epochs"])}. That leaves intrinsic conditioning / C2 as the standing explanation."
  end

  defp yes_no(true), do: "yes"
  defp yes_no(false), do: "no"
  defp yes_no(_value), do: "n/a"

  defp epoch_count_text(1), do: "1 epoch"
  defp epoch_count_text(value), do: "#{value} epochs"

  defp report_markdown(result, results_path) do
    arcs = result["arcs"]
    pooled = result["pooled"]

    [
      "# C1 slip exploration, June 2026",
      "",
      "Generated by `mix run test/fixtures/rtk/generators/rover_c1_slip_exploration_2026_06.exs`.",
      "Per-epoch JSON was emitted to `#{results_path}`.",
      "",
      "This is a single-sided exploration artifact. The RTK solver and library code were not changed.",
      "",
      "## Inputs",
      "",
      "- Phone observations: each arc's `supplemental/gnss_rinex.21o` from `/tmp/gsdc-work`.",
      "- Base observations: NOAA CORS P222 RINEX 2.11 observation files named in each oracle's provenance.",
      "- Ephemeris: BKG combined broadcast NAV (`BRDC00WRD_R_..._MN.rnx`), matching each oracle's `pos1-sateph=brdc` source.",
      "- Detector-driven segmentation is injected only by setting rover observation `ambiguity_id` values in the harness epoch builder.",
      "- The runnable key experiment caps injected detector cuts at #{@detector_max_injected_events_per_satellite} per satellite and cuts harness segments when the single-difference ambiguity count would exceed #{@detector_segment_max_sd_ambiguity_ids}; the uncapped first-arc smoke attempted 11,847 injected keys and did not finish a 60-epoch segment.",
      "",
      "## Detector availability",
      "",
      detector_availability_table(arcs),
      "",
      "GPS and Galileo expose L5/E5a on these phone files where listed. Doppler is present on the same primary signals used by the L1/B1/G1 filter path.",
      "",
      "## Detector thresholds",
      "",
      detector_threshold_table(arcs),
      "",
      "The code-phase and Doppler gates are per-arc robust thresholds: median absolute step plus #{@code_phase_sigma_multiplier} scaled MAD, with minimums of #{fmt(@code_phase_min_threshold_m)} m and #{fmt(@doppler_min_threshold_m)} m. GF and MW use the batch carrier-phase gates: #{fmt(@gf_threshold_m)} m and #{fmt(@mw_threshold_cycles)} wide-lane cycles.",
      "",
      "## Detector stats",
      "",
      detector_stats_table(arcs),
      "",
      "## Detector agreement",
      "",
      detector_agreement_table(arcs),
      "",
      "## Key experiment",
      "",
      key_experiment_table(arcs, pooled),
      "",
      "The `stock` columns are the carried-state filter without detector segmentation. The `detector` columns are the same carried-state filter with pre-segmented rover ambiguity ids. The `one-epoch` columns are the existing single-epoch comparison row.",
      "",
      "## Costs",
      "",
      costs_table(arcs, pooled),
      "",
      "## Promotion verdict",
      "",
      promotion_verdict_text(pooled),
      "",
      "## Filter options",
      "",
      option_table(),
      "",
      "## Process-noise scaling sweep",
      "",
      "This section runs hypothesis 2 against the same four arcs and harness settings, changing only `process_noise_baseline_sigma_m` over #{Enum.map_join(@process_noise_sweep_sigmas, ", ", &fmt/1)} m. `stock` is the carried-state filter without detector segmentation; `detector` is the carried-state filter with the C1 pre-segmented rover ambiguity ids.",
      "",
      process_noise_sweep_table(arcs, pooled),
      "",
      "The condition estimate is the same row-sum information-matrix estimate used above, sampled at the second carried-state epoch. Median and p95 are reported only for arcs whose cell solved every built epoch; pooled values combine only those completed arcs.",
      "",
      process_noise_sweep_verdict_text(pooled),
      ""
    ]
    |> Enum.join("\n")
  end

  defp detector_availability_table(arcs) do
    rows =
      Enum.map(arcs, fn arc ->
        detector = arc["detector"]

        [
          arc["label"],
          detector["raw_epoch_count"],
          detector["observed_satellite_count"],
          availability_systems(detector["dual_frequency"]),
          availability_systems(detector["doppler"])
        ]
        |> table_row()
      end)

    Enum.join(
      [
        "| Arc | Raw epochs | Sats | Dual-freq systems | Doppler systems |",
        "|---|---:|---:|---|---|"
      ] ++ rows,
      "\n"
    )
  end

  defp availability_systems(availability) do
    availability
    |> Enum.filter(fn {_system, values} -> values["available_observations"] > 0 end)
    |> Enum.map_join(", ", fn {system, values} ->
      "#{system} #{fmt(100.0 * values["available_fraction"])}%"
    end)
    |> case do
      "" -> "none"
      text -> text
    end
  end

  defp detector_threshold_table(arcs) do
    rows =
      Enum.map(arcs, fn arc ->
        thresholds = arc["detector"]["thresholds"]
        code_cal = thresholds["calibration"]["code_phase"]
        dop_cal = thresholds["calibration"]["doppler_phase"]

        [
          arc["label"],
          fmt(thresholds["code_phase"]),
          fmt(code_cal["median"]),
          fmt(code_cal["mad_sigma"]),
          fmt(thresholds["doppler_phase"]),
          fmt(dop_cal["median"]),
          fmt(dop_cal["mad_sigma"]),
          fmt(thresholds["geometry_free_m"]),
          fmt(thresholds["melbourne_wubbena_cycles"])
        ]
        |> table_row()
      end)

    Enum.join(
      [
        "| Arc | Code-phase gate m | Code med m | Code MAD sigma m | Doppler gate m | Dop med m | Dop MAD sigma m | GF gate m | MW gate cyc |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|"
      ] ++ rows,
      "\n"
    )
  end

  defp detector_stats_table(arcs) do
    detectors = ["code_phase", "geometry_free", "melbourne_wubbena", "doppler_phase", "rinex_lli"]

    rows =
      Enum.flat_map(arcs, fn arc ->
        Enum.map(detectors, fn detector ->
          stats = arc["detector"]["stats"][detector]

          [
            arc["label"],
            detector,
            stats["events"],
            fmt(stats["exposure_satellite_minutes"]),
            fmt(stats["events_per_satellite_minute"]),
            fmt(stats["magnitude_m_median"]),
            fmt(stats["magnitude_m_p95"]),
            fmt(stats["magnitude_cycles_median"]),
            fmt(stats["magnitude_cycles_p95"])
          ]
          |> table_row()
        end)
      end)

    Enum.join(
      [
        "| Arc | Detector | Events | Sat-min | Events/sat-min | Median m | p95 m | Median cyc | p95 cyc |",
        "|---|---|---:|---:|---:|---:|---:|---:|---:|"
      ] ++ rows,
      "\n"
    )
  end

  defp detector_agreement_table(arcs) do
    rows =
      Enum.map(arcs, fn arc ->
        agreement = arc["detector"]["agreement"]
        event_keys = agreement["event_keys"]

        [
          arc["label"],
          event_keys,
          agreement["multi_detector_keys"],
          fmt(ratio(agreement["multi_detector_keys"] * 100.0, event_keys)),
          agreement["pairwise_overlap"]["code_phase+doppler_phase"] || 0,
          agreement["pairwise_overlap"]["geometry_free+melbourne_wubbena"] || 0,
          agreement["pairwise_overlap"]["doppler_phase+geometry_free"] || 0
        ]
        |> table_row()
      end)

    Enum.join(
      [
        "| Arc | Event keys | Multi-detector keys | Multi % | Code+Dop | GF+MW | Dop+GF |",
        "|---|---:|---:|---:|---:|---:|---:|"
      ] ++ rows,
      "\n"
    )
  end

  defp key_experiment_table(arcs, pooled) do
    rows =
      Enum.map(arcs, fn arc ->
        stock = arc["stock"]
        segmented = arc["segmented"]
        comparison = arc["comparison"]

        [
          arc["label"],
          stock["longest_stable_prefix_epochs"],
          segmented["longest_stable_prefix_epochs"],
          fmt(stock["orbis"]["error_3d_median_m"]),
          fmt(stock["orbis"]["error_3d_p95_m"]),
          fmt(segmented["orbis"]["error_3d_median_m"]),
          fmt(segmented["orbis"]["error_3d_p95_m"]),
          fmt(comparison["orbis"]["error_3d_median_m"]),
          fmt(comparison["orbis"]["error_3d_p95_m"]),
          "#{segmented["solved_epoch_count"]}/#{arc["built_epoch_count"]}"
        ]
        |> table_row()
      end)

    pooled_row =
      [
        "pooled",
        "",
        pooled["promotion"]["segmented_min_stable_prefix_epochs"],
        fmt(pooled["stock"]["orbis"]["error_3d_median_m"]),
        fmt(pooled["stock"]["orbis"]["error_3d_p95_m"]),
        fmt(pooled["segmented"]["orbis"]["error_3d_median_m"]),
        fmt(pooled["segmented"]["orbis"]["error_3d_p95_m"]),
        fmt(pooled["comparison"]["orbis"]["error_3d_median_m"]),
        fmt(pooled["comparison"]["orbis"]["error_3d_p95_m"]),
        ""
      ]
      |> table_row()

    Enum.join(
      [
        "| Arc | Stock prefix | Detector prefix | Stock med m | Stock p95 m | Detector med m | Detector p95 m | One-epoch med m | One-epoch p95 m | Solved/built |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|"
      ] ++ rows ++ [pooled_row],
      "\n"
    )
  end

  defp costs_table(arcs, pooled) do
    rows =
      Enum.map(arcs, fn arc ->
        cost = arc["cost"]

        [
          arc["label"],
          cost["detector_event_keys"],
          cost["detector_eligible_event_keys"],
          cost["detector_injected_event_keys"],
          cost["split_ambiguity_ids"],
          fmt(cost["stock_sd_columns"]["max_single_difference_columns"]),
          fmt(cost["detector_sd_columns"]["max_single_difference_columns"]),
          fmt(cost["extra_max_sd_columns"]),
          cost["detector_epochs_lost"],
          cost["new_failure_mode"]
        ]
        |> table_row()
      end)

    pooled_row =
      [
        "pooled",
        pooled["cost"]["detector_event_keys"],
        pooled["cost"]["detector_eligible_event_keys"],
        pooled["cost"]["detector_injected_event_keys"],
        pooled["cost"]["split_ambiguity_ids"],
        fmt(pooled["cost"]["max_stock_sd_columns"]),
        fmt(pooled["cost"]["max_detector_sd_columns"]),
        fmt(
          nil_safe_sub(
            pooled["cost"]["max_detector_sd_columns"],
            pooled["cost"]["max_stock_sd_columns"]
          )
        ),
        pooled["cost"]["detector_epochs_lost"],
        "see per-arc rows"
      ]
      |> table_row()

    Enum.join(
      [
        "| Arc | Event keys | Eligible keys | Injected keys | Split IDs | Stock max SD cols | Detector max SD cols | Extra max cols | Detector epochs lost | New failure mode |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|---|"
      ] ++ rows ++ [pooled_row],
      "\n"
    )
  end

  defp promotion_verdict_text(pooled) do
    promotion = pooled["promotion"]
    verdict = promotion["verdict"]
    segmented_median = promotion["segmented_pooled_median_m"]
    one_epoch_median = promotion["one_epoch_pooled_median_m"]
    reference_median = promotion["reference_single_epoch_median_m"]

    if verdict == "promote_c1" do
      "Promote C1. Detector segmentation extended the minimum stable prefix from #{promotion["stock_min_stable_prefix_epochs"]} to #{promotion["segmented_min_stable_prefix_epochs"]} epochs and reached pooled median #{fmt(segmented_median)} m, beating the #{fmt(reference_median)} m reference and the measured one-epoch median #{fmt(one_epoch_median)} m."
    else
      "Do not promote C1 yet. Detector segmentation extended the minimum stable prefix from #{promotion["stock_min_stable_prefix_epochs"]} to #{promotion["segmented_min_stable_prefix_epochs"]} epochs, but pooled median was #{fmt(segmented_median)} m versus the #{fmt(reference_median)} m reference and #{fmt(one_epoch_median)} m for the measured one-epoch row."
    end
  end

  defp process_noise_sweep_table(arcs, pooled) do
    arc_labels = Enum.map(arcs, &short_arc_label/1)

    rows =
      for mode <- ["stock", "detector"], sigma_m <- @process_noise_sweep_sigmas do
        cells = Enum.map(arcs, &find_process_noise_cell(&1, mode, sigma_m))
        pooled_cell = find_pooled_process_noise_cell(pooled, mode, sigma_m)

        [
          mode,
          fmt(sigma_m),
          "#{pooled_cell["completed_arcs"]}/#{pooled_cell["total_arcs"]}",
          cells |> Enum.map_join("/", & &1["longest_stable_prefix_epochs"]),
          cells |> Enum.map_join("/", &completed_metric(&1, "completed_arc_median_m")),
          cells |> Enum.map_join("/", &completed_metric(&1, "completed_arc_p95_m")),
          pooled_completed_metric(pooled_cell),
          cells |> Enum.map_join("/", &epoch2_condition_text/1)
        ]
        |> table_row()
      end

    label_text = Enum.join(arc_labels, "/")

    Enum.join(
      [
        "| Mode | Sigma m | Complete arcs | Prefix #{label_text} | Completed median m #{label_text} | Completed p95 m #{label_text} | Pooled completed med/p95 m | Epoch-2 cond est #{label_text} |",
        "|---|---:|---:|---:|---:|---:|---:|---:|"
      ] ++ rows,
      "\n"
    )
  end

  defp short_arc_label(arc) do
    label = arc["label"]

    cond do
      String.contains?(label, "08_04") -> "SJC1-08-04"
      String.contains?(label, "svl1") -> "SVL1-08-24"
      String.contains?(label, "12_15") -> "MTV1-12-15"
      String.contains?(label, "12_28") -> "MTV1-12-28"
      true -> label
    end
  end

  defp find_process_noise_cell(arc, mode, sigma_m) do
    Enum.find(arc["process_noise_sweep"], fn cell ->
      cell["mode"] == mode and same_sigma?(cell["process_noise_baseline_sigma_m"], sigma_m)
    end)
  end

  defp find_pooled_process_noise_cell(pooled, mode, sigma_m) do
    Enum.find(pooled["process_noise_sweep"], fn cell ->
      cell["mode"] == mode and same_sigma?(cell["process_noise_baseline_sigma_m"], sigma_m)
    end)
  end

  defp same_sigma?(a, b), do: abs(a - b) < 1.0e-9

  defp completed_metric(%{"completed_arc" => true} = cell, key), do: fmt(cell[key])
  defp completed_metric(_cell, _key), do: "n/a"

  defp pooled_completed_metric(%{"completed_arcs" => 0}), do: "n/a"

  defp pooled_completed_metric(cell) do
    "#{fmt(cell["completed_pooled_median_m"])}/#{fmt(cell["completed_pooled_p95_m"])}"
  end

  defp epoch2_condition_text(%{"epoch2" => nil}), do: "n/a"

  defp epoch2_condition_text(cell) do
    cell
    |> get_in(["epoch2", "information_condition_estimate"])
    |> fmt_sci()
  end

  defp process_noise_sweep_verdict_text(pooled) do
    sweep = pooled["process_noise_sweep"]
    stock_rows = Enum.filter(sweep, &(&1["mode"] == "stock"))
    detector_rows = Enum.filter(sweep, &(&1["mode"] == "detector"))
    stock_survivors = Enum.filter(stock_rows, &(&1["min_stable_prefix_epochs"] > 1))
    detector_survivors = Enum.filter(detector_rows, &(&1["min_stable_prefix_epochs"] > 1))
    best_stock = Enum.max_by(stock_rows, & &1["min_stable_prefix_epochs"])

    if stock_survivors == [] do
      detector_sentence =
        case detector_survivors do
          [] ->
            "Detector segmentation also fails to make every arc survive epoch 2 in this sweep."

          rows ->
            sigmas = rows |> Enum.map_join(", ", &fmt(&1["process_noise_baseline_sigma_m"]))
            best_detector = best_completed_process_noise_row(rows)

            gap =
              ratio(
                best_detector["completed_pooled_median_m"],
                pooled["demo5"]["error_3d_median_m"]
              )

            "Detector segmentation survives epoch 2 for sigma #{sigmas} m, but its best completed pooled median is #{fmt(best_detector["completed_pooled_median_m"])} m versus demo5 #{fmt(pooled["demo5"]["error_3d_median_m"])} m (#{fmt_ratio(gap)}x), so survival is not accuracy recovery."
        end

      "Verdict: no. Process-noise scaling alone does not fix epoch-2 survival: the best stock row has a four-arc minimum stable prefix of #{best_stock["min_stable_prefix_epochs"]} at sigma #{fmt(best_stock["process_noise_baseline_sigma_m"])} m. Stock epoch-2 condition estimates span #{condition_range_text(stock_rows)}, which supports C2: intrinsic conditioning / missing innovation screening rather than only epoch-rate process noise. #{detector_sentence}"
    else
      best = best_completed_process_noise_row(stock_survivors)
      median_gap = ratio(best["completed_pooled_median_m"], pooled["demo5"]["error_3d_median_m"])
      p95_gap = ratio(best["completed_pooled_p95_m"], pooled["demo5"]["error_3d_p95_m"])

      "Verdict: yes. Process-noise scaling fixes epoch-2 survival without C1 segmentation at sigma #{fmt(best["process_noise_baseline_sigma_m"])} m, with completed pooled median/p95 #{fmt(best["completed_pooled_median_m"])}/#{fmt(best["completed_pooled_p95_m"])} m. Demo5 pooled median/p95 is #{fmt(pooled["demo5"]["error_3d_median_m"])}/#{fmt(pooled["demo5"]["error_3d_p95_m"])} m, leaving #{fmt_ratio(median_gap)}x median and #{fmt_ratio(p95_gap)}x p95 residual gaps. The capability candidate is epoch-rate-aware process noise guidance, likely a per-sqrt-second scaling option; remaining accuracy work still needs the next ledger item."
    end
  end

  defp best_completed_process_noise_row(rows) do
    rows
    |> Enum.filter(&(&1["completed_arcs"] > 0))
    |> case do
      [] -> Enum.min_by(rows, &(&1["completed_pooled_median_m"] || 1.0e300))
      completed -> Enum.min_by(completed, &(&1["completed_pooled_median_m"] || 1.0e300))
    end
  end

  defp condition_range_text(rows) do
    mins = rows |> Enum.map(& &1["epoch2_condition_min"]) |> Enum.reject(&is_nil/1)
    maxes = rows |> Enum.map(& &1["epoch2_condition_max"]) |> Enum.reject(&is_nil/1)

    case {mins, maxes} do
      {[], []} -> "n/a"
      _ -> "#{fmt_sci(Enum.min(mins))}-#{fmt_sci(Enum.max(maxes))}"
    end
  end

  defp fmt_ratio(nil), do: "n/a"
  defp fmt_ratio(value), do: :erlang.float_to_binary(value, decimals: 1)

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

  defp diagnosis_table(arcs) do
    rows =
      Enum.map(arcs, fn arc ->
        diagnosis = arc["diagnosis"]
        first_bad = diagnosis["first_bad_epoch"]
        previous = diagnosis["previous_epoch"] || %{}
        changes = diagnosis["changes_at_first_bad"] || %{}

        if first_bad do
          [
            arc["label"],
            diagnosis["verdict"],
            first_bad["time"],
            first_bad["segment_epoch_index"],
            fmt(previous["error_3d_m"]),
            fmt(first_bad["error_3d_m"]),
            fmt(first_bad["carried_baseline_error_3d_m"]),
            fmt(first_bad["information_condition_estimate"]),
            first_bad["sd_ambiguity_columns"],
            "#{changes["previous_hold_count"] || 0}->#{changes["hold_count"] || 0}",
            "+#{length(changes["satellites_added"] || [])}/-#{length(changes["satellites_removed"] || [])}",
            "#{fmt(first_bad["max_abs_code_residual_m"])}/#{fmt(first_bad["max_abs_phase_residual_m"])}",
            diagnosis["mechanism"]
          ]
          |> table_row()
        else
          [
            arc["label"],
            diagnosis["verdict"],
            "none",
            "",
            "",
            "",
            "",
            "",
            "",
            "",
            "",
            "",
            diagnosis["mechanism"]
          ]
          |> table_row()
        end
      end)

    Enum.join(
      [
        "| Arc | Verdict | First bad GPST | Seg idx | Prev 3D m | Bad 3D m | Carried 3D m | Cond est | SD cols | Holds | Sat +/- | Max code/phase residual m | Mechanism |",
        "|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|"
      ] ++ rows,
      "\n"
    )
  end

  defp diagnosis_summary(arcs) do
    diagnoses = Enum.map(arcs, & &1["diagnosis"])
    first_bad = Enum.map(diagnoses, & &1["first_bad_epoch"]) |> Enum.reject(&is_nil/1)
    verdicts = diagnoses |> Enum.map(& &1["verdict"]) |> Enum.uniq() |> Enum.sort()

    if verdicts == ["filter_behavior"] and length(first_bad) == length(arcs) do
      seg_indices = Enum.map(first_bad, & &1["segment_epoch_index"])
      hold_counts = Enum.map(first_bad, &(&1["hold_count"] || 0))

      if Enum.all?(seg_indices, &(&1 == 1)) and Enum.all?(hold_counts, &(&1 == 0)) do
        "A-vs-b verdict: (b) real filter behavior under the measured phone configuration. All four arcs cross the #{fmt(@divergence_threshold_m)} m threshold on the second carried-state epoch, before any integer hold is accepted; selected references are present, segmented `~ra` ambiguity ids are absent, and the residual/time sanity gates pass. Three first-bad epochs add one GPS satellite, but SVL fails with the same satellite set, so constellation churn is not required. The sequential filter completes the arcs, but the longest stable prefix is one epoch."
      else
        "A-vs-b verdict: (b) real filter behavior under the measured phone configuration. The first-bad rows occur after the sanity gates with selected references present; see the table for the triggering state changes."
      end
    else
      "A-vs-b verdict: at least one arc has a harness-bug marker; inspect the first-bad rows before treating the distributions as filter behavior."
    end
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

    comparison = pooled["comparison"]["orbis"]

    comparison_row =
      [
        "pooled per-epoch comparison only",
        "#{comparison["epochs"]}/#{pooled["demo5"]["epochs"]}",
        fmt(comparison["error_3d_median_m"]),
        fmt(pooled["demo5"]["error_3d_median_m"]),
        fmt(comparison["error_3d_p95_m"]),
        fmt(pooled["demo5"]["error_3d_p95_m"]),
        "not operating mode"
      ]
      |> table_row()

    Enum.join(
      [
        "| Arc | Epochs Orbis/demo5 | Orbis 3D median m | demo5 3D median m | Orbis 3D p95 m | demo5 3D p95 m | Bar |",
        "|---|---:|---:|---:|---:|---:|---:|"
      ] ++ rows ++ [pooled_row, comparison_row],
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

  defp min_or_nil([]), do: nil
  defp min_or_nil(values), do: Enum.min(values)

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

  defp fmt_sci(nil), do: ""

  defp fmt_sci(value) when is_float(value) do
    :io_lib.format("~.6e", [value]) |> IO.iodata_to_binary()
  end

  defp fmt_sci(value), do: to_string(value)
end

RoverC2InnovationScreening202606.main(System.argv())
