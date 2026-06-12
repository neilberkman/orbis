defmodule Orbis.GNSS.RTKRealArcTest do
  @moduledoc false

  use ExUnit.Case, async: false

  alias Orbis.GNSS.RINEX.Observations
  alias Orbis.GNSS.RTK
  alias Orbis.GNSS.SP3

  @sp3_path Path.join(__DIR__, "fixtures/sp3/GBM0MGXRAP_20201770000_01D_05M_ORB_120epoch.sp3")
  @cod_sp3_path Path.join(__DIR__, "fixtures/sp3/COD0MGXFIN_20201770000_01D_05M_ORB.SP3")
  @wtzr_obs_path Path.join(
                   __DIR__,
                   "fixtures/obs/WTZR00DEU_R_20201770000_01D_30S_MO_120epoch.rnx"
                 )
  @wtzz_obs_path Path.join(
                   __DIR__,
                   "fixtures/obs/WTZZ00DEU_R_20201770000_01D_30S_MO_120epoch.rnx"
                 )
  @rtklib_oracle_path Path.join(__DIR__, "fixtures/rtk/wtzr_wtzz_rtklib_oracle.json")
  @rtklib_precise_oracle_path Path.join(
                                __DIR__,
                                "fixtures/rtk/wtzr_wtzz_rtklib_precise_oracle.json"
                              )
  @rtklib_kinematic_oracle_path Path.join(
                                  __DIR__,
                                  "fixtures/rtk/wtzr_wtzz_kinematic_gps_rtklib_oracle.json"
                                )
  @c_m_s 299_792_458.0
  @gps_l1_hz 1_575_420_000.0
  @gps_l2_hz 1_227_600_000.0
  @gps_l1_wavelength_m @c_m_s / @gps_l1_hz

  @wtzr_marker {4_075_580.3111, 931_854.0543, 4_801_568.2808}
  @wtzz_marker {4_075_579.1913, 931_853.3696, 4_801_569.1897}

  @tag timeout: 180_000
  test "real co-located Wettzell L1 RTK solves the antenna baseline and fixes a safe partial subset" do
    sp3 = SP3.load!(@sp3_path)
    base_obs = Observations.load!(@wtzr_obs_path)
    rover_obs = Observations.load!(@wtzz_obs_path)

    base_arp = arp_position(@wtzr_marker, antenna_height_m(base_obs))
    rover_arp = arp_position(@wtzz_marker, antenna_height_m(rover_obs))
    marker_baseline = sub3(@wtzz_marker, @wtzr_marker)
    antenna_baseline = sub3(rover_arp, base_arp)
    epochs = real_gps_l1_rtk_epochs(sp3, base_obs, rover_obs, 120)

    assert length(epochs) == 120

    opts = [
      initial_baseline_m: {0.0, 0.0, 0.0},
      max_iterations: 10,
      on_cycle_slip: :split_arc,
      elevation_weighting: true,
      code_sigma_m: 2.0,
      phase_sigma_m: 0.01,
      ambiguity_wavelength_m: @gps_l1_wavelength_m,
      integer_candidate_limit: 200_000
    ]

    float_opts = Keyword.drop(opts, [:ambiguity_wavelength_m, :integer_candidate_limit])

    assert {:ok, float} = RTK.solve_float_baseline_epochs(base_arp, epochs, float_opts)

    assert {:ok, smoothed_float} =
             RTK.solve_float_baseline_epochs(
               base_arp,
               epochs,
               Keyword.put(float_opts, :code_smoothing, true)
             )

    float_antenna_error_m = position_error(float.baseline_m, antenna_baseline)
    float_marker_error_m = position_error(float.baseline_m, marker_baseline)

    assert float.metadata.n_epochs == 120
    assert float.metadata.measurement_covariance.elevation_weighting
    assert length(float.metadata.split_cycle_slip_arcs) == 4
    assert float.metadata.dropped_sats == []
    assert Enum.any?(float.used_sats, &String.starts_with?(&1, "G09"))
    assert Enum.any?(float.used_sats, &String.starts_with?(&1, "G21"))
    assert Enum.any?(float.used_sats, &String.starts_with?(&1, "G27"))
    # Epoch-local satellite use admits short post-slip fragments that the
    # previous common-satellite batch solve dropped. The baseline remains
    # centimetre-scale, while the residual gate records that this is still a
    # batch float solution rather than a production RTK filter.
    assert float.metadata.phase_rms_m < 0.03
    assert smoothed_float.metadata.code_smoothing
    assert smoothed_float.metadata.code_rms_m < float.metadata.code_rms_m * 0.5
    assert position_error(smoothed_float.baseline_m, antenna_baseline) < 0.08

    # The SSC coordinates are marker coordinates; the observations are tied to
    # the antenna reference points. Applying the RINEX antenna-height deltas is
    # the difference between a decimetre-looking residual and a real
    # centimetre-scale short-baseline float RTK result.
    assert float_marker_error_m > 0.15
    assert float_antenna_error_m < 0.08

    assert {:ok, fixed} = RTK.solve_fixed_baseline_epochs(base_arp, epochs, opts)

    assert fixed.metadata.integer_status == :not_fixed
    assert fixed.metadata.integer_ratio < 3.0
    assert fixed.metadata.integer_candidates > 0

    fixed_antenna_error_m = position_error(fixed.baseline_m, antenna_baseline)

    # The full integer set still fails the ratio test, but with receiver-specific
    # transmit-time satellite positions the corresponding float baseline is
    # already millimetre-scale on this co-located real arc.
    assert fixed_antenna_error_m < float_antenna_error_m
    assert fixed_antenna_error_m < 0.01

    assert {:ok, partial_fixed} =
             RTK.solve_fixed_baseline_epochs(
               base_arp,
               epochs,
               Keyword.merge(opts,
                 partial_ambiguity_resolution: true,
                 partial_min_ambiguities: 4
               )
             )

    partial_antenna_error_m = position_error(partial_fixed.baseline_m, antenna_baseline)

    assert partial_fixed.metadata.partial_ambiguity_resolution
    assert partial_fixed.metadata.integer_status == :fixed
    assert partial_fixed.metadata.integer_ratio > 3.0
    assert partial_fixed.metadata.partial_fixed

    assert partial_fixed.metadata.partial_fixed_ambiguities == [
             "G05",
             "G07",
             "G08",
             "G09",
             "G13",
             "G15@rover#2|ref=G30",
             "G18",
             "G27@rover#1|ref=G30",
             "G28"
           ]

    assert partial_fixed.metadata.partial_free_ambiguities != []
    assert partial_antenna_error_m < float_antenna_error_m
    assert partial_antenna_error_m < 0.06

    dual_epochs = real_gps_l1_l2_rtk_epochs(sp3, base_obs, rover_obs, 120)
    assert length(dual_epochs) == 120

    assert {:ok, wide_lane_fixed} =
             RTK.solve_widelane_fixed_baseline_epochs(base_arp, dual_epochs,
               initial_baseline_m: {0.0, 0.0, 0.0},
               max_iterations: 10,
               on_cycle_slip: :drop_satellite,
               elevation_weighting: true,
               code_sigma_m: 2.0,
               phase_sigma_m: 0.01,
               integer_candidate_limit: 200_000
             )

    assert wide_lane_fixed.wide_lane_ambiguities_cycles != nil
    assert wide_lane_fixed.metadata.integer_method == :widelane_narrowlane_lambda
    assert wide_lane_fixed.metadata.integer_status == :fixed
    assert wide_lane_fixed.metadata.integer_ratio > 3.0

    wide_lane_antenna_error_m = position_error(wide_lane_fixed.baseline_m, antenna_baseline)

    # With receiver-specific transmit-time geometry, the dual-frequency
    # wide-lane/narrow-lane path reaches the real centimetre-grade milestone as
    # a full fix rather than needing the earlier partial-AR escape hatch.
    assert wide_lane_antenna_error_m < 0.01
  end

  @tag timeout: 180_000
  test "RTKLIB two-epoch prefix fixes with receiver-specific transmit-time geometry" do
    oracle =
      @rtklib_oracle_path
      |> File.read!()
      |> Jason.decode!()

    precise_oracle =
      @rtklib_precise_oracle_path
      |> File.read!()
      |> Jason.decode!()

    sp3 = SP3.load!(@cod_sp3_path)
    base_obs = Observations.load!(@wtzr_obs_path)
    rover_obs = Observations.load!(@wtzz_obs_path)
    base_arp = arp_position(@wtzr_marker, antenna_height_m(base_obs))
    rover_arp = arp_position(@wtzz_marker, antenna_height_m(rover_obs))
    antenna_baseline = sub3(rover_arp, base_arp)
    antenna_baseline_enu = enu_map_to_tuple(oracle["truth"]["antenna_baseline_enu_m"])

    rtklib_epoch_2 = oracle["per_epoch"] |> Enum.at(1)
    rtklib_precise_epoch_2 = precise_oracle["per_epoch"] |> Enum.at(1)
    rtklib_epoch_2_baseline = enu_map_to_tuple(rtklib_epoch_2["baseline_enu_m"])

    assert oracle["reference"]["first_fixed_index"] == 1
    assert rtklib_epoch_2["fix_status"] == "fixed"
    assert rtklib_epoch_2["ratio"] >= 3.0
    assert position_error(rtklib_epoch_2_baseline, antenna_baseline_enu) < 0.006
    assert precise_oracle["reference"]["first_fixed_index"] == 1
    assert precise_oracle["broadcast_comparison"]["same_fix_status_by_epoch"]
    assert precise_oracle["broadcast_comparison"]["max_baseline_delta_m"] < 0.002
    assert rtklib_precise_epoch_2["fix_status"] == "fixed"
    assert rtklib_precise_epoch_2["ratio"] >= 3.0

    epochs =
      sp3
      |> real_gps_l1_rtk_epochs(base_obs, rover_obs, 2)

    assert length(epochs) == 2

    assert {:ok, sol} =
             RTK.solve_fixed_baseline_epochs(base_arp, epochs,
               initial_baseline_m: {0.0, 0.0, 0.0},
               max_iterations: 10,
               on_cycle_slip: :split_arc,
               elevation_weighting: true,
               elevation_mask_deg: 10.0,
               code_sigma_m: 2.0,
               phase_sigma_m: 0.01,
               ambiguity_wavelength_m: @gps_l1_wavelength_m,
               integer_candidate_limit: 200_000,
               partial_ambiguity_resolution: true,
               partial_min_ambiguities: 4
             )

    assert sol.metadata.integer_status == :fixed
    refute sol.metadata.partial_ambiguity_resolution
    refute sol.metadata.partial_fixed
    assert sol.metadata.integer_ratio >= 3.0
    assert sol.metadata.integer_method == :lambda
    assert sol.metadata.n_epochs == 2
    assert sol.metadata.elevation_mask_deg == 10.0
    assert sol.metadata.elevation_masked_sats == ["G08", "G18", "G27"]

    # The gap pinned by the earlier oracle was not LAMBDA or covariance shape:
    # the real-data builder was feeding one receive-time satellite position to
    # both receivers. With receiver-specific transmit-time positions, the same
    # two-epoch prefix now reaches the RTKLIB millimetre class instead of false
    # confidence metres away from the ARP truth.
    error_m = position_error(sol.baseline_m, antenna_baseline)
    assert error_m < 0.01

    assert {:ok, filter} =
             RTK.solve_filter_baseline_epochs(base_arp, epochs,
               initial_baseline_m: {0.0, 0.0, 0.0},
               max_iterations: 10,
               on_cycle_slip: :split_arc,
               elevation_weighting: true,
               elevation_mask_deg: 10.0,
               code_sigma_m: 2.0,
               phase_sigma_m: 0.01,
               ambiguity_wavelength_m: @gps_l1_wavelength_m,
               integer_candidate_limit: 200_000
             )

    assert filter.metadata.first_fixed_index in [0, 1]
    assert Enum.any?(filter.epochs, &(&1.integer_status == :fixed))
    assert position_error(filter.baseline_m, antenna_baseline) < 0.01

    assert {:ok, rtklib_weighted_filter} =
             RTK.solve_filter_baseline_epochs(base_arp, epochs,
               initial_baseline_m: {0.0, 0.0, 0.0},
               max_iterations: 10,
               on_cycle_slip: :split_arc,
               elevation_mask_deg: 10.0,
               stochastic_model: :rtklib,
               code_sigma_m: 0.3,
               phase_sigma_m: 0.003,
               ambiguity_wavelength_m: @gps_l1_wavelength_m,
               integer_candidate_limit: 200_000
             )

    # The RTKLIB stochastic model is still compatible with the corrected
    # geometry; it should also fix the prefix instead of reporting the old
    # covariance-matches-state-does-not gap.
    assert rtklib_weighted_filter.metadata.measurement_covariance.stochastic_model == :rtklib
    assert rtklib_weighted_filter.metadata.first_fixed_index in [0, 1]
    assert Enum.any?(rtklib_weighted_filter.epochs, &(&1.integer_status == :fixed))
    assert position_error(rtklib_weighted_filter.baseline_m, antenna_baseline) < 0.01

    full_epochs =
      sp3
      |> real_gps_l1_rtk_epochs(base_obs, rover_obs, precise_oracle["reference"]["epochs"])

    assert length(full_epochs) == precise_oracle["reference"]["epochs"]

    assert {:ok, full_filter} =
             RTK.solve_filter_baseline_epochs(base_arp, full_epochs,
               initial_baseline_m: {0.0, 0.0, 0.0},
               max_iterations: 10,
               on_cycle_slip: :split_arc,
               elevation_mask_deg: 10.0,
               stochastic_model: :rtklib,
               code_sigma_m: 0.3,
               phase_sigma_m: 0.003,
               ambiguity_wavelength_m: @gps_l1_wavelength_m,
               integer_candidate_limit: 200_000
             )

    full_fixed_epochs = Enum.count(full_filter.epochs, &(&1.integer_status == :fixed))

    # Full-arc parity gate against the precise RTKLIB fixture. The exact
    # first-fix epoch can be one epoch earlier because Orbis seeds the
    # single-difference ambiguities from phase-code differences before the first
    # update; the meaningful invariants are near-all fixed epochs and
    # millimetre-class final ARP-baseline agreement.
    assert full_filter.metadata.first_fixed_index <=
             precise_oracle["reference"]["first_fixed_index"] + 1

    assert full_fixed_epochs >= precise_oracle["reference"]["fixed_epochs"] - 1
    assert position_error(full_filter.baseline_m, antenna_baseline) < 0.01

    assert {:ok, rust_full_filter} =
             RTK.solve_filter_baseline_epochs(base_arp, full_epochs,
               initial_baseline_m: {0.0, 0.0, 0.0},
               max_iterations: 10,
               on_cycle_slip: :split_arc,
               elevation_mask_deg: 10.0,
               stochastic_model: :rtklib,
               code_sigma_m: 0.3,
               phase_sigma_m: 0.003,
               ambiguity_wavelength_m: @gps_l1_wavelength_m,
               integer_candidate_limit: 200_000,
               filter_kernel: :rust
             )

    rust_fixed_epochs = Enum.count(rust_full_filter.epochs, &(&1.integer_status == :fixed))

    assert rust_full_filter.metadata.filter_kernel == :rust
    assert length(rust_full_filter.epochs) == length(full_filter.epochs)

    assert rust_full_filter.metadata.first_fixed_index <=
             precise_oracle["reference"]["first_fixed_index"] + 1

    assert rust_fixed_epochs >= precise_oracle["reference"]["fixed_epochs"] - 1
    assert position_error(rust_full_filter.baseline_m, antenna_baseline) < 0.01
    assert_exact_position(rust_full_filter.baseline_m, full_filter.baseline_m)

    for {rust_epoch, elixir_epoch} <- Enum.zip(rust_full_filter.epochs, full_filter.epochs) do
      assert rust_epoch.index == elixir_epoch.index
      assert rust_epoch.integer_status == elixir_epoch.integer_status
      assert rust_epoch.newly_fixed_ambiguities == elixir_epoch.newly_fixed_ambiguities
      assert rust_epoch.fixed_ambiguities == elixir_epoch.fixed_ambiguities
      assert_exact_position(rust_epoch.baseline_m, elixir_epoch.baseline_m)

      if is_number(rust_epoch.integer_ratio) and is_number(elixir_epoch.integer_ratio) do
        assert rust_epoch.integer_ratio === elixir_epoch.integer_ratio
      end
    end
  end

  @tag timeout: 180_000
  test "kinematic-mode RTK filter reproduces the RTKLIB kinematic oracle on the static arc" do
    oracle = @rtklib_kinematic_oracle_path |> File.read!() |> Jason.decode!()

    sp3 = SP3.load!(@cod_sp3_path)
    base_obs = Observations.load!(@wtzr_obs_path)
    rover_obs = Observations.load!(@wtzz_obs_path)
    base_arp = arp_position(@wtzr_marker, antenna_height_m(base_obs))
    rover_arp = arp_position(@wtzz_marker, antenna_height_m(rover_obs))
    antenna_baseline = sub3(rover_arp, base_arp)

    epochs = real_gps_l1_rtk_epochs(sp3, base_obs, rover_obs, oracle["reference"]["epochs"])
    assert length(epochs) == oracle["reference"]["epochs"]

    # The oracle was generated at RTKLIB elmask 15; the Elixir harness runs at
    # mask 10 (its proven full-arc setting). Different masks select different
    # satellite sets, so this gate compares baseline-to-truth and fix fraction —
    # NOT exact per-epoch satellite counts — per the oracle-set-difference note.
    base_opts = [
      initial_baseline_m: {0.0, 0.0, 0.0},
      max_iterations: 10,
      on_cycle_slip: :split_arc,
      elevation_mask_deg: 10.0,
      stochastic_model: :rtklib,
      code_sigma_m: 0.3,
      phase_sigma_m: 0.003,
      ambiguity_wavelength_m: @gps_l1_wavelength_m,
      integer_candidate_limit: 200_000
    ]

    # Static accumulation (process noise off) is the existing default; kinematic
    # adds the between-epoch Q-inflation time update on the baseline block.
    assert {:ok, static} = RTK.solve_filter_baseline_epochs(base_arp, epochs, base_opts)

    assert {:ok, kinematic} =
             RTK.solve_filter_baseline_epochs(
               base_arp,
               epochs,
               [process_noise_baseline_sigma_m: 30.0] ++ base_opts
             )

    # The time update is genuinely exercised: loosening the baseline prior each
    # epoch moves the trajectory off the static accumulation.
    assert position_error(kinematic.baseline_m, static.baseline_m) > 0.0

    # RTKLIB kinematic fixes 119/120 on this arc; the baseline loosening can cost
    # a few fixes, so gate on the oracle's fixed count with margin.
    kin_fixed = Enum.count(kinematic.epochs, &(&1.integer_status == :fixed))

    assert kin_fixed >= oracle["reference"]["fixed_epochs"] - 5

    # Final baseline stays in RTKLIB's converged class (~mm).
    assert position_error(kinematic.baseline_m, antenna_baseline) < 0.01

    # Every fixed epoch — including epoch 0 — reports a cm-class baseline. The
    # reported baseline is the ambiguity-conditioned (fixed) solution, so the
    # first fixed epoch no longer reports its float baseline while claiming
    # `:fixed` (the cold-start false confidence that this gate first surfaced).
    fixed_errors =
      kinematic.epochs
      |> Enum.filter(&(&1.integer_status == :fixed))
      |> Enum.map(&position_error(&1.baseline_m, antenna_baseline))

    assert Enum.max(fixed_errors) < 0.02
    assert Enum.sum(fixed_errors) / length(fixed_errors) < 0.01

    # Trace gate: the Rust kernel reproduces the kinematic Elixir reference
    # bit-for-bit, including the process-noise time update (the Lambda^-1 + Q
    # covariance round-trip) and the conditioned baseline.
    assert {:ok, rust_kinematic} =
             RTK.solve_filter_baseline_epochs(
               base_arp,
               epochs,
               [process_noise_baseline_sigma_m: 30.0, filter_kernel: :rust] ++ base_opts
             )

    assert rust_kinematic.metadata.filter_kernel == :rust
    assert length(rust_kinematic.epochs) == length(kinematic.epochs)
    assert_exact_position(rust_kinematic.baseline_m, kinematic.baseline_m)

    for {rust_epoch, elixir_epoch} <- Enum.zip(rust_kinematic.epochs, kinematic.epochs) do
      assert rust_epoch.index == elixir_epoch.index
      assert rust_epoch.integer_status == elixir_epoch.integer_status
      assert rust_epoch.fixed_ambiguities == elixir_epoch.fixed_ambiguities
      assert_exact_position(rust_epoch.baseline_m, elixir_epoch.baseline_m)

      if is_number(rust_epoch.integer_ratio) and is_number(elixir_epoch.integer_ratio) do
        assert rust_epoch.integer_ratio === elixir_epoch.integer_ratio
      end
    end
  end

  defp real_gps_l1_rtk_epochs(sp3, base_obs, rover_obs, count) do
    rover_by_epoch = Map.new(Observations.epochs(rover_obs), &{&1.epoch, &1})

    base_obs
    |> Observations.epochs()
    |> Enum.take(count)
    |> Enum.flat_map(fn base_entry ->
      case Map.fetch(rover_by_epoch, base_entry.epoch) do
        {:ok, rover_entry} ->
          base_values = gps_l1_values(base_obs, base_entry.index)
          rover_values = gps_l1_values(rover_obs, rover_entry.index)

          common =
            base_values
            |> Map.keys()
            |> MapSet.new()
            |> MapSet.intersection(rover_values |> Map.keys() |> MapSet.new())
            |> MapSet.to_list()
            |> Enum.sort()

          epoch = naive_datetime(base_entry.epoch)
          positions = satellite_positions(sp3, epoch, common)

          base_positions =
            transmit_time_satellite_positions(sp3, epoch, base_values, common, :code_m)

          rover_positions =
            transmit_time_satellite_positions(sp3, epoch, rover_values, common, :code_m)

          usable =
            Enum.filter(common, fn sat ->
              Map.has_key?(positions, sat) and Map.has_key?(base_positions, sat) and
                Map.has_key?(rover_positions, sat)
            end)

          if length(usable) >= 4 do
            [
              %{
                epoch: epoch,
                satellite_positions_m: Map.take(positions, usable),
                base_satellite_positions_m: Map.take(base_positions, usable),
                rover_satellite_positions_m: Map.take(rover_positions, usable),
                base_observations: Enum.map(usable, &Map.fetch!(base_values, &1)),
                rover_observations: Enum.map(usable, &Map.fetch!(rover_values, &1))
              }
            ]
          else
            []
          end

        :error ->
          []
      end
    end)
    |> assert_receiver_position_maps()
  end

  defp gps_l1_values(obs, index) do
    {:ok, by_sat} = Observations.values(obs, index, codes: %{"G" => ["C1C", "L1C"]})

    by_sat
    |> Enum.flat_map(fn {sat, values} ->
      values_by_code = Map.new(values, &{&1.code, &1})

      with %{value: c1} when is_number(c1) <- values_by_code["C1C"],
           %{value: l1} = phase when is_number(l1) <- values_by_code["L1C"] do
        [
          {sat,
           %{
             satellite_id: sat,
             code_m: c1,
             phase_m: l1 * @gps_l1_wavelength_m,
             lli: phase.lli
           }}
        ]
      else
        _ -> []
      end
    end)
    |> Map.new()
  end

  defp real_gps_l1_l2_rtk_epochs(sp3, base_obs, rover_obs, count) do
    rover_by_epoch = Map.new(Observations.epochs(rover_obs), &{&1.epoch, &1})

    base_obs
    |> Observations.epochs()
    |> Enum.take(count)
    |> Enum.flat_map(fn base_entry ->
      case Map.fetch(rover_by_epoch, base_entry.epoch) do
        {:ok, rover_entry} ->
          base_values = gps_l1_l2_values(base_obs, base_entry.index)
          rover_values = gps_l1_l2_values(rover_obs, rover_entry.index)

          common =
            base_values
            |> Map.keys()
            |> MapSet.new()
            |> MapSet.intersection(rover_values |> Map.keys() |> MapSet.new())
            |> MapSet.to_list()
            |> Enum.sort()

          epoch = naive_datetime(base_entry.epoch)
          positions = satellite_positions(sp3, epoch, common)

          base_positions =
            transmit_time_satellite_positions(sp3, epoch, base_values, common, :p1_m)

          rover_positions =
            transmit_time_satellite_positions(sp3, epoch, rover_values, common, :p1_m)

          usable =
            Enum.filter(common, fn sat ->
              Map.has_key?(positions, sat) and Map.has_key?(base_positions, sat) and
                Map.has_key?(rover_positions, sat)
            end)

          if length(usable) >= 4 do
            [
              %{
                epoch: epoch,
                satellite_positions_m: Map.take(positions, usable),
                base_satellite_positions_m: Map.take(base_positions, usable),
                rover_satellite_positions_m: Map.take(rover_positions, usable),
                base_observations: Enum.map(usable, &Map.fetch!(base_values, &1)),
                rover_observations: Enum.map(usable, &Map.fetch!(rover_values, &1))
              }
            ]
          else
            []
          end

        :error ->
          []
      end
    end)
  end

  defp gps_l1_l2_values(obs, index) do
    {:ok, by_sat} = Observations.values(obs, index, codes: %{"G" => ["C1C", "C2W", "L1C", "L2W"]})

    by_sat
    |> Enum.flat_map(fn {sat, values} ->
      values_by_code = Map.new(values, &{&1.code, &1})

      with %{value: c1} when is_number(c1) <- values_by_code["C1C"],
           %{value: c2} when is_number(c2) <- values_by_code["C2W"],
           %{value: l1} = phase1 when is_number(l1) <- values_by_code["L1C"],
           %{value: l2} = phase2 when is_number(l2) <- values_by_code["L2W"] do
        [
          {sat,
           %{
             satellite_id: sat,
             p1_m: c1,
             p2_m: c2,
             phi1_cyc: l1,
             phi2_cyc: l2,
             f1_hz: @gps_l1_hz,
             f2_hz: @gps_l2_hz,
             lli1: phase1.lli,
             lli2: phase2.lli
           }}
        ]
      else
        _ -> []
      end
    end)
    |> Map.new()
  end

  defp satellite_positions(sp3, epoch, sats) do
    sats
    |> Enum.reduce(%{}, fn sat, acc ->
      case SP3.position(sp3, sat, epoch) do
        {:ok, %{x_m: x, y_m: y, z_m: z}} -> Map.put(acc, sat, {x, y, z})
        {:error, _reason} -> acc
      end
    end)
  end

  defp transmit_time_satellite_positions(sp3, receive_epoch, values, sats, code_key) do
    sats
    |> Enum.reduce(%{}, fn sat, acc ->
      with %{^code_key => code_m} when is_number(code_m) <- Map.get(values, sat),
           {:ok, transmit_epoch} <- transmit_epoch(receive_epoch, code_m),
           {:ok, %{x_m: x, y_m: y, z_m: z}} <- SP3.position(sp3, sat, transmit_epoch) do
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

  defp assert_receiver_position_maps(epochs) do
    # The solve path should use receiver-specific transmit-time maps when tests
    # provide them. Keep this as a light fixture sanity check instead of a
    # numeric assertion about the baseline.
    assert Enum.all?(epochs, &Map.has_key?(&1, :base_satellite_positions_m))
    assert Enum.all?(epochs, &Map.has_key?(&1, :rover_satellite_positions_m))
    epochs
  end

  defp naive_datetime({{year, month, day}, {hour, minute, second}}) do
    whole_second = trunc(second)
    microsecond = round((second - whole_second) * 1_000_000)

    NaiveDateTime.new!(
      Date.new!(year, month, day),
      Time.new!(hour, minute, whole_second, {microsecond, 6})
    )
  end

  defp arp_position(marker, antenna_h_m),
    do: add3(marker, scale3(marker, antenna_h_m / norm3(marker)))

  defp antenna_height_m(obs) do
    assert {height_m, east_m, north_m} = Observations.antenna_delta_hen(obs)
    assert east_m == 0.0
    assert north_m == 0.0
    height_m
  end

  defp add3({ax, ay, az}, {bx, by, bz}), do: {ax + bx, ay + by, az + bz}
  defp sub3({ax, ay, az}, {bx, by, bz}), do: {ax - bx, ay - by, az - bz}
  defp scale3({x, y, z}, s), do: {x * s, y * s, z * s}
  defp norm3({x, y, z}), do: :math.sqrt(x * x + y * y + z * z)

  defp enu_map_to_tuple(%{"east" => east, "north" => north, "up" => up}), do: {east, north, up}

  defp position_error(%{x_m: x, y_m: y, z_m: z}, {truth_x, truth_y, truth_z}) do
    :math.sqrt(
      (x - truth_x) * (x - truth_x) +
        (y - truth_y) * (y - truth_y) +
        (z - truth_z) * (z - truth_z)
    )
  end

  defp position_error(%{x_m: x, y_m: y, z_m: z}, %{x_m: truth_x, y_m: truth_y, z_m: truth_z}) do
    :math.sqrt(
      (x - truth_x) * (x - truth_x) +
        (y - truth_y) * (y - truth_y) +
        (z - truth_z) * (z - truth_z)
    )
  end

  defp position_error({x, y, z}, {truth_x, truth_y, truth_z}) do
    :math.sqrt(
      (x - truth_x) * (x - truth_x) +
        (y - truth_y) * (y - truth_y) +
        (z - truth_z) * (z - truth_z)
    )
  end

  defp assert_exact_position(%{x_m: x, y_m: y, z_m: z}, %{
         x_m: truth_x,
         y_m: truth_y,
         z_m: truth_z
       }) do
    assert x === truth_x
    assert y === truth_y
    assert z === truth_z
  end
end
