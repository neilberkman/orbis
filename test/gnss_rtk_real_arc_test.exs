defmodule Orbis.GNSS.RTKRealArcTest do
  @moduledoc false

  use ExUnit.Case, async: false

  alias Orbis.GNSS.RINEX.Observations
  alias Orbis.GNSS.RTK
  alias Orbis.GNSS.SP3

  @sp3_path Path.join(__DIR__, "fixtures/sp3/GBM0MGXRAP_20201770000_01D_05M_ORB_120epoch.sp3")
  @wtzr_obs_path Path.join(
                   __DIR__,
                   "fixtures/obs/WTZR00DEU_R_20201770000_01D_30S_MO_120epoch.rnx"
                 )
  @wtzz_obs_path Path.join(
                   __DIR__,
                   "fixtures/obs/WTZZ00DEU_R_20201770000_01D_30S_MO_120epoch.rnx"
                 )
  @rtklib_oracle_path Path.join(__DIR__, "fixtures/rtk/wtzr_wtzz_rtklib_oracle.json")

  @c_m_s 299_792_458.0
  @gps_l1_hz 1_575_420_000.0
  @gps_l2_hz 1_227_600_000.0
  @gps_l1_wavelength_m @c_m_s / @gps_l1_hz

  @wtzr_marker {4_075_580.3111, 931_854.0543, 4_801_568.2808}
  @wtzz_marker {4_075_579.1913, 931_853.3696, 4_801_569.1897}
  @wtzr_antenna_h_m 0.0710
  @wtzz_antenna_h_m 0.2840

  @tag timeout: 180_000
  test "real co-located Wettzell L1 RTK solves the antenna baseline and fixes a safe partial subset" do
    sp3 = SP3.load!(@sp3_path)
    base_obs = Observations.load!(@wtzr_obs_path)
    rover_obs = Observations.load!(@wtzz_obs_path)

    base_arp = arp_position(@wtzr_marker, @wtzr_antenna_h_m)
    rover_arp = arp_position(@wtzz_marker, @wtzz_antenna_h_m)
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

    # The limiting term is the code-pinned integer level on a short
    # single-frequency arc, not the carrier-phase scatter. The best integer
    # candidate is worse than the float baseline, so the integer search must refuse the
    # unsafe fix instead of reporting false centimetre confidence.
    assert fixed_antenna_error_m > float_antenna_error_m
    assert fixed_antenna_error_m < 0.2

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
             "G08",
             "G15@rover#2|ref=G30",
             "G27@rover#1|ref=G30"
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
    assert wide_lane_fixed.metadata.integer_status == :not_fixed
    assert wide_lane_fixed.metadata.integer_ratio < 3.0

    wide_lane_antenna_error_m = position_error(wide_lane_fixed.baseline_m, antenna_baseline)

    assert wide_lane_antenna_error_m < 0.2

    # Dual-frequency partial ambiguity resolution. The full narrow-lane set
    # fails the ratio test (asserted above as :not_fixed), but holding the
    # wide-lane integers fixed collapses the per-ambiguity bias for most
    # satellites. A confidence-ranked subset search therefore fixes a strictly
    # larger safe subset than the single-frequency partial (4), without
    # weakening the ratio threshold.
    assert {:ok, dual_partial} =
             RTK.solve_widelane_fixed_baseline_epochs(base_arp, dual_epochs,
               initial_baseline_m: {0.0, 0.0, 0.0},
               max_iterations: 10,
               on_cycle_slip: :drop_satellite,
               elevation_weighting: true,
               code_sigma_m: 2.0,
               phase_sigma_m: 0.01,
               integer_candidate_limit: 200_000,
               partial_ambiguity_resolution: true,
               partial_min_ambiguities: 4
             )

    dual_partial_antenna_error_m = position_error(dual_partial.baseline_m, antenna_baseline)

    assert dual_partial.metadata.integer_method == :widelane_narrowlane_lambda
    assert dual_partial.metadata.partial_ambiguity_resolution
    assert dual_partial.metadata.partial_fixed
    assert dual_partial.metadata.integer_status == :fixed
    assert dual_partial.metadata.integer_ratio > 3.0

    # Strictly larger than the single-frequency partial subset of 4.
    assert length(dual_partial.metadata.partial_fixed_ambiguities) > 4

    assert dual_partial.metadata.partial_fixed_ambiguities == [
             "G05",
             "G08",
             "G09",
             "G13",
             "G18",
             "G28"
           ]

    # G07 is the lone outlier the subset search drops (its narrow-lane float is
    # several sigma off the integer lattice), so it stays a free ambiguity.
    assert dual_partial.metadata.partial_free_ambiguities == ["G07"]

    # The accepted dual-frequency partial fix is safe: its baseline error beats
    # both the narrow-lane float and the refused full wide-lane/narrow-lane fix.
    float_dual_antenna_error_m =
      position_error(dual_partial.float_solution.baseline_m, antenna_baseline)

    assert dual_partial_antenna_error_m < float_dual_antenna_error_m
    assert dual_partial_antenna_error_m < wide_lane_antenna_error_m
    assert dual_partial_antenna_error_m < 0.06

    # The exhaustive subset search is globally bounded: it reports how many
    # candidate subsets it evaluated, and that stays under the hard cap.
    evaluated = dual_partial.metadata.partial_exhaustive_subsets_evaluated
    assert is_integer(evaluated) and evaluated > 0
    assert evaluated <= 20_000
  end

  @tag timeout: 180_000
  test "RTKLIB fixes the two-epoch prefix that the current batch partial fix gets wrong" do
    oracle =
      @rtklib_oracle_path
      |> File.read!()
      |> Jason.decode!()

    sp3 = SP3.load!(@sp3_path)
    base_obs = Observations.load!(@wtzr_obs_path)
    rover_obs = Observations.load!(@wtzz_obs_path)
    base_arp = arp_position(@wtzr_marker, @wtzr_antenna_h_m)
    rover_arp = arp_position(@wtzz_marker, @wtzz_antenna_h_m)
    antenna_baseline = sub3(rover_arp, base_arp)
    antenna_baseline_enu = enu_map_to_tuple(oracle["truth"]["antenna_baseline_enu_m"])

    rtklib_epoch_2 = oracle["per_epoch"] |> Enum.at(1)
    rtklib_epoch_2_baseline = enu_map_to_tuple(rtklib_epoch_2["baseline_enu_m"])

    assert oracle["reference"]["first_fixed_index"] == 1
    assert rtklib_epoch_2["fix_status"] == "fixed"
    assert rtklib_epoch_2["ratio"] >= 3.0
    assert position_error(rtklib_epoch_2_baseline, antenna_baseline_enu) < 0.006

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
    assert sol.metadata.partial_fixed
    assert sol.metadata.integer_ratio >= 3.0
    assert sol.metadata.integer_method == :lambda
    assert sol.metadata.n_epochs == 2
    assert sol.metadata.elevation_mask_deg == 10.0
    assert sol.metadata.elevation_masked_sats == ["G08", "G18", "G27"]

    # This pins the current gap: a batch partial-AR solve can pass the ratio
    # test on the two-epoch prefix while landing far from the same ARP truth that
    # RTKLIB fixes at epoch 2. Matching RTKLIB's 10-degree elevation mask removes
    # the low-elevation G08/G18/G27 contributors and improves the result, but the
    # sequential filter must still eliminate the remaining false confidence, not
    # merely produce a fixed status.
    error_m = position_error(sol.baseline_m, antenna_baseline)
    assert error_m > 0.3
    assert error_m < 0.5

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

    assert filter.metadata.first_fixed_index == nil
    assert Enum.all?(filter.epochs, &(&1.integer_status == :not_fixed))
    assert List.last(filter.epochs).integer_ratio < 3.0

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

    # RTKLIB's variance shape is part of the parity lane, but it is not enough
    # by itself to justify a fix on the two-epoch prefix. The remaining gap is
    # the observable/correction model, not a ratio-threshold or covariance-form
    # tweak.
    assert rtklib_weighted_filter.metadata.measurement_covariance.stochastic_model == :rtklib
    assert rtklib_weighted_filter.metadata.first_fixed_index == nil
    assert Enum.all?(rtklib_weighted_filter.epochs, &(&1.integer_status == :not_fixed))
    assert List.last(rtklib_weighted_filter.epochs).integer_ratio < 3.0
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
          usable = Enum.filter(common, &Map.has_key?(positions, &1))

          if length(usable) >= 4 do
            [
              %{
                epoch: epoch,
                satellite_positions_m: Map.take(positions, usable),
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
          usable = Enum.filter(common, &Map.has_key?(positions, &1))

          if length(usable) >= 4 do
            [
              %{
                epoch: epoch,
                satellite_positions_m: Map.take(positions, usable),
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

  defp position_error({x, y, z}, {truth_x, truth_y, truth_z}) do
    :math.sqrt(
      (x - truth_x) * (x - truth_x) +
        (y - truth_y) * (y - truth_y) +
        (z - truth_z) * (z - truth_z)
    )
  end
end
