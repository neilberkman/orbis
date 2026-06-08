defmodule Orbis.GNSS.RTKTest do
  use ExUnit.Case, async: true

  alias Orbis.GNSS.RTK

  @base {1_110_000.0, -4_840_000.0, 3_980_000.0}
  @truth_baseline {12.5, -4.25, 2.75}
  @c 299_792_458.0
  @f_l1 1_575_420_000.0
  @f_l2 1_227_600_000.0
  @l1_wavelength_m 299_792_458.0 / 1_575_420_000.0
  @l2_wavelength_m @c / @f_l2
  @narrow_lane_wavelength_m @c / (@f_l1 + @f_l2)
  @sat_positions [
    %{
      "G01" => {21_000_000.0, 14_000_000.0, 20_000_000.0},
      "G02" => {-18_000_000.0, 19_000_000.0, 18_500_000.0},
      "G03" => {15_000_000.0, -21_000_000.0, 17_000_000.0},
      "G04" => {-20_000_000.0, -12_000_000.0, 21_000_000.0},
      "G05" => {24_000_000.0, -5_000_000.0, 16_000_000.0}
    },
    %{
      "G01" => {21_020_000.0, 13_960_000.0, 20_010_000.0},
      "G02" => {-18_030_000.0, 19_020_000.0, 18_470_000.0},
      "G03" => {15_040_000.0, -20_970_000.0, 17_020_000.0},
      "G04" => {-19_980_000.0, -12_050_000.0, 21_030_000.0},
      "G05" => {23_960_000.0, -4_970_000.0, 16_050_000.0}
    },
    %{
      "G01" => {21_050_000.0, 13_910_000.0, 20_030_000.0},
      "G02" => {-18_070_000.0, 19_050_000.0, 18_430_000.0},
      "G03" => {15_090_000.0, -20_930_000.0, 17_060_000.0},
      "G04" => {-19_950_000.0, -12_090_000.0, 21_060_000.0},
      "G05" => {23_930_000.0, -4_940_000.0, 16_100_000.0}
    }
  ]
  @ambiguities %{"G01" => 0.0, "G02" => 0.42, "G03" => -0.73, "G04" => 1.12, "G05" => -0.38}
  @fixed_cycles %{"G01" => 0, "G02" => 5, "G03" => -7, "G04" => 12, "G05" => -4}
  @wide_lane_cycles %{"G01" => 0, "G02" => 3, "G03" => -5, "G04" => 8, "G05" => -2}

  describe "double_differences/3" do
    test "receiver clocks and common satellite errors cancel" do
      sats = ["G01", "G02", "G03", "G04"]
      reference = "G01"
      base_clock_m = 125.0
      rover_clock_m = -42.0

      base_ranges = %{"G01" => 20_000.0, "G02" => 21_000.0, "G03" => 22_500.0, "G04" => 23_100.0}
      rover_ranges = %{"G01" => 20_010.0, "G02" => 21_025.0, "G03" => 22_480.0, "G04" => 23_150.0}
      common_errors = %{"G01" => 3.25, "G02" => -12.0, "G03" => 8.5, "G04" => 1.0}
      base_phase_ambiguities = %{"G01" => 2.0, "G02" => -3.0, "G03" => 7.0, "G04" => 11.0}
      rover_phase_ambiguities = %{"G01" => 5.0, "G02" => 4.0, "G03" => 1.0, "G04" => 19.0}

      base =
        synth_observations(sats, base_ranges, base_clock_m, common_errors, base_phase_ambiguities)

      rover =
        synth_observations(
          sats,
          rover_ranges,
          rover_clock_m,
          common_errors,
          rover_phase_ambiguities
        )

      assert {:ok, result} =
               RTK.double_differences(base, rover, reference_satellite_id: reference)

      assert result.reference_satellite_id == reference
      assert result.dropped_sats == []

      by_sat = Map.new(result.double_differences, &{&1.satellite_id, &1})

      for sat <- sats -- [reference] do
        expected_code =
          Map.fetch!(rover_ranges, sat) - Map.fetch!(base_ranges, sat) -
            (Map.fetch!(rover_ranges, reference) - Map.fetch!(base_ranges, reference))

        expected_phase =
          expected_code +
            (Map.fetch!(rover_phase_ambiguities, sat) -
               Map.fetch!(base_phase_ambiguities, sat)) -
            (Map.fetch!(rover_phase_ambiguities, reference) -
               Map.fetch!(base_phase_ambiguities, reference))

        assert by_sat[sat].reference_satellite_id == reference
        assert by_sat[sat].code_m == expected_code
        assert by_sat[sat].phase_m == expected_phase
      end
    end

    test "selects a deterministic default reference and reports dropped satellites" do
      base = [{"G02", 210.0, 211.0}, {"G01", 100.0, 101.0}, {"G09", 900.0, 901.0}]
      rover = [{"G02", 230.0, 233.0}, {"G01", 105.0, 108.0}, {"G10", 1000.0, 1001.0}]

      assert {:ok, result} = RTK.double_differences(base, rover)

      assert result.reference_satellite_id == "G01"
      assert result.dropped_sats == ["G09", "G10"]

      assert result.double_differences == [
               %{
                 satellite_id: "G02",
                 reference_satellite_id: "G01",
                 ambiguity_id: "G02",
                 code_m: 15.0,
                 phase_m: 15.0
               }
             ]
    end

    test "reports double-difference ambiguity ids when carrier arcs are explicit" do
      base = [
        %{satellite_id: "G01", code_m: 100.0, phase_m: 101.0},
        %{satellite_id: "G02", code_m: 210.0, phase_m: 211.0}
      ]

      rover = [
        %{satellite_id: "G01", code_m: 105.0, phase_m: 108.0, ambiguity_id: "G01#2"},
        %{satellite_id: "G02", code_m: 230.0, phase_m: 233.0, ambiguity_id: "G02#2"}
      ]

      assert {:ok, result} = RTK.double_differences(base, rover, reference_satellite_id: "G01")

      assert [%{ambiguity_id: ambiguity_id}] = result.double_differences

      assert ambiguity_id == "G02#2|ref=G01#2"
    end

    test "bad inputs are tagged" do
      assert RTK.double_differences([{"G01", 1.0, 2.0}], [{"G01", 1.0, 2.0}]) ==
               {:error, {:too_few_common_satellites, 1, 2}}

      assert RTK.double_differences(
               [{"G01", 1.0, 2.0}, {"G01", 3.0, 4.0}],
               [{"G01", 1.0, 2.0}, {"G02", 3.0, 4.0}]
             ) == {:error, {:duplicate_observation, "G01"}}

      assert RTK.double_differences(
               [{"G01", 1.0, 2.0}, {"G02", 3.0, 4.0}],
               [{"G01", 1.0, 2.0}, {"G02", 3.0, 4.0}],
               reference_satellite_id: "G99"
             ) == {:error, {:reference_satellite_missing, "G99"}}

      assert RTK.double_differences([{"G01", :bad, 2.0}], [{"G01", 1.0, 2.0}]) ==
               {:error, {:invalid_base_observations, {"G01", :bad, 2.0}}}

      assert RTK.double_differences([{"G01", 1.0, 2.0}], [{"G01", 1.0, :bad}]) ==
               {:error, {:invalid_rover_observations, {"G01", 1.0, :bad}}}
    end
  end

  describe "solve_float_baseline_epochs/3" do
    test "recovers a static baseline and float DD ambiguities from a wrong seed" do
      epochs =
        @sat_positions
        |> Enum.with_index()
        |> Enum.map(fn {positions, idx} ->
          synthetic_baseline_epoch(@base, @truth_baseline, positions,
            epoch: idx,
            base_clock_m: 120.0 - 3.0 * idx,
            rover_clock_m: -45.0 + 7.0 * idx,
            common_errors_m: %{
              "G01" => 2.0 + idx,
              "G02" => -3.5,
              "G03" => 1.25 * idx,
              "G04" => -0.75,
              "G05" => 4.0
            },
            ambiguities_m: @ambiguities
          )
        end)

      assert {:ok, sol} =
               RTK.solve_float_baseline_epochs(@base, epochs,
                 reference_satellite_id: "G01",
                 initial_baseline_m: {-40.0, 35.0, 12.0}
               )

      assert sol.reference_satellite_id == "G01"
      assert sol.used_sats == ["G02", "G03", "G04", "G05"]
      assert sol.metadata.converged
      assert sol.metadata.n_epochs == 3
      assert sol.metadata.dropped_sats == []

      assert sol.metadata.measurement_covariance == %{
               model: :double_difference,
               code_sigma_m: 1.0,
               phase_sigma_m: 0.02,
               elevation_weighting: false,
               min_elevation_sin: 0.05
             }

      assert sol.metadata.ambiguity_float.order == sol.used_sats
      assert length(sol.metadata.ambiguity_float.covariance_m) == length(sol.used_sats)
      assert nonzero_off_diagonal?(sol.metadata.ambiguity_float.covariance_m)

      assert_identity(
        matmul(
          sol.metadata.ambiguity_float.covariance_m,
          sol.metadata.ambiguity_float.covariance_inverse_m
        ),
        1.0e-6
      )

      assert position_error(sol.baseline_m, @truth_baseline) < 1.0e-5

      for {sat, expected} <- Map.delete(@ambiguities, "G01") do
        assert abs(Map.fetch!(sol.ambiguities_m, sat) - expected) < 1.0e-5
      end

      for residual <- sol.residuals_m do
        assert abs(residual.code_m) < 1.0e-5
        assert abs(residual.phase_m) < 1.0e-5
      end
    end

    test "selects the highest-elevation default reference and uses epoch-local satellites" do
      [first, second | _] = @sat_positions

      epoch_a =
        synthetic_baseline_epoch(@base, @truth_baseline, first,
          ambiguities_m: @ambiguities,
          epoch: :a
        )

      epoch_b =
        synthetic_baseline_epoch(@base, @truth_baseline, Map.delete(second, "G05"),
          ambiguities_m: @ambiguities,
          epoch: :b
        )

      assert {:ok, sol} = RTK.solve_float_baseline_epochs(@base, [epoch_a, epoch_b])

      assert sol.reference_satellite_id == "G03"
      assert sol.used_sats == ["G01", "G02", "G04", "G05"]
      assert sol.metadata.dropped_sats == []
      assert sol.metadata.n_observations == 14
      assert position_error(sol.baseline_m, @truth_baseline) < 1.0e-4
    end

    test "can use elevation-dependent stochastic weighting" do
      epochs =
        @sat_positions
        |> Enum.with_index()
        |> Enum.map(fn {positions, idx} ->
          synthetic_baseline_epoch(@base, @truth_baseline, positions,
            epoch: idx,
            ambiguities_m: @ambiguities
          )
        end)

      assert {:ok, unweighted} =
               RTK.solve_float_baseline_epochs(@base, epochs,
                 reference_satellite_id: "G01",
                 initial_baseline_m: {-40.0, 35.0, 12.0}
               )

      assert {:ok, weighted} =
               RTK.solve_float_baseline_epochs(@base, epochs,
                 reference_satellite_id: "G01",
                 elevation_weighting: true,
                 initial_baseline_m: {-40.0, 35.0, 12.0}
               )

      assert weighted.metadata.measurement_covariance.elevation_weighting

      refute unweighted.metadata.ambiguity_float.covariance_m ==
               weighted.metadata.ambiguity_float.covariance_m

      assert position_error(weighted.baseline_m, @truth_baseline) < 1.0e-5
    end

    test "can Hatch-smooth code observations before forming double differences" do
      epochs =
        @sat_positions
        |> Enum.with_index()
        |> Enum.map(fn {positions, idx} ->
          synthetic_baseline_epoch(@base, @truth_baseline, positions,
            epoch: idx,
            ambiguities_m: @ambiguities
          )
        end)
        |> add_rover_code_noise(%{"G02" => [0.6, -0.3, 0.2], "G04" => [-0.5, 0.4, -0.1]})

      assert {:ok, raw} =
               RTK.solve_float_baseline_epochs(@base, epochs,
                 reference_satellite_id: "G01",
                 initial_baseline_m: {-40.0, 35.0, 12.0},
                 code_sigma_m: 100.0
               )

      assert {:ok, smoothed} =
               RTK.solve_float_baseline_epochs(@base, epochs,
                 reference_satellite_id: "G01",
                 initial_baseline_m: {-40.0, 35.0, 12.0},
                 code_sigma_m: 100.0,
                 code_smoothing: true,
                 hatch_window_cap: 2
               )

      assert smoothed.metadata.code_smoothing
      assert smoothed.metadata.code_smoothing_window_cap == 2
      assert smoothed.metadata.code_rms_m < raw.metadata.code_rms_m
    end

    test "defaults to an error on LLI cycle slips" do
      [first, second | _] = @sat_positions

      epoch_a =
        synthetic_baseline_epoch(@base, @truth_baseline, first,
          ambiguities_m: @ambiguities,
          epoch: 0
        )

      epoch_b =
        @base
        |> synthetic_baseline_epoch(@truth_baseline, second,
          ambiguities_m: @ambiguities,
          epoch: 1
        )
        |> mark_rover_lli("G02", 1)

      assert RTK.solve_float_baseline_epochs(@base, [epoch_a, epoch_b]) ==
               {:error, {:cycle_slip_detected, :rover, "G02", 1, [:lli]}}
    end

    test "can drop satellites with LLI cycle slips" do
      epochs =
        @sat_positions
        |> Enum.with_index()
        |> Enum.map(fn {positions, idx} ->
          epoch =
            synthetic_baseline_epoch(@base, @truth_baseline, positions,
              ambiguities_m: @ambiguities,
              epoch: idx
            )

          if idx == 1, do: mark_rover_lli(epoch, "G02", 1), else: epoch
        end)

      assert {:ok, sol} =
               RTK.solve_float_baseline_epochs(@base, epochs,
                 reference_satellite_id: "G01",
                 on_cycle_slip: :drop_satellite,
                 initial_baseline_m: {-40.0, 35.0, 12.0}
               )

      assert sol.used_sats == ["G03", "G04", "G05"]
      assert sol.metadata.dropped_cycle_slip_sats == ["G02"]
      assert sol.metadata.dropped_sats == ["G02"]
      assert position_error(sol.baseline_m, @truth_baseline) < 1.0e-5
    end

    test "can split ambiguity arcs at LLI cycle slips" do
      epochs =
        @sat_positions
        |> Enum.with_index()
        |> Enum.map(fn {positions, idx} ->
          ambiguities =
            if idx == 0 do
              Map.put(@ambiguities, "G02", 5.0 * @l1_wavelength_m)
            else
              Map.put(@ambiguities, "G02", 8.0 * @l1_wavelength_m)
            end

          epoch =
            synthetic_baseline_epoch(@base, @truth_baseline, positions,
              ambiguities_m: ambiguities,
              epoch: idx
            )

          if idx == 1, do: mark_rover_lli(epoch, "G02", 1), else: epoch
        end)

      assert {:ok, sol} =
               RTK.solve_float_baseline_epochs(@base, epochs,
                 reference_satellite_id: "G01",
                 on_cycle_slip: :split_arc,
                 initial_baseline_m: {-40.0, 35.0, 12.0}
               )

      g02_ids = Enum.filter(sol.used_sats, &String.contains?(&1, "G02"))
      assert length(g02_ids) == 2

      split_ids = sol.metadata.split_cycle_slip_arcs |> Enum.map(& &1.ambiguity_id) |> Enum.sort()
      assert length(split_ids) == 2

      g02_ambiguities =
        g02_ids
        |> Enum.map(&Map.fetch!(sol.ambiguities_m, &1))
        |> Enum.sort()

      assert Enum.zip(g02_ambiguities, [5.0 * @l1_wavelength_m, 8.0 * @l1_wavelength_m])
             |> Enum.all?(fn {got, expected} -> abs(got - expected) < 1.0e-5 end)

      assert position_error(sol.baseline_m, @truth_baseline) < 1.0e-5
    end

    test "bad baseline-solve inputs are tagged" do
      epoch = synthetic_baseline_epoch(@base, @truth_baseline, hd(@sat_positions))

      assert RTK.solve_float_baseline_epochs(:bad, [epoch]) == {:error, :invalid_base_position}
      assert RTK.solve_float_baseline_epochs(@base, []) == {:error, :no_epochs}
      assert RTK.solve_float_baseline_epochs(@base, :bad) == {:error, :invalid_epochs}

      assert RTK.solve_float_baseline_epochs(@base, [
               %{base_observations: [], rover_observations: []}
             ]) ==
               {:error, {:invalid_epoch_observations, 0}}

      assert RTK.solve_float_baseline_epochs(@base, [epoch], reference_satellite_id: "G99") ==
               {:error, {:reference_satellite_missing, "G99"}}

      assert RTK.solve_float_baseline_epochs(@base, [epoch], phase_sigma_m: 0.0) ==
               {:error, {:invalid_sigma, :phase_sigma_m}}

      assert RTK.solve_float_baseline_epochs(@base, [epoch], max_iterations: 0) ==
               {:error, {:invalid_option, :max_iterations}}

      assert RTK.solve_float_baseline_epochs(@base, [epoch], on_cycle_slip: :bad) ==
               {:error, {:invalid_option, :on_cycle_slip}}

      assert RTK.solve_float_baseline_epochs(@base, [epoch], elevation_weighting: :bad) ==
               {:error, {:invalid_option, :elevation_weighting}}

      assert RTK.solve_float_baseline_epochs(@base, [epoch], code_smoothing: :bad) ==
               {:error, {:invalid_option, :code_smoothing}}

      assert RTK.solve_float_baseline_epochs(@base, [epoch],
               code_smoothing: true,
               hatch_window_cap: 0
             ) == {:error, {:invalid_option, :hatch_window_cap}}
    end
  end

  describe "solve_fixed_baseline_epochs/3" do
    test "recovers a fixed static baseline and integer DD ambiguities" do
      ambiguities_m =
        Map.new(@fixed_cycles, fn {sat, cycles} -> {sat, cycles * @l1_wavelength_m} end)

      epochs =
        @sat_positions
        |> Enum.with_index()
        |> Enum.map(fn {positions, idx} ->
          synthetic_baseline_epoch(@base, @truth_baseline, positions,
            epoch: idx,
            base_clock_m: 80.0 + idx,
            rover_clock_m: -25.0 + 2.0 * idx,
            common_errors_m: %{
              "G01" => 1.5,
              "G02" => -2.0 + idx,
              "G03" => 0.25 * idx,
              "G04" => 3.0,
              "G05" => -0.5
            },
            ambiguities_m: ambiguities_m
          )
        end)

      assert {:ok, sol} =
               RTK.solve_fixed_baseline_epochs(@base, epochs,
                 reference_satellite_id: "G01",
                 ambiguity_wavelength_m: @l1_wavelength_m,
                 initial_baseline_m: {-40.0, 35.0, 12.0}
               )

      assert sol.reference_satellite_id == "G01"
      assert sol.used_sats == ["G02", "G03", "G04", "G05"]
      assert sol.metadata.integer_status == :fixed
      assert sol.metadata.integer_method == :bounded_ils
      assert sol.metadata.integer_ratio > 1.0e6
      assert sol.metadata.ambiguity_search.order == sol.used_sats
      assert position_error(sol.baseline_m, @truth_baseline) < 1.0e-5

      for sat <- sol.used_sats do
        expected_cycles = Map.fetch!(@fixed_cycles, sat)
        assert Map.fetch!(sol.fixed_ambiguities_cycles, sat) == expected_cycles

        assert abs(Map.fetch!(sol.fixed_ambiguities_m, sat) - expected_cycles * @l1_wavelength_m) <
                 1.0e-12
      end

      for residual <- sol.residuals_m do
        assert abs(residual.code_m) < 1.0e-5
        assert abs(residual.phase_m) < 1.0e-5
      end
    end

    test "fixed solve respects split ambiguity arcs" do
      epochs =
        @sat_positions
        |> Enum.with_index()
        |> Enum.map(fn {positions, idx} ->
          cycles = if idx == 0, do: 5, else: 8

          ambiguities_m =
            @fixed_cycles
            |> Map.put("G02", cycles)
            |> Map.new(fn {sat, sat_cycles} -> {sat, sat_cycles * @l1_wavelength_m} end)

          epoch =
            synthetic_baseline_epoch(@base, @truth_baseline, positions,
              epoch: idx,
              ambiguities_m: ambiguities_m
            )

          if idx == 1, do: mark_rover_lli(epoch, "G02", 1), else: epoch
        end)

      assert {:ok, sol} =
               RTK.solve_fixed_baseline_epochs(@base, epochs,
                 reference_satellite_id: "G01",
                 ambiguity_wavelength_m: @l1_wavelength_m,
                 on_cycle_slip: :split_arc,
                 initial_baseline_m: {-40.0, 35.0, 12.0}
               )

      g02_ids = Enum.filter(sol.used_sats, &String.contains?(&1, "G02"))
      assert length(g02_ids) == 2
      assert sol.metadata.integer_status == :fixed
      assert sol.metadata.dropped_cycle_slip_sats == []
      assert length(sol.metadata.split_cycle_slip_arcs) == 2

      g02_cycles =
        g02_ids
        |> Enum.map(&Map.fetch!(sol.fixed_ambiguities_cycles, &1))
        |> Enum.sort()

      assert g02_cycles == [5, 8]
      assert position_error(sol.baseline_m, @truth_baseline) < 1.0e-5
    end

    test "fixed solve supports per-ambiguity metre offsets" do
      offsets_m = %{"G02" => 0.35, "G03" => -0.22, "G04" => 0.12, "G05" => -0.41}

      ambiguities_m =
        Map.new(@fixed_cycles, fn
          {"G01", _cycles} -> {"G01", 0.0}
          {sat, cycles} -> {sat, Map.fetch!(offsets_m, sat) + cycles * @l1_wavelength_m}
        end)

      epochs =
        @sat_positions
        |> Enum.with_index()
        |> Enum.map(fn {positions, idx} ->
          synthetic_baseline_epoch(@base, @truth_baseline, positions,
            epoch: idx,
            ambiguities_m: ambiguities_m
          )
        end)

      assert {:ok, sol} =
               RTK.solve_fixed_baseline_epochs(@base, epochs,
                 reference_satellite_id: "G01",
                 ambiguity_wavelength_m: @l1_wavelength_m,
                 ambiguity_offset_m: offsets_m,
                 initial_baseline_m: {-40.0, 35.0, 12.0}
               )

      assert sol.metadata.integer_status == :fixed
      assert sol.metadata.ambiguity_offsets_m == offsets_m
      assert position_error(sol.baseline_m, @truth_baseline) < 1.0e-5

      for sat <- sol.used_sats do
        expected_cycles = Map.fetch!(@fixed_cycles, sat)
        expected_m = Map.fetch!(offsets_m, sat) + expected_cycles * @l1_wavelength_m

        assert Map.fetch!(sol.fixed_ambiguities_cycles, sat) == expected_cycles
        assert abs(Map.fetch!(sol.fixed_ambiguities_m, sat) - expected_m) < 1.0e-12
      end
    end

    test "bad fixed-baseline inputs are tagged" do
      epoch = synthetic_baseline_epoch(@base, @truth_baseline, hd(@sat_positions))

      assert RTK.solve_fixed_baseline_epochs(@base, [epoch]) ==
               {:error, :ambiguity_wavelength_required}

      assert RTK.solve_fixed_baseline_epochs(@base, [epoch], ambiguity_wavelength_m: 0.0) ==
               {:error, {:invalid_option, :ambiguity_wavelength_m}}

      assert RTK.solve_fixed_baseline_epochs(@base, [epoch],
               ambiguity_wavelength_m: %{"G02" => @l1_wavelength_m}
             ) == {:error, {:invalid_ambiguity_wavelength, "G01"}}

      assert RTK.solve_fixed_baseline_epochs(@base, [epoch],
               ambiguity_wavelength_m: @l1_wavelength_m,
               integer_ratio_threshold: -1.0
             ) == {:error, {:invalid_option, :integer_ratio_threshold}}

      assert RTK.solve_fixed_baseline_epochs(@base, [epoch],
               reference_satellite_id: "G01",
               ambiguity_wavelength_m: @l1_wavelength_m,
               ambiguity_offset_m: %{"G02" => 0.25}
             ) == {:error, {:invalid_ambiguity_offset, "G03"}}

      assert RTK.solve_fixed_baseline_epochs(@base, [epoch],
               reference_satellite_id: "G01",
               ambiguity_wavelength_m: @l1_wavelength_m,
               ambiguity_offset_m: :bad
             ) == {:error, {:invalid_option, :ambiguity_offset_m}}
    end
  end

  describe "solve_widelane_fixed_baseline_epochs/3" do
    test "fixes wide-lane then narrow-lane DD ambiguities" do
      epochs =
        @sat_positions
        |> Enum.with_index()
        |> Enum.map(fn {positions, idx} ->
          synthetic_dual_baseline_epoch(@base, @truth_baseline, positions,
            epoch: idx,
            base_clock_m: 70.0 - idx,
            rover_clock_m: -31.0 + 2.0 * idx,
            n1_cycles: @fixed_cycles,
            wide_lane_cycles: @wide_lane_cycles
          )
        end)

      assert {:ok, sol} =
               RTK.solve_widelane_fixed_baseline_epochs(@base, epochs,
                 reference_satellite_id: "G01",
                 initial_baseline_m: {-40.0, 35.0, 12.0},
                 wide_lane_tolerance_cycles: 0.01
               )

      assert sol.metadata.integer_status == :fixed
      assert sol.metadata.integer_method == :widelane_narrowlane_bounded_ils
      assert sol.metadata.wide_lane_fixed
      assert sol.wide_lane_ambiguities_cycles == Map.delete(@wide_lane_cycles, "G01")
      assert sol.metadata.wide_lane_ambiguities_cycles == sol.wide_lane_ambiguities_cycles
      assert sol.metadata.ambiguity_offsets_m == expected_narrow_lane_offsets(@wide_lane_cycles)
      assert position_error(sol.baseline_m, @truth_baseline) < 1.0e-5

      for sat <- sol.used_sats do
        assert Map.fetch!(sol.fixed_ambiguities_cycles, sat) == Map.fetch!(@fixed_cycles, sat)

        expected_m =
          narrow_lane_offset_m(Map.fetch!(@wide_lane_cycles, sat)) +
            Map.fetch!(@fixed_cycles, sat) * @narrow_lane_wavelength_m

        assert abs(Map.fetch!(sol.fixed_ambiguities_m, sat) - expected_m) < 1.0e-9
      end
    end

    test "can split dual-frequency ambiguity arcs at cycle slips" do
      epochs =
        @sat_positions
        |> Enum.with_index()
        |> Enum.map(fn {positions, idx} ->
          {n1_cycles, wide_lane_cycles} =
            if idx == 0 do
              {@fixed_cycles, @wide_lane_cycles}
            else
              {
                Map.put(@fixed_cycles, "G02", 9),
                Map.put(@wide_lane_cycles, "G02", 6)
              }
            end

          epoch =
            synthetic_dual_baseline_epoch(@base, @truth_baseline, positions,
              epoch: idx,
              n1_cycles: n1_cycles,
              wide_lane_cycles: wide_lane_cycles
            )

          if idx == 1, do: mark_dual_rover_lli(epoch, "G02", 1), else: epoch
        end)

      assert {:ok, sol} =
               RTK.solve_widelane_fixed_baseline_epochs(@base, epochs,
                 reference_satellite_id: "G01",
                 on_cycle_slip: :split_arc,
                 wide_lane_min_epochs: 1,
                 wide_lane_tolerance_cycles: 0.01,
                 initial_baseline_m: {-40.0, 35.0, 12.0}
               )

      g02_ids = Enum.filter(sol.used_sats, &String.contains?(&1, "G02"))
      assert length(g02_ids) == 2
      assert sol.metadata.integer_status == :fixed
      assert length(sol.metadata.split_cycle_slip_arcs) == 2

      assert g02_ids
             |> Enum.map(&Map.fetch!(sol.fixed_ambiguities_cycles, &1))
             |> Enum.sort() == [5, 9]

      assert g02_ids
             |> Enum.map(&Map.fetch!(sol.wide_lane_ambiguities_cycles, &1))
             |> Enum.sort() == [3, 6]

      assert position_error(sol.baseline_m, @truth_baseline) < 1.0e-5
    end

    test "split-arc dual-frequency solve skips short fragments and keeps valid fragments" do
      epochs =
        @sat_positions
        |> Enum.with_index()
        |> Enum.map(fn {positions, idx} ->
          {n1_cycles, wide_lane_cycles} =
            if idx == 0 do
              {@fixed_cycles, @wide_lane_cycles}
            else
              {
                Map.put(@fixed_cycles, "G02", 9),
                Map.put(@wide_lane_cycles, "G02", 6)
              }
            end

          epoch =
            synthetic_dual_baseline_epoch(@base, @truth_baseline, positions,
              epoch: idx,
              n1_cycles: n1_cycles,
              wide_lane_cycles: wide_lane_cycles
            )

          if idx == 1, do: mark_dual_rover_lli(epoch, "G02", 1), else: epoch
        end)

      assert {:ok, sol} =
               RTK.solve_widelane_fixed_baseline_epochs(@base, epochs,
                 reference_satellite_id: "G01",
                 on_cycle_slip: :split_arc,
                 wide_lane_min_epochs: 2,
                 wide_lane_tolerance_cycles: 0.01,
                 initial_baseline_m: {-40.0, 35.0, 12.0}
               )

      g02_ids = Enum.filter(sol.used_sats, &String.contains?(&1, "G02"))
      assert g02_ids == ["G02@rover#2|ref=G01"]
      assert sol.metadata.integer_status == :fixed
      refute Map.has_key?(sol.wide_lane_ambiguities_cycles, "G02")
      assert sol.wide_lane_ambiguities_cycles["G02@rover#2|ref=G01"] == 6
      assert position_error(sol.baseline_m, @truth_baseline) < 1.0e-5
    end

    test "bad dual-frequency inputs are tagged" do
      epoch =
        synthetic_dual_baseline_epoch(@base, @truth_baseline, hd(@sat_positions),
          n1_cycles: @fixed_cycles,
          wide_lane_cycles: @wide_lane_cycles
        )

      assert RTK.solve_widelane_fixed_baseline_epochs(@base, []) == {:error, :no_epochs}

      malformed =
        update_in(epoch.base_observations, fn observations ->
          Enum.map(observations, fn
            %{satellite_id: "G02"} = obs -> %{obs | f2_hz: @f_l1}
            obs -> obs
          end)
        end)

      assert {:error, {:wide_lane_failed, "G02", :equal_frequencies}} =
               RTK.solve_widelane_fixed_baseline_epochs(@base, [malformed],
                 reference_satellite_id: "G01"
               )

      assert RTK.solve_widelane_fixed_baseline_epochs(@base, [epoch], wide_lane_min_epochs: 0) ==
               {:error, {:invalid_option, :wide_lane_min_epochs}}
    end
  end

  defp synth_observations(sats, ranges, clock_m, errors, phase_ambiguities) do
    Enum.map(sats, fn sat ->
      code = Map.fetch!(ranges, sat) + clock_m + Map.fetch!(errors, sat)
      phase = code + Map.fetch!(phase_ambiguities, sat)
      %{satellite_id: sat, code_m: code, phase_m: phase}
    end)
  end

  defp position_error(%{x_m: x, y_m: y, z_m: z}, {tx, ty, tz}) do
    :math.sqrt((x - tx) * (x - tx) + (y - ty) * (y - ty) + (z - tz) * (z - tz))
  end

  defp nonzero_off_diagonal?(matrix) do
    matrix
    |> Enum.with_index()
    |> Enum.any?(fn {row, i} ->
      row
      |> Enum.with_index()
      |> Enum.any?(fn {value, j} -> i != j and abs(value) > 1.0e-12 end)
    end)
  end

  defp assert_identity(matrix, tol) do
    matrix
    |> Enum.with_index()
    |> Enum.each(fn {row, i} ->
      row
      |> Enum.with_index()
      |> Enum.each(fn {value, j} ->
        expected = if i == j, do: 1.0, else: 0.0
        assert abs(value - expected) < tol
      end)
    end)
  end

  defp matmul(a, b) do
    b_t = transpose(b)

    Enum.map(a, fn row ->
      Enum.map(b_t, fn col ->
        row
        |> Enum.zip(col)
        |> Enum.reduce(0.0, fn {x, y}, acc -> acc + x * y end)
      end)
    end)
  end

  defp transpose(matrix) do
    matrix
    |> Enum.zip()
    |> Enum.map(&Tuple.to_list/1)
  end

  defp mark_rover_lli(epoch, sat, lli) do
    %{epoch | rover_observations: mark_observation_lli(epoch.rover_observations, sat, lli)}
  end

  defp mark_observation_lli(observations, sat, lli) do
    Enum.map(observations, fn
      {^sat, code, phase} -> %{satellite_id: sat, code_m: code, phase_m: phase, lli: lli}
      %{satellite_id: ^sat} = obs -> Map.put(obs, :lli, lli)
      obs -> obs
    end)
  end

  defp add_rover_code_noise(epochs, noise_by_sat) do
    epochs
    |> Enum.with_index()
    |> Enum.map(fn {epoch, idx} ->
      rover =
        Enum.map(epoch.rover_observations, fn
          {sat, code, phase} ->
            {sat, code + noise_at(noise_by_sat, sat, idx), phase}

          %{satellite_id: sat, code_m: code} = obs ->
            %{obs | code_m: code + noise_at(noise_by_sat, sat, idx)}
        end)

      %{epoch | rover_observations: rover}
    end)
  end

  defp noise_at(noise_by_sat, sat, idx) do
    case Map.fetch(noise_by_sat, sat) do
      {:ok, values} -> Enum.at(values, idx, 0.0)
      :error -> 0.0
    end
  end

  defp synthetic_baseline_epoch(base, baseline, satellite_positions_m, opts \\ []) do
    base_clock_m = Keyword.get(opts, :base_clock_m, 0.0)
    rover_clock_m = Keyword.get(opts, :rover_clock_m, 0.0)
    common_errors_m = Keyword.get(opts, :common_errors_m, %{})
    ambiguities_m = Keyword.get(opts, :ambiguities_m, %{})
    epoch = Keyword.get(opts, :epoch)
    rover = add3(base, baseline)

    {base_obs, rover_obs} =
      satellite_positions_m
      |> Enum.sort_by(fn {sat, _pos} -> sat end)
      |> Enum.map(fn {sat, sat_pos} ->
        common = Map.get(common_errors_m, sat, 0.0)
        base_code = norm(sub3(sat_pos, base)) + base_clock_m + common
        rover_code = norm(sub3(sat_pos, rover)) + rover_clock_m + common

        {{sat, base_code, base_code},
         {sat, rover_code, rover_code + Map.get(ambiguities_m, sat, 0.0)}}
      end)
      |> Enum.unzip()

    %{
      epoch: epoch,
      base_observations: base_obs,
      rover_observations: rover_obs,
      satellite_positions_m: satellite_positions_m
    }
  end

  defp synthetic_dual_baseline_epoch(base, baseline, satellite_positions_m, opts) do
    base_clock_m = Keyword.get(opts, :base_clock_m, 0.0)
    rover_clock_m = Keyword.get(opts, :rover_clock_m, 0.0)
    n1_cycles = Keyword.fetch!(opts, :n1_cycles)
    wide_lane_cycles = Keyword.fetch!(opts, :wide_lane_cycles)
    epoch_idx = Keyword.get(opts, :epoch, 0)
    epoch = Keyword.get(opts, :epoch)
    rover = add3(base, baseline)

    {base_obs, rover_obs} =
      satellite_positions_m
      |> Enum.sort_by(fn {sat, _pos} -> sat end)
      |> Enum.with_index()
      |> Enum.map(fn {{sat, sat_pos}, sat_idx} ->
        base_range = norm(sub3(sat_pos, base))
        rover_range = norm(sub3(sat_pos, rover))
        common = 0.4 + 0.03 * sat_idx
        iono1_m = 1.7 + 0.04 * epoch_idx + 0.02 * sat_idx
        iono2_m = iono1_m * :math.pow(@f_l1 / @f_l2, 2)

        n1 = Map.fetch!(n1_cycles, sat)
        n2 = n1 - Map.fetch!(wide_lane_cycles, sat)

        {
          dual_observation(sat, base_range, base_clock_m, common, iono1_m, iono2_m, 0, 0),
          dual_observation(
            sat,
            rover_range,
            rover_clock_m,
            common,
            iono1_m,
            iono2_m,
            n1,
            n2
          )
        }
      end)
      |> Enum.unzip()

    %{
      epoch: epoch,
      base_observations: base_obs,
      rover_observations: rover_obs,
      satellite_positions_m: satellite_positions_m
    }
  end

  defp dual_observation(sat, range_m, clock_m, common_m, iono1_m, iono2_m, n1, n2) do
    p1 = range_m + clock_m + common_m + iono1_m
    p2 = range_m + clock_m + common_m + iono2_m
    l1 = range_m + clock_m + common_m - iono1_m + n1 * @l1_wavelength_m
    l2 = range_m + clock_m + common_m - iono2_m + n2 * @l2_wavelength_m

    %{
      satellite_id: sat,
      p1_m: p1,
      p2_m: p2,
      phi1_cyc: l1 / @l1_wavelength_m,
      phi2_cyc: l2 / @l2_wavelength_m,
      f1_hz: @f_l1,
      f2_hz: @f_l2,
      lli1: 0,
      lli2: 0
    }
  end

  defp mark_dual_rover_lli(epoch, sat, lli) do
    %{epoch | rover_observations: mark_dual_observation_lli(epoch.rover_observations, sat, lli)}
  end

  defp mark_dual_observation_lli(observations, sat, lli) do
    Enum.map(observations, fn
      %{satellite_id: ^sat} = obs -> %{obs | lli1: lli}
      obs -> obs
    end)
  end

  defp expected_narrow_lane_offsets(wide_lane_cycles) do
    wide_lane_cycles
    |> Map.delete("G01")
    |> Map.new(fn {sat, wide_lane} -> {sat, narrow_lane_offset_m(wide_lane)} end)
  end

  defp narrow_lane_offset_m(wide_lane_cycles) do
    {:ok, gamma} = Orbis.GNSS.IonosphereFree.gamma(@f_l1, @f_l2)
    (gamma - 1.0) * @l2_wavelength_m * wide_lane_cycles
  end

  defp add3({ax, ay, az}, {bx, by, bz}), do: {ax + bx, ay + by, az + bz}
  defp sub3({ax, ay, az}, {bx, by, bz}), do: {ax - bx, ay - by, az - bz}
  defp norm({x, y, z}), do: :math.sqrt(x * x + y * y + z * z)
end
