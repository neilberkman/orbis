defmodule Orbis.GNSS.RTKTest do
  use ExUnit.Case, async: true

  alias Orbis.GNSS.RTK

  @base {1_110_000.0, -4_840_000.0, 3_980_000.0}
  @truth_baseline {12.5, -4.25, 2.75}
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
               %{satellite_id: "G02", reference_satellite_id: "G01", code_m: 15.0, phase_m: 15.0}
             ]
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
               phase_sigma_m: 0.02
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

    test "selects the highest-elevation default reference and reports unusable satellites" do
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
      assert sol.used_sats == ["G01", "G02", "G04"]
      assert sol.metadata.dropped_sats == ["G05"]
      assert position_error(sol.baseline_m, @truth_baseline) < 1.0e-4
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

  defp add3({ax, ay, az}, {bx, by, bz}), do: {ax + bx, ay + by, az + bz}
  defp sub3({ax, ay, az}, {bx, by, bz}), do: {ax - bx, ay - by, az - bz}
  defp norm({x, y, z}), do: :math.sqrt(x * x + y * y + z * z)
end
