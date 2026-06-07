defmodule Orbis.GNSS.RTKTest do
  use ExUnit.Case, async: true

  alias Orbis.GNSS.RTK

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

  defp synth_observations(sats, ranges, clock_m, errors, phase_ambiguities) do
    Enum.map(sats, fn sat ->
      code = Map.fetch!(ranges, sat) + clock_m + Map.fetch!(errors, sat)
      phase = code + Map.fetch!(phase_ambiguities, sat)
      %{satellite_id: sat, code_m: code, phase_m: phase}
    end)
  end
end
