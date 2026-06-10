defmodule Orbis.GNSS.RTKRTKLIBOracleTest do
  use ExUnit.Case, async: true

  @oracle_path Path.join(__DIR__, "fixtures/rtk/wtzr_wtzz_rtklib_oracle.json")
  @precise_oracle_path Path.join(__DIR__, "fixtures/rtk/wtzr_wtzz_rtklib_precise_oracle.json")

  test "WTZR/WTZZ RTKLIB oracle fixture pins the L1+broadcast reference target" do
    oracle =
      @oracle_path
      |> File.read!()
      |> Jason.decode!()

    assert oracle["version"] == 1

    reference = oracle["reference"]
    epochs = oracle["per_epoch"]

    assert reference["label"] == "l1_brdc_fix_and_hold"
    assert reference["epochs"] == 120
    assert reference["fixed_epochs"] == 119
    assert reference["first_fixed_index"] == 1
    assert reference["first_fixed_time"] == "2020-06-25T00:00:30"
    assert reference["final_status"] == "fixed"
    assert reference["final_truth_error_m"] < 0.004

    assert length(epochs) == 120
    assert hd(epochs)["fix_status"] == "float"
    assert Enum.at(epochs, 1)["fix_status"] == "fixed"
    assert Enum.count(epochs, &(&1["q"] == 1)) == reference["fixed_epochs"]

    final = List.last(epochs)
    assert final["baseline_enu_m"] == reference["final_baseline_enu_m"]
    assert final["ratio"] == reference["final_ratio"]

    modes = Map.new(oracle["comparison_modes"], &{&1["label"], &1})

    assert modes["l1_instantaneous"]["fixed_epochs"] == 93
    assert modes["l1_instantaneous"]["first_fixed_time"] == "2020-06-25T00:02:00"
    assert modes["l1_float"]["fixed_epochs"] == 0
    assert modes["l1_l2_brdc_fix_and_hold"]["fixed_epochs"] == 120
  end

  test "WTZR/WTZZ RTKLIB precise oracle fixture pins the lowercase-SP3 precise target" do
    oracle =
      @precise_oracle_path
      |> File.read!()
      |> Jason.decode!()

    assert oracle["version"] == 1

    reference = oracle["reference"]
    epochs = oracle["per_epoch"]

    assert reference["label"] == "l1_precise_cod_sp3_grg_clk_fix_and_hold"
    assert reference["epochs"] == 120
    assert reference["fixed_epochs"] == 119
    assert reference["first_fixed_index"] == 1
    assert reference["first_fixed_time"] == "2020-06-25T00:00:30"
    assert reference["final_status"] == "fixed"
    assert reference["final_truth_error_m"] < 0.004

    assert oracle["precise_products"]["orbit"] == "COD0MGXFIN_20201770000_01D_05M_ORB.SP3"
    assert oracle["precise_products"]["clock"] == "GRG0MGXFIN_20201770000_01D_30S_CLK.CLK"
    assert oracle["precise_products"]["staged_orbit_path"] == "cod.sp3"

    assert length(epochs) == 120
    assert hd(epochs)["fix_status"] == "float"
    assert Enum.at(epochs, 1)["fix_status"] == "fixed"
    assert Enum.count(epochs, &(&1["q"] == 1)) == reference["fixed_epochs"]

    final = List.last(epochs)
    assert final["baseline_enu_m"] == reference["final_baseline_enu_m"]
    assert final["ratio"] == reference["final_ratio"]

    comparison = oracle["broadcast_comparison"]
    assert comparison["source_fixture"] == "wtzr_wtzz_rtklib_oracle.json"
    assert comparison["same_fix_status_by_epoch"]
    assert comparison["max_baseline_delta_m"] < 0.002
  end
end
