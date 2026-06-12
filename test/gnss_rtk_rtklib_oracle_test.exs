defmodule Orbis.GNSS.RTKRTKLIBOracleTest do
  use ExUnit.Case, async: true

  @oracle_path Path.join(__DIR__, "fixtures/rtk/wtzr_wtzz_rtklib_oracle.json")
  @precise_oracle_path Path.join(__DIR__, "fixtures/rtk/wtzr_wtzz_rtklib_precise_oracle.json")
  @gsdc_oracle_path Path.join(
                      __DIR__,
                      "fixtures/rtk/gsdc_2021_08_24_svl1_pixel5_p222_demo5_rtklib_oracle.json"
                    )
  @gsdc_description "RTKLIB demo5 moving-rover oracle for GSDC 2022 train/2021-08-24-US-SVL-1/GooglePixel5 against NOAA CORS P222 (G/R/E/C L1, combined, fix-and-hold, AR ratio gate 3.0). Validated fixes on this phone arc are meter-class, not cm-class; the oracle is an honest trajectory accuracy reference, not a fix-rate target."

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

  test "GSDC Pixel5/P222 demo5 oracle fixture pins the moving-rover reference" do
    oracle =
      @gsdc_oracle_path
      |> File.read!()
      |> Jason.decode!()

    assert oracle["version"] == 1
    assert oracle["description"] == @gsdc_description

    assert oracle["generator"]["rtklib"] == %{
             "program" => "rnx2rtkp",
             "version" => "EX 2.5.0",
             "commit" => "57d39e7"
           }

    assert oracle["truth"]["base_station"]["id"] == "P222"
    assert oracle["truth"]["gps_utc_offset_s"] == 18

    reference = oracle["reference"]
    epochs = oracle["per_epoch"]

    assert reference["label"] == "gsdc_svl1_pixel5_p222_grec_l1_demo5"
    assert reference["config"] == "track_a_gsdc_p222_grec_l1.conf"
    assert reference["epochs"] == 3136
    assert reference["fixed_epochs"] == 10
    assert reference["q_counts"] == %{"1" => 10, "2" => 3104, "4" => 22}
    assert reference["first_fixed_index"] == 14
    assert reference["first_fixed_time"] == "2021-08-24T20:33:14.437"
    assert reference["first_fixed_truth_time_utc"] == "2021-08-24T20:32:56.437"
    assert reference["final_status"] == "float"

    assert_in_delta reference["fix_rate"], 0.00318877551, 1.0e-12
    assert_in_delta reference["error_3d"]["median_m"], 3.9769565, 1.0e-6
    assert_in_delta reference["error_3d"]["p95_m"], 8.775371, 1.0e-6
    assert_in_delta reference["horizontal_error"]["p95_m"], 6.034325, 1.0e-6

    assert length(epochs) == reference["epochs"]
    assert hd(epochs)["time"] == "2021-08-24T20:33:00.437"
    assert hd(epochs)["truth_time_utc"] == "2021-08-24T20:32:42.437"
    assert List.last(epochs)["time"] == "2021-08-24T21:25:16.437"
    assert Enum.count(epochs, &(&1["q"] == 1)) == reference["fixed_epochs"]

    fixed = Enum.at(epochs, reference["first_fixed_index"])
    assert fixed["fix_status"] == "fixed"
    assert fixed["satellites"] >= 4
    assert fixed["ratio"] >= 3.0

    fixed_3d = status_error_values(epochs, 1, "error_3d_m")
    float_3d = status_error_values(epochs, 2, "error_3d_m")

    assert_in_delta median(fixed_3d), 3.476352, 1.0e-6
    assert_in_delta percentile(fixed_3d, 0.95), 4.562249, 1.0e-6
    assert_in_delta median(float_3d), 3.964981, 1.0e-6
    assert_in_delta percentile(float_3d, 0.95), 8.579962, 1.0e-6
  end

  @tag :local_data
  test "GSDC Pixel5/P222 demo5 oracle regenerates byte-identically from local inputs" do
    repo = Path.expand("..", __DIR__)
    generator_dir = Path.join(repo, "test/fixtures/rtk/generators")
    conf = Path.join(generator_dir, "track_a_gsdc_p222_grec_l1.conf")
    script = Path.join(generator_dir, "pos_to_oracle.py")

    rnx2rtkp = "/tmp/RTKLIB-demo5/app/consapp/rnx2rtkp/gcc/rnx2rtkp"
    work = "/tmp/gsdc-work"
    drive = Path.join(work, "train/2021-08-24-US-SVL-1/GooglePixel5")
    rover = Path.join(drive, "supplemental/gnss_rinex.21o")
    truth = Path.join(drive, "ground_truth.csv")
    base = Path.join(work, "cors/p2222360.21o")
    nav = Path.join(work, "cors/BRDC00WRD_R_20212360000_01D_MN.rnx")

    required = [rnx2rtkp, rover, truth, base, nav]

    if Enum.all?(required, &File.exists?/1) do
      tmp = Path.join(System.tmp_dir!(), "orbis-gsdc-oracle-test")
      File.rm_rf!(tmp)
      File.mkdir_p!(tmp)

      pos = Path.join(tmp, "track_a_gsdc_p222_grec_l1.pos")
      regenerated = Path.join(tmp, "gsdc_2021_08_24_svl1_pixel5_p222_demo5_rtklib_oracle.json")

      {_, 0} =
        System.cmd(
          rnx2rtkp,
          [
            "-k",
            conf,
            "-ts",
            "2021/08/24",
            "20:33:00",
            "-te",
            "2021/08/24",
            "21:25:20",
            "-o",
            pos,
            rover,
            base,
            nav
          ],
          stderr_to_stdout: true
        )

      {_, 0} =
        System.cmd(
          "python3",
          [
            script,
            pos,
            "track_a_gsdc_p222_grec_l1.conf",
            "gsdc_svl1_pixel5_p222_grec_l1_demo5",
            @gsdc_description,
            regenerated,
            "--moving-truth-csv",
            truth,
            "--truth-source",
            "train/2021-08-24-US-SVL-1/GooglePixel5/ground_truth.csv",
            "--drive",
            "train/2021-08-24-US-SVL-1/GooglePixel5",
            "--rover-source",
            "train/2021-08-24-US-SVL-1/GooglePixel5/supplemental/gnss_rinex.21o",
            "--base-source",
            "https://geodesy.noaa.gov/corsdata/rinex/2021/236/p222/p2222360.21d.gz",
            "--nav-source",
            "https://igs.bkg.bund.de/root_ftp/IGS/BRDC/2021/236/BRDC00WRD_R_20212360000_01D_MN.rnx.gz",
            "--base-station",
            "P222",
            "--base-ecef-m=-2689639.5060,-4290438.6360,3865050.9560",
            "--base-distance-km",
            "18.936",
            "--rtklib-version",
            "EX 2.5.0",
            "--rtklib-commit",
            "57d39e7"
          ],
          stderr_to_stdout: true
        )

      assert File.read!(regenerated) == File.read!(@gsdc_oracle_path)
    else
      IO.puts(
        "Skipping GSDC local-data regeneration; missing #{Enum.reject(required, &File.exists?/1) |> Enum.join(", ")}"
      )
    end
  end

  defp status_error_values(epochs, q, key) do
    epochs
    |> Enum.filter(&(&1["q"] == q))
    |> Enum.map(& &1[key])
  end

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

  defp percentile(values, pct) do
    ordered = Enum.sort(values)
    Enum.at(ordered, trunc(pct * (length(ordered) - 1)))
  end
end
