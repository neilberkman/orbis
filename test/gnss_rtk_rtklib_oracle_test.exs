defmodule Orbis.GNSS.RTKRTKLIBOracleTest do
  use ExUnit.Case, async: true

  @oracle_path Path.join(__DIR__, "fixtures/rtk/wtzr_wtzz_rtklib_oracle.json")
  @precise_oracle_path Path.join(__DIR__, "fixtures/rtk/wtzr_wtzz_rtklib_precise_oracle.json")

  @cd_oracles [
    %{
      fixture: "pasa_scoa_2026_120_l1_static_fixhold_rtklib_oracle.json",
      config: "cd_pasa_scoa_l1_static_fixhold.conf",
      pos: "cd_pasa_scoa_l1_static_fixhold.pos",
      label: "cd_pasa_scoa_2026_120_l1_static_fixhold",
      description:
        "RTKLIB 2.4.2-p13 C+D Phase 1 precise-GPS oracle for PASA00ESP rover against SCOA00FRA base on 2026-04-30 10:00-12:00 GPST (L1 static, fix-and-hold, AR ratio gate 3.0, ANTEX receiver PCV and solid earth tides enabled).",
      epochs: 240,
      fixed_epochs: 171,
      q_counts: %{"1" => 171, "2" => 69},
      fix_rate: 0.7125,
      first_fixed_index: 2,
      first_fixed_time: "2026-04-30T10:01:00",
      final_status: "fixed",
      final_ratio: 999.9,
      final_truth_error_m: 0.052353514434,
      mean_truth_error_m: 0.107036863232,
      max_truth_error_m: 0.375208123609,
      first_time: "2026-04-30T10:00:00",
      last_time: "2026-04-30T11:59:30",
      satellites_min: 4,
      satellites_max: 5
    },
    %{
      fixture: "pasa_scoa_2026_120_l1l2_static_rtklib_oracle.json",
      config: "cd_pasa_scoa_l1l2_static.conf",
      pos: "cd_pasa_scoa_l1l2_static.pos",
      label: "cd_pasa_scoa_2026_120_l1l2_static",
      description:
        "RTKLIB 2.4.2-p13 C+D Phase 1 precise-GPS oracle for PASA00ESP rover against SCOA00FRA base on 2026-04-30 10:00-12:00 GPST (dual-frequency static, continuous AR, AR ratio gate 3.0, ANTEX receiver PCV and solid earth tides enabled).",
      epochs: 240,
      fixed_epochs: 80,
      q_counts: %{"1" => 80, "2" => 160},
      fix_rate: 0.333333333333,
      first_fixed_index: 2,
      first_fixed_time: "2026-04-30T10:01:00",
      final_status: "float",
      final_ratio: 1.5,
      final_truth_error_m: 0.058098081566,
      mean_truth_error_m: 0.208126085588,
      max_truth_error_m: 0.980812363794,
      first_time: "2026-04-30T10:00:00",
      last_time: "2026-04-30T11:59:30",
      satellites_min: 4,
      satellites_max: 5
    }
  ]

  @gsdc_oracles [
    %{
      fixture: "gsdc_2021_08_04_sjc1_pixel5_p222_demo5_rtklib_oracle.json",
      config: "track_a_gsdc_2021_08_04_sjc1_p222_grec_l1.conf",
      pos: "track_a_gsdc_2021_08_04_sjc1_p222_grec_l1.pos",
      label: "gsdc_2021_08_04_sjc1_pixel5_p222_grec_l1_demo5",
      description:
        "RTKLIB demo5 moving-rover oracle for GSDC 2022 train/2021-08-04-US-SJC-1/GooglePixel5 against NOAA CORS P222 (G/R/E/C L1, combined, fix-and-hold, AR ratio gate 3.0). The 10 fixed epochs do not beat float on 3D median error, so fixed status is not a confidence target; the oracle is an honest trajectory accuracy reference.",
      drive: "train/2021-08-04-US-SJC-1/GooglePixel5",
      base_doy: "216",
      nav_doy: "216",
      start_date: "2021/08/04",
      start_time: "20:40:43",
      end_date: "2021/08/04",
      end_time: "21:06:40",
      base_distance_km: 27.403,
      truth_time_tolerance_ms: 2,
      epochs: 1554,
      fixed_epochs: 10,
      q_counts: %{"1" => 10, "2" => 1538, "4" => 6},
      first_fixed_index: 85,
      first_fixed_time: "2021-08-04T20:42:08.449",
      first_fixed_truth_time_utc: "2021-08-04T20:41:50.449",
      final_status: "float",
      fix_rate: 0.006435006435,
      error_3d_median: 4.5221765,
      error_3d_p95: 12.296687,
      horizontal_p95: 6.688258,
      first_time: "2021-08-04T20:40:43.449",
      first_truth_time_utc: "2021-08-04T20:40:25.449",
      last_time: "2021-08-04T21:06:36.450",
      fixed_3d_median: 5.044745,
      fixed_3d_p95: 7.332008,
      float_3d_median: 4.516922,
      float_3d_p95: 12.296687,
      fixed_beats_float: false
    },
    %{
      fixture: "gsdc_2021_08_24_svl1_pixel5_p222_demo5_rtklib_oracle.json",
      config: "track_a_gsdc_p222_grec_l1.conf",
      pos: "track_a_gsdc_p222_grec_l1.pos",
      label: "gsdc_svl1_pixel5_p222_grec_l1_demo5",
      description:
        "RTKLIB demo5 moving-rover oracle for GSDC 2022 train/2021-08-24-US-SVL-1/GooglePixel5 against NOAA CORS P222 (G/R/E/C L1, combined, fix-and-hold, AR ratio gate 3.0). Validated fixes on this phone arc are meter-class, not cm-class; the oracle is an honest trajectory accuracy reference, not a fix-rate target.",
      drive: "train/2021-08-24-US-SVL-1/GooglePixel5",
      base_doy: "236",
      nav_doy: "236",
      start_date: "2021/08/24",
      start_time: "20:33:00",
      end_date: "2021/08/24",
      end_time: "21:25:20",
      base_distance_km: 18.936,
      truth_time_tolerance_ms: nil,
      epochs: 3136,
      fixed_epochs: 10,
      q_counts: %{"1" => 10, "2" => 3104, "4" => 22},
      first_fixed_index: 14,
      first_fixed_time: "2021-08-24T20:33:14.437",
      first_fixed_truth_time_utc: "2021-08-24T20:32:56.437",
      final_status: "float",
      fix_rate: 0.00318877551,
      error_3d_median: 3.9769565,
      error_3d_p95: 8.775371,
      horizontal_p95: 6.034325,
      first_time: "2021-08-24T20:33:00.437",
      first_truth_time_utc: "2021-08-24T20:32:42.437",
      last_time: "2021-08-24T21:25:16.437",
      fixed_3d_median: 3.476352,
      fixed_3d_p95: 4.562249,
      float_3d_median: 3.964981,
      float_3d_p95: 8.579962,
      fixed_beats_float: true
    },
    %{
      fixture: "gsdc_2021_12_15_mtv1_pixel5_p222_demo5_rtklib_oracle.json",
      config: "track_a_gsdc_2021_12_15_mtv1_p222_grec_l1.conf",
      pos: "track_a_gsdc_2021_12_15_mtv1_p222_grec_l1.pos",
      label: "gsdc_2021_12_15_mtv1_pixel5_p222_grec_l1_demo5",
      description:
        "RTKLIB demo5 moving-rover oracle for GSDC 2022 train/2021-12-15-US-MTV-1/GooglePixel5 against NOAA CORS P222 (G/R/E/C L1, combined, fix-and-hold, AR ratio gate 3.0). This highway phone arc has only one fixed epoch; the split is underpowered and meter-class, so the oracle is an honest trajectory accuracy reference, not a fix-rate target.",
      drive: "train/2021-12-15-US-MTV-1/GooglePixel5",
      base_doy: "349",
      nav_doy: "349",
      start_date: "2021/12/15",
      start_time: "18:49:11",
      end_date: "2021/12/15",
      end_time: "19:13:40",
      base_distance_km: 13.815,
      truth_time_tolerance_ms: 2,
      epochs: 1465,
      fixed_epochs: 1,
      q_counts: %{"1" => 1, "2" => 1436, "4" => 28},
      first_fixed_index: 1312,
      first_fixed_time: "2021-12-15T19:11:05.438",
      first_fixed_truth_time_utc: "2021-12-15T19:10:47.438",
      final_status: "float",
      fix_rate: 0.000682593857,
      error_3d_median: 3.652537,
      error_3d_p95: 7.909147,
      horizontal_p95: 4.668259,
      first_time: "2021-12-15T18:49:11.438",
      first_truth_time_utc: "2021-12-15T18:48:53.438",
      last_time: "2021-12-15T19:13:37.438",
      fixed_3d_median: 3.022934,
      fixed_3d_p95: 3.022934,
      float_3d_median: 3.633223,
      float_3d_p95: 7.466983,
      fixed_beats_float: true
    },
    %{
      fixture: "gsdc_2021_12_28_mtv1_pixel5_p222_demo5_rtklib_oracle.json",
      config: "track_a_gsdc_2021_12_28_mtv1_p222_grec_l1.conf",
      pos: "track_a_gsdc_2021_12_28_mtv1_p222_grec_l1.pos",
      label: "gsdc_2021_12_28_mtv1_pixel5_p222_grec_l1_demo5",
      description:
        "RTKLIB demo5 moving-rover oracle for GSDC 2022 train/2021-12-28-US-MTV-1/GooglePixel5 against NOAA CORS P222 (G/R/E/C L1, combined, fix-and-hold, AR ratio gate 3.0). The fixed split beats float on this repeat highway route, but remains meter-class; the oracle is an honest trajectory accuracy reference, not a fix-rate target.",
      drive: "train/2021-12-28-US-MTV-1/GooglePixel5",
      base_doy: "362",
      nav_doy: "362",
      start_date: "2021/12/28",
      start_time: "20:17:25",
      end_date: "2021/12/28",
      end_time: "20:44:20",
      base_distance_km: 13.702,
      truth_time_tolerance_ms: 2,
      epochs: 1610,
      fixed_epochs: 10,
      q_counts: %{"1" => 10, "2" => 1567, "4" => 33},
      first_fixed_index: 830,
      first_fixed_time: "2021-12-28T20:31:18.438",
      first_fixed_truth_time_utc: "2021-12-28T20:31:00.437",
      final_status: "float",
      fix_rate: 0.006211180124,
      error_3d_median: 3.973879,
      error_3d_p95: 9.033375,
      horizontal_p95: 6.674816,
      first_time: "2021-12-28T20:17:25.438",
      first_truth_time_utc: "2021-12-28T20:17:07.437",
      last_time: "2021-12-28T20:44:17.438",
      fixed_3d_median: 2.642458,
      fixed_3d_p95: 3.382031,
      float_3d_median: 3.971948,
      float_3d_p95: 8.698615,
      fixed_beats_float: true
    }
  ]

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

  test "GSDC Pixel5/P222 demo5 oracle fixtures pin the moving-rover references" do
    for arc <- @gsdc_oracles do
      oracle =
        __DIR__
        |> Path.join("fixtures/rtk/#{arc.fixture}")
        |> File.read!()
        |> Jason.decode!()

      assert oracle["version"] == 1
      assert oracle["description"] == arc.description

      assert oracle["generator"]["rtklib"] == %{
               "program" => "rnx2rtkp",
               "version" => "EX 2.5.0",
               "commit" => "57d39e7"
             }

      assert oracle["truth"]["base_station"]["id"] == "P222"
      assert oracle["truth"]["gps_utc_offset_s"] == 18

      assert oracle["truth"]["base_station"]["distance_from_drive_start_km"] ==
               arc.base_distance_km

      if arc.truth_time_tolerance_ms do
        assert oracle["truth"]["time_match_tolerance_ms"] == arc.truth_time_tolerance_ms
      else
        refute Map.has_key?(oracle["truth"], "time_match_tolerance_ms")
      end

      reference = oracle["reference"]
      epochs = oracle["per_epoch"]

      assert reference["label"] == arc.label
      assert reference["config"] == arc.config
      assert reference["source_pos"] == arc.pos
      assert reference["epochs"] == arc.epochs
      assert reference["fixed_epochs"] == arc.fixed_epochs
      assert reference["q_counts"] == arc.q_counts
      assert reference["first_fixed_index"] == arc.first_fixed_index
      assert reference["first_fixed_time"] == arc.first_fixed_time
      assert reference["first_fixed_truth_time_utc"] == arc.first_fixed_truth_time_utc
      assert reference["final_status"] == arc.final_status

      assert_in_delta reference["fix_rate"], arc.fix_rate, 1.0e-12
      assert_in_delta reference["error_3d"]["median_m"], arc.error_3d_median, 1.0e-6
      assert_in_delta reference["error_3d"]["p95_m"], arc.error_3d_p95, 1.0e-6
      assert_in_delta reference["horizontal_error"]["p95_m"], arc.horizontal_p95, 1.0e-6

      assert length(epochs) == reference["epochs"]
      assert hd(epochs)["time"] == arc.first_time
      assert hd(epochs)["truth_time_utc"] == arc.first_truth_time_utc
      assert List.last(epochs)["time"] == arc.last_time
      assert Enum.count(epochs, &(&1["q"] == 1)) == reference["fixed_epochs"]

      fixed = Enum.at(epochs, reference["first_fixed_index"])
      assert fixed["fix_status"] == "fixed"
      assert fixed["satellites"] >= 4
      assert fixed["ratio"] >= 3.0

      fixed_3d = status_error_values(epochs, 1, "error_3d_m")
      float_3d = status_error_values(epochs, 2, "error_3d_m")

      assert_in_delta median(fixed_3d), arc.fixed_3d_median, 1.0e-6
      assert_in_delta percentile(fixed_3d, 0.95), arc.fixed_3d_p95, 1.0e-6
      assert_in_delta median(float_3d), arc.float_3d_median, 1.0e-6
      assert_in_delta percentile(float_3d, 0.95), arc.float_3d_p95, 1.0e-6

      assert fixed_split_beats_float?(fixed_3d, float_3d) == arc.fixed_beats_float
    end
  end

  test "PASA/SCOA C+D Phase 1 RTKLIB oracle fixtures pin the EPN references" do
    for arc <- @cd_oracles do
      oracle =
        __DIR__
        |> Path.join("fixtures/rtk/#{arc.fixture}")
        |> File.read!()
        |> Jason.decode!()

      assert oracle["version"] == "1"
      assert oracle["description"] == arc.description

      assert oracle["generator"]["rtklib"] == %{
               "program" => "rnx2rtkp",
               "version" => "v2.4.2-p13",
               "commit" => "71db0ff"
             }

      assert oracle["truth"]["frame"] ==
               "ITRF2020 ECEF metres propagated to 2026-04-30T11:00:00 GPST; ENU baseline at SCOA00FRA ARP, metres"

      assert oracle["truth"]["base_station"]["id"] == "SCOA00FRA"
      assert oracle["truth"]["rover_station"]["id"] == "PASA00ESP"
      assert oracle["truth"]["baseline_length_km"] == 21.836327792

      truth = oracle["truth"]["antenna_baseline_enu_m"]
      assert_in_delta truth["east"], -20265.52060276001, 1.0e-12
      assert_in_delta truth["north"], -8132.221127240546, 1.0e-12
      assert_in_delta truth["up"], -29.422321279917, 1.0e-12

      assert oracle["inputs"] == %{
               "rover_obs" => "test/fixtures/obs/PASA00ESP_R_20261201000_02H_30S_MO.rnx",
               "base_obs" => "test/fixtures/obs/SCOA00FRA_R_20261201000_02H_30S_MO.rnx",
               "nav" => "test/fixtures/nav/BRDC00WRD_R_20261200800_06H_MN.rnx",
               "sp3" => "test/fixtures/sp3/IGS0OPSFIN_20261200945_02H30M_15M_ORB.SP3",
               "clk" => "test/fixtures/clk/IGS0OPSFIN_2026120095930_02H01M_30S_CLK.CLK",
               "antex" => "test/fixtures/antex/igs20_pasa_scoa_gps.atx"
             }

      reference = oracle["reference"]
      epochs = oracle["per_epoch"]

      assert reference["label"] == arc.label
      assert reference["config"] == arc.config
      assert reference["source_pos"] == arc.pos
      assert reference["epochs"] == arc.epochs
      assert reference["fixed_epochs"] == arc.fixed_epochs
      assert reference["q_counts"] == arc.q_counts
      assert reference["first_fixed_index"] == arc.first_fixed_index
      assert reference["first_fixed_time"] == arc.first_fixed_time
      assert reference["final_status"] == arc.final_status
      assert reference["satellites_min"] == arc.satellites_min
      assert reference["satellites_max"] == arc.satellites_max

      assert_in_delta reference["fix_rate"], arc.fix_rate, 1.0e-12
      assert_in_delta reference["final_ratio"], arc.final_ratio, 1.0e-12
      assert_in_delta reference["final_truth_error_m"], arc.final_truth_error_m, 1.0e-12
      assert_in_delta reference["mean_truth_error_m"], arc.mean_truth_error_m, 1.0e-12
      assert_in_delta reference["max_truth_error_m"], arc.max_truth_error_m, 1.0e-12

      assert length(epochs) == arc.epochs
      assert hd(epochs)["time"] == arc.first_time
      assert List.last(epochs)["time"] == arc.last_time
      assert Enum.count(epochs, &(&1["q"] == 1)) == arc.fixed_epochs

      first_fixed = Enum.at(epochs, arc.first_fixed_index)
      assert first_fixed["fix_status"] == "fixed"
      assert first_fixed["ratio"] >= 3.0
    end
  end

  @tag :local_data
  test "GSDC Pixel5/P222 demo5 oracles regenerate byte-identically from local inputs" do
    repo = Path.expand("..", __DIR__)
    generator_dir = Path.join(repo, "test/fixtures/rtk/generators")
    script = Path.join(generator_dir, "pos_to_oracle.py")

    rnx2rtkp = "/tmp/RTKLIB-demo5/app/consapp/rnx2rtkp/gcc/rnx2rtkp"
    work = "/tmp/gsdc-work"

    for arc <- @gsdc_oracles do
      conf = Path.join(generator_dir, arc.config)
      drive = Path.join(work, arc.drive)
      rover = Path.join(drive, "supplemental/gnss_rinex.21o")
      truth = Path.join(drive, "ground_truth.csv")
      base = Path.join(work, "cors/p222#{arc.base_doy}0.21o")
      nav = Path.join(work, "cors/BRDC00WRD_R_2021#{arc.nav_doy}0000_01D_MN.rnx")

      required = [rnx2rtkp, rover, truth, base, nav]

      if Enum.all?(required, &File.exists?/1) do
        tmp = Path.join(System.tmp_dir!(), "orbis-gsdc-oracle-test/#{arc.label}")
        File.rm_rf!(tmp)
        File.mkdir_p!(tmp)

        pos = Path.join(tmp, arc.pos)
        regenerated = Path.join(tmp, arc.fixture)

        {_, 0} =
          System.cmd(
            rnx2rtkp,
            [
              "-k",
              conf,
              "-ts",
              arc.start_date,
              arc.start_time,
              "-te",
              arc.end_date,
              arc.end_time,
              "-o",
              pos,
              rover,
              base,
              nav
            ],
            stderr_to_stdout: true
          )

        args =
          [
            script,
            pos,
            arc.config,
            arc.label,
            arc.description,
            regenerated,
            "--moving-truth-csv",
            truth,
            "--truth-source",
            "#{arc.drive}/ground_truth.csv",
            "--drive",
            arc.drive,
            "--rover-source",
            "#{arc.drive}/supplemental/gnss_rinex.21o",
            "--base-source",
            "https://geodesy.noaa.gov/corsdata/rinex/2021/#{arc.base_doy}/p222/p222#{arc.base_doy}0.21d.gz",
            "--nav-source",
            "https://igs.bkg.bund.de/root_ftp/IGS/BRDC/2021/#{arc.nav_doy}/BRDC00WRD_R_2021#{arc.nav_doy}0000_01D_MN.rnx.gz",
            "--base-station",
            "P222",
            "--base-ecef-m=-2689639.5060,-4290438.6360,3865050.9560",
            "--base-distance-km",
            "#{arc.base_distance_km}",
            "--rtklib-version",
            "EX 2.5.0",
            "--rtklib-commit",
            "57d39e7"
          ] ++ truth_time_tolerance_args(arc)

        {_, 0} = System.cmd("python3", args, stderr_to_stdout: true)

        expected = Path.join(__DIR__, "fixtures/rtk/#{arc.fixture}")
        assert File.read!(regenerated) == File.read!(expected)
      else
        IO.puts(
          "Skipping #{arc.fixture} local-data regeneration; missing #{Enum.reject(required, &File.exists?/1) |> Enum.join(", ")}"
        )
      end
    end
  end

  @tag :local_data
  test "PASA/SCOA C+D Phase 1 RTKLIB oracles regenerate byte-identically" do
    repo = Path.expand("..", __DIR__)
    generator_dir = Path.join(repo, "test/fixtures/rtk/generators")
    script = Path.join(generator_dir, "pos_to_oracle.py")
    truth = Path.join(generator_dir, "cd_pasa_scoa_2026_120_truth.json")

    rnx2rtkp =
      System.get_env("RTKLIB_RNX2RTKP") ||
        "/tmp/cd-phase1-tools/RTKLIB/app/rnx2rtkp/gcc/rnx2rtkp"

    rover = Path.join(repo, "test/fixtures/obs/PASA00ESP_R_20261201000_02H_30S_MO.rnx")
    base = Path.join(repo, "test/fixtures/obs/SCOA00FRA_R_20261201000_02H_30S_MO.rnx")
    nav = Path.join(repo, "test/fixtures/nav/BRDC00WRD_R_20261200800_06H_MN.rnx")
    sp3 = Path.join(repo, "test/fixtures/sp3/IGS0OPSFIN_20261200945_02H30M_15M_ORB.SP3")
    clk = Path.join(repo, "test/fixtures/clk/IGS0OPSFIN_2026120095930_02H01M_30S_CLK.CLK")

    required = [rnx2rtkp, rover, base, nav, sp3, clk]

    if Enum.all?(required, &File.exists?/1) do
      tmp = Path.join(System.tmp_dir!(), "orbis-cd-pasa-scoa-oracle-test")
      File.rm_rf!(tmp)
      File.mkdir_p!(tmp)

      staged_sp3 = Path.join(tmp, "igs_fin.sp3")
      staged_clk = Path.join(tmp, "igs_fin.clk")
      File.cp!(sp3, staged_sp3)
      File.cp!(clk, staged_clk)

      for arc <- @cd_oracles do
        pos = Path.join(tmp, arc.pos)
        regenerated = Path.join(tmp, arc.fixture)

        {_, 0} =
          System.cmd(
            rnx2rtkp,
            ["-k", arc.config, "-o", pos, rover, base, nav, staged_sp3, staged_clk],
            cd: generator_dir,
            stderr_to_stdout: true
          )

        {_, 0} =
          System.cmd(
            "python3",
            [
              script,
              pos,
              arc.config,
              arc.label,
              arc.description,
              regenerated,
              "--static-truth-json",
              truth,
              "--rover-source",
              "test/fixtures/obs/PASA00ESP_R_20261201000_02H_30S_MO.rnx",
              "--base-source",
              "test/fixtures/obs/SCOA00FRA_R_20261201000_02H_30S_MO.rnx",
              "--nav-source",
              "test/fixtures/nav/BRDC00WRD_R_20261200800_06H_MN.rnx",
              "--sp3-source",
              "test/fixtures/sp3/IGS0OPSFIN_20261200945_02H30M_15M_ORB.SP3",
              "--clk-source",
              "test/fixtures/clk/IGS0OPSFIN_2026120095930_02H01M_30S_CLK.CLK",
              "--antex-source",
              "test/fixtures/antex/igs20_pasa_scoa_gps.atx",
              "--rtklib-version",
              "v2.4.2-p13",
              "--rtklib-commit",
              "71db0ff"
            ],
            stderr_to_stdout: true
          )

        expected = Path.join(__DIR__, "fixtures/rtk/#{arc.fixture}")
        assert File.read!(regenerated) == File.read!(expected)
      end
    else
      IO.puts(
        "Skipping PASA/SCOA local-data regeneration; missing #{Enum.reject(required, &File.exists?/1) |> Enum.join(", ")}"
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

  defp fixed_split_beats_float?(fixed_3d, float_3d) do
    median(fixed_3d) < median(float_3d) and
      percentile(fixed_3d, 0.95) < percentile(float_3d, 0.95)
  end

  defp truth_time_tolerance_args(%{truth_time_tolerance_ms: nil}), do: []

  defp truth_time_tolerance_args(%{truth_time_tolerance_ms: tolerance_ms}) do
    ["--truth-time-tolerance-ms", "#{tolerance_ms}"]
  end
end
