defmodule Orbis.GNSS.EphemerisTest do
  @moduledoc """
  Unified satellite-ephemeris sampling (`Orbis.GNSS.Ephemeris`) and the
  broadcast-vs-precise accuracy check (`Orbis.GNSS.BroadcastComparison`),
  exercised on the real 2020 DOY177 ESBC broadcast nav + GBM precise SP3
  fixtures. GPS LNAV orbits should agree with the precise product at the
  1-2 m level (IS-GPS-200 broadcast accuracy); a wild value flags a parse/eval
  regression.
  """
  use ExUnit.Case, async: true

  alias Orbis.GNSS.Broadcast
  alias Orbis.GNSS.BroadcastComparison
  alias Orbis.GNSS.Ephemeris
  alias Orbis.GNSS.SP3

  @nav_path Path.join(__DIR__, "fixtures/nav/ESBC00DNK_R_20201770000_01D_MN.rnx")
  @sp3_path Path.join(__DIR__, "fixtures/sp3/GBM0MGXRAP_20201770000_01D_05M_ORB_73epoch.sp3")

  setup_all do
    {:ok, broadcast: Broadcast.load!(@nav_path), sp3: SP3.load!(@sp3_path)}
  end

  describe "Ephemeris.sample/3" do
    test "samples a precise source into a per-satellite, per-epoch table", %{sp3: sp3} do
      window = %{from: ~N[2020-06-25 00:30:00], to: ~N[2020-06-25 01:00:00], step_s: 900}

      rows = Ephemeris.sample(sp3, ["G05", "G07"], window)

      # 2 satellites x 3 epochs (00:30, 00:45, 01:00), in sat order then epoch.
      assert length(rows) == 6
      assert Enum.map(rows, & &1.satellite_id) == ~w(G05 G05 G05 G07 G07 G07)

      assert Enum.map(rows, & &1.epoch) |> Enum.take(3) == [
               ~N[2020-06-25 00:30:00],
               ~N[2020-06-25 00:45:00],
               ~N[2020-06-25 01:00:00]
             ]

      row = hd(rows)
      assert row.status == :ok
      assert is_float(row.x_m) and is_float(row.y_m) and is_float(row.z_m)
      # GPS satellites sit at ~26 600 km geocentric radius.
      radius = :math.sqrt(row.x_m ** 2 + row.y_m ** 2 + row.z_m ** 2)
      assert radius > 25_000_000.0 and radius < 28_000_000.0
    end

    test "presents the identical surface for a broadcast source", %{broadcast: broadcast} do
      window = %{from: ~N[2020-06-25 01:00:00], to: ~N[2020-06-25 01:00:00], step_s: 900}

      [row] = Ephemeris.sample(broadcast, ["G07"], window)

      assert %Ephemeris.Row{satellite_id: "G07", status: :ok} = row
      assert is_float(row.x_m) and is_float(row.clock_s)
      radius = :math.sqrt(row.x_m ** 2 + row.y_m ** 2 + row.z_m ** 2)
      assert radius > 25_000_000.0 and radius < 28_000_000.0
    end

    test "reports an explicit gap instead of extrapolating", %{sp3: sp3} do
      # A satellite id that is not in the product, and an epoch outside coverage.
      window = %{from: ~N[2020-06-25 00:30:00], to: ~N[2020-06-25 00:30:00], step_s: 900}
      [missing_sat] = Ephemeris.sample(sp3, ["G99"], window)

      assert missing_sat.status == :no_ephemeris
      assert missing_sat.x_m == nil and missing_sat.y_m == nil
      assert missing_sat.z_m == nil and missing_sat.clock_s == nil

      out_of_window = %{from: ~N[2020-06-26 12:00:00], to: ~N[2020-06-26 12:00:00], step_s: 900}
      [out] = Ephemeris.sample(sp3, ["G05"], out_of_window)
      assert out.status == :no_ephemeris
      assert out.x_m == nil
    end

    test "rejects a non-positive step and an inverted window", %{sp3: sp3} do
      assert_raise ArgumentError, fn ->
        Ephemeris.sample(sp3, ["G05"], %{
          from: ~N[2020-06-25 00:00:00],
          to: ~N[2020-06-25 01:00:00],
          step_s: 0
        })
      end

      assert_raise ArgumentError, fn ->
        Ephemeris.sample(sp3, ["G05"], %{
          from: ~N[2020-06-25 01:00:00],
          to: ~N[2020-06-25 00:00:00],
          step_s: 900
        })
      end
    end
  end

  describe "BroadcastComparison.compare/4 on real GPS broadcast vs precise" do
    setup %{broadcast: broadcast, sp3: sp3} do
      gps =
        sp3 |> SP3.satellite_ids() |> Enum.filter(&String.starts_with?(&1, "G")) |> Enum.sort()

      window = %{from: ~N[2020-06-25 00:15:00], to: ~N[2020-06-25 05:45:00], step_s: 900}
      {:ok, report: BroadcastComparison.compare(broadcast, sp3, gps, window)}
    end

    test "GPS LNAV orbit agreement lands in the broadcast accuracy band", %{report: report} do
      overall = report.overall
      assert overall.count > 100

      # Expected GPS broadcast-vs-precise orbit RMS is ~1-2 m. The lower bound makes
      # the assertion non-tautological: a zeroed or broken eval would collapse to ~0.
      assert overall.orbit_3d_rms_m > 0.3
      assert overall.orbit_3d_rms_m < 3.0
      assert overall.orbit_3d_max_m < 6.0
    end

    test "the radial/along/cross decomposition is orthonormal and consistent", %{report: report} do
      overall = report.overall

      assert overall.radial_rms_m > 0.0
      assert overall.along_rms_m > 0.0
      assert overall.cross_rms_m > 0.0

      # RAC is an orthonormal rotation of the 3D difference, so the 3D RMS must equal
      # the quadrature sum of the component RMS values.
      quadrature =
        :math.sqrt(
          overall.radial_rms_m ** 2 + overall.along_rms_m ** 2 + overall.cross_rms_m ** 2
        )

      assert_in_delta overall.orbit_3d_rms_m, quadrature, 1.0e-6
    end

    test "clock differences are finite and reported per satellite", %{report: report} do
      overall = report.overall

      # The clock difference is the raw broadcast-minus-precise offset; it carries a
      # per-satellite datum/bias term (L1/TGD-referenced broadcast vs IF-referenced
      # precise), so it is larger than the orbit term but must stay bounded.
      assert is_float(overall.clock_rms_m)
      assert overall.clock_rms_m > 0.0 and overall.clock_rms_m < 50.0
    end

    test "per-satellite stats are populated and gaps are reported explicitly", %{report: report} do
      assert map_size(report.per_satellite) > 4

      {_sat, stats} = Enum.find(report.per_satellite, fn {_sat, s} -> s.count > 0 end)
      assert stats.orbit_3d_rms_m > 0.0

      # Some GPS satellites only broadcast valid ephemeris for part of the window, so
      # the comparison reports the skipped cells rather than silently dropping them.
      assert report.missing != []
    end
  end
end
