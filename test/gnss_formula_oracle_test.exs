defmodule Orbis.GNSS.FormulaOracleTest do
  use ExUnit.Case, async: true

  import Orbis.TestHelpers, only: [assert_ulp: 4, hex_to_float: 1]

  alias Orbis.GNSS.IonosphereFree
  alias Orbis.GNSS.Observables
  alias Orbis.GNSS.QC
  alias Orbis.GNSS.SP3

  @golden_path Path.join(__DIR__, "fixtures/orbis_gnss_formula_golden.json")
  @sp3_path Path.join(__DIR__, "fixtures/sp3/GRG0MGXFIN_20201760000_01D_15M_ORB.SP3")

  setup_all do
    {:ok, golden: @golden_path |> File.read!() |> Jason.decode!()}
  end

  describe "IonosphereFree vs Python oracle" do
    test "gamma, noise amplification, and combined ranges match the fixture", %{golden: golden} do
      for c <- golden["ionosphere_free"]["cases"] do
        f1 = h(c["f1_hz"])
        f2 = h(c["f2_hz"])
        pr1 = h(c["pr1_m"])
        pr2 = h(c["pr2_m"])

        assert {:ok, gamma} = IonosphereFree.gamma(f1, f2)
        assert_ulp(gamma, h(c["gamma"]), 0, "#{c["system"]} gamma")

        assert {:ok, amp} = IonosphereFree.noise_amplification(f1, f2)
        assert_ulp(amp, h(c["noise_amplification"]), 0, "#{c["system"]} noise amp")

        assert {:ok, pr_if} = IonosphereFree.iono_free(pr1, pr2, f1, f2)
        assert_ulp(pr_if, h(c["iono_free_m"]), 0, "#{c["system"]} iono-free")
      end
    end
  end

  describe "QC vs Python/scipy oracle" do
    test "elevation and C/N0 weighting matches the Python fixture", %{golden: golden} do
      qc = golden["qc"]

      for c <- qc["variance_cases"] do
        el = h(c["elevation_deg"])
        variance = QC.pseudorange_variance(el)
        assert_ulp(variance, h(c["variance_m2"]), 0, "variance #{el}")

        sigmas = QC.sigmas([{"G01", el}])
        weights = QC.weight_vector([{"G01", el}])
        assert_ulp(sigmas["G01"], h(c["sigma_m"]), 0, "sigma #{el}")
        assert_ulp(weights["G01"], h(c["weight"]), 0, "weight #{el}")
      end

      for c <- qc["cn0_cases"] do
        got =
          QC.pseudorange_variance(h(c["elevation_deg"]),
            model: :elevation_cn0,
            cn0: h(c["cn0_dbhz"])
          )

        assert_ulp(got, h(c["variance_m2"]), 0, "CN0 #{c["cn0_dbhz"]}")
      end
    end

    test "chi2_inv/2 matches scipy.stats.chi2.ppf", %{golden: golden} do
      for c <- golden["qc"]["chi2_cases"] do
        got = QC.chi2_inv(h(c["p"]), c["dof"])
        # Different implementation (regularized-gamma bisection vs scipy), so
        # assert numerical agreement to sub-nanometre-scale threshold precision,
        # not identical libm operation order.
        assert_in_delta got, h(c["critical"]), 1.0e-10
      end
    end
  end

  describe "Observables vs Python/scipy/pyproj oracle" do
    test "predict/5 matches the independent SP3 observable fixture", %{golden: golden} do
      sp3 = SP3.load!(@sp3_path)

      for c <- golden["observables"]["cases"] do
        opts = [light_time: c["light_time"], sagnac: c["sagnac"]]

        assert {:ok, obs} =
                 Observables.predict(sp3, c["sat"], receiver(c), ~N[2020-06-24 12:00:00], opts)

        assert_in_delta obs.geometric_range_m, h(c["geometric_range_m"]), 1.0e-4
        assert_in_delta obs.range_rate_m_s, h(c["range_rate_m_s"]), 1.0e-7
        assert_in_delta obs.doppler_hz, h(c["doppler_hz"]), 1.0e-7
        assert_in_delta obs.sat_clock_s, h(c["sat_clock_s"]), 1.0e-15
        assert_in_delta obs.elevation_deg, h(c["elevation_deg"]), 1.0e-8
        assert_in_delta obs.azimuth_deg, h(c["azimuth_deg"]), 1.0e-8

        {lx, ly, lz} = obs.los_unit
        [ex, ey, ez] = Enum.map(c["los_unit"], &h/1)
        assert_in_delta lx, ex, 5.0e-12
        assert_in_delta ly, ey, 5.0e-12
        assert_in_delta lz, ez, 5.0e-12

        assert NaiveDateTime.diff(
                 obs.transmit_time,
                 NaiveDateTime.from_iso8601!(c["transmit_time"]),
                 :microsecond
               ) == 0
      end
    end
  end

  defp h(hex), do: hex_to_float(hex)

  defp receiver(c) do
    [x, y, z] = Enum.map(c["receiver_ecef_m"], &h/1)
    {x, y, z}
  end
end
