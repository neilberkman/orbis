defmodule Orbis.GNSS.ApplicationOracleTest do
  use ExUnit.Case, async: true

  import Orbis.TestHelpers, only: [assert_ulp: 4, hex_to_float: 1]

  alias Orbis.GNSS.DGNSS
  alias Orbis.GNSS.Geometry
  alias Orbis.GNSS.Navigation.LNAV
  alias Orbis.GNSS.Navigation.LNAV.Ephemeris
  alias Orbis.GNSS.Positioning
  alias Orbis.GNSS.PrecisePositioning
  alias Orbis.GNSS.Signal.CA
  alias Orbis.GNSS.Signal.Correlator
  alias Orbis.GNSS.SolutionReport
  alias Orbis.GNSS.SP3
  alias Orbis.GNSS.Velocity

  @golden_path Path.join(__DIR__, "fixtures/orbis_gnss_application_golden.json")
  @sp3_path Path.join(__DIR__, "fixtures/sp3/GRG0MGXFIN_20201760000_01D_15M_ORB.SP3")

  setup_all do
    {:ok, golden: @golden_path |> File.read!() |> Jason.decode!(), sp3: SP3.load!(@sp3_path)}
  end

  describe "GPS C/A and correlator vs Python oracle" do
    test "C/A chips and selected correlations match the independent generator", %{golden: golden} do
      ca = golden["signal"]["ca"]

      for {prn, expected} <- ca["chips"] do
        assert {:ok, chips} = CA.code(String.to_integer(prn))
        assert Enum.take(chips, length(expected)) == expected
      end

      for c <- ca["correlations"] do
        assert {:ok, a} = CA.code(c["a"])
        assert {:ok, b} = CA.code(c["b"])
        assert CA.correlation_at(a, b, c["lag"]) == c["value"]
      end
    end

    test "replica, coherent correlation, acquisition, and loss match Python/Numpy", %{
      golden: golden
    } do
      corr = golden["signal"]["correlator"]
      rep = corr["replica_case"]

      assert {:ok, samples} =
               Correlator.replica(rep["prn"],
                 num_samples: rep["num_samples"],
                 sample_rate_hz: h(rep["sample_rate_hz"]),
                 code_phase_chips: h(rep["code_phase_chips"])
               )

      assert samples == rep["samples"]

      c = corr["correlate_case"]

      iq =
        clean_signal(
          c["prn"],
          h(c["code_phase_chips"]),
          h(c["doppler_hz"]),
          64,
          h(c["sample_rate_hz"])
        )

      assert {:ok, got} =
               Correlator.correlate(iq, c["prn"],
                 sample_rate_hz: h(c["sample_rate_hz"]),
                 doppler_hz: h(c["doppler_hz"]),
                 code_phase_chips: h(c["code_phase_chips"])
               )

      assert_ulp(got.i, h(c["i"]), 0, "correlator I")
      assert_ulp(got.q, h(c["q"]), 0, "correlator Q")
      assert_ulp(got.power, h(c["power"]), 0, "correlator power")

      a = corr["acquire_case"]
      n = a["code_phase_bins"]

      full =
        clean_signal(
          a["prn"],
          h(a["injected_code_phase_chips"]),
          h(a["injected_doppler_hz"]),
          n,
          h(a["sample_rate_hz"])
        )

      assert {:ok, acq} =
               Correlator.acquire(full, a["prn"], sample_rate_hz: h(a["sample_rate_hz"]))

      assert_ulp(acq.code_phase_chips, h(a["code_phase_chips"]), 0, "acquisition phase")
      assert_ulp(acq.doppler_hz, h(a["doppler_hz"]), 0, "acquisition doppler")
      assert_ulp(acq.peak_power, h(a["peak_power"]), 0, "acquisition peak")
      assert_ulp(acq.metric, h(a["metric"]), 4, "acquisition metric")

      for c <- corr["coherent_loss"] do
        assert_ulp(
          Correlator.coherent_loss(h(c["freq_error_hz"]), h(c["integration_time_s"])),
          h(c["loss"]),
          0,
          "coherent loss"
        )

        assert_ulp(
          Correlator.coherent_loss_db(h(c["freq_error_hz"]), h(c["integration_time_s"])),
          h(c["loss_db"]),
          1,
          "coherent loss dB"
        )
      end
    end
  end

  describe "GPS LNAV vs Python oracle" do
    test "parity vectors and encoded subframes match the independent generator", %{golden: golden} do
      lnav = golden["lnav"]

      for c <- lnav["parity_cases"] do
        assert LNAV.parity(c["data24"], c["d29_prev"], c["d30_prev"]) == c["parity"]
      end

      opts = [
        tow: lnav["options"]["tow"],
        alert: lnav["options"]["alert"],
        anti_spoof: lnav["options"]["anti_spoof"],
        integrity: lnav["options"]["integrity"],
        tlm_message: lnav["options"]["tlm_message"]
      ]

      assert {:ok, subframes} = LNAV.encode(Ephemeris.example(), opts)

      for sf <- [1, 2, 3] do
        bits = Map.fetch!(subframes, sf)
        assert bits_to_string(bits) == lnav["subframes"][Integer.to_string(sf)]

        words =
          bits
          |> Enum.chunk_every(30)
          |> Enum.map(fn word ->
            "0x" <>
              (word |> bits_to_uint() |> Integer.to_string(16) |> String.pad_leading(8, "0"))
          end)
          |> Enum.map(&String.downcase/1)

        assert words == lnav["word_hex"][Integer.to_string(sf)]
      end
    end
  end

  describe "SP3-derived application helpers vs operation-order Python oracle" do
    test "visibility and weighted DOP match the oracle fixture", %{golden: golden, sp3: sp3} do
      app = golden["sp3_application"]
      rx = vec3(app["receiver_ecef_m"])
      epoch = NaiveDateTime.from_iso8601!(app["epoch"])

      visible = Geometry.visible(sp3, rx, epoch, systems: ["G"], elevation_mask_deg: 10.0)
      expected_visible = app["visible_gps_mask10"]

      assert Enum.map(visible, & &1.satellite_id) ==
               Enum.map(expected_visible, & &1["satellite_id"])

      for {got, exp} <- Enum.zip(visible, expected_visible) do
        assert_trig_geometry_ulp(
          got.elevation_deg,
          h(exp["elevation_deg"]),
          "#{got.satellite_id} visible elevation"
        )

        assert_trig_geometry_ulp(
          got.azimuth_deg,
          h(exp["azimuth_deg"]),
          "#{got.satellite_id} visible azimuth"
        )
      end

      dop = app["dop_weighted"]

      got =
        Geometry.dop(sp3, rx, epoch,
          satellites: dop["satellites"],
          weights: :elevation,
          light_time: true
        )

      for key <- [:gdop, :pdop, :hdop, :vdop, :tdop] do
        assert_trig_geometry_ulp(
          Map.fetch!(got, key),
          h(dop[Atom.to_string(key)]),
          "#{key} weighted DOP"
        )
      end
    end

    test "receiver velocity solve matches the numpy least-squares oracle", %{
      golden: golden,
      sp3: sp3
    } do
      v = golden["sp3_application"]["velocity"]
      rx = vec3(v["receiver_ecef_m"])
      epoch = NaiveDateTime.from_iso8601!(golden["sp3_application"]["epoch"])
      observations = Enum.map(v["observations"], &{&1["sat"], h(&1["rho_dot_m_s"])})

      assert {:ok, got} = Velocity.solve(sp3, observations, epoch, rx)

      {vx, vy, vz} = got.velocity_m_s
      [ex, ey, ez] = Enum.map(v["solution_velocity_m_s"], &h/1)

      assert_ulp(vx, ex, 0, "velocity x")
      assert_ulp(vy, ey, 0, "velocity y")
      assert_ulp(vz, ez, 0, "velocity z")
      assert_ulp(got.clock_drift_s_s, h(v["solution_clock_drift_s_s"]), 0, "clock drift")
    end

    test "DGNSS corrections and corrected rover observations match the Python oracle", %{
      golden: golden,
      sp3: sp3
    } do
      d = golden["sp3_application"]["dgnss"]
      epoch = NaiveDateTime.from_iso8601!(golden["sp3_application"]["epoch"])
      base = vec3(d["base_ecef_m"])
      base_obs = Enum.map(d["base_observations"], &{&1["sat"], h(&1["pseudorange_m"])})
      rover_obs = Enum.map(d["rover_observations"], &{&1["sat"], h(&1["pseudorange_m"])})

      assert {:ok, corrections} = DGNSS.corrections(sp3, base, base_obs, epoch)

      for {sat, expected_hex} <- d["corrections_m"] do
        assert_ulp(corrections[sat], h(expected_hex), 0, "#{sat} correction")
      end

      {corrected, dropped} = DGNSS.apply(rover_obs, corrections)
      assert dropped == []

      expected_corrected =
        d["corrected_rover"]
        |> Map.new(&{&1["sat"], h(&1["pseudorange_m"])})

      for {sat, pr} <- corrected do
        assert_ulp(pr, expected_corrected[sat], 0, "#{sat} corrected rover pseudorange")
      end
    end

    test "float carrier-phase positioning matches the Python normal-equation oracle", %{
      golden: golden,
      sp3: sp3
    } do
      p = golden["sp3_application"]["precise_positioning"]
      epoch = NaiveDateTime.from_iso8601!(golden["sp3_application"]["epoch"])
      [x0, y0, z0, b0] = Enum.map(p["initial_guess"], &h/1)

      observations =
        Enum.map(p["observations"], fn o ->
          %{satellite_id: o["sat"], code_m: h(o["code_m"]), phase_m: h(o["phase_m"])}
        end)

      assert {:ok, got} =
               PrecisePositioning.solve_float(sp3, observations, epoch,
                 initial_guess: {x0, y0, z0, b0}
               )

      sol = p["solution"]
      [ex, ey, ez] = Enum.map(sol["position_m"], &h/1)

      assert_ulp(got.position.x_m, ex, 0, "precise position x")
      assert_ulp(got.position.y_m, ey, 0, "precise position y")
      assert_ulp(got.position.z_m, ez, 0, "precise position z")
      assert_ulp(got.rx_clock_m, h(sol["clock_m"]), 0, "precise receiver clock m")
      assert_ulp(got.rx_clock_s, h(sol["clock_s"]), 0, "precise receiver clock s")
      assert got.used_sats == sol["used_sats"]

      for {sat, expected} <- sol["ambiguities_m"] do
        assert_ulp(got.ambiguities_m[sat], h(expected), 0, "#{sat} float ambiguity")
      end

      for {sat, residual} <- sol["residuals_m"] do
        got_residual = got.residuals_m[sat]
        assert_ulp(got_residual.code_m, h(residual["code_m"]), 0, "#{sat} code residual")
        assert_ulp(got_residual.phase_m, h(residual["phase_m"]), 0, "#{sat} phase residual")
      end

      assert got.metadata.iterations == sol["metadata"]["iterations"]
      assert got.metadata.converged == sol["metadata"]["converged"]
      assert got.metadata.status == String.to_atom(sol["metadata"]["status"])

      assert_ulp(
        got.metadata.code_rms_m,
        h(sol["metadata"]["code_rms_m"]),
        0,
        "precise code RMS"
      )

      assert_ulp(
        got.metadata.phase_rms_m,
        h(sol["metadata"]["phase_rms_m"]),
        0,
        "precise phase RMS"
      )

      assert_ulp(
        got.metadata.weighted_rms_m,
        h(sol["metadata"]["weighted_rms_m"]),
        0,
        "precise weighted RMS"
      )
    end

    test "multi-epoch float carrier-phase positioning matches the Python oracle", %{
      golden: golden,
      sp3: sp3
    } do
      p = golden["sp3_application"]["precise_positioning_multi_epoch"]
      [x0, y0, z0, b0] = Enum.map(p["initial_guess"], &h/1)

      epoch_observations =
        Enum.map(p["epochs"], fn epoch_row ->
          %{
            epoch: NaiveDateTime.from_iso8601!(epoch_row["epoch"]),
            observations:
              Enum.map(epoch_row["observations"], fn o ->
                %{satellite_id: o["sat"], code_m: h(o["code_m"]), phase_m: h(o["phase_m"])}
              end)
          }
        end)

      assert {:ok, got} =
               PrecisePositioning.solve_float_epochs(sp3, epoch_observations,
                 initial_guess: {x0, y0, z0, b0}
               )

      sol = p["solution"]
      [ex, ey, ez] = Enum.map(sol["position_m"], &h/1)

      assert_ulp(got.position.x_m, ex, 0, "multi precise position x")
      assert_ulp(got.position.y_m, ey, 0, "multi precise position y")
      assert_ulp(got.position.z_m, ez, 0, "multi precise position z")
      assert got.used_sats == sol["used_sats"]
      assert Enum.map(got.epochs, &NaiveDateTime.to_iso8601/1) == sol["epochs"]

      for {got_clock, expected} <- Enum.zip(got.epoch_clocks, sol["epoch_clocks"]) do
        assert NaiveDateTime.to_iso8601(got_clock.epoch) == expected["epoch"]

        assert_ulp(
          got_clock.rx_clock_m,
          h(expected["clock_m"]),
          0,
          "#{expected["epoch"]} clock m"
        )

        assert_ulp(
          got_clock.rx_clock_s,
          h(expected["clock_s"]),
          0,
          "#{expected["epoch"]} clock s"
        )
      end

      for {sat, expected} <- sol["ambiguities_m"] do
        assert_ulp(got.ambiguities_m[sat], h(expected), 0, "#{sat} multi float ambiguity")
      end

      expected_residuals =
        Map.new(sol["residuals_m"], fn r -> {{r["epoch"], r["sat"]}, r} end)

      for got_residual <- got.residuals_m do
        key = {NaiveDateTime.to_iso8601(got_residual.epoch), got_residual.satellite_id}
        expected = Map.fetch!(expected_residuals, key)

        assert_ulp(
          got_residual.code_m,
          h(expected["code_m"]),
          0,
          "#{elem(key, 0)} #{elem(key, 1)} multi code residual"
        )

        assert_ulp(
          got_residual.phase_m,
          h(expected["phase_m"]),
          0,
          "#{elem(key, 0)} #{elem(key, 1)} multi phase residual"
        )
      end

      assert got.metadata.iterations == sol["metadata"]["iterations"]
      assert got.metadata.converged == sol["metadata"]["converged"]
      assert got.metadata.status == String.to_atom(sol["metadata"]["status"])
      assert got.metadata.n_epochs == sol["metadata"]["n_epochs"]
      assert got.metadata.n_observations == sol["metadata"]["n_observations"]

      assert_ulp(
        got.metadata.code_rms_m,
        h(sol["metadata"]["code_rms_m"]),
        0,
        "multi precise code RMS"
      )

      assert_ulp(
        got.metadata.phase_rms_m,
        h(sol["metadata"]["phase_rms_m"]),
        0,
        "multi precise phase RMS"
      )

      assert_ulp(
        got.metadata.weighted_rms_m,
        h(sol["metadata"]["weighted_rms_m"]),
        0,
        "multi precise weighted RMS"
      )
    end

    test "SolutionReport sky rows and residual RMS match Python-derived expectations", %{
      golden: golden,
      sp3: sp3
    } do
      app = golden["sp3_application"]
      sr = app["solution_report"]
      epoch = NaiveDateTime.from_iso8601!(app["epoch"])

      solution = %Positioning.Solution{
        position: position_map(app["receiver_ecef_m"]),
        geodetic: nil,
        rx_clock_s: 0.0,
        system_clocks_s: %{"G" => 0.0},
        dop: nil,
        residuals_m: sr["solution"]["residuals_m"],
        used_sats: sr["solution"]["used_sats"],
        rejected_sats:
          Enum.map(sr["solution"]["rejected_sats"], fn r ->
            {r["sat"], String.to_atom(r["reason"])}
          end),
        metadata: %{iterations: 1, converged: true, status: :gradient_tolerance}
      }

      assert {:ok, report} = SolutionReport.build(solution, sp3, epoch)

      assert_ulp(
        report.summary.residual_rms_m,
        h(sr["residual_rms_m"]),
        0,
        "solution residual RMS"
      )

      expected_sky = Map.new(sr["sky"], &{&1["sat"], &1})

      for row <- report.satellites do
        exp = Map.fetch!(expected_sky, row.satellite_id)

        assert_trig_geometry_ulp(
          row.elevation_deg,
          h(exp["elevation_deg"]),
          "#{row.satellite_id} report elevation"
        )

        assert_trig_geometry_ulp(
          row.azimuth_deg,
          h(exp["azimuth_deg"]),
          "#{row.satellite_id} report azimuth"
        )
      end
    end
  end

  defp assert_trig_geometry_ulp(actual, expected, label),
    do: assert_ulp(actual, expected, trig_geometry_max_ulp(), label)

  # The golden is generated in the pinned macOS/Apple-libm parity environment,
  # where these trig-derived geometry fields are required to be bit-exact. Linux
  # CI builds against glibc libm; the same operation order differs slightly.
  defp trig_geometry_max_ulp do
    case :os.type() do
      {:unix, :darwin} -> 0
      _ -> 64
    end
  end

  defp clean_signal(prn, code_phase_chips, doppler_hz, n, fs) do
    {:ok, code} =
      Correlator.replica(prn,
        num_samples: n,
        sample_rate_hz: fs,
        code_phase_chips: code_phase_chips
      )

    w = 2.0 * :math.pi() * doppler_hz / fs

    code
    |> Enum.with_index()
    |> Enum.map(fn {c, k} ->
      theta = w * k
      {c * :math.cos(theta), c * :math.sin(theta)}
    end)
  end

  defp bits_to_string(bits), do: Enum.map_join(bits, &Integer.to_string/1)

  defp bits_to_uint(bits), do: Enum.reduce(bits, 0, fn b, acc -> acc * 2 + b end)

  defp position_map(hexes) do
    [x, y, z] = Enum.map(hexes, &h/1)
    %{x_m: x, y_m: y, z_m: z}
  end

  defp vec3(hexes) do
    [x, y, z] = Enum.map(hexes, &h/1)
    {x, y, z}
  end

  defp h(nil), do: nil
  defp h(hex), do: hex_to_float(hex)
end
