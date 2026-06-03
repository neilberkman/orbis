defmodule Orbis.PointPositioningTest do
  use ExUnit.Case, async: true

  alias Orbis.PointPositioning
  alias Orbis.PointPositioning.Solution
  alias Orbis.SP3

  # End-to-end check of the Elixir -> NIF -> astrodynamics-gnss SPP path. The
  # observations, epoch parameters, atmosphere coefficients, and synthesized
  # receiver truth come from a committed parity trace fixture; the precise
  # ephemeris is the matching SP3 file. Bit-exact physics parity is asserted in
  # the crate's own test suite — here we only prove the full round trip recovers
  # the truth, so a sub-millimetre solver-agreement bound is the right bar.
  @sp3_path Path.join(__DIR__, "fixtures/sp3/GRG0MGXFIN_20201760000_01D_15M_ORB.SP3")
  @trace_path Path.join(__DIR__, "fixtures/spp_trace_L2_tropo.json")

  # The trace's epoch index 48 is 2020-06-24 12:00:00 GPST (DOY 176, noon).
  @epoch ~N[2020-06-24 12:00:00]

  # Solver-agreement bound (meters). The crate documents agreement to a few
  # nanometres; the public boundary adds only term encode/decode, so a
  # sub-millimetre bound is comfortable and proves the wiring is lossless.
  @agreement_bound_m 1.0e-3

  setup_all do
    trace = @trace_path |> File.read!() |> Jason.decode!()
    inputs = trace["fixture"]["inputs"]
    final = trace["fixture"]["final_solution"]

    observations =
      Enum.map(inputs["observations"], fn obs ->
        {obs["sat_id"], hex_to_float(obs["p_meas_m"])}
      end)

    {:ok,
     sp3: SP3.load!(@sp3_path),
     observations: observations,
     alpha: inputs["klobuchar_alpha"] |> Enum.map(&hex_to_float/1) |> List.to_tuple(),
     beta: inputs["klobuchar_beta"] |> Enum.map(&hex_to_float/1) |> List.to_tuple(),
     pressure_hpa: hex_to_float(inputs["met"]["pressure_hpa"]),
     temperature_k: hex_to_float(inputs["met"]["temperature_k"]),
     relative_humidity: hex_to_float(inputs["met"]["relative_humidity"]),
     truth_x: Enum.map(final["truth_x"], &hex_to_float/1),
     truth_rx_clock_s: hex_to_float(final["truth_rx_clock_s"])}
  end

  describe "solve/4 end-to-end" do
    test "recovers the synthesized receiver truth to sub-millimetre", ctx do
      assert {:ok, %Solution{} = sol} =
               PointPositioning.solve(ctx.sp3, ctx.observations, @epoch,
                 ionosphere: true,
                 troposphere: true,
                 klobuchar_alpha: ctx.alpha,
                 klobuchar_beta: ctx.beta,
                 pressure_hpa: ctx.pressure_hpa,
                 temperature_k: ctx.temperature_k,
                 relative_humidity: ctx.relative_humidity,
                 # The crate fixture's frozen initial guess.
                 initial_guess: {4_500_000.0, 500_000.0, 4_500_000.0, 0.0}
               )

      [tx, ty, tz, _tb] = ctx.truth_x

      assert_in_delta sol.position.x_m, tx, @agreement_bound_m
      assert_in_delta sol.position.y_m, ty, @agreement_bound_m
      assert_in_delta sol.position.z_m, tz, @agreement_bound_m

      # Clock bias agrees to the same length bound, expressed in seconds.
      c_m_s = 299_792_458.0
      assert_in_delta sol.rx_clock_s, ctx.truth_rx_clock_s, @agreement_bound_m / c_m_s

      assert sol.metadata.converged
      assert sol.metadata.ionosphere_applied
      assert sol.metadata.troposphere_applied

      # The solver's termination status is surfaced as an atom. For this fixture
      # the trust-region solve stops on the step tolerance after 7 iterations;
      # the value is deterministic, so pin it exactly rather than accepting any
      # known status.
      assert sol.metadata.status == :step_tolerance
      assert sol.metadata.iterations == 7
    end

    test "accepts a sub-second receive epoch", ctx do
      # SPP receive time is a continuous f64 second, so a fractional epoch must
      # be accepted (not rejected as a non-integer-second epoch).
      epoch = ~N[2020-06-24 12:00:00.250000]

      assert {:ok, %Solution{}} =
               PointPositioning.solve(ctx.sp3, ctx.observations, epoch,
                 ionosphere: true,
                 troposphere: true,
                 klobuchar_alpha: ctx.alpha,
                 klobuchar_beta: ctx.beta,
                 pressure_hpa: ctx.pressure_hpa,
                 temperature_k: ctx.temperature_k,
                 relative_humidity: ctx.relative_humidity,
                 initial_guess: {4_500_000.0, 500_000.0, 4_500_000.0, 0.0}
               )
    end

    test "accepts a sub-second receive epoch given as a tuple", ctx do
      # The `{{y, m, d}, {h, min, s}}` epoch form must also carry a fractional
      # second, on the same footing as a NaiveDateTime.
      epoch = {{2020, 6, 24}, {12, 0, 0.25}}

      assert {:ok, %Solution{}} =
               PointPositioning.solve(ctx.sp3, ctx.observations, epoch,
                 ionosphere: true,
                 troposphere: true,
                 klobuchar_alpha: ctx.alpha,
                 klobuchar_beta: ctx.beta,
                 pressure_hpa: ctx.pressure_hpa,
                 temperature_k: ctx.temperature_k,
                 relative_humidity: ctx.relative_humidity,
                 initial_guess: {4_500_000.0, 500_000.0, 4_500_000.0, 0.0}
               )
    end

    test "returns geodetic, DOP, residuals and used satellites", ctx do
      assert {:ok, %Solution{} = sol} =
               PointPositioning.solve(ctx.sp3, ctx.observations, @epoch,
                 ionosphere: true,
                 troposphere: true,
                 klobuchar_alpha: ctx.alpha,
                 klobuchar_beta: ctx.beta,
                 pressure_hpa: ctx.pressure_hpa,
                 temperature_k: ctx.temperature_k,
                 relative_humidity: ctx.relative_humidity,
                 initial_guess: {4_500_000.0, 500_000.0, 4_500_000.0, 0.0}
               )

      # Truth is ~45 deg N, 7 deg E, 300 m (Turin-ish).
      assert_in_delta sol.geodetic.lat_rad, :math.pi() * 45.0 / 180.0, 1.0e-6
      assert_in_delta sol.geodetic.lon_rad, :math.pi() * 7.0 / 180.0, 1.0e-6
      assert_in_delta sol.geodetic.height_m, 300.0, 1.0e-2

      assert sol.dop.pdop > 0.0
      assert is_list(sol.used_sats)
      assert length(sol.used_sats) >= 4
      assert length(sol.used_sats) == length(sol.residuals_m)
      assert Enum.all?(sol.used_sats, &is_binary/1)
      # Every observation is accounted for as either used or rejected.
      assert length(sol.used_sats) + length(sol.rejected_sats) == length(ctx.observations)

      assert Enum.all?(sol.rejected_sats, fn {sat, reason} ->
               is_binary(sat) and reason in [:no_ephemeris, :low_elevation]
             end)
    end
  end

  describe "solve/4 degenerate geometry" do
    @degenerate_sp3 Path.join(__DIR__, "fixtures/sp3/degenerate_coincident_5sat.sp3")

    test "a rank-deficient geometry returns a solution with no DOP, not a crash" do
      sp3 = SP3.load!(@degenerate_sp3)

      # All five satellites share one ECEF position, so every line of sight is
      # identical and the geometry is rank-deficient. An identical pseudorange
      # keeps the rows coincident.
      observations = for prn <- 1..5, do: {"G0#{prn}", 20_181_863.0}

      # A receive epoch inside the product's [00:00, 00:15] window.
      epoch = ~N[2020-06-24 00:03:20]

      assert {:ok, %Solution{} = sol} =
               PointPositioning.solve(sp3, observations, epoch,
                 initial_guess: {6_378_137.0, 0.0, 0.0, 0.0}
               )

      # The solver still returns a position, but a rank-deficient geometry must
      # not report a dilution of precision.
      assert sol.dop == nil
    end
  end

  describe "solve/4 from broadcast ephemeris" do
    @nav_path Path.join(__DIR__, "fixtures/nav/ESBC00DNK_R_20201770000_01D_MN.rnx")

    # GPS pseudoranges synthesized from the committed broadcast NAV product with
    # the same forward model the solver inverts, for a known receiver near the
    # ESBC station at 2020-06-25 12:00 GPST. The solve must recover that truth.
    @broadcast_truth %{x_m: 3_512_900.0, y_m: 780_500.0, z_m: 5_248_700.0}
    @broadcast_obs [
      {"G07", 24_602_022.181241553},
      {"G08", 23_676_569.520090435},
      {"G10", 23_359_996.74001386},
      {"G15", 24_308_689.12412482},
      {"G16", 20_729_337.624163955},
      {"G18", 21_218_848.782066472},
      {"G20", 21_331_195.197190672},
      {"G21", 20_769_683.82405165},
      {"G26", 22_031_046.45549123},
      {"G27", 21_170_243.258043874}
    ]

    test "recovers a known receiver from broadcast GPS pseudoranges" do
      eph = Orbis.BroadcastEphemeris.load!(@nav_path)

      assert {:ok, %Solution{} = sol} =
               PointPositioning.solve(eph, @broadcast_obs, ~N[2020-06-25 12:00:00],
                 initial_guess: {3_513_900.0, 779_500.0, 5_249_700.0, 0.0}
               )

      assert_in_delta sol.position.x_m, @broadcast_truth.x_m, 1.0e-2
      assert_in_delta sol.position.y_m, @broadcast_truth.y_m, 1.0e-2
      assert_in_delta sol.position.z_m, @broadcast_truth.z_m, 1.0e-2
      assert length(sol.used_sats) == 10
      assert sol.dop.pdop > 0.0
    end

    test "solves a mixed GPS+Galileo set together with a per-system clock" do
      eph = Orbis.BroadcastEphemeris.load!(@nav_path)
      # The 10 GPS pseudoranges plus visible Galileo sats at the same epoch,
      # synthesized with the same forward model.
      galileo = [
        {"E05", 27_038_058.346363213},
        {"E09", 25_628_329.534503363},
        {"E13", 25_860_944.73927032}
      ]

      mixed = @broadcast_obs ++ galileo

      assert {:ok, %Solution{} = sol} =
               PointPositioning.solve(eph, mixed, ~N[2020-06-25 12:00:00],
                 initial_guess: {3_513_900.0, 779_500.0, 5_249_700.0, 0.0}
               )

      # Recovers the same receiver as the GPS-only solve, now using both systems.
      assert_in_delta sol.position.x_m, @broadcast_truth.x_m, 1.0e-2
      assert_in_delta sol.position.y_m, @broadcast_truth.y_m, 1.0e-2
      assert_in_delta sol.position.z_m, @broadcast_truth.z_m, 1.0e-2

      systems = sol.used_sats |> Enum.map(&String.first/1) |> Enum.uniq() |> Enum.sort()
      assert systems == ["E", "G"], "both constellations must contribute"
      # Multi-system DOP is not yet computed.
      assert sol.dop == nil
    end

    test "a too-small broadcast observation set is rejected through the broadcast path" do
      eph = Orbis.BroadcastEphemeris.load!(@nav_path)
      few = Enum.take(@broadcast_obs, 3)

      assert {:error, {:too_few_satellites, used}} =
               PointPositioning.solve(eph, few, ~N[2020-06-25 12:00:00],
                 initial_guess: {3_513_900.0, 779_500.0, 5_249_700.0, 0.0}
               )

      assert used < 4
    end
  end

  describe "solve/4 error paths" do
    test "fewer than four observations is rejected", ctx do
      few = Enum.take(ctx.observations, 3)

      assert {:error, {:too_few_satellites, used}} =
               PointPositioning.solve(ctx.sp3, few, @epoch,
                 initial_guess: {4_500_000.0, 500_000.0, 4_500_000.0, 0.0}
               )

      assert used < 4
    end

    test "a duplicated satellite observation is rejected", ctx do
      [first | _] = ctx.observations
      dup = [first | ctx.observations]

      assert {:error, {:duplicate_observation, sat}} =
               PointPositioning.solve(ctx.sp3, dup, @epoch,
                 initial_guess: {4_500_000.0, 500_000.0, 4_500_000.0, 0.0}
               )

      assert is_binary(sat)
    end
  end

  describe "error reason mapping" do
    # The `:too_few_satellites` and `:duplicate_observation` reasons are also
    # covered end-to-end above; `:singular_geometry` and `:ephemeris_lost` are
    # defensive crate paths that real SP3 inputs do not reach, so the mapping
    # onto the public contract is exercised directly here.
    test "every advertised NIF error reason maps to its public form" do
      assert PointPositioning.map_solve_error({:error, :too_few_satellites, 3}) ==
               {:error, {:too_few_satellites, 3}}

      assert PointPositioning.map_solve_error({:error, :singular_geometry}) ==
               {:error, :singular_geometry}

      assert PointPositioning.map_solve_error({:error, :duplicate_observation, "G01"}) ==
               {:error, {:duplicate_observation, "G01"}}

      assert PointPositioning.map_solve_error({:error, :ephemeris_lost, "G07"}) ==
               {:error, {:ephemeris_lost, "G07"}}
    end

    test "an unrecognized NIF result is wrapped rather than dropped" do
      assert PointPositioning.map_solve_error(:boom) == {:error, :boom}
    end
  end

  # Decode an IEEE-754 double from its raw big-endian 8-byte hex string (the
  # fixture's bit-exact float encoding), e.g. "0x417b0050747d1762".
  defp hex_to_float("0x" <> hex) do
    bytes = hex |> String.pad_leading(16, "0") |> Base.decode16!(case: :mixed)
    <<value::float-64>> = bytes
    value
  end
end
