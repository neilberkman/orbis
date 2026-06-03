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

      # The solver's termination status is surfaced as an atom.
      assert sol.metadata.status in [
               :gradient_tolerance,
               :cost_tolerance,
               :step_tolerance,
               :max_evaluations
             ]
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

  # Decode an IEEE-754 double from its raw big-endian 8-byte hex string (the
  # fixture's bit-exact float encoding), e.g. "0x417b0050747d1762".
  defp hex_to_float("0x" <> hex) do
    bytes = hex |> String.pad_leading(16, "0") |> Base.decode16!(case: :mixed)
    <<value::float-64>> = bytes
    value
  end
end
