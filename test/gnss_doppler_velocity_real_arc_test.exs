defmodule Orbis.GNSS.DopplerVelocityRealArcTest do
  @moduledoc """
  Oracle-gated Doppler-velocity test on a real cheap-receiver arc.

  The committed `gnss_velocity_test.exs` validates `Orbis.GNSS.Velocity.solve`
  by synthetic injected-velocity self-recovery. This test pins the capability
  against a real receiver oracle: phone GPS L1 Doppler (D1C) from a Pixel 5 on a
  fast drive (2021-12-15 US-MTV-1, GSDC), recovering the receiver velocity vector
  and checking it against a decimeter-grade truth track.

  Oracle: central finite-difference of the GSDC carrier-phase post-processed
  ground-truth ECEF track (the `truth_ecef_m` field of the vendored RTKLIB demo5
  oracle), differenced over the matched GPST epoch grid. The vendored RTKLIB
  oracle JSONs are position-only and carry no Doppler and no velocity, so RTKLIB
  itself is not a usable Doppler-velocity oracle here; the truth-track finite
  difference is the defensible velocity reference.

  Method: each fixture epoch carries the GPS D1C Doppler (sign applied), the
  broadcast-code SPP receiver ECEF, the GPST epoch, and the finite-difference
  truth velocity. The test loads a thinned GPS-only broadcast NAV and runs the
  real `Velocity.solve(observable: :doppler)` so the solver computes its own
  satellite geometry from real broadcast ephemeris. The RINEX-observation parse
  and SPP are tested elsewhere and are baked into the inputs fixture.

  Pre-registered gate (see test/fixtures/rtk/generators/doppler-positioning-spec.md):
  median 3D velocity error <= 0.50 m/s, p95 <= 2.50 m/s, n >= 1000 eligible epochs.
  Tolerances are fixed assertions, never loosened to pass; out of tolerance or
  n < 1000 is a fail.
  """
  use ExUnit.Case, async: true

  alias Orbis.GNSS.Broadcast
  alias Orbis.GNSS.Velocity

  @fixture_dir Path.join(__DIR__, "fixtures/rtk")
  @nav_path Path.join(@fixture_dir, "doppler_velocity_gsdc_2021_12_15_mtv1_pixel5_l1_gps.nav")
  @inputs_path Path.join(
                 @fixture_dir,
                 "doppler_velocity_gsdc_2021_12_15_mtv1_pixel5_l1_inputs.json"
               )

  # Pre-registered tolerances.
  @median_bar_m_s 0.50
  @p95_bar_m_s 2.50
  @min_epochs 1000

  setup_all do
    nav = Broadcast.load!(@nav_path)
    inputs = @inputs_path |> File.read!() |> Jason.decode!()
    {:ok, nav: nav, inputs: inputs}
  end

  describe "Doppler receiver velocity on a real Pixel-5 arc" do
    test "recovers the 3D velocity vector within the pre-registered tolerance", ctx do
      carrier_hz = ctx.inputs["carrier_hz"]
      epochs = ctx.inputs["epochs"]

      assert is_number(carrier_hz)
      assert length(epochs) >= @min_epochs

      {errors, truth_speeds, solved, failed} =
        Enum.reduce(epochs, {[], [], 0, 0}, fn epoch, {errs, tsp, ok, bad} ->
          receiver = ecef(epoch["receiver_ecef_m"])
          truth = ecef(epoch["truth_velocity_m_s"])
          doppler = Enum.map(epoch["doppler_hz"], fn [sat, hz] -> {sat, hz} end)
          {:ok, naive} = NaiveDateTime.from_iso8601(epoch["epoch"])

          case Velocity.solve(ctx.nav, doppler, naive, receiver,
                 observable: :doppler,
                 carrier_hz: carrier_hz
               ) do
            {:ok, sol} ->
              err = norm3(sub3(sol.velocity_m_s, truth))
              {[err | errs], [norm3(truth) | tsp], ok + 1, bad}

            {:error, _reason} ->
              {errs, tsp, ok, bad + 1}
          end
        end)

      # Every eligible epoch should solve: each carries >= 4 GPS D1C sats with a
      # valid SPP position and a NAV record in window.
      assert failed == 0, "#{failed} epochs failed to solve"
      assert solved >= @min_epochs

      sorted = Enum.sort(errors)
      median = percentile(sorted, 0.5)
      p95 = percentile(sorted, 0.95)
      median_truth_speed = percentile(Enum.sort(truth_speeds), 0.5)

      # Guard the arc is the intended fast drive, not a near-static track where a
      # small velocity error would trivially pass.
      assert median_truth_speed > 5.0,
             "truth median speed #{median_truth_speed} m/s, arc is not the fast drive"

      IO.puts(
        "Doppler velocity (real Pixel-5 MTV-1): n=#{solved} " <>
          "median=#{Float.round(median, 4)} m/s p95=#{Float.round(p95, 4)} m/s " <>
          "truth_median_speed=#{Float.round(median_truth_speed, 2)} m/s"
      )

      assert median <= @median_bar_m_s,
             "median 3D velocity error #{median} m/s exceeds #{@median_bar_m_s} m/s bar"

      assert p95 <= @p95_bar_m_s,
             "p95 3D velocity error #{p95} m/s exceeds #{@p95_bar_m_s} m/s bar"
    end
  end

  defp ecef(%{"x" => x, "y" => y, "z" => z}), do: {x, y, z}

  defp sub3({ax, ay, az}, {bx, by, bz}), do: {ax - bx, ay - by, az - bz}

  defp norm3({x, y, z}), do: :math.sqrt(x * x + y * y + z * z)

  defp percentile([], _), do: nil

  defp percentile(sorted, pct) do
    n = length(sorted)
    rank = pct * (n - 1)
    lo = trunc(rank)
    hi = min(lo + 1, n - 1)
    frac = rank - lo
    Enum.at(sorted, lo) * (1.0 - frac) + Enum.at(sorted, hi) * frac
  end
end
