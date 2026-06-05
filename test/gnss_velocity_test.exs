defmodule Orbis.GnssVelocityTest do
  use ExUnit.Case, async: true

  alias Orbis.GnssGeometry
  alias Orbis.GnssObservables
  alias Orbis.GnssVelocity
  alias Orbis.SP3

  # The velocity estimate is validated by injected-velocity recovery, not by
  # self-consistency: we pick a TRUE receiver velocity and clock drift, synthesize
  # the pseudorange rate each visible GPS satellite would produce under the
  # documented model, then assert the solve recovers the truth. The forward
  # geometry (line of sight and the e.v_sat projection) comes from the same
  # precise SP3 fixture used by the positioning tests.
  @sp3_path Path.join(__DIR__, "fixtures/sp3/GRG0MGXFIN_20201760000_01D_15M_ORB.SP3")
  @epoch ~N[2020-06-24 12:00:00]

  # A receiver position with good GPS visibility for this fixture/epoch (the
  # positioning test's frozen initial guess sits in this region).
  @receiver {4_500_000.0, 500_000.0, 4_500_000.0}

  # Speed of light, m/s (matches Orbis.GnssObservables and the module under test).
  @c 299_792_458.0
  @f_l1 1_575_420_000.0

  setup_all do
    sp3 = SP3.load!(@sp3_path)

    sats =
      sp3
      |> GnssGeometry.visible(@receiver, @epoch, systems: ["G"], elevation_mask_deg: 5.0)
      |> Enum.map(& &1.satellite_id)

    # Sanity: the velocity solve needs at least four satellites.
    true = length(sats) >= 4

    {:ok, sp3: sp3, sats: sats}
  end

  # The static-receiver e.v_sat term for a satellite, from the forward model.
  defp e_dot_vsat(sp3, sat) do
    {:ok, obs} = GnssObservables.predict(sp3, sat, @receiver, @epoch)
    obs.range_rate_m_s
  end

  # Synthesize the measured pseudorange rate for `sat` under a TRUE receiver
  # velocity and clock drift, with the satellite clock drift taken as zero (the
  # estimator's default), so the model term cancels exactly:
  #   rho_dot = e.(v_sat - v_true) + c*drift_true
  #           = (e.v_sat) - e.v_true + c*drift_true.
  defp synth_rho_dot(sp3, sat, {tx, ty, tz}, drift_true) do
    {:ok, obs} = GnssObservables.predict(sp3, sat, @receiver, @epoch)
    {ex, ey, ez} = obs.los_unit
    e_dot_vtrue = ex * tx + ey * ty + ez * tz
    obs.range_rate_m_s - e_dot_vtrue + @c * drift_true
  end

  defp synth_observations(sp3, sats, v_true, drift_true) do
    Enum.map(sats, fn sat -> {sat, synth_rho_dot(sp3, sat, v_true, drift_true)} end)
  end

  describe "solve/5 injected-velocity recovery" do
    test "recovers a nonzero receiver velocity and clock drift to sub-mm/s", ctx do
      v_true = {12.0, -7.0, 3.0}
      drift_true = 1.0e-9

      observations = synth_observations(ctx.sp3, ctx.sats, v_true, drift_true)

      assert {:ok, result} = GnssVelocity.solve(ctx.sp3, observations, @epoch, @receiver)

      {vx, vy, vz} = result.velocity_m_s
      {tx, ty, tz} = v_true

      max_err =
        [abs(vx - tx), abs(vy - ty), abs(vz - tz)] |> Enum.max()

      assert max_err < 1.0e-4, "velocity max error #{max_err} m/s"
      assert abs(result.clock_drift_s_s - drift_true) < 1.0e-13
      assert result.n_satellites == length(ctx.sats)
      assert result.used_sats == ctx.sats
    end

    test "recovers a static receiver as ~zero speed", ctx do
      observations = synth_observations(ctx.sp3, ctx.sats, {0.0, 0.0, 0.0}, 0.0)

      assert {:ok, result} = GnssVelocity.solve(ctx.sp3, observations, @epoch, @receiver)

      assert result.speed_m_s < 1.0e-4
      {vx, vy, vz} = result.velocity_m_s
      assert abs(vx) < 1.0e-4 and abs(vy) < 1.0e-4 and abs(vz) < 1.0e-4
      assert abs(result.clock_drift_s_s) < 1.0e-13
    end
  end

  describe "solve/5 Doppler path" do
    test "agrees with the range-rate path", ctx do
      v_true = {12.0, -7.0, 3.0}
      drift_true = 1.0e-9

      rr_obs = synth_observations(ctx.sp3, ctx.sats, v_true, drift_true)

      doppler_obs =
        Enum.map(rr_obs, fn {sat, rho_dot} ->
          {sat, GnssVelocity.range_rate_to_doppler(rho_dot, @f_l1)}
        end)

      assert {:ok, rr} = GnssVelocity.solve(ctx.sp3, rr_obs, @epoch, @receiver)

      assert {:ok, dop} =
               GnssVelocity.solve(ctx.sp3, doppler_obs, @epoch, @receiver, observable: :doppler)

      {rx, ry, rz} = rr.velocity_m_s
      {dx, dy, dz} = dop.velocity_m_s

      assert abs(rx - dx) < 1.0e-6
      assert abs(ry - dy) < 1.0e-6
      assert abs(rz - dz) < 1.0e-6
      assert abs(rr.clock_drift_s_s - dop.clock_drift_s_s) < 1.0e-15
    end
  end

  describe "solve/5 residuals" do
    test "a consistent set has ~zero residuals", ctx do
      observations = synth_observations(ctx.sp3, ctx.sats, {12.0, -7.0, 3.0}, 1.0e-9)

      assert {:ok, result} = GnssVelocity.solve(ctx.sp3, observations, @epoch, @receiver)

      for {_sat, r} <- result.residuals_m_s do
        assert abs(r) < 1.0e-6
      end
    end

    test "perturbing one observation shows up in that satellite's residual", ctx do
      observations = synth_observations(ctx.sp3, ctx.sats, {12.0, -7.0, 3.0}, 1.0e-9)

      [{bad_sat, bad_val} | rest] = observations
      perturbed = [{bad_sat, bad_val + 2.0} | rest]

      assert {:ok, result} = GnssVelocity.solve(ctx.sp3, perturbed, @epoch, @receiver)

      perturbed_residual = abs(result.residuals_m_s[bad_sat])

      others_max =
        result.residuals_m_s
        |> Map.delete(bad_sat)
        |> Map.values()
        |> Enum.map(&abs/1)
        |> Enum.max()

      # The perturbed satellite's residual is clearly nonzero and dominates.
      assert perturbed_residual > 0.5
      assert perturbed_residual > others_max
    end
  end

  describe "solve/5 errors (no raise)" do
    test "empty observations is tagged", ctx do
      assert {:error, :no_observations} = GnssVelocity.solve(ctx.sp3, [], @epoch, @receiver)
    end

    test "fewer than four satellites is tagged", ctx do
      three = ctx.sats |> Enum.take(3) |> synth_three(ctx.sp3)

      assert {:error, {:too_few_satellites, 3, 4}} =
               GnssVelocity.solve(ctx.sp3, three, @epoch, @receiver)
    end

    test "a malformed receiver is tagged", ctx do
      observations = synth_observations(ctx.sp3, ctx.sats, {1.0, 2.0, 3.0}, 0.0)

      assert {:error, :invalid_receiver} =
               GnssVelocity.solve(ctx.sp3, observations, @epoch, {:bad, :receiver, nil})
    end

    test "a malformed observation entry is tagged", ctx do
      assert {:error, {:invalid_observation, {"G01", :nope}}} =
               GnssVelocity.solve(ctx.sp3, [{"G01", :nope}], @epoch, @receiver)
    end

    test "a duplicate satellite is tagged", ctx do
      assert {:error, {:duplicate_observation, "G01"}} =
               GnssVelocity.solve(
                 ctx.sp3,
                 [{"G01", 1.0}, {"G02", 2.0}, {"G01", 3.0}, {"G05", 4.0}],
                 @epoch,
                 @receiver
               )
    end

    test "a rank-deficient normal matrix is reported as singular (no raise)", _ctx do
      # A full-rank four-satellite geometry from a real precise ephemeris never
      # produces a rank-deficient normal matrix, so the solve's singular branch is
      # exercised at the inverse it depends on: GnssGeometry.inv4/1 returns
      # :singular for a rank-deficient matrix, which solve/5 maps to
      # {:error, :singular_geometry}. This is the same contract the positioning
      # path documents for its own non-physical singular case.
      assert :singular =
               GnssGeometry.inv4(
                 {{1.0, 0.0, 0.0, 0.0}, {0.0, 1.0, 0.0, 0.0}, {0.0, 0.0, 1.0, 0.0},
                  {0.0, 0.0, 0.0, 0.0}}
               )
    end
  end

  # Build observations for exactly three given sats (for the too-few test).
  defp synth_three(sats, sp3) do
    Enum.map(sats, fn sat -> {sat, e_dot_vsat(sp3, sat)} end)
  end
end
