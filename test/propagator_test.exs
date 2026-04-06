defmodule Orbis.PropagatorTest do
  use ExUnit.Case, async: true
  alias Orbis.Propagator
  alias Orbis.Forces.TwoBody

  describe "Two-Body RK4 propagation" do
    test "circular orbit remains circular after one step" do
      # Circular LEO at ~622km altitude (7000km radius)
      # v_circ = sqrt(mu / r) = sqrt(398600.4418 / 7000) = 7.54605 km/s
      r0 = {7000.0, 0.0, 0.0}
      v0 = {0.0, 7.5460533, 0.0}
      state = {r0, v0}
      
      forces = [&TwoBody.acceleration/2]
      
      # Small step (10s)
      {r1, v1} = Propagator.propagate_rk4(state, 10.0, forces)
      
      # Radius should still be ~7000km
      dist1 = :math.sqrt(elem(r1, 0)**2 + elem(r1, 1)**2 + elem(r1, 2)**2)
      assert_in_delta dist1, 7000.0, 1.0e-3
      
      # Speed should still be ~7.546km/s
      speed1 = :math.sqrt(elem(v1, 0)**2 + elem(v1, 1)**2 + elem(v1, 2)**2)
      assert_in_delta speed1, 7.5460533, 1.0e-6
    end

    test "propagates for half a period" do
      # r=7000km, T = 2*pi*sqrt(r^3/mu) = 5828.5s
      # Half period = 2914.25s. Position should be {-7000, 0, 0}
      r0 = {7000.0, 0.0, 0.0}
      v0 = {0.0, 7.5460533, 0.0}
      state = {r0, v0}
      forces = [&TwoBody.acceleration/2]
      
      # RK4 with 10s steps for ~2914s
      steps = 291
      dt = 10.0
      
      final_state = Enum.reduce(1..steps, state, fn _, acc ->
        Propagator.propagate_rk4(acc, dt, forces)
      end)
      
      {rf, _vf} = final_state
      assert_in_delta elem(rf, 0), -7000.0, 1.0 # 1km drift after 2900s at 10s steps is expected for RK4
    end
  end

  describe "J2 Perturbation" do
    test "nodes regress for an inclined orbit" do
      # Near sun-synchronous LEO (i=98 deg)
      # a = 7000 km, e = 0.001, i = 98 deg
      r0 = {7000.0, 0.0, 0.0}
      # v = sqrt(mu/r) * {0, cos(i), sin(i)}
      v_mag = :math.sqrt(398600.4418 / 7000.0)
      i = 98.0 * :math.pi() / 180.0
      v0 = {0.0, v_mag * :math.cos(i), v_mag * :math.sin(i)}
      state = {r0, v0}
      
      forces_kepler = [&TwoBody.acceleration/2]
      forces_j2 = [&TwoBody.acceleration/2, &Orbis.Forces.J2.acceleration/2]
      
      # Propagate for one orbit (~5800s)
      dt = 10.0
      steps = 580
      
      final_kepler = Enum.reduce(1..steps, state, fn _, acc ->
        Propagator.propagate_rk4(acc, dt, forces_kepler)
      end)
      
      final_j2 = Enum.reduce(1..steps, state, fn _, acc ->
        Propagator.propagate_rk4(acc, dt, forces_j2)
      end)
      
      {rk, _} = final_kepler
      {rj, _} = final_j2
      
      # J2 should cause a several-kilometer difference after one orbit
      diff = :math.sqrt(
        (elem(rk, 0) - elem(rj, 0))**2 +
        (elem(rk, 1) - elem(rj, 1))**2 +
        (elem(rk, 2) - elem(rj, 2))**2
      )
      
      assert diff > 5.0
      assert diff < 50.0
    end
  end

  describe "Adaptive Propagation (RK45)" do
    test "circular orbit is accurate over a full period" do
      # r=7000km, T = 5828.5s
      r0 = {7000.0, 0.0, 0.0}
      v0 = {0.0, 7.5460533, 0.0}
      state = {r0, v0}
      forces = [&TwoBody.acceleration/2]
      
      # Adaptive step for one full period
      final_state = Propagator.propagate_adaptive(state, 5828.5, forces, tolerance: 1.0e-9)
      
      {rf, _vf} = final_state
      # Should return close to the start {7000, 0, 0}
      assert_in_delta elem(rf, 0), 7000.0, 0.5
      assert_in_delta elem(rf, 1), 0.0, 0.5
    end
  end
end
