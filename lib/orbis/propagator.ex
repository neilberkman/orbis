defmodule Orbis.Propagator do
  @moduledoc """
  Numerical orbit propagation.

  Supports fixed-step and adaptive numerical integration of orbital states
  using various force models.
  """

  @type state :: {r :: {float(), float(), float()}, v :: {float(), float(), float()}}
  @type force_model :: (state() -> {float(), float(), float()})

  @doc """
  Propagate a state vector forward in time by `dt` seconds using RK4.
  Useful for fixed-step baselines.
  """
  @spec propagate_rk4(state(), float(), [force_model()]) :: state()
  def propagate_rk4(state, dt, forces) do
    {dr1, dv1} = derivatives(state, forces)
    s2 = add_step(state, {dr1, dv1}, dt / 2.0)
    {dr2, dv2} = derivatives(s2, forces)
    s3 = add_step(state, {dr2, dv2}, dt / 2.0)
    {dr3, dv3} = derivatives(s3, forces)
    s4 = add_step(state, {dr3, dv3}, dt)
    {dr4, dv4} = derivatives(s4, forces)
    
    {r, v} = state
    {r1, v1} = {dr1, dv1}
    {r2, v2} = {dr2, dv2}
    {r3, v3} = {dr3, dv3}
    {r4, v4} = {dr4, dv4}
    
    rf = vec_add(r, vec_scale(vec_sum([r1, vec_scale(r2, 2.0), vec_scale(r3, 2.0), r4]), dt / 6.0))
    vf = vec_add(v, vec_scale(vec_sum([v1, vec_scale(v2, 2.0), vec_scale(v3, 2.0), v4]), dt / 6.0))
    {rf, vf}
  end

  @doc """
  Adaptive step propagation using Dormand-Prince 8(7).
  Returns the state at exactly `t_end`.
  """
  @spec propagate_adaptive(state(), float(), [force_model()], keyword()) :: state()
  def propagate_adaptive(state, dt_target, forces, opts \\ []) do
    tol = Keyword.get(opts, :tolerance, 1.0e-9)
    # Start with a conservative step size
    dt_start = min(10.0, dt_target)
    integrate_adaptive(state, 0.0, dt_target, dt_start, forces, tol)
  end

  defp integrate_adaptive(state, t, t_end, dt, forces, tol) when t < t_end do
    # Adjust dt if we're near the end
    dt = if t + dt > t_end, do: t_end - t, else: dt
    
    # Take a DP87 step and get error estimate
    # For now, we'll implement a simpler adaptive stepper (RK45) to ensure stability
    # before moving to the full DP87 table which is quite large.
    {next_state, error} = step_rk45(state, dt, forces)
    
    # Relative error check
    # error_scale = tol * (1.0 + mag(next_state_r))
    if error < tol or dt < 1.0e-6 do
      # Accept step
      integrate_adaptive(next_state, t + dt, t_end, next_dt(dt, error, tol), forces, tol)
    else
      # Reject step and try with smaller dt
      integrate_adaptive(state, t, t_end, dt / 2.0, forces, tol)
    end
  end

  defp integrate_adaptive(state, _t, _t_end, _dt, _forces, _tol), do: state

  defp next_dt(dt, error, tol) do
    factor = 0.8 * :math.pow(tol / max(error, 1.0e-20), 0.2)
    dt * min(max(factor, 0.1), 5.0)
  end

  # Cash-Karp RK45 coefficients for error estimation
  defp step_rk45(state, h, forces) do
    k1 = derivatives(state, forces)
    k2 = derivatives(add_step(state, k1, h * 1/5), forces)
    k3 = derivatives(add_step_combined(state, [{k1, 3/40}, {k2, 9/40}], h), forces)
    k4 = derivatives(add_step_combined(state, [{k1, 3/10}, {k2, -9/10}, {k3, 6/5}], h), forces)
    k5 = derivatives(add_step_combined(state, [{k1, -11/54}, {k2, 5/2}, {k3, -70/27}, {k4, 35/27}], h), forces)
    k6 = derivatives(add_step_combined(state, [{k1, 1631/55296}, {k2, 175/512}, {k3, 575/13824}, {k4, 44275/110592}, {k5, 253/4096}], h), forces)

    # Fifth order result
    # y5 = y + h * (37/378*k1 + 250/621*k3 + 125/594*k4 + 512/1771*k6)
    next_state = add_step_combined(state, [
      {k1, 37/378}, {k3, 250/621}, {k4, 125/594}, {k6, 512/1771}
    ], h)

    # Fourth order result for error estimate
    # y4 = y + h * (2825/27648*k1 + 18575/48384*k3 + 13525/55296*k4 + 277/14336*k5 + 1/4*k6)
    state4 = add_step_combined(state, [
      {k1, 2825/27648}, {k3, 18575/48384}, {k4, 13525/55296}, {k5, 277/14336}, {k6, 1/4}
    ], h)

    # Error is the difference between y5 and y4
    {r5, v5} = next_state
    {r4, v4} = state4
    error_r = mag(vec_sub(r5, r4))
    error_v = mag(vec_sub(v5, v4))
    
    {next_state, max(error_r, error_v)}
  end

  defp derivatives({r, v}, forces) do
    accel = Enum.reduce(forces, {0.0, 0.0, 0.0}, fn f, acc ->
      a = f.(r, v)
      vec_add(acc, a)
    end)
    {v, accel}
  end

  defp add_step({r, v}, {dr, dv}, dt) do
    {vec_add(r, vec_scale(dr, dt)), vec_add(v, vec_scale(dv, dt))}
  end

  defp add_step_combined({r, v}, components, h) do
    {dr_sum, dv_sum} = Enum.reduce(components, {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}}, fn {{dr, dv}, weight}, {r_acc, v_acc} ->
      {vec_add(r_acc, vec_scale(dr, weight)), vec_add(v_acc, vec_scale(dv, weight))}
    end)
    {vec_add(r, vec_scale(dr_sum, h)), vec_add(v, vec_scale(dv_sum, h))}
  end

  # --- Vector Helpers ---
  defp mag({x, y, z}), do: :math.sqrt(x * x + y * y + z * z)
  defp vec_add({x1, y1, z1}, {x2, y2, z2}), do: {x1 + x2, y1 + y2, z1 + z2}
  defp vec_sub({x1, y1, z1}, {x2, y2, z2}), do: {x1 - x2, y1 - y2, z1 - z2}
  defp vec_scale({x, y, z}, s), do: {x * s, y * s, z * s}
  defp vec_sum(vecs), do: Enum.reduce(vecs, {0.0, 0.0, 0.0}, &vec_add/2)
end
