defmodule Orbis.GNSS.Velocity do
  @moduledoc """
  Recover a receiver's velocity and clock drift from one epoch of Doppler or
  pseudorange-rate measurements against a precise (SP3) ephemeris source.

  This is the velocity counterpart to `Orbis.GNSS.Positioning`: where single-point
  positioning inverts pseudoranges for position and clock bias, this module inverts
  pseudorange rates (or equivalently Doppler shifts) for the receiver velocity and
  clock drift, given a known receiver position. It reuses the forward model in
  `Orbis.GNSS.Observables` for the satellite geometry and velocity, and the explicit
  4x4 inverse in `Orbis.GNSS.Geometry.inv4/1` for the normal-equation solution.

  ## Observation model (standard Doppler / range-rate least squares)

  For satellite `i` with line-of-sight unit vector `e_i` (receiver -> satellite,
  ECEF), satellite velocity `v_sat_i`, receiver velocity `v_rx`, the measured
  pseudorange rate is

      rho_dot_i = e_i . (v_sat_i - v_rx) + c * (rx_clock_drift - sat_clock_drift_i)

  Collecting the four unknowns into `x = [v_rx_x, v_rx_y, v_rx_z, c * rx_clock_drift]`
  (the clock term carried in length-rate units, m/s), the model rearranges to a
  linear system

      e_i . v_sat_i - rho_dot_i - c * sat_clock_drift_i = e_i . v_rx - c * rx_clock_drift

  which is one row `H_i . x = y_i` with

      H_i = [-e_ix, -e_iy, -e_iz, 1]
      y_i = rho_dot_i - (e_i . v_sat_i) + c * sat_clock_drift_i

  Stacking the visible satellites and solving the (unweighted) normal equations
  gives `x = (H^T H)^-1 H^T y`, from which the receiver velocity is `{x0, x1, x2}`
  and the clock drift is `x3 / c` (seconds per second).

  ## How the geometry terms are obtained

  Every quantity comes from `Orbis.GNSS.Observables.predict/5` evaluated at the
  known receiver position, which is treated as static for the forward model:

    * `los_unit` is `e_i`, the receiver -> satellite ECEF unit vector;
    * `range_rate_m_s` for a static receiver is exactly `e_i . v_sat_i` (the
      `los . (v_sat - v_rx)` projection with `v_rx = 0`), already expressed in the
      consistent receive-epoch ECEF frame with light-time and Sagnac corrections
      applied. The satellite velocity is finite-differenced inside `predict/5`, so
      it is never accessed directly here.

  ## Satellite clock drift

  The per-satellite term `c * sat_clock_drift_i` is a known correction folded into
  `y_i`. By default it is taken as **zero**: the SP3 forward model exposes the
  satellite clock offset but not its time derivative, and a constant common
  satellite-clock-drift bias is absorbed by the estimated receiver clock drift,
  exactly as a common clock offset is absorbed by the receiver clock bias in
  single-point positioning. Callers who have per-satellite drift estimates may
  supply them through `opts[:sat_clock_drift]` (a `%{sat => drift_s_s}` map or a
  one-argument function of the satellite id); the same value is then subtracted
  from `y_i`.

  ## Doppler observations

  A Doppler shift relates to the pseudorange rate by `doppler_hz = -rho_dot * f / c`
  (the convention in `Orbis.GNSS.Observables`), so with `opts[:observable] = :doppler`
  each measurement is converted via `rho_dot = -doppler_hz * c / f` before the solve,
  with the carrier `f` taken from `opts[:carrier_hz]` (default L1, 1575.42 MHz).

  ## Result map

      %{
        velocity_m_s:    {vx, vy, vz},   # receiver velocity, ECEF m/s
        speed_m_s:       float(),        # norm of velocity_m_s
        clock_drift_s_s: float(),        # receiver clock drift, s/s
        residuals_m_s:   %{sat => r},    # post-fit range-rate residuals, m/s
        used_sats:       [sat],          # satellites contributing, build order
        n_satellites:    non_neg_integer()
      }
  """

  alias Orbis.GNSS.Core.Constants
  alias Orbis.GNSS.Core.Types
  alias Orbis.GNSS.Geometry
  alias Orbis.GNSS.Observables
  alias Orbis.GNSS.SP3

  @type vec3 :: {float(), float(), float()}
  @type receiver :: vec3() | %{x_m: number(), y_m: number(), z_m: number()}
  @type observation :: {String.t(), number()}

  @type result :: %{
          velocity_m_s: vec3(),
          speed_m_s: float(),
          clock_drift_s_s: float(),
          residuals_m_s: %{String.t() => float()},
          used_sats: [String.t()],
          n_satellites: non_neg_integer()
        }

  @doc """
  Solve for the receiver velocity and clock drift at one receive epoch.

  `source` is a loaded `Orbis.GNSS.SP3` precise product. `observations` is a list of
  `{satellite_id, value}` pairs: pseudorange rates in m/s by default, or Doppler
  shifts in Hz when `opts[:observable]` is `:doppler`. `epoch` is the receive
  epoch (a `NaiveDateTime` in the product's time scale). `receiver_position` is
  the known receiver ECEF position in metres, as `{x_m, y_m, z_m}` or
  `%{x_m: _, y_m: _, z_m: _}` (typically a prior position solve).

  ## Options

    * `:observable` - `:range_rate` (default) or `:doppler`.
    * `:carrier_hz` - carrier frequency for the Doppler conversion, default the
      L1 carrier `1575.42 MHz`.
    * `:sat_clock_drift` - per-satellite clock drift in s/s, as a `%{sat => drift}`
      map or a one-argument function of the satellite id; default zero for every
      satellite (absorbed into the receiver clock drift; see the module docs).
    * `:light_time` - apply the light-time correction in the forward model,
      default `true`.
    * `:sagnac` - apply the Sagnac / Earth-rotation correction, default `true`.

  A satellite whose forward-model prediction fails (e.g. no ephemeris at this
  epoch) is dropped from the used set rather than failing the whole solve.

  Returns `{:ok, result}` or `{:error, reason}`, where `reason` is one of
  `:no_observations` (empty list), `{:too_few_satellites, used, 4}` (fewer than
  four usable satellites), `:singular_geometry` (rank-deficient normal matrix),
  `:invalid_receiver` (malformed receiver), `{:invalid_observation, entry}`
  (malformed observation entry), or `{:duplicate_observation, sat}` (a satellite
  appears more than once). Never raises.
  """
  @spec solve(SP3.t(), [observation()], NaiveDateTime.t(), receiver(), keyword()) ::
          {:ok, result()} | {:error, term()}
  def solve(source, observations, epoch, receiver_position, opts \\ [])

  def solve(%SP3{} = sp3, observations, %NaiveDateTime{} = epoch, receiver_position, opts)
      when is_list(observations) do
    observable = Keyword.get(opts, :observable, :range_rate)
    carrier_hz = Keyword.get(opts, :carrier_hz, Constants.gps_l1_hz()) * 1.0
    sat_drift_fun = sat_clock_drift_fun(Keyword.get(opts, :sat_clock_drift))

    predict_opts = [
      light_time: Keyword.get(opts, :light_time, true),
      sagnac: Keyword.get(opts, :sagnac, true)
    ]

    with :ok <- ensure_nonempty(observations),
         {:ok, receiver} <- Types.normalize_ecef(receiver_position),
         {:ok, normalized} <- normalize_observations(observations, observable, carrier_hz),
         {:ok, rows} <-
           build_rows(sp3, normalized, epoch, receiver, sat_drift_fun, predict_opts),
         :ok <- ensure_enough(rows),
         {:ok, x} <- solve_normal_equations(rows) do
      {:ok, assemble(x, rows)}
    end
  end

  def solve(%SP3{}, observations, %NaiveDateTime{}, _receiver, _opts)
      when not is_list(observations), do: {:error, :no_observations}

  @doc """
  Convert a Doppler shift in Hz to a pseudorange rate in m/s.

  Uses `rho_dot = -doppler_hz * c / carrier_hz`, the inverse of the Doppler
  convention in `Orbis.GNSS.Observables` (`doppler_hz = -rho_dot * carrier_hz / c`).
  `carrier_hz` defaults to the L1 carrier.
  """
  @spec doppler_to_range_rate(number(), number()) :: float()
  def doppler_to_range_rate(doppler_hz, carrier_hz \\ Constants.gps_l1_hz()) do
    -doppler_hz * Constants.speed_of_light_m_s() / carrier_hz
  end

  @doc """
  Convert a pseudorange rate in m/s to a Doppler shift in Hz.

  Uses `doppler_hz = -rho_dot * carrier_hz / c`, matching `Orbis.GNSS.Observables`.
  `carrier_hz` defaults to the L1 carrier.
  """
  @spec range_rate_to_doppler(number(), number()) :: float()
  def range_rate_to_doppler(rho_dot_m_s, carrier_hz \\ Constants.gps_l1_hz()) do
    -rho_dot_m_s * carrier_hz / Constants.speed_of_light_m_s()
  end

  # --- observation normalization -------------------------------------------

  defp ensure_nonempty([]), do: {:error, :no_observations}
  defp ensure_nonempty(_), do: :ok

  # Validate shape, convert Doppler -> range rate, and reject duplicate sats.
  # Produces a list of {sat, rho_dot_m_s} in input order.
  defp normalize_observations(observations, observable, carrier_hz) do
    Enum.reduce_while(observations, {:ok, [], MapSet.new()}, fn entry, {:ok, acc, seen} ->
      case normalize_one(entry, observable, carrier_hz) do
        {:ok, {sat, _} = pair} ->
          if MapSet.member?(seen, sat) do
            {:halt, {:error, {:duplicate_observation, sat}}}
          else
            {:cont, {:ok, [pair | acc], MapSet.put(seen, sat)}}
          end

        {:error, _} = err ->
          {:halt, err}
      end
    end)
    |> case do
      {:ok, acc, _seen} -> {:ok, Enum.reverse(acc)}
      {:error, _} = err -> err
    end
  end

  defp normalize_one({sat, value}, :range_rate, _carrier_hz)
       when is_binary(sat) and is_number(value), do: {:ok, {sat, value * 1.0}}

  defp normalize_one({sat, value}, :doppler, carrier_hz) when is_binary(sat) and is_number(value),
    do: {:ok, {sat, doppler_to_range_rate(value * 1.0, carrier_hz)}}

  defp normalize_one(entry, _observable, _carrier_hz), do: {:error, {:invalid_observation, entry}}

  # --- design-row construction ---------------------------------------------

  # For each usable satellite produce a row %{sat, h, y, rho_dot}, where h is the
  # design row {-ex, -ey, -ez, 1.0} and y = rho_dot - e.v_sat + c*sat_drift. A
  # satellite whose prediction fails is dropped; a malformed receiver fails fast.
  defp build_rows(sp3, normalized, epoch, receiver, sat_drift_fun, predict_opts) do
    Enum.reduce_while(normalized, {:ok, []}, fn {sat, rho_dot}, {:ok, acc} ->
      case Observables.predict(sp3, sat, receiver, epoch, predict_opts) do
        {:ok, obs} ->
          {ex, ey, ez} = obs.los_unit
          e_dot_vsat = obs.range_rate_m_s
          sat_drift = sat_drift_fun.(sat)
          y = rho_dot - e_dot_vsat + Constants.speed_of_light_m_s() * sat_drift

          row = %{sat: sat, h: {-ex, -ey, -ez, 1.0}, y: y}
          {:cont, {:ok, [row | acc]}}

        {:error, _reason} ->
          # No usable geometry for this satellite at this epoch: drop it.
          {:cont, {:ok, acc}}
      end
    end)
    |> case do
      {:ok, acc} -> {:ok, Enum.reverse(acc)}
      {:error, _} = err -> err
    end
  end

  defp ensure_enough(rows) when length(rows) >= 4, do: :ok
  defp ensure_enough(rows), do: {:error, {:too_few_satellites, length(rows), 4}}

  # --- normal equations -----------------------------------------------------

  # Accumulate H^T H (4x4) and H^T y (4-vector), invert via inv4, x = inv * H^T y.
  defp solve_normal_equations(rows) do
    {ata, aty} =
      Enum.reduce(rows, {zero4x4(), {0.0, 0.0, 0.0, 0.0}}, fn %{h: h, y: y}, {ata, aty} ->
        {accumulate_ata(ata, h), accumulate_aty(aty, h, y)}
      end)

    case Geometry.inv4(ata) do
      {:ok, inv} -> {:ok, mat4_vec4(inv, aty)}
      :singular -> {:error, :singular_geometry}
    end
  end

  defp zero4x4 do
    row = {0.0, 0.0, 0.0, 0.0}
    {row, row, row, row}
  end

  # H^T H += h h^T (outer product of the design row with itself).
  defp accumulate_ata(ata, h) do
    hv = Tuple.to_list(h)

    rows =
      for {hi, i} <- Enum.with_index(hv) do
        existing = elem(ata, i) |> Tuple.to_list()

        for {hj, j} <- Enum.with_index(hv) do
          Enum.at(existing, j) + hi * hj
        end
        |> List.to_tuple()
      end

    List.to_tuple(rows)
  end

  # H^T y += h * y.
  defp accumulate_aty({a0, a1, a2, a3}, {h0, h1, h2, h3}, y) do
    {a0 + h0 * y, a1 + h1 * y, a2 + h2 * y, a3 + h3 * y}
  end

  defp mat4_vec4(m, {v0, v1, v2, v3}) do
    {r0, r1, r2, r3} = m

    {dot4(r0, v0, v1, v2, v3), dot4(r1, v0, v1, v2, v3), dot4(r2, v0, v1, v2, v3),
     dot4(r3, v0, v1, v2, v3)}
  end

  defp dot4({m0, m1, m2, m3}, v0, v1, v2, v3), do: m0 * v0 + m1 * v1 + m2 * v2 + m3 * v3

  # --- result assembly ------------------------------------------------------

  defp assemble({vx, vy, vz, c_drift}, rows) do
    x = {vx, vy, vz, c_drift}

    residuals =
      Map.new(rows, fn %{sat: sat, h: h, y: y} ->
        {sat, y - hx(h, x)}
      end)

    %{
      velocity_m_s: {vx, vy, vz},
      speed_m_s: :math.sqrt(vx * vx + vy * vy + vz * vz),
      clock_drift_s_s: c_drift / Constants.speed_of_light_m_s(),
      residuals_m_s: residuals,
      used_sats: Enum.map(rows, & &1.sat),
      n_satellites: length(rows)
    }
  end

  defp hx({h0, h1, h2, h3}, {x0, x1, x2, x3}), do: h0 * x0 + h1 * x1 + h2 * x2 + h3 * x3

  # --- satellite clock drift lookup ----------------------------------------

  defp sat_clock_drift_fun(nil), do: fn _sat -> 0.0 end

  defp sat_clock_drift_fun(map) when is_map(map), do: fn sat -> (map[sat] || 0.0) * 1.0 end

  defp sat_clock_drift_fun(fun) when is_function(fun, 1), do: fn sat -> fun.(sat) * 1.0 end
end
