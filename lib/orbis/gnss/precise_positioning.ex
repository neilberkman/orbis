defmodule Orbis.GNSS.PrecisePositioning do
  @moduledoc """
  Float-ambiguity carrier-phase positioning primitives.

  This is the first precise-positioning layer above the code and carrier-phase
  combinations in `Orbis.GNSS.IonosphereFree` / `Orbis.GNSS.CarrierPhase`. It
  solves one SP3-backed epoch from dual-frequency ionosphere-free code and phase
  observations:

      P_IF_i = rho_i(x) + b - c * dt_sat_i
      L_IF_i = rho_i(x) + b - c * dt_sat_i + N_i

  where `x` is the receiver ECEF position, `b` is the receiver clock in metres,
  and `N_i` is one float carrier-phase ambiguity per satellite, also in metres.
  The state is linearized and iterated over `[x, y, z, b, N_1, N_2, ...]`.

  This is deliberately **float ambiguity only**: no integer fixing, no LAMBDA,
  no double differencing, and no stochastic PPP process model. The output exposes
  the float ambiguities and residuals so later RTK/PPP layers can build on them.

  ## Observation shape

  Observations may be maps or tuples:

      %{satellite_id: "G05", code_m: 24_000_000.0, phase_m: 24_012_345.0}
      {"G05", 24_000_000.0, 24_012_345.0}

  `code_m` and `phase_m` should normally be ionosphere-free combinations. Use
  `Orbis.GNSS.IonosphereFree.iono_free/4` and
  `Orbis.GNSS.IonosphereFree.iono_free_phase_cycles/4` to form them from raw
  dual-frequency RINEX observations.
  """

  alias Orbis.GNSS.Core.Constants
  alias Orbis.GNSS.Observables
  alias Orbis.GNSS.Positioning
  alias Orbis.GNSS.SP3

  @default_max_iterations 8
  @default_position_tolerance_m 1.0e-4
  @default_clock_tolerance_m 1.0e-4
  @default_code_sigma_m 1.0
  @default_phase_sigma_m 0.01
  @pivot_epsilon 1.0e-12

  defmodule Solution do
    @moduledoc """
    Float-ambiguity phase positioning solution for one epoch.
    """

    @enforce_keys [
      :position,
      :rx_clock_s,
      :rx_clock_m,
      :ambiguities_m,
      :residuals_m,
      :used_sats,
      :metadata
    ]
    defstruct [
      :position,
      :rx_clock_s,
      :rx_clock_m,
      :ambiguities_m,
      :residuals_m,
      :used_sats,
      :metadata
    ]

    @type position :: %{x_m: float(), y_m: float(), z_m: float()}
    @type residual :: %{code_m: float(), phase_m: float()}

    @type t :: %__MODULE__{
            position: position(),
            rx_clock_s: float(),
            rx_clock_m: float(),
            ambiguities_m: %{String.t() => float()},
            residuals_m: %{String.t() => residual()},
            used_sats: [String.t()],
            metadata: %{
              iterations: pos_integer(),
              converged: boolean(),
              status: :position_tolerance | :max_iterations,
              code_rms_m: float(),
              phase_rms_m: float(),
              weighted_rms_m: float()
            }
          }
  end

  @typedoc "A dual-frequency ionosphere-free code/phase observation."
  @type observation ::
          %{satellite_id: String.t(), code_m: number(), phase_m: number()}
          | {String.t(), number(), number()}

  @typedoc "A receiver ECEF position in metres."
  @type receiver ::
          {number(), number(), number()} | %{x_m: number(), y_m: number(), z_m: number()}

  @doc """
  Solve a float-ambiguity carrier-phase position for one SP3-backed epoch.

  `source` is a loaded `Orbis.GNSS.SP3` product. `observations` is a list of
  ionosphere-free code/phase pairs for one epoch. `epoch` is interpreted in the
  SP3 product's time scale.

  ## Options

    * `:initial_guess` - `{x_m, y_m, z_m, clock_m}`. If omitted, the code
      observations are first passed through `Orbis.GNSS.Positioning.solve/4`
      with ionosphere/troposphere disabled, and that code-only solution seeds
      the float solve.
    * `:spp_initial_guess` - code-only SPP seed used only when `:initial_guess`
      is omitted (default `{0, 0, 0, 0}`).
    * `:code_sigma_m` - code row standard deviation (default `1.0` m).
    * `:phase_sigma_m` - phase row standard deviation (default `0.01` m).
    * `:max_iterations` - maximum nonlinear iterations (default `8`).
    * `:position_tolerance_m` - position-update convergence threshold
      (default `1.0e-4` m).
    * `:clock_tolerance_m` - receiver-clock update threshold (default
      `1.0e-4` m).

  Returns `{:ok, %Solution{}}` or `{:error, reason}`. Reasons include
  `:no_observations`, `{:too_few_satellites, used, 4}`,
  `{:duplicate_observation, sat}`, `{:invalid_observation, entry}`,
  `:invalid_initial_guess`, `{:invalid_sigma, key}`, `{:invalid_option, key}`,
  `{:code_seed_failed, reason}`, `{:no_ephemeris, sat, reason}`, and
  `:singular_geometry`. If the iteration limit is reached after a valid solve
  step, the function returns a solution with `metadata.converged == false` and
  `metadata.status == :max_iterations` so callers can inspect the residuals and
  decide whether to reject it.
  """
  @spec solve_float(SP3.t(), [observation()], NaiveDateTime.t(), keyword()) ::
          {:ok, Solution.t()} | {:error, term()}
  def solve_float(source, observations, epoch, opts \\ [])

  def solve_float(%SP3{} = sp3, observations, %NaiveDateTime{} = epoch, opts)
      when is_list(observations) do
    with :ok <- ensure_nonempty(observations),
         {:ok, obs} <- normalize_observations(observations),
         :ok <- ensure_enough(obs),
         {:ok, weights} <- weights(opts),
         {:ok, solve_opts} <- solve_options(opts),
         {:ok, state} <- initial_state(sp3, obs, epoch, opts) do
      iterate(sp3, obs, epoch, state, weights, solve_opts, 1)
    end
  end

  def solve_float(%SP3{}, observations, %NaiveDateTime{}, _opts) when not is_list(observations),
    do: {:error, :no_observations}

  # --- input normalization -------------------------------------------------

  defp ensure_nonempty([]), do: {:error, :no_observations}
  defp ensure_nonempty(_), do: :ok

  defp normalize_observations(observations) do
    Enum.reduce_while(observations, {:ok, [], MapSet.new()}, fn entry, {:ok, acc, seen} ->
      case normalize_one(entry) do
        {:ok, %{satellite_id: sat} = obs} ->
          if MapSet.member?(seen, sat) do
            {:halt, {:error, {:duplicate_observation, sat}}}
          else
            {:cont, {:ok, [obs | acc], MapSet.put(seen, sat)}}
          end

        {:error, _} = err ->
          {:halt, err}
      end
    end)
    |> case do
      {:ok, acc, _seen} ->
        sorted = acc |> Enum.reverse() |> Enum.sort_by(& &1.satellite_id)
        {:ok, sorted}

      {:error, _} = err ->
        err
    end
  end

  defp normalize_one(%{satellite_id: sat, code_m: code, phase_m: phase} = obs)
       when is_binary(sat) and is_number(code) and is_number(phase) do
    {:ok, %{satellite_id: sat, code_m: code / 1.0, phase_m: phase / 1.0, raw: obs}}
  end

  defp normalize_one({sat, code, phase} = obs)
       when is_binary(sat) and is_number(code) and is_number(phase) do
    {:ok, %{satellite_id: sat, code_m: code / 1.0, phase_m: phase / 1.0, raw: obs}}
  end

  defp normalize_one(entry), do: {:error, {:invalid_observation, entry}}

  defp ensure_enough(obs) when length(obs) >= 4, do: :ok
  defp ensure_enough(obs), do: {:error, {:too_few_satellites, length(obs), 4}}

  defp weights(opts) do
    code_sigma = Keyword.get(opts, :code_sigma_m, @default_code_sigma_m)
    phase_sigma = Keyword.get(opts, :phase_sigma_m, @default_phase_sigma_m)

    cond do
      not is_number(code_sigma) or code_sigma <= 0.0 ->
        {:error, {:invalid_sigma, :code_sigma_m}}

      not is_number(phase_sigma) or phase_sigma <= 0.0 ->
        {:error, {:invalid_sigma, :phase_sigma_m}}

      true ->
        {:ok, %{code: 1.0 / code_sigma, phase: 1.0 / phase_sigma}}
    end
  end

  defp solve_options(opts) do
    max_iterations = Keyword.get(opts, :max_iterations, @default_max_iterations)
    pos_tol = Keyword.get(opts, :position_tolerance_m, @default_position_tolerance_m)
    clock_tol = Keyword.get(opts, :clock_tolerance_m, @default_clock_tolerance_m)

    cond do
      not is_integer(max_iterations) or max_iterations < 1 ->
        {:error, {:invalid_option, :max_iterations}}

      not is_number(pos_tol) or pos_tol < 0.0 ->
        {:error, {:invalid_option, :position_tolerance_m}}

      not is_number(clock_tol) or clock_tol < 0.0 ->
        {:error, {:invalid_option, :clock_tolerance_m}}

      true ->
        {:ok,
         %{
           max_iterations: max_iterations,
           position_tolerance_m: pos_tol / 1.0,
           clock_tolerance_m: clock_tol / 1.0
         }}
    end
  end

  # --- initialization ------------------------------------------------------

  defp initial_state(sp3, obs, epoch, opts) do
    case Keyword.fetch(opts, :initial_guess) do
      {:ok, guess} ->
        with {:ok, {x, y, z, clock_m}} <- normalize_guess(guess) do
          {:ok, state_from_guess(obs, {x, y, z}, clock_m)}
        end

      :error ->
        spp_seed(sp3, obs, epoch, opts)
    end
  end

  defp normalize_guess({x, y, z, clock_m})
       when is_number(x) and is_number(y) and is_number(z) and is_number(clock_m),
       do: {:ok, {x / 1.0, y / 1.0, z / 1.0, clock_m / 1.0}}

  defp normalize_guess(_guess), do: {:error, :invalid_initial_guess}

  defp state_from_guess(obs, position, clock_m) do
    ambiguities =
      Map.new(obs, fn o ->
        {o.satellite_id, o.phase_m - o.code_m}
      end)

    %{position: position, clock_m: clock_m, ambiguities: ambiguities}
  end

  defp spp_seed(sp3, obs, epoch, opts) do
    observations = Enum.map(obs, &{&1.satellite_id, &1.code_m})
    spp_initial = Keyword.get(opts, :spp_initial_guess, {0.0, 0.0, 0.0, 0.0})

    case Positioning.solve(sp3, observations, epoch,
           ionosphere: false,
           troposphere: false,
           initial_guess: spp_initial,
           with_geodetic: false
         ) do
      {:ok, sol} ->
        pos = {sol.position.x_m, sol.position.y_m, sol.position.z_m}
        state = state_from_guess(obs, pos, sol.rx_clock_s * Constants.speed_of_light_m_s())
        {:ok, state}

      {:error, reason} ->
        {:error, {:code_seed_failed, reason}}
    end
  end

  # --- nonlinear solve -----------------------------------------------------

  defp iterate(sp3, obs, epoch, state, weights, opts, iter) do
    with {:ok, rows} <- build_rows(sp3, obs, epoch, state, weights),
         {:ok, dx} <- solve_normal_equations(rows, length(obs)) do
      next = apply_delta(state, obs, dx)
      {pos_step, clock_step} = step_norms(dx)

      cond do
        pos_step <= opts.position_tolerance_m and abs(clock_step) <= opts.clock_tolerance_m ->
          finalize(sp3, obs, epoch, next, weights, iter, true, :position_tolerance)

        iter >= opts.max_iterations ->
          finalize(sp3, obs, epoch, next, weights, iter, false, :max_iterations)

        true ->
          iterate(sp3, obs, epoch, next, weights, opts, iter + 1)
      end
    end
  end

  defp build_rows(sp3, obs, epoch, state, weights) do
    {rx, ry, rz} = state.position

    Enum.reduce_while(obs, {:ok, []}, fn o, {:ok, acc} ->
      case Observables.predict(sp3, o.satellite_id, {rx, ry, rz}, epoch) do
        {:ok, pred} ->
          {ex, ey, ez} = pred.los_unit
          sat_clock_m = Constants.speed_of_light_m_s() * (pred.sat_clock_s || 0.0)
          model_code = pred.geometric_range_m + state.clock_m - sat_clock_m
          ambiguity = Map.fetch!(state.ambiguities, o.satellite_id)

          code_row = %{
            kind: :code,
            sat: o.satellite_id,
            h: design_row({-ex, -ey, -ez, 1.0}, nil, obs),
            y: o.code_m - model_code,
            weight: weights.code
          }

          phase_row = %{
            kind: :phase,
            sat: o.satellite_id,
            h: design_row({-ex, -ey, -ez, 1.0}, o.satellite_id, obs),
            y: o.phase_m - (model_code + ambiguity),
            weight: weights.phase
          }

          {:cont, {:ok, [phase_row, code_row | acc]}}

        {:error, reason} ->
          {:halt, {:error, {:no_ephemeris, o.satellite_id, reason}}}
      end
    end)
    |> case do
      {:ok, rows} -> {:ok, Enum.reverse(rows)}
      {:error, _} = err -> err
    end
  end

  defp design_row({dx, dy, dz, dc}, ambiguity_sat, obs) do
    ambiguity_cols =
      Enum.map(obs, fn o ->
        if o.satellite_id == ambiguity_sat, do: 1.0, else: 0.0
      end)

    [dx, dy, dz, dc | ambiguity_cols]
  end

  defp apply_delta(state, obs, dx) do
    {rx, ry, rz} = state.position
    [dx0, dx1, dx2, dclock | ambiguity_deltas] = dx

    ambiguities =
      obs
      |> Enum.zip(ambiguity_deltas)
      |> Map.new(fn {o, d} ->
        {o.satellite_id, Map.fetch!(state.ambiguities, o.satellite_id) + d}
      end)

    %{
      position: {rx + dx0, ry + dx1, rz + dx2},
      clock_m: state.clock_m + dclock,
      ambiguities: ambiguities
    }
  end

  defp step_norms([dx, dy, dz, dclock | _]) do
    {:math.sqrt(dx * dx + dy * dy + dz * dz), dclock}
  end

  # --- dynamic least squares ----------------------------------------------

  defp solve_normal_equations(rows, n_obs) do
    n = 4 + n_obs

    {ata, aty} =
      Enum.reduce(rows, {zero_matrix(n), zero_vector(n)}, fn row, {ata, aty} ->
        h = Enum.map(row.h, &(&1 * row.weight))
        y = row.y * row.weight
        {accumulate_ata(ata, h), accumulate_aty(aty, h, y)}
      end)

    solve_linear(ata, aty)
  end

  defp zero_matrix(n), do: for(_ <- 1..n, do: zero_vector(n))
  defp zero_vector(n), do: for(_ <- 1..n, do: 0.0)

  defp accumulate_ata(ata, h) do
    Enum.with_index(ata)
    |> Enum.map(fn {row, i} ->
      hi = Enum.at(h, i)

      row
      |> Enum.with_index()
      |> Enum.map(fn {aij, j} -> aij + hi * Enum.at(h, j) end)
    end)
  end

  defp accumulate_aty(aty, h, y) do
    aty
    |> Enum.zip(h)
    |> Enum.map(fn {acc, hi} -> acc + hi * y end)
  end

  # Gaussian elimination with partial pivoting, returning :singular for
  # rank-deficient normal equations.
  defp solve_linear(a, b) do
    augmented = Enum.zip_with(a, b, fn row, bi -> row ++ [bi] end)

    case eliminate(augmented, 0, length(b)) do
      :singular -> {:error, :singular_geometry}
      upper -> {:ok, back_substitute(upper)}
    end
  end

  defp eliminate(rows, col, n) when col >= n, do: rows

  defp eliminate(rows, col, n) do
    {pivot_row, pivot_abs} =
      rows
      |> Enum.with_index()
      |> Enum.drop(col)
      |> Enum.map(fn {row, idx} -> {idx, abs(Enum.at(row, col))} end)
      |> Enum.max_by(fn {_idx, value} -> value end)

    if pivot_abs <= @pivot_epsilon do
      :singular
    else
      rows = swap_rows(rows, col, pivot_row)
      pivot = Enum.at(rows, col)
      pivot_value = Enum.at(pivot, col)

      rows =
        rows
        |> Enum.with_index()
        |> Enum.map(fn {row, idx} ->
          if idx <= col do
            row
          else
            factor = Enum.at(row, col) / pivot_value

            row
            |> Enum.zip(pivot)
            |> Enum.map(fn {rij, pij} -> rij - factor * pij end)
          end
        end)

      eliminate(rows, col + 1, n)
    end
  end

  defp swap_rows(rows, i, i), do: rows

  defp swap_rows(rows, i, j) do
    ri = Enum.at(rows, i)
    rj = Enum.at(rows, j)

    rows
    |> List.replace_at(i, rj)
    |> List.replace_at(j, ri)
  end

  defp back_substitute(rows) do
    n = length(rows)

    solved =
      (n - 1)..0//-1
      |> Enum.reduce(%{}, fn i, solved ->
        row = Enum.at(rows, i)

        known =
          if i == n - 1 do
            0.0
          else
            Enum.reduce((i + 1)..(n - 1), 0.0, fn j, acc ->
              acc + Enum.at(row, j) * Map.fetch!(solved, j)
            end)
          end

        xi = (Enum.at(row, n) - known) / Enum.at(row, i)
        Map.put(solved, i, xi)
      end)

    Enum.map(0..(n - 1), &Map.fetch!(solved, &1))
  end

  # --- result assembly -----------------------------------------------------

  defp finalize(sp3, obs, epoch, state, weights, iterations, converged, status) do
    with {:ok, residual_rows} <- residual_rows(sp3, obs, epoch, state) do
      code = residual_rows |> Enum.map(& &1.code_m)
      phase = residual_rows |> Enum.map(& &1.phase_m)

      {x, y, z} = state.position

      {:ok,
       %Solution{
         position: %{x_m: x, y_m: y, z_m: z},
         rx_clock_s: state.clock_m / Constants.speed_of_light_m_s(),
         rx_clock_m: state.clock_m,
         ambiguities_m: state.ambiguities,
         residuals_m:
           Map.new(residual_rows, fn r -> {r.sat, %{code_m: r.code_m, phase_m: r.phase_m}} end),
         used_sats: Enum.map(obs, & &1.satellite_id),
         metadata: %{
           iterations: iterations,
           converged: converged,
           status: status,
           code_rms_m: rms(code),
           phase_rms_m: rms(phase),
           weighted_rms_m: weighted_rms(residual_rows, weights)
         }
       }}
    end
  end

  defp residual_rows(sp3, obs, epoch, state) do
    {rx, ry, rz} = state.position

    Enum.reduce_while(obs, {:ok, []}, fn o, {:ok, acc} ->
      case Observables.predict(sp3, o.satellite_id, {rx, ry, rz}, epoch) do
        {:ok, pred} ->
          sat_clock_m = Constants.speed_of_light_m_s() * (pred.sat_clock_s || 0.0)
          model_code = pred.geometric_range_m + state.clock_m - sat_clock_m
          ambiguity = Map.fetch!(state.ambiguities, o.satellite_id)

          row = %{
            sat: o.satellite_id,
            code_m: o.code_m - model_code,
            phase_m: o.phase_m - (model_code + ambiguity)
          }

          {:cont, {:ok, [row | acc]}}

        {:error, reason} ->
          {:halt, {:error, {:no_ephemeris, o.satellite_id, reason}}}
      end
    end)
    |> case do
      {:ok, rows} -> {:ok, Enum.reverse(rows)}
      {:error, _} = err -> err
    end
  end

  defp rms([]), do: 0.0

  defp rms(xs) do
    :math.sqrt(Enum.reduce(xs, 0.0, fn x, acc -> acc + x * x end) / length(xs))
  end

  defp weighted_rms(rows, weights) do
    values =
      Enum.flat_map(rows, fn r ->
        [r.code_m * weights.code, r.phase_m * weights.phase]
      end)

    rms(values)
  end
end
