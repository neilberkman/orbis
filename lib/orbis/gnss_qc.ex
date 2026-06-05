defmodule Orbis.GnssQC do
  @moduledoc """
  Measurement-quality control for single-point positioning: a generic
  "is this measurement set self-consistent?" layer over
  `Orbis.PointPositioning`.

  Three standard receiver-autonomous integrity tools are provided, all from the
  textbook GNSS integrity literature:

    * **Elevation- (and optionally C/N0-) dependent measurement weighting.** A
      low-elevation or low-carrier-to-noise observation carries more noise, so it
      should be down-weighted in the solve and in the fault test. The model here
      is the standard RTKLIB-style elevation form (see
      `pseudorange_variance/2`).

    * **Residual-based RAIM (receiver autonomous integrity monitoring).** A
      chi-square goodness-of-fit test on the post-fit pseudorange residuals: the
      (optionally weighted) sum of squared residuals is compared against the
      chi-square critical value for the redundancy of the geometry. A statistic
      above the threshold flags an inconsistent measurement set (a likely
      fault). See `raim/2`.

    * **Fault detection and exclusion (FDE).** The standard leave-one-out
      exclusion loop: when RAIM detects a fault, the satellite with the largest
      normalized residual is removed and the position is re-solved, repeating
      until the set is self-consistent or too few satellites remain. See
      `fde/4`.

  All math is standard practice; no positioning math is duplicated here — the
  re-solves go back through `Orbis.PointPositioning.solve/4`.

  ## Degrees of freedom

  The position estimate has three position coordinates plus **one receiver clock
  per distinct GNSS system** (a mixed-constellation set carries one clock column per
  constellation). The number of estimated states is therefore
  `n_states = 3 + n_systems`, and the redundancy available to the fault test is

      dof = n_used - n_states = n_used - (3 + n_systems)

  where `n_systems` is the count of distinct system letters (the leading
  character of each satellite id, e.g. `"G"`, `"E"`) in the used satellites.
  With `dof <= 0` the geometry has no redundancy and no fault test is possible.
  """

  alias Orbis.PointPositioning
  alias Orbis.PointPositioning.Solution

  # Default elevation-weighting coefficients (metres). RTKLIB uses a = b = 0.3 m
  # for the carrier-phase/code error model; the same values are a sensible code
  # default here.
  @default_a 0.3
  @default_b 0.3

  # Default false-alarm probability for the chi-square fault-detection test.
  @default_p_fa 1.0e-3

  @typedoc "A `{satellite_id, elevation_deg}` or `{satellite_id, elevation_deg, cn0_dbhz}` entry."
  @type weight_entry ::
          {String.t(), number()} | {String.t(), number(), number()}

  @typedoc """
  The result of `raim/2`.

  `fault_detected?` is `true` when the test statistic exceeds the chi-square
  threshold. `test_statistic` is the (optionally weighted) sum of squared
  residuals; `threshold` is the chi-square critical value at the requested
  false-alarm probability for `dof` degrees of freedom (or `nil` when the
  geometry is not testable). `testable?` is `false` exactly when `dof <= 0`.
  `normalized_residuals` maps each used satellite to its standardized residual
  `r_i / sigma_i`, and `worst_sat` is the satellite with the largest-magnitude
  normalized residual (the leave-one-out exclusion candidate), or `nil` when
  there are no residuals.
  """
  @type raim_result :: %{
          fault_detected?: boolean(),
          test_statistic: float(),
          threshold: float() | nil,
          dof: integer(),
          testable?: boolean(),
          normalized_residuals: %{String.t() => float()},
          worst_sat: String.t() | nil
        }

  # ===========================================================================
  # Part 1 — measurement variance model
  # ===========================================================================

  @doc """
  Pseudorange measurement variance (m^2) from satellite elevation, using the
  standard elevation-dependent weighting model.

  The chosen form is the RTKLIB-style **variance** model

      sigma^2(el) = a^2 + b^2 / sin^2(el)

  with `el` the topocentric elevation. As `el -> 90 deg` the variance tends to
  `a^2 + b^2` (the zenith floor); as `el -> 0` it grows without bound, so a
  low-elevation observation is heavily down-weighted. The model is therefore
  monotonically decreasing in elevation: the lowest satellite always has the
  largest variance.

  ## Options

    * `:a` - zenith term, metres (default `0.3`)
    * `:b` - elevation-scaled term, metres (default `0.3`)
    * `:model` - `:elevation` (default) or `:elevation_cn0`
    * `:cn0` - carrier-to-noise density `C/N0` in dB-Hz, required for
      `:elevation_cn0`

  With `model: :elevation_cn0` a standard carrier-to-noise term is **added** to
  the elevation variance:

      sigma^2 = a^2 + b^2 / sin^2(el) + c0 * 10^(-cn0/10)

  A higher `C/N0` (a stronger signal) yields a smaller added variance, matching
  standard C/N0-based measurement weighting. `c0` is a reference scale
  (`:cn0_scale`, default `1.0` m^2).

  Returns the variance as a float, or `{:error, :invalid_elevation}` when
  `elevation_deg <= 0` (a satellite at or below the horizon has no usable
  weighting), and `{:error, :missing_cn0}` when the C/N0 model is selected
  without a `:cn0` value.
  """
  @spec pseudorange_variance(number(), keyword()) ::
          float() | {:error, :invalid_elevation | :missing_cn0}
  def pseudorange_variance(elevation_deg, opts \\ [])

  def pseudorange_variance(elevation_deg, _opts) when elevation_deg <= 0,
    do: {:error, :invalid_elevation}

  def pseudorange_variance(elevation_deg, opts) do
    a = Keyword.get(opts, :a, @default_a)
    b = Keyword.get(opts, :b, @default_b)
    model = Keyword.get(opts, :model, :elevation)

    sin_el = :math.sin(elevation_deg * :math.pi() / 180.0)
    elevation_var = a * a + b * b / (sin_el * sin_el)

    case model do
      :elevation ->
        elevation_var

      :elevation_cn0 ->
        case Keyword.get(opts, :cn0) do
          nil ->
            {:error, :missing_cn0}

          cn0 ->
            scale = Keyword.get(opts, :cn0_scale, 1.0)
            elevation_var + scale * :math.pow(10.0, -cn0 / 10.0)
        end
    end
  end

  @doc """
  Build a `satellite => sigma_m` map (the per-observation standard deviation in
  metres) for a list of weight entries.

  Each entry is `{sat, elevation_deg}` or, for the C/N0 model,
  `{sat, elevation_deg, cn0_dbhz}`. Options are forwarded to
  `pseudorange_variance/2`. An entry with an invalid elevation is dropped from
  the result.
  """
  @spec sigmas([weight_entry()], keyword()) :: %{String.t() => float()}
  def sigmas(entries, opts \\ []) when is_list(entries) do
    entries
    |> Enum.map(&entry_variance(&1, opts))
    |> Enum.reject(fn {_sat, var} -> match?({:error, _}, var) end)
    |> Map.new(fn {sat, var} -> {sat, :math.sqrt(var)} end)
  end

  @doc """
  Build a `satellite => weight` map, where each weight is `1 / sigma^2`, the
  inverse-variance weight used by a weighted least-squares solve and by the
  weighted RAIM test.

  Entries and options are as for `sigmas/2`. An entry with an invalid elevation
  is dropped.
  """
  @spec weight_vector([weight_entry()], keyword()) :: %{String.t() => float()}
  def weight_vector(entries, opts \\ []) when is_list(entries) do
    entries
    |> Enum.map(&entry_variance(&1, opts))
    |> Enum.reject(fn {_sat, var} -> match?({:error, _}, var) end)
    |> Map.new(fn {sat, var} -> {sat, 1.0 / var} end)
  end

  defp entry_variance({sat, el}, opts), do: {sat, pseudorange_variance(el, opts)}

  defp entry_variance({sat, el, cn0}, opts) do
    opts = opts |> Keyword.put(:model, :elevation_cn0) |> Keyword.put(:cn0, cn0)
    {sat, pseudorange_variance(el, opts)}
  end

  # ===========================================================================
  # Part 2 — residual-based RAIM (chi-square fault detection)
  # ===========================================================================

  @doc """
  Residual-based RAIM: a chi-square goodness-of-fit test on a solution's
  post-fit pseudorange residuals.

  Given an `Orbis.PointPositioning.Solution`, this computes the (optionally
  weighted) sum of squared residuals

      T = sum_i r_i^2 / sigma_i^2        (weighted)
      T = sum_i r_i^2                    (unit weights, the default)

  which under the no-fault hypothesis is distributed chi-square with `dof`
  degrees of freedom (see the moduledoc for `dof = n_used - (3 + n_systems)`).
  A fault is declared when `T` exceeds the chi-square critical value
  `chi2_inv(1 - p_fa, dof)` for the requested false-alarm probability.

  The per-satellite **normalized residual** `r_i / sigma_i` is the standardized
  residual; `worst_sat` is the satellite with the largest-magnitude normalized
  residual, i.e. the standard largest-normalized-residual exclusion candidate.

  ## Options

    * `:p_fa` - false-alarm probability, a number strictly between 0 and 1
      (default `1.0e-3`); any other value raises `ArgumentError`
    * `:weights` - either `:unit` (default; all sigma = 1) or a
      `%{sat => weight}` map of positive inverse-variance weights `1 / sigma_i^2`
      as built by `weight_vector/2`; an unspecified satellite defaults to unit
      weight. A non-positive or non-numeric weight raises `ArgumentError`
    * `:n_systems` - override the distinct-system count; by default it is
      derived from the leading character of each used satellite id

  When `dof <= 0` the geometry carries no redundancy: the result reports
  `fault_detected?: false`, `testable?: false`, `threshold: nil`, and the
  computed `dof` without raising.
  """
  @spec raim(Solution.t(), keyword()) :: raim_result()
  def raim(%Solution{} = solution, opts \\ []) do
    p_fa = Keyword.get(opts, :p_fa, @default_p_fa)
    weights_opt = Keyword.get(opts, :weights, :unit)

    validate_p_fa!(p_fa)
    validate_weights!(weights_opt)

    used = solution.used_sats
    residuals = solution.residuals_m
    pairs = Enum.zip(used, residuals)

    n_used = length(used)
    n_systems = Keyword.get(opts, :n_systems) || distinct_systems(used)
    n_states = 3 + n_systems
    dof = n_used - n_states

    {test_statistic, normalized} = statistic_and_normalized(pairs, weights_opt)

    worst_sat =
      case normalized do
        map when map_size(map) == 0 ->
          nil

        map ->
          {sat, _v} = Enum.max_by(map, fn {_sat, v} -> abs(v) end)
          sat
      end

    if dof <= 0 do
      %{
        fault_detected?: false,
        test_statistic: test_statistic,
        threshold: nil,
        dof: dof,
        testable?: false,
        normalized_residuals: normalized,
        worst_sat: worst_sat
      }
    else
      threshold = chi2_inv(1.0 - p_fa, dof)

      %{
        fault_detected?: test_statistic > threshold,
        test_statistic: test_statistic,
        threshold: threshold,
        dof: dof,
        testable?: true,
        normalized_residuals: normalized,
        worst_sat: worst_sat
      }
    end
  end

  defp statistic_and_normalized(pairs, weights_opt) do
    Enum.reduce(pairs, {0.0, %{}}, fn {sat, r}, {acc_t, acc_n} ->
      weight = weight_for(sat, weights_opt)
      # weight = 1 / sigma^2, so sigma = 1 / sqrt(weight).
      contribution = r * r * weight
      normalized = r * :math.sqrt(weight)
      {acc_t + contribution, Map.put(acc_n, sat, normalized)}
    end)
  end

  defp weight_for(_sat, :unit), do: 1.0
  defp weight_for(sat, weights) when is_map(weights), do: Map.get(weights, sat, 1.0)

  # The chi-square threshold maps p_fa through an inverse-normal that is only
  # defined on the open interval (0, 1); reject the endpoints (and any
  # non-number) up front rather than crashing inside the quantile.
  defp validate_p_fa!(p) when is_number(p) and p > 0.0 and p < 1.0, do: :ok

  defp validate_p_fa!(p) do
    raise ArgumentError,
          "raim :p_fa must be a number strictly between 0 and 1, got: #{inspect(p)}"
  end

  # Weights are inverse variances; a normalized residual takes sqrt(weight), so a
  # non-positive or non-numeric weight is rejected rather than producing a NaN or
  # an arithmetic error.
  defp validate_weights!(:unit), do: :ok

  defp validate_weights!(weights) when is_map(weights) do
    if Enum.all?(weights, fn {_sat, w} -> is_number(w) and w > 0.0 end) do
      :ok
    else
      raise ArgumentError, "raim :weights must all be positive numbers, got: #{inspect(weights)}"
    end
  end

  defp validate_weights!(other) do
    raise ArgumentError,
          "raim :weights must be :unit or a %{sat => weight} map, got: #{inspect(other)}"
  end

  defp distinct_systems(used) do
    used
    |> Enum.map(&String.first/1)
    |> Enum.uniq()
    |> length()
  end

  # ===========================================================================
  # Part 3 — fault detection and exclusion (FDE)
  # ===========================================================================

  @doc """
  Fault detection and exclusion: solve, run RAIM, and if a fault is detected,
  exclude the worst satellite and re-solve, repeating until the measurement set
  is self-consistent or too few satellites remain.

  This is the standard leave-one-out FDE loop. At each step the satellite with
  the largest normalized residual (`worst_sat` from `raim/2`) is dropped from
  the observation list and `Orbis.PointPositioning.solve/4` is re-run on the
  remainder. The loop stops when RAIM passes, when there is no redundancy left
  to test (`dof <= 0`), or when the iteration bound is reached.

  ## Options

  All `Orbis.PointPositioning.solve/4` options are forwarded to every re-solve
  (e.g. `:initial_guess`, `:ionosphere`, `:troposphere`). Additionally:

    * `:p_fa` - false-alarm probability for the RAIM test (default `1.0e-3`)
    * `:weights` - forwarded to `raim/2` (default `:unit`)
    * `:max_iterations` - maximum number of exclusions to attempt (default
      `length(observations) - 4`, never less than `0`)

  Returns `{:ok, %{solution: Solution.t(), excluded: [{sat, :raim_excluded}],
  iterations: non_neg_integer()}}` — `excluded` is the list of removed
  satellites in exclusion order (empty for a clean set), and `iterations` is the
  number of exclusion steps performed. Returns `{:error, reason}` if any solve
  fails (the `solve/4` reason is propagated, e.g.
  `{:too_few_satellites, used, required}`).
  """
  @spec fde(term(), [PointPositioning.observation()], PointPositioning.epoch(), keyword()) ::
          {:ok,
           %{
             solution: Solution.t(),
             excluded: [{String.t(), :raim_excluded}],
             iterations: non_neg_integer()
           }}
          | {:error, term()}
  def fde(source, observations, epoch, opts \\ []) when is_list(observations) do
    max_iterations =
      Keyword.get(opts, :max_iterations, max(length(observations) - 4, 0))

    raim_opts = Keyword.take(opts, [:p_fa, :weights, :n_systems])
    solve_opts = Keyword.drop(opts, [:p_fa, :weights, :n_systems, :max_iterations])

    fde_loop(source, observations, epoch, solve_opts, raim_opts, max_iterations, [], 0)
  end

  defp fde_loop(
         source,
         observations,
         epoch,
         solve_opts,
         raim_opts,
         max_iterations,
         excluded,
         iter
       ) do
    case PointPositioning.solve(source, observations, epoch, solve_opts) do
      {:ok, solution} ->
        result = raim(solution, raim_opts)

        cond do
          not result.fault_detected? ->
            {:ok, %{solution: solution, excluded: Enum.reverse(excluded), iterations: iter}}

          iter >= max_iterations or is_nil(result.worst_sat) ->
            # No more exclusions permitted/possible: return the best solution we
            # have, still carrying whatever was excluded so far.
            {:ok, %{solution: solution, excluded: Enum.reverse(excluded), iterations: iter}}

          true ->
            worst = result.worst_sat
            remaining = Enum.reject(observations, fn {sat, _pr} -> sat == worst end)

            fde_loop(
              source,
              remaining,
              epoch,
              solve_opts,
              raim_opts,
              max_iterations,
              [{worst, :raim_excluded} | excluded],
              iter + 1
            )
        end

      {:error, reason} ->
        {:error, reason}
    end
  end

  # ===========================================================================
  # Chi-square inverse CDF
  # ===========================================================================

  @doc """
  Chi-square inverse CDF (quantile): the value `x` such that
  `P(X <= x) = p` for a chi-square distribution with `k` degrees of freedom.

  Computed via the **Wilson-Hilferty approximation**

      x ~= k * (1 - 2/(9k) + z_p * sqrt(2/(9k)))^3

  where `z_p` is the standard-normal quantile for cumulative probability `p`
  (obtained from the Acklam rational inverse-normal approximation). For
  fault-detection use, `p = 1 - p_fa`.

  In the deep upper tail at small `k` (the regime the fault test uses, e.g.
  `p = 0.999` for `p_fa = 1e-3`) the Wilson-Hilferty form is accurate to a few
  percent and is a consistently **conservative** over-estimate of the true
  critical value (which is the safe direction for fault detection: it errs
  toward fewer false alarms). The published reference critical values at
  `p = 0.999`, alongside what this function actually returns, are:

      dof:        1       2       3       4       5
      reference: 10.828  13.816  16.266  18.467  20.515   (chi-square tables)
      chi2_inv:  11.157  14.133  16.551  18.724  20.751   (this function)

  i.e. an over-estimate of ~3.0% at `dof = 1` shrinking to ~1.2% at `dof = 5`.

  The result is clamped to be non-negative: a chi-square quantile is always
  `>= 0`, but the cubed Wilson-Hilferty base can go negative for strongly
  negative `z_p` (deep *lower* tail at small `k`), which the documented
  fault-detection usage (`p` near `1`) never reaches.
  """
  @spec chi2_inv(float(), pos_integer()) :: float()
  def chi2_inv(p, k) when k >= 1 do
    z = inv_normal_cdf(p)
    t = 2.0 / (9.0 * k)
    max(0.0, k * :math.pow(1.0 - t + z * :math.sqrt(t), 3.0))
  end

  # Acklam's rational approximation to the inverse standard-normal CDF.
  # Maximum relative error ~1.15e-9 over the full range.
  @a1 -3.969683028665376e1
  @a2 2.209460984245205e2
  @a3 -2.759285104469687e2
  @a4 1.383577518672690e2
  @a5 -3.066479806614716e1
  @a6 2.506628277459239e0

  @b1 -5.447609879822406e1
  @b2 1.615858368580409e2
  @b3 -1.556989798598866e2
  @b4 6.680131188771972e1
  @b5 -1.328068155288572e1

  @c1 -7.784894002430293e-3
  @c2 -3.223964580411365e-1
  @c3 -2.400758277161838e0
  @c4 -2.549732539343734e0
  @c5 4.374664141464968e0
  @c6 2.938163982698783e0

  @d1 7.784695709041462e-3
  @d2 3.224671290700398e-1
  @d3 2.445134137142996e0
  @d4 3.754408661907416e0

  @p_low 0.02425
  @p_high 1.0 - 0.02425

  defp inv_normal_cdf(p) when p > 0.0 and p < @p_low do
    q = :math.sqrt(-2.0 * :math.log(p))

    (((((@c1 * q + @c2) * q + @c3) * q + @c4) * q + @c5) * q + @c6) /
      ((((@d1 * q + @d2) * q + @d3) * q + @d4) * q + 1.0)
  end

  defp inv_normal_cdf(p) when p >= @p_low and p <= @p_high do
    q = p - 0.5
    r = q * q

    (((((@a1 * r + @a2) * r + @a3) * r + @a4) * r + @a5) * r + @a6) * q /
      (((((@b1 * r + @b2) * r + @b3) * r + @b4) * r + @b5) * r + 1.0)
  end

  defp inv_normal_cdf(p) when p > @p_high and p < 1.0 do
    q = :math.sqrt(-2.0 * :math.log(1.0 - p))

    -(((((@c1 * q + @c2) * q + @c3) * q + @c4) * q + @c5) * q + @c6) /
      ((((@d1 * q + @d2) * q + @d3) * q + @d4) * q + 1.0)
  end
end
