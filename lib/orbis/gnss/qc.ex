defmodule Orbis.GNSS.QC do
  @moduledoc """
  Measurement-quality control for single-point positioning: a generic
  "is this measurement set self-consistent?" layer over
  `Orbis.GNSS.Positioning`.

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
  re-solves go back through `Orbis.GNSS.Positioning.solve/4`.

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

  alias Orbis.GNSS.Positioning
  alias Orbis.GNSS.Positioning.Solution

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

  Given an `Orbis.GNSS.Positioning.Solution`, this computes the (optionally
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
  the observation list and `Orbis.GNSS.Positioning.solve/4` is re-run on the
  remainder. The loop stops when RAIM passes, when there is no redundancy left
  to test (`dof <= 0`), or when the iteration bound is reached.

  ## Options

  All `Orbis.GNSS.Positioning.solve/4` options are forwarded to every re-solve
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
  @spec fde(term(), [Positioning.observation()], Positioning.epoch(), keyword()) ::
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
    case Positioning.solve(source, observations, epoch, solve_opts) do
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
  `p` must be strictly between 0 and 1, and `k` must be a positive integer.

  The chi-square CDF is the regularized lower incomplete gamma function
  `P(k/2, x/2)`. This function inverts that CDF by bracketing and bisection,
  evaluating `P(a, x)` with the standard series / continued-fraction split from
  Numerical Recipes. It is dependency-free but checked against scipy's
  `scipy.stats.chi2.ppf` oracle in the test fixture.
  """
  @spec chi2_inv(float(), pos_integer()) :: float()
  def chi2_inv(p, k) when is_number(p) and p > 0.0 and p < 1.0 and is_integer(k) and k >= 1 do
    a = 0.5 * k
    hi0 = max(k + 10.0 * :math.sqrt(2.0 * k), 1.0)
    hi = chi2_bracket_hi(p, a, hi0)
    chi2_bisect(p, a, 0.0, hi, 0)
  end

  def chi2_inv(p, k) do
    raise ArgumentError,
          "chi2_inv probability must be strictly between 0 and 1 and dof must be a positive integer, got p=#{inspect(p)}, dof=#{inspect(k)}"
  end

  defp chi2_bracket_hi(p, a, hi) do
    if chi2_cdf(hi, a) >= p do
      hi
    else
      chi2_bracket_hi(p, a, hi * 2.0)
    end
  end

  defp chi2_bisect(_p, _a, lo, hi, 120), do: 0.5 * (lo + hi)

  defp chi2_bisect(p, a, lo, hi, iter) do
    mid = 0.5 * (lo + hi)

    if chi2_cdf(mid, a) < p do
      chi2_bisect(p, a, mid, hi, iter + 1)
    else
      chi2_bisect(p, a, lo, mid, iter + 1)
    end
  end

  defp chi2_cdf(x, a), do: regularized_gamma_p(a, 0.5 * x)

  @gamma_eps 1.0e-15
  @gamma_fpmin 1.0e-300
  @gamma_itmax 1_000

  defp regularized_gamma_p(_a, x) when x <= 0.0, do: 0.0

  defp regularized_gamma_p(a, x) when x < a + 1.0 do
    gln = log_gamma(a)
    {sum, _del, _ap} = gamma_series(a, x, 1.0 / a, 1.0 / a, a, 1)
    sum * :math.exp(-x + a * :math.log(x) - gln)
  end

  defp regularized_gamma_p(a, x) do
    gln = log_gamma(a)
    q = gamma_continued_fraction(a, x) * :math.exp(-x + a * :math.log(x) - gln)
    1.0 - q
  end

  defp gamma_series(_a, _x, sum, del, ap, n) when n > @gamma_itmax, do: {sum, del, ap}

  defp gamma_series(a, x, sum, del, ap, n) do
    ap = ap + 1.0
    del = del * x / ap
    sum = sum + del

    if abs(del) < abs(sum) * @gamma_eps do
      {sum, del, ap}
    else
      gamma_series(a, x, sum, del, ap, n + 1)
    end
  end

  defp gamma_continued_fraction(a, x) do
    b = x + 1.0 - a
    c = 1.0 / @gamma_fpmin
    d = safe_denominator(b)
    d = 1.0 / d
    gamma_cf_iter(a, x, b, c, d, d, 1)
  end

  defp gamma_cf_iter(_a, _x, _b, _c, _d, h, n) when n > @gamma_itmax, do: h

  defp gamma_cf_iter(a, x, b, c, d, h, n) do
    an = -n * (n - a)
    b = b + 2.0
    d = safe_denominator(an * d + b)
    c = safe_denominator(b + an / c)
    d = 1.0 / d
    delta = d * c
    h = h * delta

    if abs(delta - 1.0) < @gamma_eps do
      h
    else
      gamma_cf_iter(a, x, b, c, d, h, n + 1)
    end
  end

  defp safe_denominator(x) when abs(x) < @gamma_fpmin, do: @gamma_fpmin
  defp safe_denominator(x), do: x

  # Lanczos log-gamma approximation, g=7, coefficients from Numerical Recipes /
  # common CPython-compatible implementations. The chi-square gate uses positive
  # half-integer/integer `a = k/2`, so the reflection branch is only here for
  # completeness.
  @lanczos [
    0.9999999999998099,
    676.5203681218851,
    -1259.1392167224028,
    771.3234287776531,
    -176.6150291621406,
    12.507343278686905,
    -0.13857109526572012,
    9.984369578019572e-6,
    1.5056327351493116e-7
  ]
  @sqrt_2pi 2.5066282746310002

  defp log_gamma(z) when z < 0.5 do
    :math.log(:math.pi()) - :math.log(:math.sin(:math.pi() * z)) - log_gamma(1.0 - z)
  end

  defp log_gamma(z) do
    z = z - 1.0
    [head | rest] = @lanczos

    x =
      rest
      |> Enum.with_index(1)
      |> Enum.reduce(head, fn {coef, i}, acc -> acc + coef / (z + i) end)

    t = z + 7.5
    :math.log(@sqrt_2pi) + (z + 0.5) * :math.log(t) - t + :math.log(x)
  end
end
