defmodule Orbis.IonosphereFree do
  @moduledoc """
  The dual-frequency ionosphere-free linear combination of pseudoranges.

  The first-order ionospheric group delay on a GNSS pseudorange is dispersive: to
  first order it scales as `1 / f^2`, so a signal at carrier `f` is delayed by
  `I(f) = K / f^2` for a slant-path constant `K` proportional to the total
  electron content (TEC). Measuring the same range on two carriers `f1` and `f2`
  therefore gives two observations that share the geometry but carry different
  ionospheric delays, and a fixed linear combination of the two cancels the
  `1 / f^2` term exactly:

      PR_IF = (f1^2 * PR1 - f2^2 * PR2) / (f1^2 - f2^2)

  Writing `gamma = f1^2 / (f1^2 - f2^2)` this is the affine combination

      PR_IF = gamma * PR1 - (gamma - 1) * PR2

  Substituting `PR_i = R + K / f_i^2` (a true range `R` plus the first-order
  ionospheric delay on band `i`) the `K` terms cancel and `PR_IF = R`. A position
  solve fed these combined pseudoranges therefore needs no ionosphere model — call
  `Orbis.PointPositioning.solve/4` with `ionosphere: false` (the troposphere term,
  which is non-dispersive and does not cancel, still applies).

  ## Frequency table

  The standard carrier frequencies used here, in hertz:

  | System  | Band  | Frequency (MHz) |
  |---------|-------|-----------------|
  | GPS     | L1    | 1575.42         |
  | GPS     | L2    | 1227.60         |
  | Galileo | E1    | 1575.42         |
  | Galileo | E5a   | 1176.45         |
  | BeiDou  | B1I   | 1561.098        |
  | BeiDou  | B3I   | 1268.52         |

  The default combination pair per system is GPS `{:l1, :l2}`, Galileo
  `{:e1, :e5a}`, BeiDou `{:b1i, :b3i}`.

  ## Noise amplification

  The combination is not free: because it is a weighted difference of two noisy
  observations, uncorrelated band noise of equal standard deviation `sigma` is
  amplified to `sigma * sqrt(gamma^2 + (gamma - 1)^2)`. For GPS L1/L2 this factor
  is about `2.978`; for Galileo E1/E5a about `2.588`. See `noise_amplification/2`.

  ## Non-goals

  This module builds only the first-order ionosphere-free pseudorange
  combination. It deliberately does not implement carrier-phase combinations,
  the Melbourne-Wubbena / wide-lane / geometry-free combinations, ambiguity
  resolution, the second-order ionospheric term, or triple-frequency
  combinations.
  """

  alias Orbis.RinexObs

  @typedoc "A pseudorange observation `{satellite_id, range_m}`."
  @type observation :: {String.t(), float()}

  @typedoc "A reason a satellite was dropped from the paired set."
  @type drop_reason :: :missing_band1 | :missing_band2 | :unknown_system

  # Standard carrier frequencies in hertz. GPS L1 and Galileo E1 share the
  # 1575.42 MHz carrier; they remain distinct named bands per system.
  @frequencies %{
    "G" => %{l1: 1_575_420_000.0, l2: 1_227_600_000.0},
    "E" => %{e1: 1_575_420_000.0, e5a: 1_176_450_000.0},
    "C" => %{b1i: 1_561_098_000.0, b3i: 1_268_520_000.0}
  }

  # The default band pair to combine, per system.
  @default_pairs %{
    "G" => {:l1, :l2},
    "E" => {:e1, :e5a},
    "C" => {:b1i, :b3i}
  }

  @doc """
  The full carrier-frequency table as `%{system => %{band => f_hz}}`.
  """
  @spec frequencies() :: %{String.t() => %{atom() => float()}}
  def frequencies, do: @frequencies

  @doc """
  The carrier frequency in hertz for `system` (`"G"`, `"E"`, `"C"`) and `band`.

  Returns `{:ok, f_hz}` or `{:error, {:unknown_band, system, band}}`.
  """
  @spec frequency(String.t(), atom()) ::
          {:ok, float()} | {:error, {:unknown_band, String.t(), atom()}}
  def frequency(system, band) when is_binary(system) and is_atom(band) do
    case @frequencies do
      %{^system => %{^band => f}} -> {:ok, f}
      _ -> {:error, {:unknown_band, system, band}}
    end
  end

  @doc """
  The standard combination band pair for `system`.

  Returns `{:ok, {band1, band2}}` or `{:error, {:unknown_system, system}}`.
  """
  @spec default_pair(String.t()) ::
          {:ok, {atom(), atom()}} | {:error, {:unknown_system, String.t()}}
  def default_pair(system) when is_binary(system) do
    case @default_pairs do
      %{^system => pair} -> {:ok, pair}
      _ -> {:error, {:unknown_system, system}}
    end
  end

  @doc """
  The ionosphere-free combination coefficient `gamma = f1^2 / (f1^2 - f2^2)`.

  Returns `{:error, :equal_frequencies}` when `f1 == f2` (the combination is
  undefined; the denominator vanishes). Never raises.
  """
  @spec gamma(float(), float()) :: {:ok, float()} | {:error, :equal_frequencies}
  def gamma(f1, f2) when is_number(f1) and is_number(f2) do
    if f1 == f2 do
      {:error, :equal_frequencies}
    else
      f1sq = f1 * f1
      f2sq = f2 * f2
      {:ok, f1sq / (f1sq - f2sq)}
    end
  end

  @doc """
  The noise-amplification factor `sqrt(gamma^2 + (gamma - 1)^2)`.

  This is the factor by which uncorrelated equal-variance band noise is amplified
  into the combined pseudorange. About `2.978` for GPS L1/L2 and `2.588` for
  Galileo E1/E5a. Returns `{:error, :equal_frequencies}` when `f1 == f2`. Never
  raises.
  """
  @spec noise_amplification(float(), float()) :: {:ok, float()} | {:error, :equal_frequencies}
  def noise_amplification(f1, f2) when is_number(f1) and is_number(f2) do
    case gamma(f1, f2) do
      {:ok, g} -> {:ok, :math.sqrt(g * g + (g - 1.0) * (g - 1.0))}
      {:error, _} = err -> err
    end
  end

  @doc """
  The ionosphere-free pseudorange from two carrier-band pseudoranges.

      PR_IF = (f1^2 * pr1 - f2^2 * pr2) / (f1^2 - f2^2)
            = gamma * pr1 - (gamma - 1) * pr2

  `pr1`/`pr2` are in metres on carriers `f1`/`f2` (hertz). Returns `{:ok, pr_if}`,
  or `{:error, :equal_frequencies}` when `f1 == f2`. Never raises.
  """
  @spec iono_free(float(), float(), float(), float()) ::
          {:ok, float()} | {:error, :equal_frequencies}
  def iono_free(pr1, pr2, f1, f2)
      when is_number(pr1) and is_number(pr2) and is_number(f1) and is_number(f2) do
    case gamma(f1, f2) do
      {:ok, g} -> {:ok, g * pr1 - (g - 1.0) * pr2}
      {:error, _} = err -> err
    end
  end

  @doc """
  Combine two per-satellite pseudorange bands into ionosphere-free pseudoranges.

  `band1` and `band2` are `[{satellite_id, range_m}]` lists for the two carriers.
  Satellites are paired by id, and each pair is combined with the frequency pair
  for that satellite's system (the leading letter of the id), so a mixed
  GPS+Galileo+BeiDou set uses each system's own carriers.

  Returns `{combined, dropped}` where `combined` is the ascending-by-id list of
  `{satellite_id, pr_if}` and `dropped` reports every satellite that could not be
  combined as `{satellite_id, reason}`:

    * `:missing_band2` — present in `band1` only;
    * `:missing_band1` — present in `band2` only;
    * `:unknown_system` — the system letter has no known frequency pair.

  Empty input yields `{[], []}`. Never raises.

  ## Options

    * `:pairs` — override the band pair per system, e.g.
      `pairs: %{"G" => {:l1, :l2}}`. A system without an override uses its
      standard default pair.
  """
  @spec iono_free_pseudoranges([observation()], [observation()], keyword()) ::
          {[observation()], [{String.t(), drop_reason()}]}
  def iono_free_pseudoranges(band1, band2, opts \\ []) when is_list(band1) and is_list(band2) do
    overrides = Keyword.get(opts, :pairs, %{})
    m1 = Map.new(band1)
    m2 = Map.new(band2)

    ids1 = MapSet.new(Map.keys(m1))
    ids2 = MapSet.new(Map.keys(m2))

    both = MapSet.intersection(ids1, ids2)
    only1 = MapSet.difference(ids1, ids2)
    only2 = MapSet.difference(ids2, ids1)

    {combined, dropped_unknown} =
      both
      |> Enum.sort()
      |> Enum.reduce({[], []}, fn sat, {acc, drops} ->
        case combine_one(sat, Map.fetch!(m1, sat), Map.fetch!(m2, sat), overrides) do
          {:ok, pr_if} -> {[{sat, pr_if} | acc], drops}
          {:error, reason} -> {acc, [{sat, reason} | drops]}
        end
      end)

    dropped =
      dropped_unknown ++
        Enum.map(Enum.sort(only1), fn sat -> {sat, :missing_band2} end) ++
        Enum.map(Enum.sort(only2), fn sat -> {sat, :missing_band1} end)

    {Enum.reverse(combined), Enum.sort(dropped)}
  end

  @doc """
  Convenience: pull two carrier bands from a parsed observation handle and combine
  them into ionosphere-free pseudoranges for one epoch.

  Calls `Orbis.RinexObs.pseudoranges/3` twice — once for each band's code
  preference — then `iono_free_pseudoranges/3`. `epoch` is an epoch index or
  `{{y, mo, d}, {h, mi, s}}` tuple, exactly as `RinexObs.pseudoranges/3` accepts.

  Returns `{:ok, {combined, dropped}}` or `{:error, reason}` (propagated from
  either extraction).

  ## Options

    * `:codes` — a map `%{system => {band1_codes, band2_codes}}` of the
      observation codes to extract for each band, e.g.
      `%{"G" => {["C1C"], ["C2W", "C2L"]}, "E" => {["C1C"], ["C5Q"]}}`. When
      omitted, the standard band-1/band-2 codes for the systems present in the
      file's `observation_codes/1` are used (GPS C1C/C2W, Galileo C1C/C5Q,
      BeiDou C2I/C6I), restricted to those whose codes the file actually carries.
    * `:pairs` — forwarded to `iono_free_pseudoranges/3`.
  """
  @spec iono_free_from_obs(RinexObs.t(), non_neg_integer() | tuple(), keyword()) ::
          {:ok, {[observation()], [{String.t(), drop_reason()}]}} | {:error, term()}
  def iono_free_from_obs(%RinexObs{} = obs, epoch, opts \\ []) do
    code_map = Keyword.get(opts, :codes) || default_obs_codes(obs)

    band1_codes = Map.new(code_map, fn {sys, {b1, _b2}} -> {sys, b1} end)
    band2_codes = Map.new(code_map, fn {sys, {_b1, b2}} -> {sys, b2} end)

    with {:ok, band1} <- RinexObs.pseudoranges(obs, epoch, codes: band1_codes),
         {:ok, band2} <- RinexObs.pseudoranges(obs, epoch, codes: band2_codes) do
      {:ok, iono_free_pseudoranges(band1, band2, Keyword.take(opts, [:pairs]))}
    end
  end

  # --- helpers -------------------------------------------------------------

  defp combine_one(sat, pr1, pr2, overrides) do
    system = String.first(sat)

    with {:ok, {b1, b2}} <- pair_for(system, overrides),
         {:ok, f1} <- frequency(system, b1),
         {:ok, f2} <- frequency(system, b2),
         {:ok, pr_if} <- iono_free(pr1, pr2, f1, f2) do
      {:ok, pr_if}
    else
      _ -> {:error, :unknown_system}
    end
  end

  defp pair_for(system, overrides) do
    case overrides do
      %{^system => pair} -> {:ok, pair}
      _ -> default_pair(system)
    end
  end

  # The standard band-1/band-2 codes per system, kept only for systems whose
  # codes the file actually carries.
  @standard_obs_codes %{
    "G" => {["C1C"], ["C2W", "C2L"]},
    "E" => {["C1C"], ["C5Q"]},
    "C" => {["C2I"], ["C6I"]}
  }

  defp default_obs_codes(%RinexObs{} = obs) do
    available = RinexObs.observation_codes(obs)

    @standard_obs_codes
    |> Enum.filter(fn {sys, {b1, b2}} ->
      file_codes = Map.get(available, sys, [])
      Enum.any?(b1, &(&1 in file_codes)) and Enum.any?(b2, &(&1 in file_codes))
    end)
    |> Map.new()
  end
end
