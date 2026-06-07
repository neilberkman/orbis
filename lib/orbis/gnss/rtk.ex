defmodule Orbis.GNSS.RTK do
  @moduledoc """
  RTK-facing carrier/code double-difference primitives.

  A base receiver and a rover receiver observing the same satellites have
  receiver-clock terms that differ by station but are common to every satellite.
  A *single difference* subtracts base from rover for the same satellite; a
  *double difference* subtracts a reference satellite's single difference:

      DD_s = (rover_s - base_s) - (rover_ref - base_ref)

  The receiver clocks cancel in the second subtraction. Satellite-clock,
  ephemeris, and short-baseline atmosphere errors that are common between base
  and rover also cancel in the receiver single difference. The remaining
  carrier-phase double differences are the measurement surface used by RTK
  baseline estimation and integer ambiguity fixing.

  This module deliberately stops at the normalized double-difference
  observations; baseline solving and integer fixing live in later RTK layers.

  ## Example

      iex> base = [
      ...>   {"G01", 20_100.0, 20_103.0},
      ...>   {"G02", 21_105.0, 21_110.0}
      ...> ]
      iex> rover = [
      ...>   {"G01", 20_040.0, 20_044.0},
      ...>   {"G02", 21_060.0, 21_066.0}
      ...> ]
      iex> {:ok, result} = Orbis.GNSS.RTK.double_differences(base, rover, reference_satellite_id: "G01")
      iex> result.double_differences
      [%{satellite_id: "G02", reference_satellite_id: "G01", code_m: 15.0, phase_m: 17.0}]
  """

  @typedoc "Code and carrier-phase observation in metres."
  @type observation ::
          %{satellite_id: String.t(), code_m: number(), phase_m: number()}
          | {String.t(), number(), number()}

  @typedoc "One non-reference satellite's double-difference observation."
  @type double_difference :: %{
          satellite_id: String.t(),
          reference_satellite_id: String.t(),
          code_m: float(),
          phase_m: float()
        }

  @typedoc "Double-difference result with deterministic satellite ordering."
  @type result :: %{
          reference_satellite_id: String.t(),
          double_differences: [double_difference()],
          dropped_sats: [String.t()]
        }

  @doc """
  Build code and carrier-phase double differences from base and rover observations.

  Observations can be maps with `:satellite_id`, `:code_m`, and `:phase_m`, or
  `{satellite_id, code_m, phase_m}` tuples. Satellites are paired by id; any
  satellite not present at both receivers is reported in `:dropped_sats`.

  Options:

    * `:reference_satellite_id` - reference satellite for the second
      difference. When omitted, the lexicographically first common satellite is
      selected deterministically.

  Returns `{:ok, result}` or a tagged error. At least two common satellites are
  required so one non-reference double difference can be produced.
  """
  @spec double_differences([observation()], [observation()], keyword()) ::
          {:ok, result()} | {:error, term()}
  def double_differences(base_observations, rover_observations, opts \\ [])

  def double_differences(base_observations, rover_observations, opts)
      when is_list(base_observations) and is_list(rover_observations) do
    with {:ok, base} <- normalize_observations(base_observations, :invalid_base_observations),
         {:ok, rover} <- normalize_observations(rover_observations, :invalid_rover_observations),
         {:ok, common, dropped} <- common_observations(base, rover),
         {:ok, reference_sat} <- reference_satellite(common, opts) do
      ref_base = Map.fetch!(base, reference_sat)
      ref_rover = Map.fetch!(rover, reference_sat)
      ref_code_sd = ref_rover.code_m - ref_base.code_m
      ref_phase_sd = ref_rover.phase_m - ref_base.phase_m

      dds =
        common
        |> Enum.reject(&(&1 == reference_sat))
        |> Enum.map(fn sat ->
          base_obs = Map.fetch!(base, sat)
          rover_obs = Map.fetch!(rover, sat)

          %{
            satellite_id: sat,
            reference_satellite_id: reference_sat,
            code_m: rover_obs.code_m - base_obs.code_m - ref_code_sd,
            phase_m: rover_obs.phase_m - base_obs.phase_m - ref_phase_sd
          }
        end)

      {:ok,
       %{
         reference_satellite_id: reference_sat,
         double_differences: dds,
         dropped_sats: dropped
       }}
    end
  end

  def double_differences(_base_observations, _rover_observations, _opts),
    do: {:error, :invalid_observations}

  defp normalize_observations(observations, error_tag) do
    observations
    |> Enum.reduce_while({:ok, %{}}, fn observation, {:ok, acc} ->
      case normalize_observation(observation) do
        {:ok, %{satellite_id: sat} = obs} ->
          if Map.has_key?(acc, sat) do
            {:halt, {:error, {:duplicate_observation, sat}}}
          else
            {:cont, {:ok, Map.put(acc, sat, obs)}}
          end

        {:error, _} ->
          {:halt, {:error, {error_tag, observation}}}
      end
    end)
  end

  defp normalize_observation(%{satellite_id: sat, code_m: code, phase_m: phase})
       when is_binary(sat) and is_number(code) and is_number(phase) do
    {:ok, %{satellite_id: sat, code_m: code / 1.0, phase_m: phase / 1.0}}
  end

  defp normalize_observation({sat, code, phase})
       when is_binary(sat) and is_number(code) and is_number(phase) do
    {:ok, %{satellite_id: sat, code_m: code / 1.0, phase_m: phase / 1.0}}
  end

  defp normalize_observation(_observation), do: {:error, :invalid_observation}

  defp common_observations(base, rover) do
    base_sats = base |> Map.keys() |> MapSet.new()
    rover_sats = rover |> Map.keys() |> MapSet.new()

    common =
      base_sats
      |> MapSet.intersection(rover_sats)
      |> MapSet.to_list()
      |> Enum.sort()

    if length(common) < 2 do
      {:error, {:too_few_common_satellites, length(common), 2}}
    else
      dropped =
        base_sats
        |> MapSet.union(rover_sats)
        |> MapSet.difference(MapSet.new(common))
        |> MapSet.to_list()
        |> Enum.sort()

      {:ok, common, dropped}
    end
  end

  defp reference_satellite(common, opts) do
    case Keyword.get(opts, :reference_satellite_id) do
      nil ->
        {:ok, hd(common)}

      sat when is_binary(sat) ->
        if sat in common do
          {:ok, sat}
        else
          {:error, {:reference_satellite_missing, sat}}
        end

      _other ->
        {:error, {:invalid_option, :reference_satellite_id}}
    end
  end
end
