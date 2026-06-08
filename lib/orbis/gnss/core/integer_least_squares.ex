defmodule Orbis.GNSS.Core.IntegerLeastSquares do
  @moduledoc false

  import Orbis.GNSS.Core.LinearAlgebra, only: [invert_matrix: 1]

  @spec search(%{String.t() => number()}, [[number()]], map()) ::
          {:ok, %{String.t() => integer()}, map()} | {:error, term()}
  def search(float_cycles_by_id, covariance_cycles, opts)
      when is_map(float_cycles_by_id) and is_list(covariance_cycles) do
    ids = float_cycles_by_id |> Map.keys() |> Enum.sort()
    q_cycles = symmetrize_matrix(covariance_cycles)

    case invert_matrix(q_cycles) do
      {:ok, q_inv} ->
        q_inv = symmetrize_matrix(q_inv)

        with {:ok, candidates, evaluated} <-
               bounded_integer_search(float_cycles_by_id, q_inv, opts),
             {:ok, fixed_cycles, fixed_meta} <-
               integer_search_result(candidates, evaluated, opts) do
          {:ok, fixed_cycles,
           Map.put(fixed_meta, :ambiguity_search, %{
             order: ids,
             float_cycles: float_cycles_by_id,
             covariance_cycles: q_cycles,
             covariance_inverse_cycles: q_inv
           })}
        end

      {:error, _} = err ->
        err
    end
  end

  defp integer_search_result(candidates, evaluated, opts) do
    [{best_score, fixed_cycles} | rest] = candidates

    second_best_score =
      case rest do
        [{score, _} | _] -> score
        [] -> nil
      end

    ratio = integer_ratio(best_score, second_best_score)
    status = if integer_ratio_pass?(ratio, opts.ratio_threshold), do: :fixed, else: :not_fixed

    {:ok, fixed_cycles,
     %{
       integer_status: status,
       integer_method: :bounded_ils,
       integer_ratio: ratio,
       integer_best_score: best_score,
       integer_second_best_score: second_best_score,
       integer_candidates: evaluated
     }}
  end

  defp symmetrize_matrix(matrix) do
    n = length(matrix)

    for i <- 0..(n - 1) do
      for j <- 0..(n - 1) do
        ((matrix |> Enum.at(i) |> Enum.at(j)) + (matrix |> Enum.at(j) |> Enum.at(i))) / 2.0
      end
    end
  end

  defp bounded_integer_search(float_cycles_by_id, q_inv, opts) do
    ids = float_cycles_by_id |> Map.keys() |> Enum.sort()
    float_cycles = Enum.map(ids, &Map.fetch!(float_cycles_by_id, &1))

    ranges =
      Enum.map(float_cycles, fn float ->
        center = round(float)
        integers_near(float, center - opts.radius_cycles, center + opts.radius_cycles)
      end)

    case bounded_integer_recurse(
           ids,
           float_cycles,
           q_inv,
           ranges,
           [],
           [],
           0,
           opts.candidate_limit
         ) do
      {:ok, [], evaluated} ->
        {:error, {:no_integer_candidates, evaluated}}

      {:ok, candidates, evaluated} ->
        {:ok, bounded_candidates_to_maps(ids, candidates), evaluated}

      {:error, _} = err ->
        err
    end
  end

  defp bounded_integer_recurse(
         _ids,
         float_cycles,
         q_inv,
         [],
         cycles,
         candidates,
         evaluated,
         limit
       ) do
    evaluated = evaluated + 1

    if evaluated > limit do
      {:error, {:too_many_integer_candidates, evaluated, limit}}
    else
      score = quadratic_integer_score(float_cycles, cycles, q_inv)
      candidates = integer_top_two([{score, cycles} | candidates])
      {:ok, candidates, evaluated}
    end
  end

  defp bounded_integer_recurse(
         ids,
         float_cycles,
         q_inv,
         [range | rest],
         cycles,
         candidates,
         evaluated,
         limit
       ) do
    Enum.reduce_while(range, {:ok, candidates, evaluated}, fn value,
                                                              {:ok, acc_candidates, acc_evaluated} ->
      case bounded_integer_recurse(
             ids,
             float_cycles,
             q_inv,
             rest,
             cycles ++ [value],
             acc_candidates,
             acc_evaluated,
             limit
           ) do
        {:ok, next_candidates, next_evaluated} ->
          {:cont, {:ok, next_candidates, next_evaluated}}

        {:error, _} = err ->
          {:halt, err}
      end
    end)
  end

  defp bounded_candidates_to_maps(ids, candidates) do
    Enum.map(candidates, fn {score, cycles} ->
      {score, Map.new(Enum.zip(ids, cycles))}
    end)
  end

  defp quadratic_integer_score(float_cycles, fixed_cycles, q_inv) do
    deltas =
      fixed_cycles
      |> Enum.zip(float_cycles)
      |> Enum.map(fn {z, a} -> a - z end)

    n = length(deltas)

    Enum.reduce(0..(n - 1), 0.0, fn i, acc ->
      Enum.reduce(0..(n - 1), acc, fn j, inner ->
        inner + Enum.at(deltas, i) * (q_inv |> Enum.at(i) |> Enum.at(j)) * Enum.at(deltas, j)
      end)
    end)
  end

  defp integer_top_two(candidates) do
    candidates
    |> Enum.sort_by(fn {score, cycles} -> {score, cycles} end)
    |> Enum.take(2)
  end

  defp integers_near(center, low, high) do
    low..high
    |> Enum.to_list()
    |> Enum.sort_by(fn value -> {abs(value - center), value} end)
  end

  # A ratio test requires both a best and a runner-up candidate. If the search
  # bound yields only one lattice point, the fix is not independently validated;
  # report a non-passing ratio rather than treating the missing runner-up as
  # infinite confidence.
  defp integer_ratio(_best_score, nil), do: 0.0

  defp integer_ratio(best_score, second_best_score) do
    cond do
      best_score == 0.0 and second_best_score > 0.0 -> :infinity
      best_score == 0.0 -> 0.0
      true -> second_best_score / best_score
    end
  end

  defp integer_ratio_pass?(:infinity, _threshold), do: true
  defp integer_ratio_pass?(ratio, threshold), do: ratio >= threshold
end
