defmodule Orbis.GNSS.Core.IntegerLeastSquares do
  @moduledoc false

  import Orbis.GNSS.Core.LinearAlgebra,
    only: [
      identity_matrix: 1,
      invert_matrix: 1,
      matmul: 2,
      matvec_transpose: 2,
      solve_linear: 2,
      transpose: 1
    ]

  @spec search(%{String.t() => number()}, [[number()]], map()) ::
          {:ok, %{String.t() => integer()}, map()} | {:error, term()}
  def search(float_cycles_by_id, covariance_cycles, opts)
      when is_map(float_cycles_by_id) and is_list(covariance_cycles) do
    ids = float_cycles_by_id |> Map.keys() |> Enum.sort()
    q_cycles = symmetrize_matrix(covariance_cycles)

    case invert_matrix(q_cycles) do
      {:ok, q_inv} ->
        q_inv = symmetrize_matrix(q_inv)

        case lambda_decorrelate(q_cycles) do
          {:ok, q_decorrelated, z_transform} ->
            with {:ok, candidates, evaluated} <-
                   lambda_sphere_search(
                     float_cycles_by_id,
                     q_inv,
                     q_decorrelated,
                     z_transform,
                     opts
                   ),
                 {:ok, fixed_cycles, fixed_meta} <-
                   lambda_search_result(candidates, evaluated, opts) do
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

      {:error, _} = err ->
        err
    end
  end

  defp lambda_search_result(candidates, evaluated, opts) do
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
       integer_method: :lambda,
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

  defp lambda_sphere_search(float_cycles_by_id, q_inv, q_decorrelated, z_transform, opts) do
    ids = float_cycles_by_id |> Map.keys() |> Enum.sort()
    float_cycles = Enum.map(ids, &Map.fetch!(float_cycles_by_id, &1))
    bound = lambda_initial_bound(float_cycles, q_inv, opts.radius_cycles)
    decorrelated_float_cycles = matvec_transpose(z_transform, float_cycles)

    with {:ok, lower, diagonal} <- ldl_decompose(q_decorrelated) do
      case lambda_recurse(
             lower,
             diagonal,
             decorrelated_float_cycles,
             0,
             [],
             empty_solution(length(float_cycles)),
             0.0,
             bound,
             [],
             0,
             opts.candidate_limit
           ) do
        {:error, _} = err ->
          err

        {:ok, [], evaluated, _bound} ->
          {:error, {:no_integer_candidates, evaluated}}

        {:ok, candidates, evaluated, _bound} ->
          sorted =
            candidates
            |> Enum.map(fn {_score, decorrelated_cycles} ->
              original_cycles = original_integer_cycles(z_transform, decorrelated_cycles)

              {
                quadratic_integer_score(float_cycles, original_cycles, q_inv),
                Map.new(Enum.zip(ids, original_cycles))
              }
            end)
            |> Enum.sort_by(fn {score, fixed_cycles} ->
              {score, Enum.map(ids, &Map.fetch!(fixed_cycles, &1))}
            end)

          {:ok, sorted, evaluated}
      end
    end
  end

  defp lambda_decorrelate(q_cycles) do
    n = length(q_cycles)
    max_steps = max(100, n * n * 20)

    lambda_reduce(q_cycles, identity_matrix(n), n - 2, n - 2, 0, max_steps)
  end

  defp lambda_reduce(q, z, j, _k, _steps, _max_steps) when j < 0, do: {:ok, q, z}

  defp lambda_reduce(_q, _z, _j, _k, steps, max_steps) when steps > max_steps,
    do: {:error, :singular_geometry}

  defp lambda_reduce(q, z, j, k, steps, max_steps) do
    n = length(q)

    with {:ok, q, z} <- maybe_gauss_reduce(q, z, j, k),
         {:ok, lower, diagonal} <- ldl_decompose(q) do
      lij = lower |> Enum.at(j + 1) |> Enum.at(j)
      delta = Enum.at(diagonal, j) + lij * lij * Enum.at(diagonal, j + 1)

      if delta + 1.0e-12 < Enum.at(diagonal, j + 1) do
        t = permutation_transform(n, j, j + 1)
        q = transform_covariance(q, t)
        z = matmul(z, t)
        lambda_reduce(q, z, n - 2, j, steps + 1, max_steps)
      else
        lambda_reduce(q, z, j - 1, k, steps + 1, max_steps)
      end
    end
  end

  defp maybe_gauss_reduce(q, z, j, k) when j > k, do: {:ok, q, z}

  defp maybe_gauss_reduce(q, z, j, _k) do
    n = length(q)

    Enum.reduce_while((j + 1)..(n - 1), {:ok, q, z}, fn i, {:ok, q_acc, z_acc} ->
      case ldl_decompose(q_acc) do
        {:ok, lower, _diagonal} ->
          mu = lower |> Enum.at(i) |> Enum.at(j) |> round()

          if mu == 0 do
            {:cont, {:ok, q_acc, z_acc}}
          else
            t = integer_gauss_transform(n, i, j, mu)
            {:cont, {:ok, transform_covariance(q_acc, t), matmul(z_acc, t)}}
          end

        {:error, _} = err ->
          {:halt, err}
      end
    end)
  end

  defp integer_gauss_transform(n, i, j, mu) do
    identity_matrix(n)
    |> List.update_at(i, fn row -> List.replace_at(row, j, -mu / 1.0) end)
  end

  defp permutation_transform(n, i, j) do
    identity_matrix(n)
    |> swap_columns(i, j)
  end

  defp transform_covariance(q, transform) do
    transform
    |> transpose()
    |> matmul(matmul(q, transform))
    |> symmetrize_matrix()
  end

  defp original_integer_cycles(z_transform, decorrelated_cycles) do
    {:ok, original} =
      solve_linear(transpose(z_transform), Enum.map(decorrelated_cycles, &(&1 / 1.0)))

    Enum.map(original, &round/1)
  end

  defp lambda_initial_bound(float_cycles, q_inv, radius) do
    n = length(float_cycles)

    rounded = Enum.map(float_cycles, &round/1)

    alternatives =
      if radius == 0 do
        []
      else
        for idx <- 0..(n - 1),
            step <- 1..radius,
            sign <- [-1, 1] do
          List.update_at(rounded, idx, &(&1 + sign * step))
        end
      end

    scores =
      [rounded | alternatives]
      |> Enum.map(&quadratic_integer_score(float_cycles, &1, q_inv))
      |> Enum.sort()

    bound =
      case scores do
        [_best, second | _] -> second * (1.0 + 1.0e-12)
        [best] -> best + 1.0
      end

    if bound > 0.0, do: bound, else: 1.0
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

  defp empty_solution(n), do: List.duplicate(nil, n)

  defp lambda_recurse(
         _lower,
         _diagonal,
         _float_cycles,
         k,
         _y_prefix,
         z,
         dist,
         bound,
         candidates,
         evaluated,
         limit
       )
       when k == length(z) do
    evaluated = evaluated + 1

    if evaluated > limit do
      {:error, {:too_many_integer_candidates, evaluated, limit}}
    else
      candidates = lambda_top_two([{dist, z} | candidates])
      bound = lambda_live_bound(candidates, bound)
      {:ok, candidates, evaluated, bound}
    end
  end

  defp lambda_recurse(
         lower,
         diagonal,
         float_cycles,
         k,
         y_prefix,
         z,
         dist,
         bound,
         candidates,
         evaluated,
         limit
       ) do
    dk = Enum.at(diagonal, k)
    center = Enum.at(float_cycles, k) + lambda_forward_offset(lower, y_prefix, k)
    remaining = bound - dist

    if remaining < 0.0 do
      {:ok, candidates, evaluated, bound}
    else
      span = :math.sqrt(remaining * dk)
      low = Float.floor(center - span) |> trunc()
      high = Float.ceil(center + span) |> trunc()

      integers_near(center, low, high)
      |> Enum.reduce_while({:ok, candidates, evaluated, bound}, fn value,
                                                                   {:ok, acc_candidates,
                                                                    acc_evaluated, acc_bound} ->
        term = center - value
        next_dist = dist + term * term / dk

        if next_dist > acc_bound do
          {:cont, {:ok, acc_candidates, acc_evaluated, acc_bound}}
        else
          next_z = List.replace_at(z, k, value)
          next_y_prefix = y_prefix ++ [value - center]

          case lambda_recurse(
                 lower,
                 diagonal,
                 float_cycles,
                 k + 1,
                 next_y_prefix,
                 next_z,
                 next_dist,
                 acc_bound,
                 acc_candidates,
                 acc_evaluated,
                 limit
               ) do
            {:ok, next_candidates, next_evaluated, next_bound} ->
              {:cont, {:ok, next_candidates, next_evaluated, next_bound}}

            {:error, _} = err ->
              {:halt, err}
          end
        end
      end)
    end
  end

  defp lambda_top_two(candidates) do
    candidates
    |> Enum.sort_by(fn {score, cycles} -> {score, cycles} end)
    |> Enum.take(2)
  end

  defp lambda_live_bound([{_best_score, _best_cycles}, {second_score, _second_cycles}], _bound),
    do: second_score

  defp lambda_live_bound(_candidates, bound), do: bound

  defp lambda_forward_offset(_lower, [], _k), do: 0.0

  defp lambda_forward_offset(lower, y_prefix, k) do
    row = Enum.at(lower, k)

    y_prefix
    |> Enum.with_index()
    |> Enum.reduce(0.0, fn {y_j, j}, acc ->
      acc + Enum.at(row, j) * y_j
    end)
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

  defp swap_columns(matrix, i, j) do
    Enum.map(matrix, fn row ->
      vi = Enum.at(row, i)
      vj = Enum.at(row, j)

      row
      |> List.replace_at(i, vj)
      |> List.replace_at(j, vi)
    end)
  end

  defp ldl_decompose(a) do
    n = length(a)

    0..(n - 1)
    |> Enum.reduce_while({:ok, [], []}, fn i, {:ok, rows, diagonal} ->
      case ldl_row(a, rows, diagonal, i, n) do
        {:ok, row, d_i} -> {:cont, {:ok, rows ++ [row], diagonal ++ [d_i]}}
        {:error, _} = err -> {:halt, err}
      end
    end)
    |> case do
      {:ok, rows, diagonal} -> {:ok, rows, diagonal}
      {:error, _} = err -> err
    end
  end

  defp ldl_row(a, rows, diagonal, i, n) do
    l_values =
      if i == 0 do
        []
      else
        Enum.reduce(0..(i - 1), [], fn j, acc ->
          sum = (a |> Enum.at(i) |> Enum.at(j)) - ldl_cross_sum(rows, diagonal, acc, j)
          d_j = Enum.at(diagonal, j)
          acc ++ [sum / d_j]
        end)
      end

    d_i = (a |> Enum.at(i) |> Enum.at(i)) - ldl_diag_sum(l_values, diagonal)

    if not is_number(d_i) or d_i <= 0.0 do
      {:error, :singular_geometry}
    else
      {:ok, l_values ++ [1.0] ++ List.duplicate(0.0, n - i - 1), d_i}
    end
  end

  defp ldl_cross_sum(_rows, _diagonal, _current, 0), do: 0.0

  defp ldl_cross_sum(rows, diagonal, current, j) do
    Enum.reduce(0..(j - 1), 0.0, fn k, acc ->
      lik = Enum.at(current, k)
      ljk = rows |> Enum.at(j) |> Enum.at(k)
      acc + lik * ljk * Enum.at(diagonal, k)
    end)
  end

  defp ldl_diag_sum([], _diagonal), do: 0.0

  defp ldl_diag_sum(l_values, diagonal) do
    l_values
    |> Enum.zip(diagonal)
    |> Enum.reduce(0.0, fn {lik, d_k}, acc -> acc + lik * lik * d_k end)
  end
end
