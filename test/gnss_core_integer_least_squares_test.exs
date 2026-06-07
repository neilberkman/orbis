defmodule Orbis.GNSS.Core.IntegerLeastSquaresTest do
  use ExUnit.Case, async: true

  alias Orbis.GNSS.Core.IntegerLeastSquares
  alias Orbis.GNSS.Core.LinearAlgebra

  test "search matches an independent brute-force ILS scan for a non-rounding optimum" do
    float_cycles = %{"A" => 0.49, "B" => -0.49}

    covariance = [
      [1.0, 0.9],
      [0.9, 1.0]
    ]

    opts = %{radius_cycles: 1, ratio_threshold: 3.0, candidate_limit: 100}

    assert {:ok, fixed, meta} = IntegerLeastSquares.search(float_cycles, covariance, opts)
    assert fixed != coordinate_round(float_cycles)
    assert {fixed, meta.integer_best_score} == brute_force_best(float_cycles, covariance, 1)
    assert meta.integer_status == :not_fixed
    assert meta.ambiguity_search.order == ["A", "B"]
  end

  defp coordinate_round(float_cycles) do
    Map.new(float_cycles, fn {id, value} -> {id, round(value)} end)
  end

  defp brute_force_best(float_cycles, covariance, radius) do
    ids = float_cycles |> Map.keys() |> Enum.sort()
    floats = Enum.map(ids, &Map.fetch!(float_cycles, &1))
    rounded = Enum.map(floats, &round/1)
    {:ok, q_inv} = LinearAlgebra.invert_matrix(covariance)

    candidates =
      for a <- (Enum.at(rounded, 0) - radius)..(Enum.at(rounded, 0) + radius),
          b <- (Enum.at(rounded, 1) - radius)..(Enum.at(rounded, 1) + radius) do
        cycles = [a, b]
        score = quadratic_score(floats, cycles, q_inv)
        {score, Map.new(Enum.zip(ids, cycles))}
      end

    {score, fixed} =
      Enum.min_by(candidates, fn {score, fixed} ->
        {score, Enum.map(ids, &Map.fetch!(fixed, &1))}
      end)

    {fixed, score}
  end

  defp quadratic_score(float_cycles, fixed_cycles, q_inv) do
    deltas =
      fixed_cycles
      |> Enum.zip(float_cycles)
      |> Enum.map(fn {z, a} -> a - z end)

    Enum.reduce(0..1, 0.0, fn i, acc ->
      Enum.reduce(0..1, acc, fn j, inner ->
        inner + Enum.at(deltas, i) * (q_inv |> Enum.at(i) |> Enum.at(j)) * Enum.at(deltas, j)
      end)
    end)
  end
end
