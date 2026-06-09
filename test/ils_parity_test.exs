defmodule Orbis.GNSS.Core.IlsParityTest do
  @moduledoc """
  Proves the Rust bounded-ILS kernel (`IntegerLeastSquares.bounded_search/3`, via the NIF)
  is bit-identical to the pure-Elixir reference (`reference_search/3`): same fixed
  integers, status, ratio, scores, candidate count, and symmetrized covariance /
  inverse. Same IEEE-754 ops in the same order (no FMA) → exact equality, not a
  tolerance.
  """
  use ExUnit.Case, async: true

  alias Orbis.GNSS.Core.IntegerLeastSquares, as: ILS

  @opts %{radius_cycles: 1, candidate_limit: 200_000, ratio_threshold: 3.0}

  defp floats_by_id(floats) do
    floats
    |> Enum.with_index()
    |> Map.new(fn {f, i} -> {"A#{String.pad_leading(Integer.to_string(i), 2, "0")}", f} end)
  end

  # Symmetric, diagonally-dominant (=> SPD => invertible) covariance with a
  # tunable off-diagonal correlation structure.
  defp dd_cov(n, off) do
    for i <- 0..(n - 1) do
      for j <- 0..(n - 1) do
        if i == j, do: 1.0, else: off / (1.0 + abs(i - j))
      end
    end
  end

  defp assert_parity(floats, cov, opts) do
    fbi = floats_by_id(floats)

    case {ILS.reference_search(fbi, cov, opts), ILS.bounded_search(fbi, cov, opts)} do
      {{:ok, rf, rm}, {:ok, nf, nm}} ->
        assert rf == nf, "fixed integers differ"
        assert rm.integer_status == nm.integer_status
        assert rm.integer_ratio == nm.integer_ratio
        assert rm.integer_best_score == nm.integer_best_score
        assert rm.integer_second_best_score == nm.integer_second_best_score
        assert rm.integer_candidates == nm.integer_candidates

        assert rm.ambiguity_search.covariance_cycles ==
                 nm.ambiguity_search.covariance_cycles

        assert rm.ambiguity_search.covariance_inverse_cycles ==
                 nm.ambiguity_search.covariance_inverse_cycles

      {ref_err, nif_err} ->
        assert ref_err == nif_err, "error parity: #{inspect({ref_err, nif_err})}"
    end
  end

  test "bit-identical across dimensions, correlation structures, and float offsets" do
    for n <- 1..6,
        off <- [0.0, 0.2, -0.15, 0.3],
        offset <- [0.02, 0.5, 0.27, 0.48, -0.49] do
      cov = dd_cov(n, off)

      floats =
        for i <- 0..(n - 1) do
          i + 1 + if rem(i, 2) == 0, do: offset, else: -offset
        end

      assert_parity(floats, cov, @opts)
    end
  end

  test "error-path parity: singular covariance" do
    singular = [[1.0, 2.0], [2.0, 4.0]]
    assert_parity([0.1, 0.2], singular, @opts)
  end

  test "error-path parity: candidate-limit exceeded" do
    cov = dd_cov(3, 0.0)
    tight = %{@opts | candidate_limit: 5}
    assert_parity([0.0, 0.0, 0.0], cov, tight)
  end
end
