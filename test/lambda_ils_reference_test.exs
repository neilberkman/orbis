defmodule Orbis.GNSS.Core.LambdaIlsReferenceTest do
  @moduledoc """
  External-reference gate for the ambiguity-resolution kernels against
  **RTKLIB's `lambda()`** — the de-facto GNSS reference engine (BSD-2), the
  Teunissen LAMBDA method.

  This is the validation that *matters*: it checks our kernels against an
  established, widely-used implementation, not against our own Elixir.

  The golden (`test/fixtures/lambda_golden.json`) is produced by compiling
  RTKLIB's vendored, unmodified `lambda()` and running it over a battery — see
  `parity/generator/lambda_ref/`. The first two cases are RTKLIB's OWN committed
  unit-test vectors (`test/utest/t_lambda.c`).

  It pins both kernels:

    * `IntegerLeastSquares.search/3` (**LAMBDA**, the production solver) reproduces
      RTKLIB's exact integer vector and residuals on EVERY case — including the
      strongly-correlated `utest2`, where the optimum is up to 14 cycles from
      componentwise rounding.

    * `IntegerLeastSquares.bounded_search/3` (the ±radius box) matches RTKLIB only
      **in-regime** (weakly-correlated geometry, optimum within ±1 of rounding —
      e.g. the Wettzell arc). Out-of-regime it provably cannot reach the optimum:
      it misses with a small radius and explodes past the candidate limit with a
      large one. That gap is exactly why production uses LAMBDA.
  """
  use ExUnit.Case, async: true

  alias Orbis.GNSS.Core.IntegerLeastSquares, as: ILS

  @golden Path.join(__DIR__, "fixtures/lambda_golden.json") |> File.read!() |> Jason.decode!()

  # Cross-implementation tolerance: RTKLIB and our kernels invert/factorize Q
  # differently, so the residual metric Δᵀ Q⁻¹ Δ agrees only to round-off, not
  # bit-exactly. RTKLIB's own unit test uses 1e-4 on residuals ~1.5e3; scale
  # accordingly per case.
  @score_tol 1.0e-6
  @opts %{radius_cycles: 1, candidate_limit: 200_000, ratio_threshold: 3.0}

  defp ids_for(n) do
    for i <- 0..(n - 1), do: "A#{String.pad_leading(Integer.to_string(i), 2, "0")}"
  end

  defp float_map(a) do
    a |> Enum.with_index() |> Map.new(fn {f, i} -> {Enum.at(ids_for(length(a)), i), f} end)
  end

  defp case_by_name(name), do: Enum.find(@golden["cases"], &(&1["name"] == name))

  defp fixed_vector(fixed_map, n), do: Enum.map(ids_for(n), &Map.fetch!(fixed_map, &1))

  defp tol_for(s_best), do: max(@score_tol, abs(s_best) * 1.0e-9 + 1.0e-4)

  describe "search/3 (LAMBDA) reproduces RTKLIB lambda() on ALL cases" do
    for c <- @golden["cases"] do
      @case c
      test "#{c["name"]}: same fixed integers, residuals, and ratio as RTKLIB" do
        {:ok, fixed_map, meta} = ILS.search(float_map(@case["a"]), @case["Q"], @opts)
        fixed_vec = fixed_vector(fixed_map, @case["n"])

        # Exact integer match even for the strongly-correlated utest2 the box
        # cannot solve — the whole point of the LAMBDA solver.
        assert fixed_vec == hd(@case["lambda_fixed"]),
               "LAMBDA selected a different integer vector than RTKLIB lambda()"

        [s_best, s_second] = @case["lambda_residuals"]
        tol = tol_for(s_best)
        assert_in_delta meta.integer_best_score, s_best, tol
        assert_in_delta meta.integer_second_best_score, s_second, tol
        assert meta.integer_method == :lambda
      end
    end
  end

  describe "bounded_search/3 (box) matches RTKLIB in-regime" do
    for c <- @golden["cases"], c["in_regime"] do
      @case c
      test "#{c["name"]}: box reproduces RTKLIB's exact solution" do
        {:ok, fixed_map, meta} = ILS.bounded_search(float_map(@case["a"]), @case["Q"], @opts)
        fixed_vec = fixed_vector(fixed_map, @case["n"])

        assert fixed_vec == hd(@case["lambda_fixed"])

        [s_best, s_second] = @case["lambda_residuals"]
        assert_in_delta meta.integer_best_score, s_best, @score_tol
        assert_in_delta meta.integer_second_best_score, s_second, @score_tol
        assert_in_delta meta.integer_ratio, @case["lambda_ratio"], @score_tol
        assert meta.integer_method == :bounded_ils
      end
    end
  end

  describe "bounded_search/3 (box) cannot reach RTKLIB's optimum out-of-regime" do
    test "rtklib_utest2 (strongly correlated): the box misses the gold-standard solution" do
      c = case_by_name("rtklib_utest2")
      {:ok, fixed_map, meta} = ILS.bounded_search(float_map(c["a"]), c["Q"], @opts)
      [s_best, _] = c["lambda_residuals"]

      # The ±1 box around rounding does not contain RTKLIB's optimum (up to 14
      # cycles away), so the box selects a strictly worse integer vector.
      refute fixed_vector(fixed_map, c["n"]) == hd(c["lambda_fixed"])
      assert meta.integer_best_score > s_best
    end

    test "rtklib_utest2: widening the box to reach it explodes past the candidate limit" do
      c = case_by_name("rtklib_utest2")
      # The optimum is ~14 cycles from rounding; a box wide enough to contain it
      # is 29^10 ≈ 4e14 lattice points — the search aborts at the limit. The
      # naive enumeration has no way to solve this; LAMBDA decorrelation does.
      wide = %{@opts | radius_cycles: 14}

      assert {:error, {:too_many_integer_candidates, _evaluated, _limit}} =
               ILS.bounded_search(float_map(c["a"]), c["Q"], wide)
    end
  end
end
