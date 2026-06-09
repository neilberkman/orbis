defmodule Orbis.GNSS.Core.LambdaIlsReferenceTest do
  @moduledoc """
  External-reference gate for the bounded integer-least-squares kernel
  (`IntegerLeastSquares.search/3`) against **RTKLIB's `lambda()`** — the de-facto
  GNSS reference engine (BSD-2), the Teunissen LAMBDA method.

  This is the validation that *matters*: it checks our kernel against an
  established, widely-used implementation, not against our own Elixir (that is
  what `ils_parity_test.exs` does, and that is only a no-regression check).

  The golden (`test/fixtures/lambda_golden.json`) is produced by compiling
  RTKLIB's vendored, unmodified `lambda()` and running it over a battery — see
  `parity/generator/lambda_ref/`. The first two cases are RTKLIB's OWN committed
  unit-test vectors (`test/utest/t_lambda.c`).

  What it establishes:

    * **In-regime** (weakly-correlated geometry, where the ILS solution lies
      within ±1 cycle of componentwise rounding — the regime the kernel is
      actually used in, e.g. the Wettzell arc): our bounded ±radius search finds
      the EXACT same integer vector RTKLIB does, with matching residuals/ratio.

    * **Out-of-regime** (strongly-correlated geometry — RTKLIB's `utest2`, where
      the ILS optimum is up to 14 cycles from rounding): the naive box search
      provably CANNOT reproduce the gold-standard solution — it either misses the
      optimum (small radius) or explodes combinatorially (large radius). This is
      exactly the gap that LAMBDA's Z-transform decorrelation closes, and is the
      documented motivation for the parked decorrelation work.
  """
  use ExUnit.Case, async: true

  alias Orbis.GNSS.Core.IntegerLeastSquares, as: ILS

  @golden Path.join(__DIR__, "fixtures/lambda_golden.json") |> File.read!() |> Jason.decode!()

  # Cross-implementation tolerance: RTKLIB inverts Q via its own LU (ludcmp);
  # we invert via Gaussian elimination with partial pivoting. The residual
  # metric Δᵀ Q⁻¹ Δ therefore agrees only to floating-point round-off, not
  # bit-exactly. RTKLIB's own unit test uses 1e-4 on residuals ~3.5e3; we hold
  # the tighter 1e-6 (our largest residual is ~3.5).
  @score_tol 1.0e-6
  @opts %{radius_cycles: 1, candidate_limit: 200_000, ratio_threshold: 3.0}

  defp ids_for(n) do
    for i <- 0..(n - 1), do: "A#{String.pad_leading(Integer.to_string(i), 2, "0")}"
  end

  defp float_map(a) do
    a |> Enum.with_index() |> Map.new(fn {f, i} -> {Enum.at(ids_for(length(a)), i), f} end)
  end

  defp case_by_name(name), do: Enum.find(@golden["cases"], &(&1["name"] == name))

  defp run(c) do
    ids = ids_for(c["n"])
    {:ok, fixed_map, meta} = ILS.search(float_map(c["a"]), c["Q"], @opts)
    fixed_vec = Enum.map(ids, &Map.fetch!(fixed_map, &1))
    {fixed_vec, meta}
  end

  describe "in-regime: orbis bounded-ILS matches RTKLIB lambda() exactly" do
    for c <- @golden["cases"], c["in_regime"] do
      @case c
      test "#{c["name"]}: same fixed integers, residuals, and ratio as RTKLIB" do
        {fixed_vec, meta} = run(@case)

        # The integer vector is the hard assertion — exact, no tolerance.
        assert fixed_vec == @case["lambda_fixed"] |> hd(),
               "orbis selected a different integer vector than RTKLIB lambda()"

        [s_best, s_second] = @case["lambda_residuals"]
        assert_in_delta meta.integer_best_score, s_best, @score_tol
        assert_in_delta meta.integer_second_best_score, s_second, @score_tol
        assert_in_delta meta.integer_ratio, @case["lambda_ratio"], @score_tol
      end
    end
  end

  describe "out-of-regime: the naive box search cannot reach RTKLIB's optimum" do
    test "rtklib_utest2 (strongly correlated): orbis misses the gold-standard solution" do
      c = case_by_name("rtklib_utest2")
      {fixed_vec, meta} = run(c)
      rtklib_best = hd(c["lambda_fixed"])
      [s_best, _] = c["lambda_residuals"]

      # The ±1 box around rounding does not contain RTKLIB's optimum (which is
      # up to 14 cycles away), so we select a strictly worse integer vector.
      refute fixed_vec == rtklib_best
      assert meta.integer_best_score > s_best
    end

    test "rtklib_utest2: widening the box to reach it explodes past the candidate limit" do
      c = case_by_name("rtklib_utest2")
      # The optimum is ~14 cycles from rounding; a box wide enough to contain it
      # is 29^10 ≈ 4e14 lattice points — the search aborts at the limit. The
      # naive enumeration has no way to solve this case; LAMBDA decorrelation does.
      wide = %{@opts | radius_cycles: 14}

      assert {:error, {:too_many_integer_candidates, _evaluated, _limit}} =
               ILS.search(float_map(c["a"]), c["Q"], wide)
    end
  end
end
