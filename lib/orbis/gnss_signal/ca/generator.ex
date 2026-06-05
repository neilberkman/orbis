defmodule Orbis.GnssSignal.CA.Generator do
  @moduledoc false

  # Internal C/A code generator. Clocks the two IS-GPS-200 LFSRs (G1 and G2)
  # to produce one 1023-chip Gold code for a given G2 phase-select stage pair.
  # Stages are 1-based (1..10); the register output is taken from stage 10.

  import Bitwise

  @code_length 1023

  @doc """
  Generate the bipolar (`±1`) C/A code for a G2 stage pair `{tap_a, tap_b}`.

  Binary `0` maps to `+1` and binary `1` maps to `-1`.
  """
  @spec bipolar_code({1..10, 1..10}) :: [integer()]
  def bipolar_code({tap_a, tap_b}) do
    {tap_a, tap_b}
    |> raw_code()
    |> Enum.map(fn bit -> 1 - 2 * bit end)
  end

  @doc """
  Generate the raw binary (`0`/`1`) C/A code for a G2 stage pair.
  """
  @spec raw_code({1..10, 1..10}) :: [integer()]
  def raw_code({tap_a, tap_b}) do
    g1_init = List.to_tuple(List.duplicate(1, 10))
    g2_init = g1_init

    {chips, _g1, _g2} =
      Enum.reduce(1..@code_length, {[], g1_init, g2_init}, fn _i, {acc, g1, g2} ->
        g1_out = elem(g1, 9)
        g2i = bxor(elem(g2, tap_a - 1), elem(g2, tap_b - 1))
        chip = bxor(g1_out, g2i)

        {[chip | acc], step_g1(g1), step_g2(g2)}
      end)

    Enum.reverse(chips)
  end

  # G1 feedback: stages 3 and 10 (1-based) -> indices 2 and 9.
  defp step_g1(g1) do
    feedback = bxor(elem(g1, 2), elem(g1, 9))
    shift(g1, feedback)
  end

  # G2 feedback: stages 2,3,6,8,9,10 -> indices 1,2,5,7,8,9.
  defp step_g2(g2) do
    feedback =
      elem(g2, 1)
      |> bxor(elem(g2, 2))
      |> bxor(elem(g2, 5))
      |> bxor(elem(g2, 7))
      |> bxor(elem(g2, 8))
      |> bxor(elem(g2, 9))

    shift(g2, feedback)
  end

  # Shift right: new stage 1 (index 0) is the feedback bit; stages 2..10 take
  # the previous stages 1..9.
  defp shift(reg, feedback) do
    {feedback, elem(reg, 0), elem(reg, 1), elem(reg, 2), elem(reg, 3), elem(reg, 4), elem(reg, 5),
     elem(reg, 6), elem(reg, 7), elem(reg, 8)}
  end
end
