defmodule Orbis.Forces.TwoBody do
  @moduledoc """
  Standard Keplerian Two-Body gravity force model.
  """

  # km³/s²
  @mu 398_600.4418

  @doc """
  Compute the acceleration due to two-body gravity.
  """
  @spec acceleration({float(), float(), float()}, {float(), float(), float()}) ::
          {float(), float(), float()}
  def acceleration({x, y, z}, _v) do
    r2 = x * x + y * y + z * z
    r = :math.sqrt(r2)
    r3 = r2 * r

    f = -@mu / r3

    {x * f, y * f, z * f}
  end
end
