defmodule Orbis.Forces.J2 do
  @moduledoc """
  Earth's oblateness (J2) perturbation force model.
  """

  # km³/s²
  @mu 398_600.4418
  # km (WGS84)
  @re 6378.137
  @j2 1.08262668e-3

  @doc """
  Compute the acceleration due to J2 perturbation in ECI.
  Ref: Vallado (4th ed), Eq 8-30.
  """
  @spec acceleration({float(), float(), float()}, {float(), float(), float()}) ::
          {float(), float(), float()}
  def acceleration({x, y, z}, _v) do
    r2 = x * x + y * y + z * z
    r = :math.sqrt(r2)

    # Common factor
    # a_j2 = -3/2 * J2 * (mu/r^2) * (Re/r)^2 * [...]
    f = 1.5 * @j2 * (@mu / r2) * :math.pow(@re / r, 2)

    z_r2 = z * z / r2

    # Acceleration components
    ax = f * (x / r) * (5.0 * z_r2 - 1.0)
    ay = f * (y / r) * (5.0 * z_r2 - 1.0)
    az = f * (z / r) * (5.0 * z_r2 - 3.0)

    {ax, ay, az}
  end
end
