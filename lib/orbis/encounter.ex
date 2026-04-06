defmodule Orbis.Encounter do
  @moduledoc """
  Encounter geometry helpers for conjunction assessment.

  This module builds the relative frame used by collision-probability
  calculations and projects states and covariance into the encounter plane.

  It is the geometry layer beneath `Orbis.Collision`.
  """

  alias Orbis.Encounter.Frame

  @type vec3 :: {float(), float(), float()}

  @doc """
  Build an orthonormal encounter frame from two objects' states.

  Returns `{:ok, %Frame{}}` or `{:error, reason}`.
  """
  @spec frame(vec3(), vec3(), vec3(), vec3()) :: {:ok, Frame.t()} | {:error, String.t()}
  def frame(r1, v1, r2, v2) do
    dr = sub(r2, r1)
    dv = sub(v2, v1)

    _miss = mag(dr)
    rel_speed = mag(dv)

    if rel_speed < 1.0e-12 do
      {:error, "zero relative velocity"}
    else
      y_hat = normalize(dv)

      # dr_orthogonal = dr - (dr . y_hat) * y_hat
      dr_dot_y = dot(dr, y_hat)
      dr_ortho = sub(dr, scale_vec(y_hat, dr_dot_y))
      miss_ortho = mag(dr_ortho)

      x_hat =
        if miss_ortho < 1.0e-12 do
          # Objects on collision course. Find any orthogonal vector to y_hat.
          {yx, _yy, _yz} = y_hat
          trial = if abs(yx) < 0.9, do: {1.0, 0.0, 0.0}, else: {0.0, 1.0, 0.0}
          normalize(cross(y_hat, trial))
        else
          normalize(dr_ortho)
        end

      z_hat = cross(y_hat, x_hat)

      {:ok,
       %Frame{
         x_hat: x_hat,
         y_hat: y_hat,
         z_hat: z_hat,
         relative_position_km: dr,
         relative_velocity_km_s: dv,
         miss_km: miss_ortho,
         relative_speed_km_s: rel_speed
       }}
    end
  end

  @doc """
  Project a 3x3 ECI covariance matrix into the 2D encounter plane (x, z).
  """
  @spec encounter_plane_covariance(Frame.t(), [[float()]]) :: [[float()]]
  def encounter_plane_covariance(frame, cov_eci_3x3) do
    # Rotation matrix R: ECI -> Encounter
    # Row 1: x_hat
    # Row 2: z_hat
    r = [
      [elem(frame.x_hat, 0), elem(frame.x_hat, 1), elem(frame.x_hat, 2)],
      [elem(frame.z_hat, 0), elem(frame.z_hat, 1), elem(frame.z_hat, 2)]
    ]

    rt = transpose(r)

    # C_encounter = R * C_eci * R'
    mat_mul(mat_mul(r, cov_eci_3x3), rt)
  end

  # --- Internal Helpers ---

  defp mag({x, y, z}), do: :math.sqrt(x * x + y * y + z * z)
  defp sub({x1, y1, z1}, {x2, y2, z2}), do: {x1 - x2, y1 - y2, z1 - z2}
  defp dot({x1, y1, z1}, {x2, y2, z2}), do: x1 * x2 + y1 * y2 + z1 * z2
  defp scale_vec({x, y, z}, s), do: {x * s, y * s, z * s}
  defp normalize(v), do: {elem(v, 0) / mag(v), elem(v, 1) / mag(v), elem(v, 2) / mag(v)}

  defp cross({ax, ay, az}, {bx, by, bz}),
    do: {ay * bz - az * by, az * bx - ax * bz, ax * by - ay * bx}

  defp transpose(m), do: m |> Enum.zip() |> Enum.map(&Tuple.to_list/1)

  defp mat_mul(a, b) do
    bt = transpose(b)

    for row <- a,
        do:
          for(col <- bt, do: Enum.zip(row, col) |> Enum.map(fn {x, y} -> x * y end) |> Enum.sum())
  end
end
