defmodule Orbis.Collision do
  @moduledoc """
  Collision probability calculation.

  Computes the probability of collision (Pc) between two objects at
  the time of closest approach (TCA) given their states and covariance
  matrices. Implements the 2D Foster method (Foster 1992, Alfano 2005).

  The algorithm:
  1. Combine position covariances in ECI
  2. Construct the encounter plane from relative velocity
  3. Project combined covariance onto 2D encounter plane
  4. Eigendecompose to get principal axes
  5. Integrate 2D Gaussian over the hard-body radius circle

  ## Examples

      result = Orbis.Collision.probability(%{
        r1: {378.39559, 4305.721887, 5752.767554},
        v1: {2.360800244, 5.580331936, -4.322349039},
        cov1: [...],  # 3x3 position covariance in km²
        r2: {374.5180598, 4307.560983, 5751.130418},
        v2: {-5.388125081, -3.946827739, 3.322820358},
        cov2: [...],  # 3x3 position covariance in km²
        hard_body_radius_km: 0.020
      })

      result.pc        # collision probability
      result.miss_km   # miss distance in km
  """

  @doc """
  Compute collision probability from two objects' states and covariances.

  All positions in km, velocities in km/s, covariances in km².
  Hard body radius in km.

  Uses the equal-area-square approximation of the 2D Foster circular Pc,
  which requires no numerical integration and matches the circular method
  to high precision for typical conjunction geometries.

  ## Input map keys

    * `:r1`, `:r2` — position tuples `{x, y, z}` in km (ECI)
    * `:v1`, `:v2` — velocity tuples `{vx, vy, vz}` in km/s (ECI)
    * `:cov1`, `:cov2` — 3x3 position covariance as list of lists in km²
    * `:hard_body_radius_km` — combined hard body radius in km

  ## Returns

  A map with:
    * `:pc` — collision probability
    * `:miss_km` — miss distance in km
    * `:relative_speed_km_s` — relative speed at TCA
    * `:method` — `:foster_2d_equal_area`
  """
  @spec probability(map()) :: map()
  def probability(%{
        r1: {r1x, r1y, r1z},
        v1: {v1x, v1y, v1z},
        cov1: cov1,
        r2: {r2x, r2y, r2z},
        v2: {v2x, v2y, v2z},
        cov2: cov2,
        hard_body_radius_km: hbr
      }) do
    # Relative state
    dr = {r1x - r2x, r1y - r2y, r1z - r2z}
    dv = {v1x - v2x, v1y - v2y, v1z - v2z}

    miss = mag(dr)
    rel_speed = mag(dv)

    # Guard: zero relative velocity means no well-defined encounter plane
    if rel_speed < 1.0e-30 do
      %{
        pc: 0.0,
        miss_km: miss,
        relative_speed_km_s: 0.0,
        sigma_x_km: 0.0,
        sigma_z_km: 0.0,
        method: :foster_2d_equal_area
      }
    else
      probability_2d(dr, dv, miss, rel_speed, cov1, cov2, hbr)
    end
  end

  defp probability_2d(dr, dv, miss, rel_speed, cov1, cov2, hbr) do
    # Step 1: Combine position covariances
    c_combined = mat_add(cov1, cov2)

    # Step 2: Encounter frame (x = in-plane cross-track, y = along velocity, z = normal)
    y_hat = normalize(dv)
    h = cross(dr, dv)

    # Guard: if dr and dv are parallel (or dr is zero), pick an arbitrary normal
    z_hat =
      if mag(h) < 1.0e-30 do
        # Find a vector not parallel to dv
        {yx, _, _} = y_hat
        trial = if abs(yx) < 0.9, do: {1.0, 0.0, 0.0}, else: {0.0, 1.0, 0.0}
        normalize(cross(y_hat, trial))
      else
        normalize(h)
      end

    x_hat = cross(y_hat, z_hat)

    # Rotation matrix: ECI → encounter frame
    t_mat = [
      [elem(x_hat, 0), elem(x_hat, 1), elem(x_hat, 2)],
      [elem(y_hat, 0), elem(y_hat, 1), elem(y_hat, 2)],
      [elem(z_hat, 0), elem(z_hat, 1), elem(z_hat, 2)]
    ]

    # Step 3: Rotate covariance to encounter frame, extract 2D (x, z)
    c_enc = mat_mul(mat_mul(t_mat, c_combined), transpose(t_mat))

    c_2d = [
      [at(c_enc, 0, 0), at(c_enc, 0, 2)],
      [at(c_enc, 2, 0), at(c_enc, 2, 2)]
    ]

    # Step 4: Eigendecompose 2x2
    {l1, l2, v1_eig} = eigen_2x2(c_2d)

    # Guard against non-positive eigenvalues
    min_var = 1.0e-4 * hbr * (1.0e-4 * hbr)
    l1 = max(l1, min_var)
    l2 = max(l2, min_var)

    sigma_x = :math.sqrt(l1)
    sigma_z = :math.sqrt(l2)

    # Step 5: Miss distance in principal axes
    # Project relative position onto encounter plane principal axes
    dr_enc = mat_vec(t_mat, dr)
    miss_enc = {elem(dr_enc, 0), elem(dr_enc, 2)}

    # Rotate to principal axes
    xm = abs(elem(miss_enc, 0) * elem(v1_eig, 0) + elem(miss_enc, 1) * elem(v1_eig, 1))
    zm = abs(-elem(miss_enc, 0) * elem(v1_eig, 1) + elem(miss_enc, 1) * elem(v1_eig, 0))

    # Step 6: Equal-area square Pc (no numerical integration needed)
    hsq = :math.sqrt(:math.pi() / 4.0) * hbr
    sqrt2 = :math.sqrt(2.0)

    pc =
      0.25 *
        (erf((xm + hsq) / (sqrt2 * sigma_x)) - erf((xm - hsq) / (sqrt2 * sigma_x))) *
        (erf((zm + hsq) / (sqrt2 * sigma_z)) - erf((zm - hsq) / (sqrt2 * sigma_z)))

    %{
      pc: pc,
      miss_km: miss,
      relative_speed_km_s: rel_speed,
      sigma_x_km: sigma_x,
      sigma_z_km: sigma_z,
      method: :foster_2d_equal_area
    }
  end

  # --- Vector math ---

  defp mag({x, y, z}), do: :math.sqrt(x * x + y * y + z * z)

  defp normalize(v) do
    m = mag(v)
    {elem(v, 0) / m, elem(v, 1) / m, elem(v, 2) / m}
  end

  defp cross({ax, ay, az}, {bx, by, bz}) do
    {ay * bz - az * by, az * bx - ax * bz, ax * by - ay * bx}
  end

  # --- Matrix math (3x3 as list of lists) ---

  defp mat_add(a, b) do
    for {ar, br} <- Enum.zip(a, b) do
      for {av, bv} <- Enum.zip(ar, br), do: av + bv
    end
  end

  defp transpose(m) do
    m |> Enum.zip() |> Enum.map(&Tuple.to_list/1)
  end

  defp mat_mul(a, b) do
    bt = transpose(b)

    for row <- a do
      for col <- bt do
        Enum.zip(row, col) |> Enum.map(fn {x, y} -> x * y end) |> Enum.sum()
      end
    end
  end

  defp mat_vec(m, {x, y, z}) do
    [r0, r1, r2] = m

    {
      Enum.at(r0, 0) * x + Enum.at(r0, 1) * y + Enum.at(r0, 2) * z,
      Enum.at(r1, 0) * x + Enum.at(r1, 1) * y + Enum.at(r1, 2) * z,
      Enum.at(r2, 0) * x + Enum.at(r2, 1) * y + Enum.at(r2, 2) * z
    }
  end

  defp at(m, i, j), do: Enum.at(m, i) |> Enum.at(j)

  # --- 2x2 eigendecomposition ---

  defp eigen_2x2([[a, b], [c, d]]) do
    trace = a + d
    det = a * d - b * c
    disc = :math.sqrt(max(trace * trace / 4.0 - det, 0.0))

    l1 = trace / 2.0 + disc
    l2 = trace / 2.0 - disc

    # Eigenvector for l1: from (A - λI)v = 0
    # Row 1: (a-λ)x + by = 0  →  direction (b, λ-a)
    # Row 2: cx + (d-λ)y = 0  →  direction (λ-d, c)
    {vx, vy} =
      if abs(b) > 1.0e-30 do
        normalize_2d({l1 - d, c})
      else
        if abs(c) > 1.0e-30 do
          normalize_2d({c, l1 - a})
        else
          {1.0, 0.0}
        end
      end

    {l1, l2, {vx, vy}}
  end

  defp normalize_2d({x, y}) do
    m = :math.sqrt(x * x + y * y)
    if m > 1.0e-30, do: {x / m, y / m}, else: {1.0, 0.0}
  end

  defp erf(x), do: :math.erf(x)
end
