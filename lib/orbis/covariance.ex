defmodule Orbis.Covariance do
  @moduledoc """
  Covariance matrix helpers for conjunction and orbit analysis.

  Supports covariance extraction, frame transforms, scaling, combination,
  and validation checks such as positive semidefiniteness.
  """

  @type mat3 :: [[float()]]
  @type vec3 :: {float(), float(), float()}

  @doc """
  Add two 3x3 matrices.
  """
  @spec add(mat3(), mat3()) :: mat3()
  def add(a, b) do
    for {ra, rb} <- Enum.zip(a, b),
        do: for({va, vb} <- Enum.zip(ra, rb), do: va + vb)
  end

  @doc """
  Transpose a 3x3 matrix.
  """
  @spec transpose(mat3()) :: mat3()
  def transpose(m), do: m |> Enum.zip() |> Enum.map(&Tuple.to_list/1)

  @doc """
  Matrix multiplication (3x3).
  """
  @spec mat_mul(mat3(), mat3()) :: mat3()
  def mat_mul(a, b) do
    bt = transpose(b)

    for row <- a,
        do:
          for(col <- bt, do: Enum.zip(row, col) |> Enum.map(fn {x, y} -> x * y end) |> Enum.sum())
  end

  @doc """
  Scale a matrix by a scalar.
  """
  @spec scale(mat3(), float()) :: mat3()
  def scale(m, s), do: for(row <- m, do: for(v <- row, do: v * s))

  @doc """
  Transform a 3x3 RTN covariance matrix to ECI.
  """
  @spec rtn_to_eci(mat3(), vec3(), vec3()) :: {:ok, mat3()} | {:error, String.t()}
  def rtn_to_eci(cov_rtn, r, v) do
    if valid_matrix?(cov_rtn) do
      case build_rtn_rotation(r, v) do
        {:ok, rot} ->
          rot_t = transpose(rot)
          {:ok, mat_mul(mat_mul(rot, cov_rtn), rot_t)}

        error ->
          error
      end
    else
      {:error, "invalid covariance matrix: must be a 3x3 numeric matrix"}
    end
  end

  @doc """
  Extract a 3x3 position covariance matrix from a 6-element lower triangle (RTN).
  Expected order: CR_R (0,0), CT_R (1,0), CT_T (1,1), CN_R (2,0), CN_T (2,1), CN_N (2,2).
  """
  @spec extract_pos_cov([float()]) :: {:ok, mat3()} | {:error, String.t()}
  def extract_pos_cov([cr_r, ct_r, ct_t, cn_r, cn_t, cn_n | _]) do
    res = [
      [cr_r, ct_r, cn_r],
      [ct_r, ct_t, cn_t],
      [cn_r, cn_t, cn_n]
    ]

    if Enum.all?(List.flatten(res), &is_number/1) do
      {:ok, res}
    else
      {:error, "non-numeric values in covariance list"}
    end
  end

  def extract_pos_cov(_), do: {:error, "invalid covariance list length"}

  @doc """
  Validate that the input is a 3x3 numeric matrix.
  """
  @spec valid_matrix?(any()) :: boolean()
  def valid_matrix?(m) when is_list(m) and length(m) == 3 do
    Enum.all?(m, fn row ->
      is_list(row) and length(row) == 3 and Enum.all?(row, &is_number/1)
    end)
  end

  def valid_matrix?(_), do: false

  @doc """
  Check if a 3x3 matrix is symmetric and positive semidefinite (PSD).

  A symmetric 3x3 matrix is PSD if all its principal minors are non-negative.
  """
  @spec positive_semidefinite?(mat3()) :: boolean()
  def positive_semidefinite?(m) do
    if valid_matrix?(m) and symmetric?(m) do
      # 1st order principal minors (diagonal elements)
      m11 = at(m, 0, 0)
      m22 = at(m, 1, 1)
      m33 = at(m, 2, 2)

      # 2nd order principal minors
      m12 = at(m, 0, 1)
      m13 = at(m, 0, 2)
      m23 = at(m, 1, 2)

      det12 = m11 * m22 - m12 * m12
      det13 = m11 * m33 - m13 * m13
      det23 = m22 * m33 - m23 * m23

      # 3rd order principal minor (determinant)
      det123 = det3x3(m)

      # All principal minors must be >= 0 (Sylvester's criterion for PSD)
      m11 >= -1.0e-15 and m22 >= -1.0e-15 and m33 >= -1.0e-15 and
        det12 >= -1.0e-12 and det13 >= -1.0e-12 and det23 >= -1.0e-12 and
        det123 >= -1.0e-12
    else
      false
    end
  end

  @doc """
  Check if a matrix is symmetric.
  """
  def symmetric?(m) do
    if valid_matrix?(m) do
      # Use small epsilon for symmetry due to possible float precision
      abs(at(m, 0, 1) - at(m, 1, 0)) < 1.0e-12 and
        abs(at(m, 0, 2) - at(m, 2, 0)) < 1.0e-12 and
        abs(at(m, 1, 2) - at(m, 2, 1)) < 1.0e-12
    else
      false
    end
  end

  # --- Internal Helpers ---

  defp at(m, i, j), do: Enum.at(Enum.at(m, i), j)

  defp det3x3(m) do
    a = at(m, 0, 0)
    b = at(m, 0, 1)
    c = at(m, 0, 2)
    d = at(m, 1, 0)
    e = at(m, 1, 1)
    f = at(m, 1, 2)
    g = at(m, 2, 0)
    h = at(m, 2, 1)
    i = at(m, 2, 2)

    a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g)
  end

  defp build_rtn_rotation(r, v) do
    r_mag = mag(r)

    if r_mag < 1.0e-30 do
      {:error, "zero position vector"}
    else
      r_hat = normalize(r)
      h = cross(r, v)
      h_mag = mag(h)

      if h_mag < 1.0e-30 do
        {:error, "position and velocity are parallel"}
      else
        n_hat = normalize(h)
        t_hat = cross(n_hat, r_hat)

        # Rotation matrix: RTN → ECI (columns are R, T, N unit vectors)
        {:ok,
         [
           [elem(r_hat, 0), elem(t_hat, 0), elem(n_hat, 0)],
           [elem(r_hat, 1), elem(t_hat, 1), elem(n_hat, 1)],
           [elem(r_hat, 2), elem(t_hat, 2), elem(n_hat, 2)]
         ]}
      end
    end
  end

  defp mag({x, y, z}), do: :math.sqrt(x * x + y * y + z * z)
  defp normalize(v), do: {elem(v, 0) / mag(v), elem(v, 1) / mag(v), elem(v, 2) / mag(v)}

  defp cross({ax, ay, az}, {bx, by, bz}),
    do: {ay * bz - az * by, az * bx - ax * bz, ax * by - ay * bx}
end
