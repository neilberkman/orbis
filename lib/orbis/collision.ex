defmodule Orbis.Collision do
  @moduledoc """
  Collision probability calculation for close approaches.

  Computes `Pc` from relative state, covariance, and hard-body radius
  using encounter-plane methods such as Foster 2D.

  This module is intended for operational conjunction screening and
  standards-based workflows such as CCSDS CDM ingestion.
  """

  alias Orbis.Collision.Result
  alias Orbis.Covariance
  alias Orbis.Encounter

  @type params :: %{
          r1: {float(), float(), float()},
          v1: {float(), float(), float()},
          cov1: [[float()]],
          r2: {float(), float(), float()},
          v2: {float(), float(), float()},
          cov2: [[float()]],
          hard_body_radius_km: float()
        }

  @doc """
  Compute collision probability from two objects' states and covariances.

  All positions in km, velocities in km/s, covariances in km².

  ## Options
    * `:method` — `:equal_area` (default) or `:numerical`

  Returns `{:ok, %Result{}}` or `{:error, reason}`.
  """
  @spec probability(params(), keyword()) :: {:ok, Result.t()} | {:error, String.t()}
  def probability(params, opts \\ []) do
    %{cov1: cov1, cov2: cov2} = params

    cond do
      not Covariance.valid_matrix?(cov1) ->
        {:error, "cov1 is not a 3x3 numeric matrix"}

      not Covariance.valid_matrix?(cov2) ->
        {:error, "cov2 is not a 3x3 numeric matrix"}

      not Covariance.positive_semidefinite?(cov1) ->
        {:error, "cov1 is not positive semidefinite"}

      not Covariance.positive_semidefinite?(cov2) ->
        {:error, "cov2 is not positive semidefinite"}

      true ->
        method = Keyword.get(opts, :method, :equal_area)

        case method do
          :equal_area -> probability_equal_area(params)
          :numerical -> probability_numerical(params)
          _ -> {:error, "unsupported method: #{method}"}
        end
    end
  end

  @doc """
  Compute Pc using the 2D Foster method with equal-area square approximation.
  """
  @spec probability_equal_area(params()) :: {:ok, Result.t()} | {:error, String.t()}
  def probability_equal_area(params) do
    %{r1: r1, v1: v1, cov1: cov1, r2: r2, v2: v2, cov2: cov2, hard_body_radius_km: hbr} = params

    case Encounter.frame(r1, v1, r2, v2) do
      {:ok, frame} ->
        c_combined = Covariance.add(cov1, cov2)
        c_encounter = Encounter.encounter_plane_covariance(frame, c_combined)
        {:ok, compute_foster_2d_equal_area(frame, c_encounter, hbr)}

      {:error, reason} ->
        {:error, reason}
    end
  end

  @doc """
  Compute Pc using the 2D Foster method with numerical integration over the circle.
  """
  @spec probability_numerical(params()) :: {:ok, Result.t()} | {:error, String.t()}
  def probability_numerical(params) do
    %{r1: r1, v1: v1, cov1: cov1, r2: r2, v2: v2, cov2: cov2, hard_body_radius_km: hbr} = params

    case Encounter.frame(r1, v1, r2, v2) do
      {:ok, frame} ->
        c_combined = Covariance.add(cov1, cov2)
        c_encounter = Encounter.encounter_plane_covariance(frame, c_combined)
        {:ok, compute_foster_2d_numerical(frame, c_encounter, hbr)}

      {:error, reason} ->
        {:error, reason}
    end
  end

  # --- Encounter-plane methods ---

  defp compute_foster_2d_equal_area(frame, c_enc, hbr) do
    {sigma_x, sigma_z, xm, zm, _theta} = principal_components(frame, c_enc)

    # Equal-area square Pc
    hsq = :math.sqrt(:math.pi() / 4.0) * hbr
    sqrt2 = :math.sqrt(2.0)

    pc =
      if sigma_x > 1.0e-12 and sigma_z > 1.0e-12 do
        0.25 *
          (:math.erf((xm + hsq) / (sqrt2 * sigma_x)) - :math.erf((xm - hsq) / (sqrt2 * sigma_x))) *
          (:math.erf((zm + hsq) / (sqrt2 * sigma_z)) - :math.erf((zm - hsq) / (sqrt2 * sigma_z)))
      else
        0.0
      end

    %Result{
      pc: max(0.0, pc),
      miss_km: frame.miss_km,
      relative_speed_km_s: frame.relative_speed_km_s,
      sigma_x_km: sigma_x,
      sigma_z_km: sigma_z,
      method: :foster_2d_equal_area
    }
  end

  defp compute_foster_2d_numerical(frame, c_enc, hbr) do
    {sigma_x, sigma_z, xm, zm, _theta} = principal_components(frame, c_enc)

    # Numerical integration over the circle using a simple polar grid (20x20)
    # This is a baseline implementation; higher precision would use adaptive quadrature.
    pc =
      if sigma_x > 1.0e-12 and sigma_z > 1.0e-12 do
        steps_r = 20
        steps_theta = 40
        dr = hbr / steps_r
        dtheta = 2.0 * :math.pi() / steps_theta

        # Integral sum: f(r, theta) * r * dr * dtheta
        for i <- 0..(steps_r - 1),
            j <- 0..(steps_theta - 1),
            reduce: 0.0 do
          acc ->
            r = (i + 0.5) * dr
            theta = (j + 0.5) * dtheta
            x = r * :math.cos(theta)
            z = r * :math.sin(theta)

            # Evaluate Gaussian centered at (xm, zm)
            arg =
              :math.pow(x - xm, 2) / (2.0 * sigma_x * sigma_x) +
                :math.pow(z - zm, 2) / (2.0 * sigma_z * sigma_z)

            f = :math.exp(-arg) / (2.0 * :math.pi() * sigma_x * sigma_z)
            acc + f * r * dr * dtheta
        end
      else
        0.0
      end

    %Result{
      pc: max(0.0, pc),
      miss_km: frame.miss_km,
      relative_speed_km_s: frame.relative_speed_km_s,
      sigma_x_km: sigma_x,
      sigma_z_km: sigma_z,
      method: :foster_2d_numerical
    }
  end

  defp principal_components(frame, c_enc) do
    [[a, b], [c, d]] = c_enc

    trace = a + d
    det = a * d - b * c
    disc = :math.sqrt(max(0.0, trace * trace / 4.0 - det))

    l1 = trace / 2.0 + disc
    l2 = trace / 2.0 - disc

    sigma_x = :math.sqrt(max(0.0, l1))
    sigma_z = :math.sqrt(max(0.0, l2))

    {vx, vy} =
      cond do
        abs(b) > 1.0e-30 -> normalize_2d({l1 - d, c})
        abs(c) > 1.0e-30 -> normalize_2d({c, l1 - a})
        true -> {1.0, 0.0}
      end

    theta = :math.atan2(vy, vx)
    xm = frame.miss_km * :math.cos(theta)
    zm = -frame.miss_km * :math.sin(theta)

    {sigma_x, sigma_z, xm, zm, theta}
  end

  defp normalize_2d({x, y}) do
    m = :math.sqrt(x * x + y * y)
    if m > 1.0e-30, do: {x / m, y / m}, else: {1.0, 0.0}
  end
end
