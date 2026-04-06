defmodule Orbis.CollisionTest do
  @moduledoc """
  Collision probability tests.

  Reference: NASA CARA Analysis Tools (Omitron test case).
  """
  use ExUnit.Case

  # NASA CARA Omitron test case — states in ECI km/km/s, covariances in km²
  @omitron_params %{
    r1: {378.39559, 4305.721887, 5752.767554},
    v1: {2.360800244, 5.580331936, -4.322349039},
    cov1: [
      [44.5757544811362, 81.6751751052616, -67.8687662707124],
      [81.6751751052616, 158.453402956163, -128.616921644857],
      [-67.8687662707124, -128.616921644858, 105.490542562701]
    ],
    r2: {374.5180598, 4307.560983, 5751.130418},
    v2: {-5.388125081, -3.946827739, 3.322820358},
    cov2: [
      [2.31067077720423, 1.69905293875632, -1.4170164577661],
      [1.69905293875632, 1.24957388457206, -1.04174164279599],
      [-1.4170164577661, -1.04174164279599, 0.869260558223714]
    ],
    hard_body_radius_km: 0.020
  }

  describe "probability/1" do
    test "Omitron test case: equal-area Pc matches CARA reference" do
      result = Orbis.Collision.probability(@omitron_params)

      # CARA reference: equal-area square Pc = 2.70601573490111e-05
      assert_in_delta result.pc, 2.70601573490111e-05, 1.0e-09
      assert result.miss_km > 0 and result.miss_km < 10
      assert result.method == :foster_2d_equal_area
    end

    test "zero miss distance with tight covariance gives high Pc" do
      # Small covariance (10m σ) with 15m HBR and zero miss = high Pc
      # 10m = 0.01 km, squared
      sigma_km2 = 0.01 * 0.01

      result =
        Orbis.Collision.probability(%{
          r1: {7000.0, 0.0, 0.0},
          v1: {0.0, 7.5, 0.0},
          cov1: [[sigma_km2, 0.0, 0.0], [0.0, sigma_km2, 0.0], [0.0, 0.0, sigma_km2]],
          r2: {7000.0, 0.0, 0.0},
          v2: {0.0, -7.5, 0.0},
          cov2: [[sigma_km2, 0.0, 0.0], [0.0, sigma_km2, 0.0], [0.0, 0.0, sigma_km2]],
          hard_body_radius_km: 0.015
        })

      assert result.pc > 0.1
    end

    test "large miss distance gives near-zero Pc" do
      result =
        Orbis.Collision.probability(%{
          r1: {7000.0, 0.0, 0.0},
          v1: {0.0, 7.5, 0.0},
          cov1: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
          r2: {7100.0, 0.0, 0.0},
          v2: {0.0, -7.5, 0.0},
          cov2: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
          hard_body_radius_km: 0.015
        })

      assert result.pc < 1.0e-10
    end

    test "larger HBR increases Pc" do
      small = Orbis.Collision.probability(Map.put(@omitron_params, :hard_body_radius_km, 0.010))
      large = Orbis.Collision.probability(Map.put(@omitron_params, :hard_body_radius_km, 0.040))

      assert large.pc > small.pc
    end

    test "zero relative velocity returns Pc = 0 without crashing" do
      result =
        Orbis.Collision.probability(%{
          r1: {7000.0, 0.0, 0.0},
          v1: {0.0, 7.5, 0.0},
          cov1: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
          r2: {7000.01, 0.0, 0.0},
          v2: {0.0, 7.5, 0.0},
          cov2: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
          hard_body_radius_km: 0.015
        })

      assert result.pc == 0.0
      assert result.relative_speed_km_s == 0.0
    end
  end
end
