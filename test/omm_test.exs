defmodule Orbis.OMMTest do
  use ExUnit.Case

  @fixtures_dir Path.join(__DIR__, "fixtures/celestrak")

  setup do
    omms = Path.join(@fixtures_dir, "stations.json") |> File.read!() |> Jason.decode!()
    iss_omm = Enum.find(omms, &(&1["NORAD_CAT_ID"] == 25544))
    %{omms: omms, iss_omm: iss_omm}
  end

  describe "from_omm/1" do
    test "parses ISS OMM record", %{iss_omm: omm} do
      {:ok, tle} = Orbis.Format.OMM.parse(omm)
      assert tle.catalog_number == "25544"
      assert tle.inclination_deg > 51.0 and tle.inclination_deg < 52.0
      assert tle.eccentricity > 0.0 and tle.eccentricity < 0.01
      assert tle.mean_motion > 15.0 and tle.mean_motion < 16.0
      assert tle.object_name == "ISS (ZARYA)"
    end

    test "parses all station OMMs", %{omms: omms} do
      results = Enum.map(omms, &Orbis.Format.OMM.parse/1)
      ok_count = Enum.count(results, &match?({:ok, _}, &1))
      assert ok_count == length(omms)
    end
  end

  describe "propagation from OMM" do
    test "OMM-sourced TLE propagates correctly", %{iss_omm: omm} do
      {:ok, tle} = Orbis.Format.OMM.parse(omm)
      {:ok, teme} = Orbis.propagate(tle, tle.epoch)

      {x, y, z} = teme.position
      radius = :math.sqrt(x * x + y * y + z * z)
      assert radius > 6500 and radius < 7200
    end
  end

  describe "to_omm/1" do
    test "round-trips through OMM", %{iss_omm: original} do
      {:ok, tle} = Orbis.Format.OMM.parse(original)
      omm = Orbis.Format.OMM.encode(tle)

      assert omm["NORAD_CAT_ID"] == 25544
      assert_in_delta omm["INCLINATION"], original["INCLINATION"], 1.0e-10
      assert_in_delta omm["ECCENTRICITY"], original["ECCENTRICITY"], 1.0e-10
      assert_in_delta omm["MEAN_MOTION"], original["MEAN_MOTION"], 1.0e-10
    end
  end
end
