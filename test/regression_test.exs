defmodule Orbis.RegressionTest do
  @moduledoc """
  Regression tests for bugs found during code review.
  """
  use ExUnit.Case

  describe "OMM timezone handling" do
    test "offset timestamp converts to UTC" do
      omm = %{
        "EPOCH" => "2026-04-05T13:16:46+05:00",
        "NORAD_CAT_ID" => 1,
        "INCLINATION" => 0.0,
        "RA_OF_ASC_NODE" => 0.0,
        "ECCENTRICITY" => 0.0,
        "ARG_OF_PERICENTER" => 0.0,
        "MEAN_ANOMALY" => 0.0,
        "MEAN_MOTION" => 1.0
      }

      {:ok, el} = Orbis.Format.OMM.parse(omm)
      # +05:00 means 13:16 local = 08:16 UTC
      assert el.epoch.hour == 8
      assert el.epoch.minute == 16
    end

    test "Z timestamp parses as UTC" do
      omm = %{
        "EPOCH" => "2026-04-05T13:16:46Z",
        "NORAD_CAT_ID" => 1,
        "INCLINATION" => 0.0,
        "RA_OF_ASC_NODE" => 0.0,
        "ECCENTRICITY" => 0.0,
        "ARG_OF_PERICENTER" => 0.0,
        "MEAN_ANOMALY" => 0.0,
        "MEAN_MOTION" => 1.0
      }

      {:ok, el} = Orbis.Format.OMM.parse(omm)
      assert el.epoch.hour == 13
    end

    test "bare timestamp treated as UTC" do
      omm = %{
        "EPOCH" => "2026-04-05T13:16:46.804800",
        "NORAD_CAT_ID" => 1,
        "INCLINATION" => 0.0,
        "RA_OF_ASC_NODE" => 0.0,
        "ECCENTRICITY" => 0.0,
        "ARG_OF_PERICENTER" => 0.0,
        "MEAN_ANOMALY" => 0.0,
        "MEAN_MOTION" => 1.0
      }

      {:ok, el} = Orbis.Format.OMM.parse(omm)
      assert el.epoch.hour == 13
    end
  end

  describe "SGP4 malformed elements" do
    test "nil epoch returns error" do
      el = %Orbis.Elements{epoch: nil, catalog_number: "0"}
      assert {:error, _} = Orbis.SGP4.propagate(el, ~U[2024-01-01 00:00:00Z])
    end

    test "nil catalog_number returns error" do
      el = %Orbis.Elements{epoch: ~U[2024-01-01 00:00:00Z], catalog_number: nil}
      assert {:error, _} = Orbis.SGP4.propagate(el, ~U[2024-01-01 00:00:00Z])
    end

    test "empty struct returns error, not crash" do
      el = %Orbis.Elements{}
      assert {:error, _} = Orbis.SGP4.propagate(el, ~U[2024-01-01 00:00:00Z])
    end
  end

  describe "TLE encode overflow" do
    test "raises on catalog_number > 5 chars" do
      el = %Orbis.Elements{
        catalog_number: "123456",
        epoch: ~U[2024-01-01 00:00:00Z],
        mean_motion: 15.0
      }

      assert_raise ArgumentError, ~r/catalog_number/, fn ->
        Orbis.Format.TLE.encode(el)
      end
    end

    test "raises on rev_number > 99999" do
      el = %Orbis.Elements{
        catalog_number: "25544",
        epoch: ~U[2024-01-01 00:00:00Z],
        mean_motion: 15.0,
        rev_number: 100_000
      }

      assert_raise ArgumentError, ~r/rev_number/, fn ->
        Orbis.Format.TLE.encode(el)
      end
    end
  end

  describe "TLE parse catalog_number trimming" do
    test "low NORAD IDs are trimmed" do
      {:ok, el} =
        Orbis.Format.TLE.parse(
          "1 00005U 58002B   00179.78495062  .00000023  00000-0  28098-4 0  4753",
          "2 00005  34.2682 348.7242 1859667 331.7664  19.3264 10.82419157413667"
        )

      assert el.catalog_number == "00005"
    end
  end

  describe "Passes.predict step_seconds guard" do
    test "raises on step_seconds: 0" do
      {:ok, el} =
        Orbis.Format.TLE.parse(
          "1 25544U 98067A   18184.80969102  .00001614  00000-0  31745-4 0  9993",
          "2 25544  51.6414 295.8524 0003435 262.6267 204.2868 15.54005638121106"
        )

      station = %{latitude: 51.5074, longitude: -0.1278, altitude_m: 11.0}

      assert_raise ArgumentError, ~r/step_seconds/, fn ->
        Orbis.Passes.predict(
          el,
          station,
          ~U[2024-12-19 00:00:00Z],
          ~U[2024-12-19 12:00:00Z],
          step_seconds: 0
        )
      end
    end
  end
end
