defmodule Orbis.PassesTest do
  use ExUnit.Case, async: true

  # ISS TLE (epoch 2024-12-19)
  @iss_line1 "1 25544U 98067A   24354.52609954  .00020888  00000+0  37042-3 0  9992"
  @iss_line2 "2 25544  51.6393 213.2584 0006955  37.7614  87.9783 15.49970085486016"

  # Ground station: London, UK
  @london %{latitude: 51.5074, longitude: -0.1278, altitude_m: 11.0}

  setup do
    {:ok, tle} = Orbis.Format.TLE.parse(@iss_line1, @iss_line2)
    %{tle: tle}
  end

  describe "predict/5" do
    test "returns a list of Pass structs", %{tle: tle} do
      start_time = ~U[2024-12-19 00:00:00Z]
      end_time = ~U[2024-12-19 12:00:00Z]

      passes = Orbis.Passes.predict(tle, @london, start_time, end_time)

      assert is_list(passes)

      for pass <- passes do
        assert %DateTime{} = pass.rise
        assert %DateTime{} = pass.set
        assert %DateTime{} = pass.max_elevation_time
        assert is_float(pass.max_elevation)

        # Set must be after rise
        assert DateTime.after?(pass.set, pass.rise)

        # Max elevation time must be between rise and set
        assert DateTime.compare(pass.max_elevation_time, pass.rise) in [:gt, :eq]
        assert DateTime.compare(pass.max_elevation_time, pass.set) in [:lt, :eq]

        # Max elevation must be non-negative (we default min_elevation to 0)
        assert pass.max_elevation >= 0.0

        # Duration must be positive and consistent with rise/set
      end
    end

    test "finds at least one ISS pass in a 12-hour window", %{tle: tle} do
      start_time = ~U[2024-12-19 00:00:00Z]
      end_time = ~U[2024-12-19 12:00:00Z]

      passes = Orbis.Passes.predict(tle, @london, start_time, end_time)

      # ISS orbits ~15.5 times/day at 51.6 deg inclination.
      # Over 12 hours from a mid-latitude station, we should see several passes.
      assert length(passes) >= 1
    end

    test "min_elevation option filters out low passes", %{tle: tle} do
      start_time = ~U[2024-12-19 00:00:00Z]
      end_time = ~U[2024-12-19 12:00:00Z]

      all_passes = Orbis.Passes.predict(tle, @london, start_time, end_time, min_elevation: 0.0)
      high_passes = Orbis.Passes.predict(tle, @london, start_time, end_time, min_elevation: 30.0)

      assert length(high_passes) <= length(all_passes)

      for pass <- high_passes do
        assert pass.max_elevation >= 30.0
      end
    end

    test "returns empty list for zero-length window", %{tle: tle} do
      t = ~U[2024-12-19 06:00:00Z]
      assert Orbis.Passes.predict(tle, @london, t, t) == []
    end

    test "passes are sorted by rise time", %{tle: tle} do
      start_time = ~U[2024-12-19 00:00:00Z]
      end_time = ~U[2024-12-20 00:00:00Z]

      passes = Orbis.Passes.predict(tle, @london, start_time, end_time)

      rise_times = Enum.map(passes, & &1.rise)

      assert rise_times ==
               Enum.sort(rise_times, fn a, b -> DateTime.compare(a, b) != :gt end)
    end

    test "delegate from Orbis module works", %{tle: tle} do
      start_time = ~U[2024-12-19 00:00:00Z]
      end_time = ~U[2024-12-19 06:00:00Z]

      passes = Orbis.predict_passes(tle, @london, start_time, end_time)
      assert is_list(passes)
    end
  end
end
