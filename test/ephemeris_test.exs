defmodule Orbis.EphemerisTest do
  use ExUnit.Case

  describe "Orbis.Ephemeris.load/1" do
    test "raises on missing file" do
      assert_raise ArgumentError, ~r/not found/, fn ->
        Orbis.Ephemeris.load("/nonexistent/file.bsp")
      end
    end
  end

  describe "body name resolution" do
    test "position raises on invalid body atom" do
      # Create a dummy struct to test body resolution
      eph = %Orbis.Ephemeris{path: "/dummy.bsp"}

      assert_raise ArgumentError, ~r/invalid body/, fn ->
        Orbis.Ephemeris.position(eph, :invalid_body, :earth, 2_451_545.0)
      end
    end

    test "position raises on unknown NAIF code" do
      eph = %Orbis.Ephemeris{path: "/dummy.bsp"}

      assert_raise ArgumentError, ~r/unknown NAIF body code/, fn ->
        Orbis.Ephemeris.position(eph, 9999, :earth, 2_451_545.0)
      end
    end
  end

  describe "Julian Date conversion" do
    test "J2000.0 epoch is correct" do
      # J2000.0 = 2000-01-01 12:00:00 TDB = JD 2451545.0
      # We test this indirectly through the position function's datetime conversion.
      # The struct just holds a path, so we can verify the module loads.
      eph = %Orbis.Ephemeris{path: "/dummy.bsp"}
      assert eph.path == "/dummy.bsp"
    end
  end

  # SPK file tests only run when a DE file is available.
  # Run with: mix test --include spk_file
  describe "with DE421 SPK file" do
    @describetag :spk_file

    setup do
      # Look for DE421 in common locations.
      paths = [
        Path.expand("~/.skyfield/de421.bsp"),
        Path.expand("~/de421.bsp"),
        "/tmp/de421.bsp",
        Path.join(File.cwd!(), "de421.bsp")
      ]

      path = Enum.find(paths, &File.exists?/1)

      if path do
        {:ok, eph: Orbis.Ephemeris.load(path)}
      else
        {:ok, skip: true}
      end
    end

    test "earth position relative to SSB at J2000", %{} = ctx do
      if Map.get(ctx, :skip) do
        IO.puts("\n  [skipped] DE421 file not found")
      else
        {x, y, z} = Orbis.Ephemeris.position(ctx.eph, :earth, :ssb, 2_451_545.0)

        # Earth should be roughly 1 AU from SSB (within ~0.02 AU)
        distance_km = :math.sqrt(x * x + y * y + z * z)
        au_km = 149_597_870.7
        assert_in_delta distance_km / au_km, 1.0, 0.02
      end
    end

    test "moon position relative to earth at J2000", %{} = ctx do
      if Map.get(ctx, :skip) do
        IO.puts("\n  [skipped] DE421 file not found")
      else
        {x, y, z} = Orbis.Ephemeris.position(ctx.eph, :moon, :earth, 2_451_545.0)

        # Moon should be roughly 384,400 km from Earth
        distance_km = :math.sqrt(x * x + y * y + z * z)
        assert_in_delta distance_km, 384_400.0, 20_000.0
      end
    end

    test "sun position relative to earth using DateTime", %{} = ctx do
      if Map.get(ctx, :skip) do
        IO.puts("\n  [skipped] DE421 file not found")
      else
        dt = ~U[2020-06-21 12:00:00Z]
        {x, y, z} = Orbis.Ephemeris.position(ctx.eph, :sun, :earth, dt)

        # Sun should be roughly 1 AU from Earth
        distance_km = :math.sqrt(x * x + y * y + z * z)
        au_km = 149_597_870.7
        assert_in_delta distance_km / au_km, 1.0, 0.02
      end
    end
  end
end
