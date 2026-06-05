defmodule Orbis.GNSS.SP3Test do
  use ExUnit.Case, async: true

  alias Orbis.GNSS.SP3

  # Minimal but standards-shaped SP3-c position+clock file with two GPS sats and
  # two epochs (mirrors the astrodynamics-gnss parser fixture). G01 has a clock
  # estimate at both epochs; G02's second epoch is a missing-orbit (all-zero)
  # record.
  @sp3c """
  #cP2020  6 24  0  0  0.00000000       2 ORBIT IGS14 FIT  TST
  ## 2111 432000.00000000   900.00000000 59024 0.0000000000000
  +    2   G01G02  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
  ++         0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
  %c G  cc GPS ccc cccc cccc cccc cccc ccccc ccccc ccccc ccccc
  %c cc cc ccc ccc cccc cccc cccc cccc ccccc ccccc ccccc ccccc
  %f  1.2500000  1.025000000  0.00000000000  0.000000000000000
  %f  0.0000000  0.000000000  0.00000000000  0.000000000000000
  %i    0    0    0    0      0      0      0      0         0
  %i    0    0    0    0      0      0      0      0         0
  /* TEST SP3-c FIXTURE
  *  2020  6 24  0  0  0.00000000
  PG01  15000.000000 -20000.000000   5000.000000    123.456789
  PG02  -1234.567890   2345.678901  -3456.789012 999999.999999
  *  2020  6 24  0 15  0.00000000
  PG01  15100.000000 -20100.000000   5100.000000   -987.654321
  PG02      0.000000      0.000000      0.000000    100.000000
  EOF
  """

  describe "parse/1" do
    test "parses an in-memory SP3-c buffer into a handle" do
      assert {:ok, %SP3{} = sp3} = SP3.parse(@sp3c)
      assert is_reference(sp3.handle)
      # Header time scale comes from the %c descriptor (GPS -> GPST).
      assert sp3.time_scale == "GPST"
    end

    test "returns an error tuple on a malformed buffer" do
      assert {:error, _reason} = SP3.parse("not an sp3 file\n")
    end
  end

  describe "load/1 and load!/1" do
    setup do
      path =
        Path.join(System.tmp_dir!(), "orbis_sp3_test_#{System.unique_integer([:positive])}.sp3")

      File.write!(path, @sp3c)
      on_exit(fn -> File.rm(path) end)
      {:ok, path: path}
    end

    test "load/1 returns {:ok, handle} for a real file", %{path: path} do
      assert {:ok, %SP3{}} = SP3.load(path)
    end

    test "load/1 returns {:error, _} for a missing path" do
      assert {:error, :enoent} = SP3.load("/no/such/sp3/file.sp3")
    end

    test "load!/1 returns the handle", %{path: path} do
      assert %SP3{} = SP3.load!(path)
    end

    test "load!/1 raises on a missing path" do
      assert_raise ArgumentError, fn -> SP3.load!("/no/such/sp3/file.sp3") end
    end
  end

  describe "position/3" do
    setup do
      {:ok, sp3} = SP3.parse(@sp3c)
      {:ok, sp3: sp3}
    end

    test "evaluates at a node epoch, returning ITRF meters + clock seconds", %{sp3: sp3} do
      # First epoch is a spline node, so the value equals the parsed record,
      # converted km -> m and clock us -> s.
      assert {:ok, state} = SP3.position(sp3, "G01", ~N[2020-06-24 00:00:00])
      assert_in_delta state.x_m, 15_000_000.0, 1.0e-3
      assert_in_delta state.y_m, -20_000_000.0, 1.0e-3
      assert_in_delta state.z_m, 5_000_000.0, 1.0e-3
      assert_in_delta state.clock_s, 123.456789e-6, 1.0e-15
    end

    test "accepts an erlang datetime tuple", %{sp3: sp3} do
      assert {:ok, state} = SP3.position(sp3, "G01", {{2020, 6, 24}, {0, 0, 0}})
      assert_in_delta state.x_m, 15_000_000.0, 1.0e-3
    end

    test "errors for an unknown satellite", %{sp3: sp3} do
      assert {:error, _reason} = SP3.position(sp3, "G31", ~N[2020-06-24 00:00:00])
    end

    test "errors for a malformed satellite token", %{sp3: sp3} do
      assert {:error, {:bad_sat_id, _}} = SP3.position(sp3, "GXX", ~N[2020-06-24 00:00:00])
    end
  end
end
