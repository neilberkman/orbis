defmodule Orbis.CCSDS.CDMTest do
  @moduledoc """
  CDM KVN parser tests.

  Primary fixture: CCSDS 508.0-B-1 CDMExample2 (JSPOC, SATELLITE A vs FENGYUN 1C DEB).
  Reference Pc = 4.835e-05, miss distance = 715 m, relative speed = 14762 m/s.
  """
  use ExUnit.Case

  alias Orbis.CCSDS.CDM

  @fixture_path "test/fixtures/cdm/ccsds_example2.kvn"

  setup_all do
    kvn = File.read!(@fixture_path)
    {:ok, cdm} = CDM.parse(kvn)
    %{cdm: cdm, kvn: kvn}
  end

  describe "parse/1" do
    test "parses header fields", %{cdm: cdm} do
      assert cdm.originator == "JSPOC"
      assert cdm.message_id == "201113719185"
      assert cdm.creation_date == ~U[2010-03-12 22:31:12.000Z]
    end

    test "parses TCA", %{cdm: cdm} do
      assert cdm.tca == ~U[2010-03-13 22:37:52.618Z]
    end

    test "parses relative metadata", %{cdm: cdm} do
      assert_in_delta cdm.miss_distance_m, 715.0, 0.1
      assert_in_delta cdm.relative_speed_m_s, 14762.0, 0.1
    end

    test "parses collision probability", %{cdm: cdm} do
      assert_in_delta cdm.collision_probability, 4.835e-05, 1.0e-08
      assert cdm.collision_probability_method == "FOSTER-1992"
    end

    test "parses object1 metadata", %{cdm: cdm} do
      obj = cdm.object1
      assert obj.object_designator == "12345"
      assert obj.object_name == "SATELLITE A"
      assert obj.international_designator == "1997-030E"
      assert obj.object_type == "PAYLOAD"
      assert obj.ref_frame == "EME2000"
    end

    test "parses object1 state vector", %{cdm: cdm} do
      {r, v} = cdm.object1.state

      assert_in_delta elem(r, 0), 2570.097065, 1.0e-6
      assert_in_delta elem(r, 1), 2244.654904, 1.0e-6
      assert_in_delta elem(r, 2), 6281.497978, 1.0e-6

      assert_in_delta elem(v, 0), 4.418769571, 1.0e-9
      assert_in_delta elem(v, 1), 4.833547743, 1.0e-9
      assert_in_delta elem(v, 2), -3.526774282, 1.0e-9
    end

    test "parses object1 RTN covariance", %{cdm: cdm} do
      [cr_r, ct_r, ct_t, cn_r, cn_t, cn_n] = cdm.object1.covariance_rtn

      assert_in_delta cr_r, 4.142e+01, 1.0e-2
      assert_in_delta ct_r, -8.579e+00, 1.0e-3
      assert_in_delta ct_t, 2.533e+03, 1.0e-1
      assert_in_delta cn_r, -2.313e+01, 1.0e-2
      assert_in_delta cn_t, 1.336e+01, 1.0e-2
      assert_in_delta cn_n, 7.098e+01, 1.0e-2
    end

    test "parses object2 metadata", %{cdm: cdm} do
      obj = cdm.object2
      assert obj.object_designator == "30337"
      assert obj.object_name == "FENGYUN 1C DEB"
      assert obj.object_type == "DEBRIS"
    end

    test "parses object2 state vector", %{cdm: cdm} do
      {r, v} = cdm.object2.state

      assert_in_delta elem(r, 0), 2569.540800, 1.0e-6
      assert_in_delta elem(r, 1), 2245.093614, 1.0e-6
      assert_in_delta elem(r, 2), 6281.599946, 1.0e-6

      assert_in_delta elem(v, 0), -2.888612500, 1.0e-9
      assert_in_delta elem(v, 1), -6.007247516, 1.0e-9
      assert_in_delta elem(v, 2), 3.328770172, 1.0e-9
    end

    test "parses object2 RTN covariance", %{cdm: cdm} do
      [cr_r, ct_r, ct_t, cn_r, cn_t, cn_n] = cdm.object2.covariance_rtn

      assert_in_delta cr_r, 1.337e+03, 1.0e-1
      assert_in_delta ct_r, -4.806e+04, 1.0
      assert_in_delta ct_t, 2.492e+06, 1.0
      assert_in_delta cn_r, -3.298e+01, 1.0e-2
      assert_in_delta cn_t, -7.5888e+02, 1.0e-1
      assert_in_delta cn_n, 7.105e+01, 1.0e-2
    end

    test "returns error for missing TCA" do
      assert {:error, _} = CDM.parse("CREATION_DATE = 2024-01-01T00:00:00.000\n")
    end

    test "returns error for garbage input" do
      assert {:error, _} = CDM.parse("not a CDM at all")
    end

    test "returns error for incomplete object blocks" do
      kvn = """
      CCSDS_CDM_VERS = 1.0
      CREATION_DATE = 2024-01-01T00:00:00.000
      ORIGINATOR = TEST
      MESSAGE_ID = INCOMPLETE
      TCA = 2024-01-01T12:00:00.000
      OBJECT = OBJECT1
      OBJECT_DESIGNATOR = 00001
      X = 7000.0 [km]
      OBJECT = OBJECT2
      OBJECT_DESIGNATOR = 00002
      """

      assert {:error, "incomplete state vector"} = CDM.parse(kvn)
    end

    test "HBR is nil when not in COMMENT", %{cdm: cdm} do
      assert cdm.hard_body_radius_m == nil
    end

    test "round-trips through encode/1", %{cdm: cdm} do
      kvn = CDM.encode(cdm)
      {:ok, cdm2} = CDM.parse(kvn)

      assert cdm2.message_id == cdm.message_id
      assert cdm2.tca == cdm.tca
      assert_in_delta cdm2.miss_distance_m, cdm.miss_distance_m, 1.0e-6
      assert cdm2.object1.object_designator == cdm.object1.object_designator
      assert cdm2.object1.state == cdm.object1.state
    end

    test "parses HBR from COMMENT" do
      kvn = """
      CCSDS_CDM_VERS = 1.0
      CREATION_DATE = 2024-01-01T00:00:00.000
      ORIGINATOR = TEST
      MESSAGE_ID = HBR_TEST
      COMMENT HBR = 15.5
      TCA = 2024-01-01T12:00:00.000
      MISS_DISTANCE = 100 [m]
      COLLISION_PROBABILITY = 1.0E-06
      COLLISION_PROBABILITY_METHOD = FOSTER-1992
      OBJECT = OBJECT1
      OBJECT_DESIGNATOR = 00001
      X = 7000.0 [km]
      Y = 0.0 [km]
      Z = 0.0 [km]
      X_DOT = 0.0 [km/s]
      Y_DOT = 7.5 [km/s]
      Z_DOT = 0.0 [km/s]
      CR_R = 1.0 [m**2]
      CT_R = 0.0 [m**2]
      CT_T = 1.0 [m**2]
      CN_R = 0.0 [m**2]
      CN_T = 0.0 [m**2]
      CN_N = 1.0 [m**2]
      OBJECT = OBJECT2
      OBJECT_DESIGNATOR = 00002
      X = 7000.1 [km]
      Y = 0.0 [km]
      Z = 0.0 [km]
      X_DOT = 0.0 [km/s]
      Y_DOT = -7.5 [km/s]
      Z_DOT = 0.0 [km/s]
      CR_R = 1.0 [m**2]
      CT_R = 0.0 [m**2]
      CT_T = 1.0 [m**2]
      CN_R = 0.0 [m**2]
      CN_T = 0.0 [m**2]
      CN_N = 1.0 [m**2]
      """

      {:ok, cdm} = CDM.parse(kvn)
      assert_in_delta cdm.hard_body_radius_m, 15.5, 0.01
    end
  end

  describe "to_collision_params/1" do
    test "produces valid input for Collision.probability/1", %{cdm: cdm} do
      # Override HBR since fixture doesn't have one in COMMENT
      cdm = %{cdm | hard_body_radius_m: 20.0}
      params = CDM.to_collision_params(cdm)

      assert is_tuple(params.r1)
      assert is_tuple(params.v1)
      assert is_tuple(params.r2)
      assert is_tuple(params.v2)
      assert is_list(params.cov1)
      assert is_list(params.cov2)
      assert is_float(params.hard_body_radius_km)

      # Covariances should be 3x3 lists
      assert length(params.cov1) == 3
      assert length(hd(params.cov1)) == 3

      # HBR should be converted from m to km
      assert_in_delta params.hard_body_radius_km, 0.020, 1.0e-6
    end

    test "default HBR when CDM has none", %{cdm: cdm} do
      params = CDM.to_collision_params(cdm)
      # Default is 15m = 0.015 km
      assert_in_delta params.hard_body_radius_km, 0.015, 1.0e-6
    end

    test "CDM -> Pc pipeline produces physically reasonable Pc", %{cdm: cdm} do
      # CDM states Pc = 4.835e-05 (JSPOC FOSTER-1992)
      cdm = %{cdm | hard_body_radius_m: 20.0}
      params = CDM.to_collision_params(cdm)
      {:ok, result} = Orbis.Collision.probability(params)

      # IO.inspect(result.pc, label: "CDM Pc")
      assert result.pc > 1.0e-9 and result.pc < 1.0e-4
      assert_in_delta result.miss_km, 0.716, 0.001
      assert result.method == :foster_2d_equal_area
    end
  end
end
