defmodule Orbis.GNSS.ConstellationTest do
  # Not async: the live-fetch error test toggles app env to simulate the optional
  # CelesTrak Req dependency being absent.
  use ExUnit.Case, async: false

  alias Orbis.GNSS.Constellation
  alias Orbis.GNSS.Constellation.Record
  alias Orbis.GNSS.SP3

  @fixtures Path.join(__DIR__, "fixtures/gnss_constellation")

  @sp3 """
  #cP2020  6 24  0  0  0.00000000       1 ORBIT IGS14 FIT  TST
  ## 2111 432000.00000000   900.00000000 59024 0.0000000000000
  +    2   G03G32  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
  ++         0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
  %c G  cc GPS ccc cccc cccc cccc cccc ccccc ccccc ccccc ccccc
  %c cc cc ccc ccc cccc cccc cccc cccc ccccc ccccc ccccc ccccc
  %f  1.2500000  1.025000000  0.00000000000  0.000000000000000
  %f  0.0000000  0.000000000  0.00000000000  0.000000000000000
  %i    0    0    0    0      0      0      0      0         0
  %i    0    0    0    0      0      0      0      0         0
  /* TEST SP3-c FIXTURE
  *  2020  6 24  0  0  0.00000000
  PG03  15000.000000 -20000.000000   5000.000000    123.456789
  PG32  -1234.567890   2345.678901  -3456.789012    100.000000
  EOF
  """

  defp celestrak_omms do
    @fixtures
    |> Path.join("gps_ops_sample.json")
    |> File.read!()
    |> Jason.decode!()
  end

  defp navcen_html do
    @fixtures
    |> Path.join("navcen_gps_sample.html")
    |> File.read!()
  end

  defp merged_records do
    {:ok, records} = Constellation.from_celestrak_omm(celestrak_omms())
    {:ok, statuses} = Constellation.parse_navcen_html(navcen_html())
    Constellation.merge_navcen(records, statuses)
  end

  describe "fetch_gps/1" do
    test "surfaces the typed CelesTrak dependency error" do
      Application.put_env(:orbis, :celestrak_req_available, false)
      on_exit(fn -> Application.delete_env(:orbis, :celestrak_req_available) end)

      assert {:error, :req_not_available} = Constellation.fetch_gps()
    end
  end

  describe "from_celestrak_omm/1" do
    test "normalizes GPS PRN/NORAD records from gps-ops OMM JSON" do
      assert {:ok, records} = Constellation.from_celestrak_omm(celestrak_omms())
      assert Enum.map(records, & &1.prn) == [3, 5, 13, 19]

      prn3 = Enum.find(records, &(&1.prn == 3))
      assert %Record{} = prn3
      assert prn3.system == :gps
      assert prn3.svn == nil
      assert prn3.norad_id == 40294
      assert prn3.sp3_id == "G03"
      assert prn3.active?
      assert prn3.usable?
      assert prn3.source.celestrak.group == "gps-ops"
      assert prn3.source.celestrak.block_type == "IIF"
    end

    test "rejects a gps-ops record without a PRN in OBJECT_NAME" do
      bad = [%{"OBJECT_NAME" => "GPS WITHOUT PRN", "NORAD_CAT_ID" => 1}]

      assert {:error, {:bad_celestrak_record, {:missing_prn, "GPS WITHOUT PRN"}, _}} =
               Constellation.from_celestrak_omm(bad)
    end
  end

  describe "parse_navcen_html/1 and merge_navcen/2" do
    test "parses SVN and active NANU status rows" do
      assert {:ok, statuses} = Constellation.parse_navcen_html(navcen_html())
      assert Enum.map(statuses, &{&1.prn, &1.svn}) == [{3, 69}, {5, 50}, {13, 43}, {19, 59}]

      prn19 = Enum.find(statuses, &(&1.prn == 19))
      assert prn19.active_nanu?
      refute prn19.usable?
      assert prn19.nanu_type == "UNUSABLE"
      assert prn19.source.navcen.active_nanu?
    end

    test "overlays NAVCEN SVN and usability on CelesTrak identity records" do
      records = merged_records()

      assert Enum.map(records, &{&1.prn, &1.svn, &1.usable?}) == [
               {3, 69, true},
               {5, 50, true},
               {13, nil, true},
               {19, 59, false}
             ]

      prn3 = Enum.find(records, &(&1.prn == 3))
      assert prn3.norad_id == 40294
      assert prn3.source.navcen.nanu_type == "FCSTSUMM"

      prn13 = Enum.find(records, &(&1.prn == 13))
      assert prn13.norad_id == 68791
      assert prn13.source.celestrak.block_type == "III"
      assert prn13.source.navcen_conflict.svn == 43
      assert prn13.source.navcen_conflict.block_type == "IIR"
      refute Map.has_key?(prn13.source, :navcen)
    end
  end

  describe "CSV export" do
    test "exports the compact mapping CSV with active=false for unusable rows" do
      assert Constellation.to_csv(merged_records()) ==
               """
               prn,norad_cat_id,active,sp3_id
               3,40294,true,G03
               5,35752,true,G05
               13,68791,true,G13
               19,28190,false,G19
               """
    end
  end

  describe "validation" do
    test "reports duplicate PRNs, duplicate NORAD ids, and inactive/unusable PRNs" do
      records = [
        %Record{
          system: :gps,
          prn: 3,
          svn: 69,
          norad_id: 40294,
          sp3_id: "G03",
          active?: true,
          usable?: true,
          source: %{}
        },
        %Record{
          system: :gps,
          prn: 3,
          svn: 70,
          norad_id: 40294,
          sp3_id: "G03",
          active?: false,
          usable?: true,
          source: %{}
        }
      ]

      report = Constellation.validate(records)
      assert report.duplicate_prns == [3]
      assert report.duplicate_norad_ids == [40294]
      assert report.inactive_unusable_prns == [3]
      refute Constellation.valid?(report)
    end

    test "compares active usable catalog ids against a loaded SP3 product" do
      {:ok, sp3} = SP3.parse(@sp3)
      assert SP3.satellite_ids(sp3) == ["G03", "G32"]

      report = Constellation.validate_sp3(merged_records(), sp3)
      assert report.missing_sp3_ids == ["G05", "G13"]
      assert report.extra_sp3_ids == ["G32"]
      assert report.inactive_unusable_prns == [19]
      refute Constellation.valid?(report)
    end

    test "accepts a plain SP3 id list for validation" do
      report = Constellation.validate_sp3(merged_records(), ["G03", "G05", "G13"])
      assert report.missing_sp3_ids == []
      assert report.extra_sp3_ids == []
      assert report.inactive_unusable_prns == [19]
    end
  end
end
