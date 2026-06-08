defmodule Orbis.GNSS.SP3MergeTest do
  use ExUnit.Case, async: true

  alias Orbis.GNSS.SP3

  # Build a single-epoch SP3-c buffer from explicit
  # `{satellite_token, [x_km, y_km, z_km], clock_us | nil}` records, so each test
  # controls which satellites a "center" reports and where. Mirrors the crate's
  # `sp3_records` test helper.
  defp sp3_records(records) do
    n = length(records)

    sats =
      Enum.map_join(records, "", fn {sat, _, _} -> sat end) <>
        String.duplicate("  0", 17 - n)

    header = [
      "#cP2020  6 24  0  0  0.00000000       1 ORBIT IGS14 FIT  TST",
      "## 2111 432000.00000000   900.00000000 59024 0.0000000000000",
      "+   #{String.pad_leading(Integer.to_string(n), 2)}   #{sats}",
      "++         0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0",
      "%c G  cc GPS ccc cccc cccc cccc cccc ccccc ccccc ccccc ccccc",
      "%c cc cc ccc ccc cccc cccc cccc cccc ccccc ccccc ccccc ccccc",
      "%f  1.2500000  1.025000000  0.00000000000  0.000000000000000",
      "%f  0.0000000  0.000000000  0.00000000000  0.000000000000000",
      "%i    0    0    0    0      0      0      0      0         0",
      "%i    0    0    0    0      0      0      0      0         0",
      "/* TEST SP3-c FIXTURE",
      "*  2020  6 24  0  0  0.00000000"
    ]

    recs =
      Enum.map(records, fn {sat, [x, y, z], clk} ->
        c = clk || 999_999.999999
        "P" <> sat <> fmt(x) <> fmt(y) <> fmt(z) <> fmt(c)
      end)

    {:ok, sp3} = SP3.parse(Enum.join(header ++ recs ++ ["EOF", ""], "\n"))
    sp3
  end

  defp fmt(v), do: :io_lib.format(~c"~14.6f", [v]) |> IO.iodata_to_binary()

  describe "merge/2" do
    test "union coverage: merged product covers a satellite a center is missing" do
      a =
        sp3_records([
          {"G01", [15000.0, -20000.0, 5000.0], 100.0},
          {"G02", [16000.0, -21000.0, 6000.0], 200.0},
          {"G03", [17000.0, -22000.0, 7000.0], 300.0}
        ])

      b =
        sp3_records([
          {"G01", [15000.0, -20000.0, 5000.0], 100.0},
          {"G02", [16000.0, -21000.0, 6000.0], 200.0}
        ])

      assert {:ok, merged, report} = SP3.merge([a, b])

      ids = SP3.satellite_ids(merged)
      assert "G03" in ids, "merged output must cover G03 from the center that has it"
      assert Enum.sort(ids) == ["G01", "G02", "G03"]

      assert report.quarantined == []
      # G03 had a single source (index 0) -> carried through, recorded.
      assert [%{satellite: "G03", sources: [0]}] = report.single_source
    end

    test "quarantines a satellite all centers disagree on" do
      # Three centers, mutually beyond the default 0.5 m tolerance on G01.
      a = sp3_records([{"G01", [15000.000, -20000.0, 5000.0], 100.0}])
      b = sp3_records([{"G01", [15000.010, -20000.0, 5000.0], 100.0}])
      c = sp3_records([{"G01", [15000.020, -20000.0, 5000.0], 100.0}])

      assert {:ok, merged, report} = SP3.merge([a, b, c])

      refute "G01" in SP3.satellite_ids(merged),
             "no consensus -> G01 omitted, not averaged across disagreeing centers"

      assert [%{satellite: "G01"}] = report.quarantined
    end

    test "rejects an outlier and combines the agreeing centers" do
      # A and B agree on G01; C is 10 m off in X.
      a = sp3_records([{"G01", [15000.000, -20000.0, 5000.0], 100.0}])
      b = sp3_records([{"G01", [15000.000, -20000.0, 5000.0], 100.0}])
      c = sp3_records([{"G01", [15000.010, -20000.0, 5000.0], 100.0}])

      assert {:ok, merged, report} = SP3.merge([a, b, c])

      assert "G01" in SP3.satellite_ids(merged)
      assert [%{satellite: "G01", sources: [2]}] = report.position_outliers
      assert report.quarantined == []
    end
  end

  describe "clock_reference_offset/3 and align_clock_reference/3" do
    defp shifted_pair do
      pos = [
        {"G01", [15000.0, -20000.0, 5000.0]},
        {"G02", [16000.0, -21000.0, 6000.0]},
        {"G03", [17000.0, -22000.0, 7000.0]}
      ]

      a = sp3_records(Enum.map(pos, fn {s, p} -> {s, p, 100.0} end))
      # `b`'s clocks all run +50 us (= 5e-5 s) ahead of `a`'s.
      b = sp3_records(Enum.map(pos, fn {s, p} -> {s, p, 150.0} end))
      {a, b}
    end

    test "clock_reference_offset recovers a uniform datum shift" do
      {a, b} = shifted_pair()

      assert [offset] = SP3.clock_reference_offset(a, b, min_common: 3)
      assert offset.satellites == 3
      assert_in_delta offset.offset_s, 5.0e-5, 1.0e-12
    end

    test "align_clock_reference removes the datum (residual offset ~ 0)" do
      {a, b} = shifted_pair()

      assert {:ok, aligned} = SP3.align_clock_reference(a, b, min_common: 3)
      assert [residual] = SP3.clock_reference_offset(a, aligned, min_common: 3)
      assert_in_delta residual.offset_s, 0.0, 1.0e-12
    end
  end
end
