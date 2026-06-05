defmodule Orbis.GNSS.CarrierPhaseOracleTest do
  use ExUnit.Case, async: true

  import Orbis.TestHelpers, only: [assert_ulp: 4, hex_to_float: 1]

  alias Orbis.GNSS.CarrierPhase
  alias Orbis.GNSS.RINEX.Observations

  @obs_path Path.join(__DIR__, "fixtures/obs/ESBC00DNK_R_20201770000_01D_30S_MO_trim.crx")
  @golden_path Path.join(__DIR__, "fixtures/carrier_phase_golden.json")

  setup_all do
    golden = @golden_path |> File.read!() |> Jason.decode!()
    {:ok, obs: Observations.load!(@obs_path), golden: golden}
  end

  describe "RINEX observation values vs georinex oracle" do
    test "values/3 matches georinex for selected multi-system raw observations", %{
      obs: obs,
      golden: golden
    } do
      for row <- golden["rinex_observations"]["rows"] do
        system = String.first(row["sat"])

        {:ok, values} =
          Observations.values(obs, row["epoch_index"], codes: %{system => [row["code"]]})

        actual =
          values
          |> Map.fetch!(row["sat"])
          |> Enum.find(&(&1.code == row["code"]))

        assert actual.kind == String.to_existing_atom(row["kind"])
        assert actual.units == String.to_existing_atom(row["units"])
        assert actual.lli == row["lli"]
        assert actual.ssi == row["ssi"]
        assert_float_or_nil(actual.value, row["value"], "#{row["sat"]} #{row["code"]}")
      end
    end

    test "phases/3 wavelengths and metre values match the oracle fixture", %{
      obs: obs,
      golden: golden
    } do
      for row <- golden["rinex_observations"]["phases"] do
        system = String.first(row["sat"])

        {:ok, values} =
          Observations.phases(obs, row["epoch_index"], codes: %{system => [row["code"]]})

        actual =
          values
          |> Map.fetch!(row["sat"])
          |> Enum.find(&(&1.code == row["code"]))

        assert actual.lli == row["lli"]
        assert actual.ssi == row["ssi"]

        assert_float_or_nil(
          actual.value_cycles,
          row["value_cycles"],
          "#{row["sat"]} #{row["code"]} cycles"
        )

        assert_float_or_nil(
          actual.wavelength_m,
          row["wavelength_m"],
          "#{row["sat"]} #{row["code"]} lambda"
        )

        assert_float_or_nil(actual.value_m, row["value_m"], "#{row["sat"]} #{row["code"]} metres")
      end
    end
  end

  describe "CarrierPhase combinations vs Python oracle" do
    test "scalar combinations match the Python fixture", %{golden: golden} do
      scalars = golden["carrier_phase"]["scalar_cases"]
      f1 = h(golden["carrier_phase"]["constants"]["f_l1_hz"])
      f2 = h(golden["carrier_phase"]["constants"]["f_l2_hz"])

      assert_ulp(
        CarrierPhase.geometry_free(100.0, 60.0),
        h(scalars["geometry_free_m"]),
        0,
        "GF scalar"
      )

      assert {:ok, lambda_wl} = CarrierPhase.wide_lane_wavelength(f1, f2)
      assert_ulp(lambda_wl, h(scalars["wide_lane_wavelength_m"]), 0, "wide-lane wavelength")

      assert {:ok, p_nl} = CarrierPhase.narrow_lane_code(10.0, 12.0, f1, f2)
      assert_ulp(p_nl, h(scalars["narrow_lane_code_m"]), 0, "narrow-lane code")

      assert {:ok, mw} = CarrierPhase.melbourne_wubbena(5.0, 3.0, 10.0, 12.0, f1, f2)
      assert_ulp(mw, h(scalars["melbourne_wubbena_m"]), 0, "Melbourne-Wubbena")
    end

    test "cycle-slip classification matches the Python arc oracle", %{golden: golden} do
      arc = decode_arc(golden["carrier_phase"]["arc"])
      actual = CarrierPhase.detect_cycle_slips(arc)
      expected = golden["carrier_phase"]["detect_cycle_slips"]

      for {a, e} <- Enum.zip(actual, expected) do
        assert a.epoch == e["epoch"]
        assert a.slip == e["slip"]
        assert Enum.map(a.reasons, &Atom.to_string/1) == e["reasons"]
        assert a.skipped == e["skipped"]
        assert_float_or_nil(a.gf, e["gf"], "GF epoch #{e["epoch"]}")
        assert_float_or_nil(a.mw, e["mw"], "MW epoch #{e["epoch"]}")
      end
    end

    test "Hatch-smoothed code matches the Python arc oracle", %{golden: golden} do
      arc = decode_arc(golden["carrier_phase"]["arc"])
      actual = CarrierPhase.smooth_code(arc)
      expected = golden["carrier_phase"]["smooth_code"]

      for {a, e} <- Enum.zip(actual, expected) do
        assert a.epoch == e["epoch"]
        assert a.window == e["window"]
        assert a.reset == e["reset"]
        assert_float_or_nil(a.p_smooth, e["p_smooth"], "Hatch epoch #{e["epoch"]}")
      end
    end
  end

  defp decode_arc(rows) do
    Enum.map(rows, fn row ->
      %{
        epoch: row["epoch"],
        phi1: h(row["phi1"]),
        phi2: h(row["phi2"]),
        p1: h(row["p1"]),
        p2: h(row["p2"]),
        lli1: row["lli1"],
        lli2: row["lli2"],
        f1: h(row["f1"]),
        f2: h(row["f2"])
      }
    end)
  end

  defp assert_float_or_nil(actual, nil, _label), do: assert(actual == nil)

  defp assert_float_or_nil(actual, expected_hex, label) do
    assert is_number(actual), "#{label}: expected numeric actual, got #{inspect(actual)}"
    assert_ulp(actual, h(expected_hex), 0, label)
  end

  defp h(nil), do: nil
  defp h(hex), do: hex_to_float(hex)
end
