defmodule Orbis.GNSS.PreciseRealArcTest do
  @moduledoc false

  use ExUnit.Case, async: false

  alias Orbis.GNSS.IonosphereFree
  alias Orbis.GNSS.PrecisePositioning
  alias Orbis.GNSS.RINEX.Observations
  alias Orbis.GNSS.SP3

  @sp3_path Path.join(__DIR__, "fixtures/sp3/GBM0MGXRAP_20201770000_01D_05M_ORB_120epoch.sp3")
  @obs_path Path.join(__DIR__, "fixtures/obs/ESBC00DNK_R_20201770000_01D_30S_MO_120epoch.rnx")

  test "a real multi-epoch ionosphere-free arc improves with troposphere correction" do
    sp3 = SP3.load!(@sp3_path)
    obs = Observations.load!(@obs_path)
    truth = Observations.approx_position(obs)
    {x0, y0, z0} = truth

    epoch_observations = real_gps_iono_free_arc(obs, 120)
    assert length(epoch_observations) == 120

    opts = [
      initial_guess: {x0 + 100.0, y0 - 100.0, z0 + 100.0, 0.0},
      max_iterations: 8
    ]

    assert {:ok, uncorrected} =
             PrecisePositioning.solve_float_epochs(sp3, epoch_observations, opts)

    assert {:ok, corrected} =
             PrecisePositioning.solve_float_epochs(
               sp3,
               epoch_observations,
               Keyword.put(opts, :troposphere, true)
             )

    assert {:ok, ztd_estimated} =
             PrecisePositioning.solve_float_epochs(
               sp3,
               epoch_observations,
               Keyword.merge(opts, troposphere: true, estimate_ztd: true)
             )

    uncorrected_error_m = position_error(uncorrected.position, truth)
    corrected_error_m = position_error(corrected.position, truth)

    assert uncorrected.metadata.n_epochs == 120
    assert corrected.metadata.n_epochs == 120
    assert corrected.metadata.n_observations == 1282
    assert corrected.metadata.troposphere_applied
    assert ztd_estimated.metadata.troposphere_applied
    assert ztd_estimated.metadata.ztd_estimated
    assert ztd_estimated.ztd_residual_m > 0.0
    assert ztd_estimated.ztd_residual_m < 1.0

    # This is a real ESBC00DNK arc, not a synthetic zero-troposphere fixture. The
    # a-priori Saastamoinen/Niell slant correction should remove the dominant
    # bias, but the one-hour float arc is still not a cm-grade PPP gate.
    assert uncorrected_error_m > 20.0
    assert corrected_error_m < 5.0
    assert uncorrected_error_m - corrected_error_m > 18.0
    assert ztd_estimated.metadata.weighted_rms_m < corrected.metadata.weighted_rms_m
    assert ztd_estimated.metadata.phase_rms_m < corrected.metadata.phase_rms_m
  end

  defp real_gps_iono_free_arc(obs, count) do
    {:ok, f1} = IonosphereFree.frequency("G", :l1)
    {:ok, f2} = IonosphereFree.frequency("G", :l2)

    obs
    |> Observations.epochs()
    |> Enum.take(count)
    |> Enum.map(fn entry ->
      rows = epoch_rows(obs, entry.index, f1, f2)

      if length(rows) < 6 do
        raise "fixture epoch #{inspect(entry.epoch)} has only #{length(rows)} complete GPS L1/L2 code+phase rows"
      end

      %{epoch: naive_datetime(entry.epoch), observations: rows}
    end)
  end

  defp epoch_rows(obs, index, f1, f2) do
    {:ok, by_sat} =
      Observations.values(obs, index, codes: %{"G" => ["C1C", "C2W", "L1C", "L2W"]})

    by_sat
    |> Enum.flat_map(fn {sat, values} ->
      values_by_code = Map.new(values, &{&1.code, &1.value})

      with c1 when is_number(c1) <- values_by_code["C1C"],
           c2 when is_number(c2) <- values_by_code["C2W"],
           l1 when is_number(l1) <- values_by_code["L1C"],
           l2 when is_number(l2) <- values_by_code["L2W"],
           {:ok, code_m} <- IonosphereFree.iono_free(c1, c2, f1, f2),
           {:ok, phase_m} <- IonosphereFree.iono_free_phase_cycles(l1, l2, f1, f2) do
        [%{satellite_id: sat, code_m: code_m, phase_m: phase_m}]
      else
        _ -> []
      end
    end)
    |> Enum.sort_by(& &1.satellite_id)
  end

  defp naive_datetime({{year, month, day}, {hour, minute, second}}) do
    whole_second = trunc(second)
    microsecond = round((second - whole_second) * 1_000_000)

    NaiveDateTime.new!(
      Date.new!(year, month, day),
      Time.new!(hour, minute, whole_second, {microsecond, 6})
    )
  end

  defp position_error(%{x_m: x, y_m: y, z_m: z}, {truth_x, truth_y, truth_z}) do
    :math.sqrt(
      (x - truth_x) * (x - truth_x) +
        (y - truth_y) * (y - truth_y) +
        (z - truth_z) * (z - truth_z)
    )
  end
end
