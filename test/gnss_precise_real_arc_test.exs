defmodule Orbis.GNSS.PreciseRealArcTest do
  @moduledoc false

  use ExUnit.Case, async: false

  alias Orbis.GNSS.IonosphereFree
  alias Orbis.GNSS.PrecisePositioning
  alias Orbis.GNSS.RINEX.Observations
  alias Orbis.GNSS.SP3

  @sp3_path Path.join(__DIR__, "fixtures/sp3/GBM0MGXRAP_20201770000_01D_05M_ORB_120epoch.sp3")
  @obs_path Path.join(__DIR__, "fixtures/obs/ESBC00DNK_R_20201770000_01D_30S_MO_120epoch.rnx")

  @tag timeout: 180_000
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

    assert {:ok, elevation_weighted} =
             PrecisePositioning.solve_float_epochs(
               sp3,
               epoch_observations,
               Keyword.put(opts, :elevation_weighting, true)
             )

    elevation_weighted_error_m = position_error(elevation_weighted.position, truth)

    # Elevation weighting is a stochastic-model improvement, not an integer-fix
    # breakthrough. On this no-troposphere real arc it materially improves the
    # float position and weighted RMS by down-weighting low-elevation rows.
    assert elevation_weighted_error_m < uncorrected_error_m / 3.0
    assert elevation_weighted.metadata.weighted_rms_m < uncorrected.metadata.weighted_rms_m / 5.0
  end

  @tag timeout: 180_000
  test "real noisy narrow-lane fixing returns a not-fixed candidate set" do
    sp3 = SP3.load!(@sp3_path)
    obs = Observations.load!(@obs_path)
    {x0, y0, z0} = Observations.approx_position(obs)

    # G21 has a detected slip on this arc; exclude it so this test exercises the
    # noisy narrow-lane integer-search path rather than the slip hard-abort.
    epoch_observations = real_gps_iono_free_arc(obs, 120, exclude: ["G21"])

    {:ok, f1} = IonosphereFree.frequency("G", :l1)
    {:ok, f2} = IonosphereFree.frequency("G", :l2)

    assert {:ok, fixed} =
             PrecisePositioning.solve_fixed_epochs(sp3, epoch_observations,
               initial_guess: {x0 + 100.0, y0 - 100.0, z0 + 100.0, 0.0},
               max_iterations: 8,
               troposphere: true,
               ambiguity_wavelength_m: 299_792_458.0 / (f1 + f2)
             )

    assert fixed.metadata.integer_status == :not_fixed
    assert fixed.metadata.integer_ratio < 3.0
    assert fixed.metadata.integer_candidates == 2
    assert position_error(fixed.position, {x0, y0, z0}) < 6.0
  end

  @tag timeout: 180_000
  test "real elevation-weighted narrow-lane fixing still refuses an unsafe fix" do
    sp3 = SP3.load!(@sp3_path)
    obs = Observations.load!(@obs_path)
    truth = Observations.approx_position(obs)
    {x0, y0, z0} = truth

    # Keep G21 excluded so this test measures noisy narrow-lane separability
    # rather than cycle-slip handling.
    epoch_observations = real_gps_iono_free_arc(obs, 120, exclude: ["G21"])

    {:ok, f1} = IonosphereFree.frequency("G", :l1)
    {:ok, f2} = IonosphereFree.frequency("G", :l2)

    opts = [
      initial_guess: {x0 + 100.0, y0 - 100.0, z0 + 100.0, 0.0},
      max_iterations: 8,
      ambiguity_wavelength_m: 299_792_458.0 / (f1 + f2)
    ]

    assert {:ok, unweighted} =
             PrecisePositioning.solve_fixed_epochs(sp3, epoch_observations, opts)

    assert {:ok, weighted} =
             PrecisePositioning.solve_fixed_epochs(
               sp3,
               epoch_observations,
               Keyword.put(opts, :elevation_weighting, true)
             )

    assert unweighted.metadata.integer_status == :not_fixed
    assert weighted.metadata.integer_status == :not_fixed
    assert weighted.metadata.integer_ratio < 3.0
    assert weighted.metadata.integer_candidates == 2

    assert position_error(weighted.position, truth) <
             position_error(unweighted.position, truth) / 2.0
  end

  @tag timeout: 180_000
  test "real slipped arcs can be split before the narrow-lane search" do
    sp3 = SP3.load!(@sp3_path)
    obs = Observations.load!(@obs_path)
    {x0, y0, z0} = Observations.approx_position(obs)
    dual_epoch_observations = real_gps_dual_frequency_arc(obs, 30)

    assert {:error, {:cycle_slip_detected, "G21", _epoch, _reasons}} =
             PrecisePositioning.solve_widelane_fixed_epochs(sp3, dual_epoch_observations,
               initial_guess: {x0 + 100.0, y0 - 100.0, z0 + 100.0, 0.0},
               max_iterations: 5,
               troposphere: true
             )

    assert {:ok, fixed} =
             PrecisePositioning.solve_widelane_fixed_epochs(sp3, dual_epoch_observations,
               initial_guess: {x0 + 100.0, y0 - 100.0, z0 + 100.0, 0.0},
               max_iterations: 5,
               troposphere: true,
               on_cycle_slip: :split_arc,
               wide_lane_tolerance_cycles: 2.0,
               integer_candidate_limit: 2_000_000
             )

    assert fixed.metadata.integer_status == :not_fixed
    assert fixed.metadata.integer_ratio < 3.0
    assert fixed.metadata.integer_candidates == 2
    assert length(fixed.metadata.split_cycle_slip_arcs) == 2
    assert position_error(fixed.position, {x0, y0, z0}) < 9.0
  end

  defp real_gps_iono_free_arc(obs, count, opts \\ []) do
    {:ok, f1} = IonosphereFree.frequency("G", :l1)
    {:ok, f2} = IonosphereFree.frequency("G", :l2)
    excluded = opts |> Keyword.get(:exclude, []) |> MapSet.new()

    obs
    |> Observations.epochs()
    |> Enum.take(count)
    |> Enum.map(fn entry ->
      rows = epoch_rows(obs, entry.index, f1, f2, excluded)

      if length(rows) < 6 do
        raise "fixture epoch #{inspect(entry.epoch)} has only #{length(rows)} complete GPS L1/L2 code+phase rows"
      end

      %{epoch: naive_datetime(entry.epoch), observations: rows}
    end)
  end

  defp real_gps_dual_frequency_arc(obs, count) do
    {:ok, f1} = IonosphereFree.frequency("G", :l1)
    {:ok, f2} = IonosphereFree.frequency("G", :l2)

    obs
    |> Observations.epochs()
    |> Enum.take(count)
    |> Enum.map(fn entry ->
      rows = dual_epoch_rows(obs, entry.index, f1, f2)

      if length(rows) < 6 do
        raise "fixture epoch #{inspect(entry.epoch)} has only #{length(rows)} complete GPS L1/L2 code+phase rows"
      end

      %{epoch: naive_datetime(entry.epoch), observations: rows}
    end)
  end

  defp epoch_rows(obs, index, f1, f2, excluded) do
    {:ok, by_sat} =
      Observations.values(obs, index, codes: %{"G" => ["C1C", "C2W", "L1C", "L2W"]})

    by_sat
    |> Enum.flat_map(fn {sat, values} ->
      values_by_code = Map.new(values, &{&1.code, &1.value})

      with c1 when is_number(c1) <- values_by_code["C1C"],
           c2 when is_number(c2) <- values_by_code["C2W"],
           l1 when is_number(l1) <- values_by_code["L1C"],
           l2 when is_number(l2) <- values_by_code["L2W"],
           false <- MapSet.member?(excluded, sat),
           {:ok, code_m} <- IonosphereFree.iono_free(c1, c2, f1, f2),
           {:ok, phase_m} <- IonosphereFree.iono_free_phase_cycles(l1, l2, f1, f2) do
        [%{satellite_id: sat, code_m: code_m, phase_m: phase_m}]
      else
        _ -> []
      end
    end)
    |> Enum.sort_by(& &1.satellite_id)
  end

  defp dual_epoch_rows(obs, index, f1, f2) do
    {:ok, by_sat} =
      Observations.values(obs, index, codes: %{"G" => ["C1C", "C2W", "L1C", "L2W"]})

    by_sat
    |> Enum.flat_map(fn {sat, values} ->
      values_by_code = Map.new(values, &{&1.code, &1.value})

      with c1 when is_number(c1) <- values_by_code["C1C"],
           c2 when is_number(c2) <- values_by_code["C2W"],
           l1 when is_number(l1) <- values_by_code["L1C"],
           l2 when is_number(l2) <- values_by_code["L2W"] do
        [
          %{
            satellite_id: sat,
            p1_m: c1,
            p2_m: c2,
            phi1_cyc: l1,
            phi2_cyc: l2,
            f1_hz: f1,
            f2_hz: f2
          }
        ]
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
