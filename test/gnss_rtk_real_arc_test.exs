defmodule Orbis.GNSS.RTKRealArcTest do
  @moduledoc false

  use ExUnit.Case, async: false

  alias Orbis.GNSS.RINEX.Observations
  alias Orbis.GNSS.RTK
  alias Orbis.GNSS.SP3

  @sp3_path Path.join(__DIR__, "fixtures/sp3/GBM0MGXRAP_20201770000_01D_05M_ORB_120epoch.sp3")
  @wtzr_obs_path Path.join(
                   __DIR__,
                   "fixtures/obs/WTZR00DEU_R_20201770000_01D_30S_MO_120epoch.rnx"
                 )
  @wtzz_obs_path Path.join(
                   __DIR__,
                   "fixtures/obs/WTZZ00DEU_R_20201770000_01D_30S_MO_120epoch.rnx"
                 )

  @c_m_s 299_792_458.0
  @gps_l1_hz 1_575_420_000.0
  @gps_l1_wavelength_m @c_m_s / @gps_l1_hz

  @wtzr_marker {4_075_580.3111, 931_854.0543, 4_801_568.2808}
  @wtzz_marker {4_075_579.1913, 931_853.3696, 4_801_569.1897}
  @wtzr_antenna_h_m 0.0710
  @wtzz_antenna_h_m 0.2840

  @tag timeout: 180_000
  test "real co-located Wettzell L1 RTK solves the antenna baseline and refuses an unsafe fix" do
    sp3 = SP3.load!(@sp3_path)
    base_obs = Observations.load!(@wtzr_obs_path)
    rover_obs = Observations.load!(@wtzz_obs_path)

    base_arp = arp_position(@wtzr_marker, @wtzr_antenna_h_m)
    rover_arp = arp_position(@wtzz_marker, @wtzz_antenna_h_m)
    marker_baseline = sub3(@wtzz_marker, @wtzr_marker)
    antenna_baseline = sub3(rover_arp, base_arp)
    epochs = real_gps_l1_rtk_epochs(sp3, base_obs, rover_obs, 120)

    assert length(epochs) == 120

    opts = [
      initial_baseline_m: {0.0, 0.0, 0.0},
      max_iterations: 10,
      on_cycle_slip: :split_arc,
      elevation_weighting: true,
      code_sigma_m: 2.0,
      phase_sigma_m: 0.01,
      ambiguity_wavelength_m: @gps_l1_wavelength_m,
      integer_candidate_limit: 200_000
    ]

    assert {:ok, float} = RTK.solve_float_baseline_epochs(base_arp, epochs, opts)

    float_antenna_error_m = position_error(float.baseline_m, antenna_baseline)
    float_marker_error_m = position_error(float.baseline_m, marker_baseline)

    assert float.metadata.n_epochs == 120
    assert float.metadata.measurement_covariance.elevation_weighting
    assert length(float.metadata.split_cycle_slip_arcs) == 4
    assert float.metadata.dropped_sats == ["G09", "G21", "G27"]
    assert float.metadata.phase_rms_m < 0.01

    # The SSC coordinates are marker coordinates; the observations are tied to
    # the antenna reference points. Applying the RINEX antenna-height deltas is
    # the difference between a decimetre-looking residual and a real
    # centimetre-scale short-baseline float RTK result.
    assert float_marker_error_m > 0.15
    assert float_antenna_error_m < 0.08

    assert {:ok, fixed} = RTK.solve_fixed_baseline_epochs(base_arp, epochs, opts)

    assert fixed.metadata.integer_status == :not_fixed
    assert fixed.metadata.integer_ratio < 3.0
    assert fixed.metadata.integer_candidates > 0

    fixed_antenna_error_m = position_error(fixed.baseline_m, antenna_baseline)

    # Real L1 phase on this one-hour arc is not separable enough for an integer
    # fix; the important behavior is that LAMBDA refuses the unsafe fix instead
    # of reporting centimetre confidence from noisy data.
    assert fixed_antenna_error_m < 0.2
  end

  defp real_gps_l1_rtk_epochs(sp3, base_obs, rover_obs, count) do
    rover_by_epoch = Map.new(Observations.epochs(rover_obs), &{&1.epoch, &1})

    base_obs
    |> Observations.epochs()
    |> Enum.take(count)
    |> Enum.flat_map(fn base_entry ->
      case Map.fetch(rover_by_epoch, base_entry.epoch) do
        {:ok, rover_entry} ->
          base_values = gps_l1_values(base_obs, base_entry.index)
          rover_values = gps_l1_values(rover_obs, rover_entry.index)

          common =
            base_values
            |> Map.keys()
            |> MapSet.new()
            |> MapSet.intersection(rover_values |> Map.keys() |> MapSet.new())
            |> MapSet.to_list()
            |> Enum.sort()

          epoch = naive_datetime(base_entry.epoch)
          positions = satellite_positions(sp3, epoch, common)
          usable = Enum.filter(common, &Map.has_key?(positions, &1))

          if length(usable) >= 4 do
            [
              %{
                epoch: epoch,
                satellite_positions_m: Map.take(positions, usable),
                base_observations: Enum.map(usable, &Map.fetch!(base_values, &1)),
                rover_observations: Enum.map(usable, &Map.fetch!(rover_values, &1))
              }
            ]
          else
            []
          end

        :error ->
          []
      end
    end)
  end

  defp gps_l1_values(obs, index) do
    {:ok, by_sat} = Observations.values(obs, index, codes: %{"G" => ["C1C", "L1C"]})

    by_sat
    |> Enum.flat_map(fn {sat, values} ->
      values_by_code = Map.new(values, &{&1.code, &1})

      with %{value: c1} when is_number(c1) <- values_by_code["C1C"],
           %{value: l1} = phase when is_number(l1) <- values_by_code["L1C"] do
        [
          {sat,
           %{
             satellite_id: sat,
             code_m: c1,
             phase_m: l1 * @gps_l1_wavelength_m,
             lli: phase.lli
           }}
        ]
      else
        _ -> []
      end
    end)
    |> Map.new()
  end

  defp satellite_positions(sp3, epoch, sats) do
    sats
    |> Enum.reduce(%{}, fn sat, acc ->
      case SP3.position(sp3, sat, epoch) do
        {:ok, %{x_m: x, y_m: y, z_m: z}} -> Map.put(acc, sat, {x, y, z})
        {:error, _reason} -> acc
      end
    end)
  end

  defp naive_datetime({{year, month, day}, {hour, minute, second}}) do
    whole_second = trunc(second)
    microsecond = round((second - whole_second) * 1_000_000)

    NaiveDateTime.new!(
      Date.new!(year, month, day),
      Time.new!(hour, minute, whole_second, {microsecond, 6})
    )
  end

  defp arp_position(marker, antenna_h_m),
    do: add3(marker, scale3(marker, antenna_h_m / norm3(marker)))

  defp add3({ax, ay, az}, {bx, by, bz}), do: {ax + bx, ay + by, az + bz}
  defp sub3({ax, ay, az}, {bx, by, bz}), do: {ax - bx, ay - by, az - bz}
  defp scale3({x, y, z}, s), do: {x * s, y * s, z * s}
  defp norm3({x, y, z}), do: :math.sqrt(x * x + y * y + z * z)

  defp position_error(%{x_m: x, y_m: y, z_m: z}, {truth_x, truth_y, truth_z}) do
    :math.sqrt(
      (x - truth_x) * (x - truth_x) +
        (y - truth_y) * (y - truth_y) +
        (z - truth_z) * (z - truth_z)
    )
  end
end
