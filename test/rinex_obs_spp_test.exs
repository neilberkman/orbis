defmodule Orbis.RinexObsSppTest do
  @moduledoc """
  End-to-end single-point positioning from a real station observation file: load
  the CRINEX, extract single-frequency pseudoranges, solve against the matching
  broadcast navigation product, and recover the receiver's surveyed position.

  This proves the whole last mile — CRINEX decode, RINEX 3 observation parse,
  pseudorange extraction, and the solve — on real data for the ESBC00DNK station
  (Esbjerg, Denmark) at 2020-06-25 00:00 GPST.
  """
  use ExUnit.Case, async: true

  alias Orbis.BroadcastEphemeris
  alias Orbis.PointPositioning
  alias Orbis.RinexObs

  @obs_path Path.join(__DIR__, "fixtures/obs/ESBC00DNK_R_20201770000_01D_30S_MO_trim.crx")
  @nav_path Path.join(__DIR__, "fixtures/nav/ESBC00DNK_R_20201770000_01D_MN.rnx")

  # GPS broadcast Klobuchar coefficients from the committed NAV header (GPSA/GPSB).
  @gps_alpha {4.6566e-09, 1.4901e-08, -5.9605e-08, -1.1921e-07}
  @gps_beta {8.1920e+04, 9.8304e+04, -6.5536e+04, -5.2429e+05}

  test "recovers the surveyed station position from real GPS observations" do
    obs = RinexObs.load!(@obs_path)
    eph = BroadcastEphemeris.load!(@nav_path)

    {truth_x, truth_y, truth_z} = RinexObs.approx_position(obs)

    [%{index: index, epoch: epoch} | _] = RinexObs.epochs(obs)

    # Default GPS code (C1C), single frequency.
    {:ok, prs} = RinexObs.pseudoranges(obs, index, codes: %{"G" => ["C1C"]})

    # An over-determined GPS-only set at epoch 0.
    assert length(prs) >= 6

    # Seed with a coarse a-priori position ~45 km off the surveyed truth (not the
    # answer), so the test demonstrates convergence to the surveyed position
    # rather than assuming it. The solver freezes its elevation mask and weights
    # at the initial geometry, so the seed must be a near-surface point.
    coarse_guess = {truth_x + 30_000.0, truth_y - 20_000.0, truth_z + 25_000.0, 0.0}

    assert {:ok, sol} =
             PointPositioning.solve(eph, prs, epoch,
               ionosphere: true,
               troposphere: true,
               klobuchar_alpha: @gps_alpha,
               klobuchar_beta: @gps_beta,
               initial_guess: coarse_guess
             )

    assert sol.metadata.converged

    # Single-frequency broadcast SPP is metre-level; assert each axis is within a
    # few metres of the surveyed APPROX POSITION XYZ.
    assert_in_delta sol.position.x_m, truth_x, 5.0
    assert_in_delta sol.position.y_m, truth_y, 5.0
    assert_in_delta sol.position.z_m, truth_z, 5.0

    # The 3D position error is metre-level (the real solve lands within ~3 m).
    err =
      :math.sqrt(
        :math.pow(sol.position.x_m - truth_x, 2) +
          :math.pow(sol.position.y_m - truth_y, 2) +
          :math.pow(sol.position.z_m - truth_z, 2)
      )

    assert err < 5.0
  end

  test "pseudoranges/3 accepts an epoch tuple as well as an index" do
    obs = RinexObs.load!(@obs_path)
    [%{index: index, epoch: epoch} | _] = RinexObs.epochs(obs)

    {:ok, by_index} = RinexObs.pseudoranges(obs, index, codes: %{"G" => ["C1C"]})
    {:ok, by_tuple} = RinexObs.pseudoranges(obs, epoch, codes: %{"G" => ["C1C"]})

    assert by_index == by_tuple
  end

  test "pseudoranges/3 reports an out-of-range epoch index" do
    obs = RinexObs.load!(@obs_path)
    assert {:error, :epoch_out_of_range} = RinexObs.pseudoranges(obs, 9_999)
  end
end
