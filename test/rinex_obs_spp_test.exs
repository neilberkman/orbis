defmodule Orbis.GNSS.RINEX.ObservationsSppTest do
  @moduledoc """
  End-to-end single-point positioning from a real station observation file: load
  the CRINEX, extract single-frequency pseudoranges, solve against the matching
  broadcast navigation product, and recover the receiver's surveyed position.

  This proves the whole last mile — CRINEX decode, RINEX 3 observation parse,
  pseudorange extraction, and the solve — on real data for the ESBC00DNK station
  (Esbjerg, Denmark) at 2020-06-25 00:00 GPST.
  """
  use ExUnit.Case, async: true

  alias Orbis.GNSS.Broadcast
  alias Orbis.GNSS.Positioning
  alias Orbis.GNSS.RINEX.Observations

  @obs_path Path.join(__DIR__, "fixtures/obs/ESBC00DNK_R_20201770000_01D_30S_MO_trim.crx")
  @nav_path Path.join(__DIR__, "fixtures/nav/ESBC00DNK_R_20201770000_01D_MN.rnx")

  # GPS broadcast Klobuchar coefficients from the committed NAV header (GPSA/GPSB).
  @gps_alpha {4.6566e-09, 1.4901e-08, -5.9605e-08, -1.1921e-07}
  @gps_beta {8.1920e+04, 9.8304e+04, -6.5536e+04, -5.2429e+05}

  test "recovers the surveyed station position from real GPS observations" do
    obs = Observations.load!(@obs_path)
    eph = Broadcast.load!(@nav_path)

    {truth_x, truth_y, truth_z} = Observations.approx_position(obs)

    [%{index: index, epoch: epoch} | _] = Observations.epochs(obs)

    # Default GPS code (C1C), single frequency.
    {:ok, prs} = Observations.pseudoranges(obs, index, codes: %{"G" => ["C1C"]})

    # An over-determined GPS-only set at epoch 0.
    assert length(prs) >= 6

    # Seed with a coarse a-priori position ~45 km off the surveyed truth (not the
    # answer), so the test demonstrates convergence to the surveyed position
    # rather than assuming it. The solver freezes its elevation mask and weights
    # at the initial geometry, so the seed must be a near-surface point.
    coarse_guess = {truth_x + 30_000.0, truth_y - 20_000.0, truth_z + 25_000.0, 0.0}

    assert {:ok, sol} =
             Positioning.solve(eph, prs, epoch,
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

  describe "coarse cold-start (:coarse_search)" do
    # The convergence-basin capability: recover the surveyed position from a
    # degraded or absent prior, with no hardcoded seed. The truth is the same
    # RINEX APPROX POSITION the first test trusts; the declared tolerance is the
    # pre-registered 5 m bar (single-frequency SPP floor here is about 2 m).
    setup do
      obs = Observations.load!(@obs_path)
      eph = Broadcast.load!(@nav_path)
      {truth_x, truth_y, truth_z} = Observations.approx_position(obs)
      [%{index: index, epoch: epoch} | _] = Observations.epochs(obs)
      {:ok, prs} = Observations.pseudoranges(obs, index, codes: %{"G" => ["C1C"]})
      assert length(prs) >= 6

      base_opts = [
        ionosphere: true,
        troposphere: true,
        klobuchar_alpha: @gps_alpha,
        klobuchar_beta: @gps_beta
      ]

      {:ok,
       eph: eph, prs: prs, epoch: epoch, truth: {truth_x, truth_y, truth_z}, base_opts: base_opts}
    end

    defp err_3d(%{position: p}, {tx, ty, tz}) do
      :math.sqrt((p.x_m - tx) ** 2 + (p.y_m - ty) ** 2 + (p.z_m - tz) ** 2)
    end

    test "converges from the earth-center default prior with no hardcoded seed", ctx do
      # Without :coarse_search the {0,0,0,0} default returns the seed flagged
      # converged at ~6.36e6 m (the kernel step-tolerance test fires on iteration
      # 0). The coarse search must instead land on the true fix.
      assert {:ok, baseline} =
               Positioning.solve(ctx.eph, ctx.prs, ctx.epoch, ctx.base_opts)

      assert err_3d(baseline, ctx.truth) > 1.0e6

      assert {:ok, sol} =
               Positioning.solve(
                 ctx.eph,
                 ctx.prs,
                 ctx.epoch,
                 Keyword.put(ctx.base_opts, :coarse_search, true)
               )

      assert sol.metadata.converged
      assert length(sol.used_sats) >= 5
      assert err_3d(sol, ctx.truth) <= 5.0
    end

    test "converges from an antipodal prior that starves the single solve", ctx do
      {tx, ty, tz} = ctx.truth
      antipodal = {-tx, -ty, -tz, 0.0}

      # The single solve from the antipodal seed drops every satellite below the
      # frozen horizon and returns too-few-satellites.
      assert {:error, {:too_few_satellites, 0, 4}} =
               Positioning.solve(
                 ctx.eph,
                 ctx.prs,
                 ctx.epoch,
                 Keyword.put(ctx.base_opts, :initial_guess, antipodal)
               )

      # The coarse search recovers the true fix regardless of the bad prior.
      assert {:ok, sol} =
               Positioning.solve(
                 ctx.eph,
                 ctx.prs,
                 ctx.epoch,
                 ctx.base_opts
                 |> Keyword.put(:initial_guess, antipodal)
                 |> Keyword.put(:coarse_search, true)
               )

      assert sol.metadata.converged
      assert err_3d(sol, ctx.truth) <= 5.0
    end

    test "an explicit seed count is honored and meets the bar", ctx do
      for n <- [12, 24] do
        assert {:ok, sol} =
                 Positioning.solve(
                   ctx.eph,
                   ctx.prs,
                   ctx.epoch,
                   Keyword.put(ctx.base_opts, :coarse_search, seeds: n)
                 )

        assert sol.metadata.converged
        assert length(sol.used_sats) >= 5
        assert err_3d(sol, ctx.truth) <= 5.0
      end
    end

    test "coarse_search off is identical to the plain single solve", ctx do
      good_prior =
        Keyword.put(ctx.base_opts, :initial_guess, {
          elem(ctx.truth, 0) + 30_000.0,
          elem(ctx.truth, 1) - 20_000.0,
          elem(ctx.truth, 2) + 25_000.0,
          0.0
        })

      assert {:ok, plain} = Positioning.solve(ctx.eph, ctx.prs, ctx.epoch, good_prior)

      assert {:ok, off} =
               Positioning.solve(
                 ctx.eph,
                 ctx.prs,
                 ctx.epoch,
                 Keyword.put(good_prior, :coarse_search, nil)
               )

      # Byte-for-byte identical: the off path runs exactly one solve from the
      # caller's seed, with no scorer or extra seeds in the way.
      assert off.position == plain.position
      assert off.rx_clock_s == plain.rx_clock_s
      assert off.residuals_m == plain.residuals_m
      assert off.used_sats == plain.used_sats
    end

    test "an invalid :coarse_search value is rejected", ctx do
      assert_raise ArgumentError, fn ->
        Positioning.solve(
          ctx.eph,
          ctx.prs,
          ctx.epoch,
          Keyword.put(ctx.base_opts, :coarse_search, 0)
        )
      end
    end
  end

  test "pseudoranges/3 accepts an epoch tuple as well as an index" do
    obs = Observations.load!(@obs_path)
    [%{index: index, epoch: epoch} | _] = Observations.epochs(obs)

    {:ok, by_index} = Observations.pseudoranges(obs, index, codes: %{"G" => ["C1C"]})
    {:ok, by_tuple} = Observations.pseudoranges(obs, epoch, codes: %{"G" => ["C1C"]})

    assert by_index == by_tuple
  end

  test "pseudoranges/3 reports an out-of-range epoch index" do
    obs = Observations.load!(@obs_path)
    assert {:error, :epoch_out_of_range} = Observations.pseudoranges(obs, 9_999)
  end

  describe "values/3 (raw multi-code observations)" do
    test "exposes pseudorange, carrier-phase, Doppler, and signal-strength values" do
      obs = Observations.load!(@obs_path)
      {:ok, by_sat} = Observations.values(obs, 0)

      g05 = Map.fetch!(by_sat, "G05")
      codes = Enum.map(g05, & &1.code)
      # The fixture carries L1/L2/L5 code + phase + Doppler + SNR for G05.
      assert "C1C" in codes and "L1C" in codes and "L2W" in codes and "D1C" in codes

      l1c = Enum.find(g05, &(&1.code == "L1C"))
      assert l1c.kind == :carrier_phase
      assert l1c.units == :cycles
      assert is_float(l1c.value) and l1c.value > 0.0

      c1c = Enum.find(g05, &(&1.code == "C1C"))
      assert c1c.kind == :pseudorange and c1c.units == :meters

      assert Enum.find(g05, &(&1.code == "D1C")).kind == :doppler
      assert Enum.find(g05, &(&1.code == "S1C")).kind == :signal_strength
    end

    test "the :codes option scopes the systems and codes that cross the boundary" do
      obs = Observations.load!(@obs_path)

      {:ok, all} = Observations.values(obs, 0)

      # Restrict to GPS L1C/L2W only.
      {:ok, gps_two} = Observations.values(obs, 0, codes: %{"G" => ["L1C", "L2W"]})
      assert Enum.all?(Map.keys(gps_two), &String.starts_with?(&1, "G"))
      assert Enum.sort(Enum.map(Map.fetch!(gps_two, "G05"), & &1.code)) == ["L1C", "L2W"]
      # Non-GPS systems present in the unfiltered result are excluded.
      assert map_size(gps_two) < map_size(all)

      # A system mapped to [] keeps all of that system's codes (GPS-only).
      {:ok, gps_all} = Observations.values(obs, 0, codes: %{"G" => []})
      assert Enum.all?(Map.keys(gps_all), &String.starts_with?(&1, "G"))
      assert length(Map.fetch!(gps_all, "G05")) == length(Map.fetch!(all, "G05"))
    end

    test "index and epoch-tuple selection agree; out-of-range is tagged" do
      obs = Observations.load!(@obs_path)
      [%{index: index, epoch: epoch} | _] = Observations.epochs(obs)

      assert Observations.values(obs, index) == Observations.values(obs, epoch)
      assert {:error, :epoch_out_of_range} = Observations.values(obs, 9_999)
    end
  end

  describe "glonass_slots/1" do
    test "exposes the RINEX GLONASS FDMA channel map" do
      obs = Observations.load!(@obs_path)
      slots = Observations.glonass_slots(obs)

      assert map_size(slots) == 23
      assert slots["R01"] == 1
      assert slots["R10"] == -7
      assert slots["R21"] == 4
    end
  end

  describe "phases/3 (carrier phase with wavelength)" do
    test "returns cycles plus metres for GPS, and the L1/L2 geometry-free offset is small" do
      obs = Observations.load!(@obs_path)
      {:ok, by_sat} = Observations.phases(obs, 0)

      g05 = Map.fetch!(by_sat, "G05")
      assert Enum.all?(g05, &(&1.code |> String.starts_with?("L")))

      l1 = Enum.find(g05, &(&1.code == "L1C"))
      l2 = Enum.find(g05, &(&1.code == "L2W"))

      # GPS L1 wavelength is ~0.190294 m; L2 ~0.244210 m.
      assert_in_delta l1.wavelength_m, 0.190_294, 1.0e-5
      assert_in_delta l2.wavelength_m, 0.244_210, 1.0e-5

      # value_m = cycles * wavelength; both phases track the same geometric range,
      # so L1 - L2 in metres is the small (sub-metre) geometry-free combination,
      # not a full pseudorange-scale difference.
      assert abs(l1.value_m - l2.value_m) < 5.0
    end

    test "returns channel-dependent GLONASS G1/G2 wavelengths from the header slot map" do
      obs = Observations.load!(@obs_path)
      {:ok, by_sat} = Observations.phases(obs, 0, codes: %{"R" => ["L1C", "L2C"]})

      r01 = Map.fetch!(by_sat, "R01")
      l1 = Enum.find(r01, &(&1.code == "L1C"))
      l2 = Enum.find(r01, &(&1.code == "L2C"))

      # R01 has frequency channel +1 in the fixture's GLONASS SLOT / FRQ # map.
      f1 = 1_602_000_000.0 + 562_500.0
      f2 = 1_246_000_000.0 + 437_500.0

      assert_in_delta l1.wavelength_m, 299_792_458.0 / f1, 1.0e-15
      assert_in_delta l2.wavelength_m, 299_792_458.0 / f2, 1.0e-15
      assert_in_delta l1.value_m, l1.value_cycles * l1.wavelength_m, 1.0e-8
      assert_in_delta l2.value_m, l2.value_cycles * l2.wavelength_m, 1.0e-8
    end
  end
end
