defmodule Orbis.GNSS.CarrierPhaseTest do
  use ExUnit.Case, async: true

  alias Orbis.GNSS.CarrierPhase
  alias Orbis.GNSS.Observables
  alias Orbis.GNSS.RINEX.Observations
  alias Orbis.GNSS.SP3

  @grg Path.join(__DIR__, "fixtures/sp3/GRG0MGXFIN_20201760000_01D_15M_ORB.SP3")
  @rx {3_512_900.0, 780_500.0, 5_248_700.0}
  @epoch ~N[2020-06-24 12:00:00]

  @c 299_792_458.0
  @f_l1 1_575_420_000.0
  @f_l2 1_227_600_000.0

  # Constant per-band integer ambiguities for the synthetic phase.
  @amb1 123_456
  @amb2 234_567

  # Number of epochs and cadence for the synthetic arc.
  @n_epochs 25
  @cadence_s 30
  # Slip injected from this (0-based) index onward.
  @slip_index 12

  setup do
    sp3 = SP3.load!(@grg)

    # Highest GPS satellite, chosen empirically (same approach as the
    # observables test) so it stays visible across the short arc.
    results = Observables.predict_all(sp3, @rx, @epoch)

    {high_id, _} =
      results
      |> Enum.filter(fn {id, r} -> String.starts_with?(id, "G") and match?({:ok, _}, r) end)
      |> Enum.map(fn {id, {:ok, obs}} -> {id, obs} end)
      |> Enum.max_by(fn {_id, obs} -> obs.elevation_deg end)

    {:ok, sp3: sp3, sat: high_id}
  end

  # Build a clean dual-frequency arc for one satellite from the SP3 geometry.
  #
  # phi_i(k) = (range(k) + amb_i*lambda_i + iono_i(k)) / lambda_i   [cycles]
  # p_i(k)   = range(k) - iono_i(k) + code_noise                    [metres]
  #
  # Ionospheric delay on phase is an advance (opposite sign to the code group
  # delay); we use a slow linear ramp on band 1 and scale band 2 by (f1/f2)^2
  # (the 1/f^2 dispersion). Code carries seeded Gaussian noise; phase is clean.
  defp clean_arc(sp3, sat, opts \\ []) do
    rx = Keyword.get(opts, :rx, @rx)
    code_sigma = Keyword.get(opts, :code_sigma, 0.3)
    lambda1 = @c / @f_l1
    lambda2 = @c / @f_l2
    iono2_scale = :math.pow(@f_l1 / @f_l2, 2)

    for k <- 0..(@n_epochs - 1) do
      epoch = NaiveDateTime.add(@epoch, k * @cadence_s, :second)
      {:ok, obs} = Observables.predict(sp3, sat, rx, epoch)
      range = obs.geometric_range_m

      # Slow ionospheric ramp, a few cm over the arc.
      iono1 = 0.02 + 0.001 * k
      iono2 = iono1 * iono2_scale

      phi1 = (range + @amb1 * lambda1 + iono1) / lambda1
      phi2 = (range + @amb2 * lambda2 + iono2) / lambda2

      noise1 = code_sigma * :rand.normal()
      noise2 = code_sigma * :rand.normal()

      %{
        epoch: epoch,
        true_range: range,
        phi1: phi1,
        phi2: phi2,
        p1: range - iono1 + noise1,
        p2: range - iono2 + noise2,
        lli1: 0,
        lli2: 0,
        f1: @f_l1,
        f2: @f_l2
      }
    end
  end

  defp seed, do: :rand.seed(:exsss, {101, 202, 303})

  describe "wide_lane_wavelength/2" do
    test "GPS L1/L2 wide-lane wavelength is about 0.8619 m" do
      assert {:ok, lambda_wl} = CarrierPhase.wide_lane_wavelength(@f_l1, @f_l2)
      assert_in_delta lambda_wl, 0.8619, 1.0e-3
    end

    test "equal frequencies returns a tagged error" do
      assert CarrierPhase.wide_lane_wavelength(@f_l1, @f_l1) == {:error, :equal_frequencies}
    end
  end

  describe "scalar combinations" do
    test "geometry_free is the difference of the two phases in metres" do
      assert CarrierPhase.geometry_free(100.0, 60.0) == 40.0
    end

    test "narrow_lane_code is the frequency-weighted code mean" do
      assert {:ok, p_nl} = CarrierPhase.narrow_lane_code(100.0, 200.0, @f_l1, @f_l2)
      expected = (@f_l1 * 100.0 + @f_l2 * 200.0) / (@f_l1 + @f_l2)
      assert_in_delta p_nl, expected, 1.0e-9
    end

    test "melbourne_wubbena equals wide-lane phase minus narrow-lane code" do
      {:ok, lambda_wl} = CarrierPhase.wide_lane_wavelength(@f_l1, @f_l2)
      {:ok, p_nl} = CarrierPhase.narrow_lane_code(10.0, 12.0, @f_l1, @f_l2)
      expected = lambda_wl * (5.0 - 3.0) - p_nl
      assert {:ok, mw} = CarrierPhase.melbourne_wubbena(5.0, 3.0, 10.0, 12.0, @f_l1, @f_l2)
      assert_in_delta mw, expected, 1.0e-9
    end

    test "melbourne_wubbena with equal frequencies is a tagged error" do
      assert CarrierPhase.melbourne_wubbena(1.0, 2.0, 3.0, 4.0, @f_l1, @f_l1) ==
               {:error, :equal_frequencies}
    end
  end

  describe "detect_cycle_slips/2 geometry-free (STRONG)" do
    test "a clean arc flags no geometry-free slip", %{sp3: sp3, sat: sat} do
      seed()
      arc = clean_arc(sp3, sat)
      results = CarrierPhase.detect_cycle_slips(arc)

      gf_flagged = Enum.filter(results, &(:geometry_free in &1.reasons))

      assert gf_flagged == [],
             "clean arc flagged GF at #{inspect(Enum.map(gf_flagged, & &1.epoch))}"

      # Clean-arc geometry-free variation across the arc (documented bound).
      gfs = results |> Enum.map(& &1.gf) |> Enum.reject(&is_nil/1)
      gf_var = Enum.max(gfs) - Enum.min(gfs)
      assert gf_var < 0.05, "clean GF variation #{gf_var} m"
    end

    test "an integer cycle slip on phi1 is flagged at exactly that epoch via geometry_free",
         %{sp3: sp3, sat: sat} do
      seed()
      arc = clean_arc(sp3, sat)
      slip_cycles = 5

      slipped =
        arc
        |> Enum.with_index()
        |> Enum.map(fn {ep, idx} ->
          if idx >= @slip_index, do: %{ep | phi1: ep.phi1 + slip_cycles}, else: ep
        end)

      results = CarrierPhase.detect_cycle_slips(slipped)

      flagged_idx =
        results
        |> Enum.with_index()
        |> Enum.filter(fn {r, _i} -> :geometry_free in r.reasons end)
        |> Enum.map(fn {_r, i} -> i end)

      assert flagged_idx == [@slip_index],
             "expected GF slip only at #{@slip_index}, got #{inspect(flagged_idx)}"
    end
  end

  describe "detect_cycle_slips/2 Melbourne-Wubbena (STRONG)" do
    test "clean arc: MW variation is within a few wide-lane cycles", %{sp3: sp3, sat: sat} do
      seed()
      arc = clean_arc(sp3, sat)
      results = CarrierPhase.detect_cycle_slips(arc)

      mws = results |> Enum.map(& &1.mw) |> Enum.reject(&is_nil/1)
      mw_var = Enum.max(mws) - Enum.min(mws)
      {:ok, lambda_wl} = CarrierPhase.wide_lane_wavelength(@f_l1, @f_l2)
      mw_var_cycles = mw_var / lambda_wl

      # Code-noise-driven scatter; bound well under the slip gate.
      assert mw_var < 1.5, "clean MW variation #{mw_var} m (#{mw_var_cycles} WL cycles)"
      assert Enum.all?(results, &(:melbourne_wubbena not in &1.reasons))
    end

    test "a wide-lane slip shifts MW by ~n*lambda_WL and is flagged via melbourne_wubbena",
         %{sp3: sp3, sat: sat} do
      seed()
      arc = clean_arc(sp3, sat)
      n_wl = 8
      {:ok, lambda_wl} = CarrierPhase.wide_lane_wavelength(@f_l1, @f_l2)

      # phi_WL = phi1 - phi2; bumping phi1 by n_wl cycles shifts the wide-lane
      # phase (hence MW) by n_wl * lambda_WL metres.
      slipped =
        arc
        |> Enum.with_index()
        |> Enum.map(fn {ep, idx} ->
          if idx >= @slip_index, do: %{ep | phi1: ep.phi1 + n_wl}, else: ep
        end)

      results = CarrierPhase.detect_cycle_slips(slipped)

      mw_before = Enum.at(results, @slip_index - 1).mw
      mw_at = Enum.at(results, @slip_index).mw
      step = mw_at - mw_before

      assert_in_delta step, n_wl * lambda_wl, 1.0, "MW step #{step} vs #{n_wl * lambda_wl}"
      assert :melbourne_wubbena in Enum.at(results, @slip_index).reasons
    end

    test "MW is geometry-free: a different receiver geometry leaves MW unchanged up to noise",
         %{sp3: sp3, sat: sat} do
      seed()
      arc_a = clean_arc(sp3, sat)
      seed()
      arc_b = clean_arc(sp3, sat, rx: {3_400_000.0, 900_000.0, 5_300_000.0})

      res_a = CarrierPhase.detect_cycle_slips(arc_a)
      res_b = CarrierPhase.detect_cycle_slips(arc_b)

      # Same ambiguities + identical seeded code noise => MW (geometry-free)
      # series agree closely despite very different ranges.
      pairs =
        Enum.zip(res_a, res_b)
        |> Enum.reject(fn {a, b} -> is_nil(a.mw) or is_nil(b.mw) end)

      max_diff = pairs |> Enum.map(fn {a, b} -> abs(a.mw - b.mw) end) |> Enum.max()
      assert max_diff < 1.0e-6, "MW differed by #{max_diff} m across geometries"
    end
  end

  describe "detect_cycle_slips/2 LLI" do
    test "an LLI loss-of-lock bit flags :lli even with continuous GF/MW",
         %{sp3: sp3, sat: sat} do
      seed()
      arc = clean_arc(sp3, sat)

      lli_idx = 10

      tagged =
        arc
        |> Enum.with_index()
        |> Enum.map(fn {ep, idx} ->
          if idx == lli_idx, do: %{ep | lli1: 1}, else: ep
        end)

      results = CarrierPhase.detect_cycle_slips(tagged)
      target = Enum.at(results, lli_idx)

      assert :lli in target.reasons
      refute :geometry_free in target.reasons
      refute :melbourne_wubbena in target.reasons
      assert target.slip

      # No other epoch flags :lli.
      others = results |> List.delete_at(lli_idx)
      assert Enum.all?(others, &(:lli not in &1.reasons))
    end
  end

  describe "detect_cycle_slips/2 degenerate / GLONASS" do
    test "an empty arc returns empty" do
      assert CarrierPhase.detect_cycle_slips([]) == []
    end

    test "a single-epoch arc flags no GF/MW slip" do
      arc = [%{epoch: 0, phi1: 1.0, phi2: 1.0, p1: 1.0, p2: 1.0, f1: @f_l1, f2: @f_l2}]
      [r] = CarrierPhase.detect_cycle_slips(arc)
      refute r.slip
      assert r.reasons == []
    end

    test "a GLONASS epoch (nil frequency) is skipped, never raised" do
      arc = [
        %{epoch: 0, phi1: 1.0, phi2: 1.0, p1: 1.0, p2: 1.0, f1: nil, f2: nil, lli1: 0, lli2: 0},
        %{epoch: 1, phi1: 2.0, phi2: 2.0, p1: 2.0, p2: 2.0, f1: nil, f2: nil, lli1: 0, lli2: 0}
      ]

      results = CarrierPhase.detect_cycle_slips(arc)
      assert Enum.all?(results, & &1.skipped)
      assert Enum.all?(results, &(&1.slip == false))
      assert Enum.all?(results, &(&1.gf == nil and &1.mw == nil))
    end

    test "a real GLONASS satellite from the SP3 product is skipped, never raised",
         %{sp3: sp3} do
      glonass =
        sp3 |> SP3.satellite_ids() |> Enum.find(&String.starts_with?(&1, "R"))

      assert glonass != nil, "expected a GLONASS satellite in the GRG product"

      # band_frequency_hz returns nil for GLONASS, so the arc is all skipped.
      f1 = Observations.band_frequency_hz(String.first(glonass), "1")
      assert f1 == nil

      arc = [%{epoch: @epoch, phi1: 1.0, phi2: 1.0, p1: 1.0, p2: 1.0, f1: f1, f2: f1}]
      assert [%{skipped: true, slip: false}] = CarrierPhase.detect_cycle_slips(arc)
    end
  end

  describe "smooth_code/2 (Hatch)" do
    test "smoothed code scatter about the true range is much smaller than raw",
         %{sp3: sp3, sat: sat} do
      seed()
      arc = clean_arc(sp3, sat)
      smoothed = CarrierPhase.smooth_code(arc)

      # Compare scatter only after the window has grown (skip the first few).
      warmup = 5

      raw_resid =
        arc
        |> Enum.drop(warmup)
        |> Enum.map(fn ep -> ep.p1 - ep.true_range end)

      smooth_resid =
        Enum.zip(arc, smoothed)
        |> Enum.drop(warmup)
        |> Enum.map(fn {ep, s} -> s.p_smooth - ep.true_range end)

      raw_std = std(raw_resid)
      smooth_std = std(smooth_resid)

      assert smooth_std < 0.4 * raw_std,
             "smoothed std #{smooth_std} not << raw std #{raw_std}"
    end

    test "a cycle slip / LLI resets the filter (no smoothing across the slip)",
         %{sp3: sp3, sat: sat} do
      seed()
      arc = clean_arc(sp3, sat)
      reset_idx = 14

      tagged =
        arc
        |> Enum.with_index()
        |> Enum.map(fn {ep, idx} ->
          if idx == reset_idx, do: %{ep | lli1: 1}, else: ep
        end)

      smoothed = CarrierPhase.smooth_code(tagged)
      at = Enum.at(smoothed, reset_idx)
      tagged_ep = Enum.at(tagged, reset_idx)

      assert at.reset, "expected a reset at #{reset_idx}"
      assert at.window == 1, "expected window reset to 1, got #{at.window}"
      assert_in_delta at.p_smooth, tagged_ep.p1, 1.0e-9
    end

    test "an empty arc returns empty" do
      assert CarrierPhase.smooth_code([]) == []
    end

    test "an epoch missing band-1 code yields a nil smoothed value" do
      arc = [
        %{epoch: 0, phi1: 1.0, p1: 100.0, lli1: 0, f1: @f_l1, f2: @f_l2, phi2: 1.0, p2: 100.0},
        %{epoch: 1, phi1: 2.0, p1: nil, lli1: 0, f1: @f_l1, f2: @f_l2, phi2: 2.0, p2: nil}
      ]

      [_a, b] = CarrierPhase.smooth_code(arc)
      assert b.p_smooth == nil
      assert b.window == 0
    end
  end

  defp std([]), do: 0.0

  defp std(xs) do
    n = length(xs)
    mean = Enum.sum(xs) / n
    var = Enum.reduce(xs, 0.0, fn x, acc -> acc + (x - mean) * (x - mean) end) / n
    :math.sqrt(var)
  end
end
