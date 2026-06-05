defmodule Orbis.GNSS.QCTest do
  use ExUnit.Case, async: true

  alias Orbis.GNSS.Observables
  alias Orbis.GNSS.Positioning
  alias Orbis.GNSS.QC
  alias Orbis.GNSS.SP3

  @grg Path.join(__DIR__, "fixtures/sp3/GRG0MGXFIN_20201760000_01D_15M_ORB.SP3")

  # Known ground receiver in ITRF/ECEF metres (near IGS ESBC, Copenhagen), as
  # used by the observables/positioning tests.
  @rx {3_512_900.0, 780_500.0, 5_248_700.0}

  # Interior epoch of the SP3 span (2020-06-24 00:00 -> 23:45 GPST).
  @epoch ~N[2020-06-24 12:00:00]

  @c 299_792_458.0

  # Chosen clean receiver clock bias (~30 km of range); GPS-only -> one clock.
  @rx_bias_s 1.0e-4

  # An initial guess a few km off truth with a zero clock seed.
  @initial_guess {3_513_900.0, 779_500.0, 5_249_700.0, 0.0}

  @solve_opts [initial_guess: @initial_guess]

  setup_all do
    sp3 = SP3.load!(@grg)

    visible =
      sp3
      |> Observables.predict_all(@rx, @epoch)
      |> Enum.filter(fn {id, r} -> String.starts_with?(id, "G") and match?({:ok, _}, r) end)
      |> Enum.map(fn {id, {:ok, obs}} -> {id, obs} end)
      |> Enum.filter(fn {_id, obs} -> obs.elevation_deg > 0.0 end)
      |> Enum.sort_by(fn {_id, obs} -> -obs.elevation_deg end)
      |> Enum.take(7)

    # Clean pseudorange set: P = geometric_range + c*(rx_bias - sat_clock).
    clean_obs =
      Enum.map(visible, fn {id, obs} ->
        {id, obs.geometric_range_m + @c * (@rx_bias_s - obs.sat_clock_s)}
      end)

    {:ok, sp3: sp3, visible: visible, clean_obs: clean_obs}
  end

  defp position_error(solution) do
    {x, y, z} = @rx
    p = solution.position
    :math.sqrt(:math.pow(p.x_m - x, 2) + :math.pow(p.y_m - y, 2) + :math.pow(p.z_m - z, 2))
  end

  describe "pseudorange_variance/2 (elevation-dependent weighting)" do
    test "monotonically decreases as elevation rises" do
      vars = Enum.map([5, 15, 30, 60, 90], &QC.pseudorange_variance/1)

      assert vars == Enum.sort(vars, :desc)
      assert Enum.uniq(vars) == vars
    end

    test "matches sigma^2 = a^2 + b^2 / sin^2(el) at sample elevations" do
      a = 0.3
      b = 0.3

      for el <- [30.0, 90.0] do
        expected = a * a + b * b / :math.pow(:math.sin(el * :math.pi() / 180.0), 2)
        assert_in_delta QC.pseudorange_variance(el), expected, 1.0e-12
      end

      # At zenith sin(el) = 1, so variance = a^2 + b^2 = 0.18.
      assert_in_delta QC.pseudorange_variance(90.0), 0.18, 1.0e-9
    end

    test "C/N0 variant returns smaller variance for a higher C/N0" do
      strong = QC.pseudorange_variance(30.0, model: :elevation_cn0, cn0: 50.0)
      weak = QC.pseudorange_variance(30.0, model: :elevation_cn0, cn0: 30.0)

      assert is_float(strong) and is_float(weak)
      assert strong < weak
    end

    test "invalid (non-positive) elevation is a tagged error" do
      assert QC.pseudorange_variance(0.0) == {:error, :invalid_elevation}
      assert QC.pseudorange_variance(-5.0) == {:error, :invalid_elevation}
    end

    test "C/N0 model without a cn0 value is a tagged error" do
      assert QC.pseudorange_variance(30.0, model: :elevation_cn0) ==
               {:error, :missing_cn0}
    end

    test "sigmas/2 and weight_vector/2 are consistent and drop invalid entries" do
      entries = [{"G01", 90.0}, {"G02", 30.0}, {"G03", -1.0}]
      sigmas = QC.sigmas(entries)
      weights = QC.weight_vector(entries)

      refute Map.has_key?(sigmas, "G03")
      refute Map.has_key?(weights, "G03")

      for {sat, sigma} <- sigmas do
        assert_in_delta weights[sat], 1.0 / (sigma * sigma), 1.0e-12
      end
    end
  end

  describe "chi2_inv/2 (chi-square threshold)" do
    test "matches published 99.9th-percentile critical values for dof 1..5" do
      # Standard chi-square distribution critical values at the 0.999 quantile.
      published = %{1 => 10.828, 2 => 13.816, 3 => 16.266, 4 => 18.467, 5 => 20.515}

      for {dof, ref} <- published do
        got = QC.chi2_inv(0.999, dof)
        # Wilson-Hilferty in the far upper tail (p = 0.999) at small dof is a
        # slight, consistently conservative *over*-estimate of the true critical
        # value (~0.24-0.33 here), shrinking as dof grows. A conservative
        # threshold is the safe direction for fault detection (fewer false
        # alarms). Bound the absolute error and confirm the conservative sign.
        assert_in_delta got, ref, 0.35
        assert got >= ref
      end
    end
  end

  describe "raim/2 on a clean SP3-synthesized set" do
    test "clean solve recovers truth and RAIM passes", ctx do
      assert {:ok, sol} = Positioning.solve(ctx.sp3, ctx.clean_obs, @epoch, @solve_opts)

      assert position_error(sol) < 1.0e-2

      result = QC.raim(sol)

      assert result.fault_detected? == false
      assert result.testable? == true
      # GPS-only: n_systems = 1, n_states = 4.
      assert result.dof == length(sol.used_sats) - 4
      assert result.test_statistic < result.threshold
    end
  end

  describe "raim/2 fault detection" do
    setup ctx do
      # Inject the fault on a known-localizing satellite. Detection
      # (fault_detected?, T > threshold) and the FDE exclusion of the
      # largest-normalized-residual satellite are robust for any choice, but the
      # `worst_sat == biased_sat` equality is geometry-dependent: with unit
      # weights least squares can spread a single bias across the residual
      # vector so the largest post-fit residual lands on a neighbour rather than
      # the faulted satellite. Which satellites localize is a non-monotonic
      # function of the full geometry (not simply elevation), so this test fixes
      # on one index empirically confirmed to localize for this fixture/epoch;
      # it is not a claim that every satellite would.
      biased_sat = elem(Enum.at(ctx.clean_obs, 4), 0)

      faulted_obs =
        Enum.map(ctx.clean_obs, fn {sat, pr} ->
          if sat == biased_sat, do: {sat, pr + 200.0}, else: {sat, pr}
        end)

      {:ok, biased_sat: biased_sat, faulted_obs: faulted_obs}
    end

    test "a +200 m bias on one satellite is detected and is the worst sat", ctx do
      assert {:ok, sol} = Positioning.solve(ctx.sp3, ctx.faulted_obs, @epoch, @solve_opts)

      result = QC.raim(sol)

      assert result.fault_detected? == true
      assert result.test_statistic > result.threshold
      assert result.worst_sat == ctx.biased_sat
    end

    test "FDE excludes exactly the biased satellite and recovers the position", ctx do
      # The faulted solve's own error, for the recovery comparison.
      {:ok, faulted_sol} = Positioning.solve(ctx.sp3, ctx.faulted_obs, @epoch, @solve_opts)
      faulted_error = position_error(faulted_sol)

      assert {:ok, fde} = QC.fde(ctx.sp3, ctx.faulted_obs, @epoch, @solve_opts)

      assert fde.excluded == [{ctx.biased_sat, :raim_excluded}]
      assert fde.iterations == 1

      recovered_error = position_error(fde.solution)
      assert recovered_error < 1.0e-2
      assert recovered_error < faulted_error

      # The cleaned solution passes RAIM.
      assert QC.raim(fde.solution).fault_detected? == false
    end
  end

  describe "fde/4 on a clean set" do
    test "excludes nothing and converges immediately", ctx do
      assert {:ok, fde} = QC.fde(ctx.sp3, ctx.clean_obs, @epoch, @solve_opts)

      assert fde.excluded == []
      assert fde.iterations == 0
      assert position_error(fde.solution) < 1.0e-2
    end
  end

  describe "degenerate geometry" do
    test "dof <= 0 -> RAIM reports a non-testable result without raising", ctx do
      # Exactly four GPS sats: n_used == n_states (4), so dof == 0.
      four_obs = Enum.take(ctx.clean_obs, 4)
      assert {:ok, sol} = Positioning.solve(ctx.sp3, four_obs, @epoch, @solve_opts)
      assert length(sol.used_sats) == 4

      result = QC.raim(sol)

      assert result.testable? == false
      assert result.fault_detected? == false
      assert result.dof <= 0
      assert result.threshold == nil
    end

    test "fde/4 with too few satellites returns a tagged error", ctx do
      three_obs = Enum.take(ctx.clean_obs, 3)

      assert {:error, {:too_few_satellites, _used, _required}} =
               QC.fde(ctx.sp3, three_obs, @epoch, @solve_opts)
    end
  end

  describe "raim/2 option validation" do
    test "an out-of-range p_fa raises ArgumentError, not an obscure math error", ctx do
      assert {:ok, sol} = Positioning.solve(ctx.sp3, ctx.clean_obs, @epoch, @solve_opts)

      assert_raise ArgumentError, fn -> QC.raim(sol, p_fa: 0.0) end
      assert_raise ArgumentError, fn -> QC.raim(sol, p_fa: 1.0) end
      assert_raise ArgumentError, fn -> QC.raim(sol, p_fa: -0.1) end
    end

    test "a non-positive custom weight raises ArgumentError", ctx do
      assert {:ok, sol} = Positioning.solve(ctx.sp3, ctx.clean_obs, @epoch, @solve_opts)
      bad_weights = Map.new(sol.used_sats, fn s -> {s, -1.0} end)

      assert_raise ArgumentError, fn -> QC.raim(sol, weights: bad_weights) end
    end
  end
end
