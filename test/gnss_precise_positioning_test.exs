defmodule Orbis.GNSS.PrecisePositioningTest do
  @moduledoc false

  use ExUnit.Case, async: true

  alias Orbis.GNSS.Observables
  alias Orbis.GNSS.PrecisePositioning
  alias Orbis.GNSS.PrecisePositioning.Solution
  alias Orbis.GNSS.SP3

  @sp3_path Path.join(__DIR__, "fixtures/sp3/GRG0MGXFIN_20201760000_01D_15M_ORB.SP3")
  @epoch ~N[2020-06-24 12:00:00]
  @epoch2 ~N[2020-06-24 12:15:00]
  @epoch3 ~N[2020-06-24 12:30:00]
  @epochs [@epoch, @epoch2, @epoch3]
  @truth {3_512_900.0, 780_500.0, 5_248_700.0}
  @clock_m 12.5
  @epoch_clocks_m [12.5, -8.25, 4.0]
  @c 299_792_458.0

  setup_all do
    sp3 = SP3.load!(@sp3_path)

    sats =
      sp3
      |> SP3.satellite_ids()
      |> Enum.filter(&String.starts_with?(&1, "G"))
      |> Enum.flat_map(fn sat ->
        case Observables.predict(sp3, sat, @truth, @epoch) do
          {:ok, obs} when obs.elevation_deg > 10.0 -> [{sat, obs}]
          _ -> []
        end
      end)
      |> Enum.take(8)

    multi_sats =
      sp3
      |> SP3.satellite_ids()
      |> Enum.filter(&String.starts_with?(&1, "G"))
      |> Enum.flat_map(fn sat ->
        with [_ | _] = predictions <-
               Enum.map(@epochs, fn epoch ->
                 case Observables.predict(sp3, sat, @truth, epoch) do
                   {:ok, obs} when obs.elevation_deg > 10.0 -> {epoch, obs}
                   _ -> nil
                 end
               end),
             false <- Enum.any?(predictions, &is_nil/1) do
          [{sat, predictions}]
        else
          _ -> []
        end
      end)
      |> Enum.take(8)

    true = length(sats) >= 6
    true = length(multi_sats) >= 6

    {:ok,
     sp3: sp3,
     sats: sats,
     observations: synth_observations(sats),
     multi_sats: multi_sats,
     epoch_observations: synth_epoch_observations(multi_sats)}
  end

  describe "solve_float/4" do
    test "recovers position, receiver clock, and float ambiguities from exact code+phase", ctx do
      # Known-truth oracle: the observations are synthesized from the SP3 forward
      # model with hidden receiver position, clock, and per-satellite float
      # ambiguities. The estimator must recover those hidden values from a
      # separate initial guess; it is not checked against its own output.
      assert {:ok, %Solution{} = sol} =
               PrecisePositioning.solve_float(ctx.sp3, ctx.observations, @epoch,
                 initial_guess: {3_513_400.0, 780_100.0, 5_249_000.0, -20.0}
               )

      assert position_error(sol.position, @truth) < 1.0e-3
      assert abs(sol.rx_clock_m - @clock_m) < 1.0e-4
      assert abs(sol.rx_clock_s - @clock_m / @c) < 1.0e-13

      for {sat, expected} <- true_ambiguities(ctx.sats) do
        assert abs(sol.ambiguities_m[sat] - expected) < 1.0e-4
      end

      for {_sat, residual} <- sol.residuals_m do
        assert abs(residual.code_m) < 1.0e-4
        assert abs(residual.phase_m) < 1.0e-4
      end

      assert sol.used_sats == Enum.map(ctx.sats, &elem(&1, 0))
      assert sol.metadata.converged
      assert sol.metadata.status == :position_tolerance
      assert sol.metadata.code_rms_m < 1.0e-4
      assert sol.metadata.phase_rms_m < 1.0e-4
    end

    test "can seed itself from the code-only SPP solution", ctx do
      assert {:ok, %Solution{} = sol} =
               PrecisePositioning.solve_float(ctx.sp3, ctx.observations, @epoch,
                 spp_initial_guess: {3_513_400.0, 780_100.0, 5_249_000.0, 0.0}
               )

      assert position_error(sol.position, @truth) < 1.0e-3
      assert abs(sol.rx_clock_m - @clock_m) < 1.0e-4
    end

    test "a phase fault appears in the phase residuals, not the code residuals", ctx do
      [{bad_sat, _obs} | _] = ctx.sats

      faulted =
        Enum.map(ctx.observations, fn
          %{satellite_id: ^bad_sat} = obs -> %{obs | phase_m: obs.phase_m + 1.25}
          obs -> obs
        end)

      assert {:ok, %Solution{} = sol} =
               PrecisePositioning.solve_float(ctx.sp3, faulted, @epoch,
                 initial_guess: {3_513_400.0, 780_100.0, 5_249_000.0, -20.0}
               )

      assert abs(sol.residuals_m[bad_sat].code_m) < 1.0e-3

      # A float ambiguity per satellite absorbs a constant phase offset on a
      # one-epoch solve, so the residual stays small and the ambiguity records
      # the faulted offset.
      assert abs(sol.residuals_m[bad_sat].phase_m) < 1.0e-3

      assert abs(
               sol.ambiguities_m[bad_sat] -
                 (Map.fetch!(true_ambiguities(ctx.sats), bad_sat) + 1.25)
             ) < 1.0e-3
    end
  end

  describe "solve_float_epochs/3" do
    test "recovers static position, per-epoch clocks, and constant ambiguities", ctx do
      assert {:ok, sol} =
               PrecisePositioning.solve_float_epochs(ctx.sp3, ctx.epoch_observations,
                 initial_guess: {3_513_400.0, 780_100.0, 5_249_000.0, -20.0}
               )

      assert position_error(sol.position, @truth) < 1.0e-3
      assert sol.epochs == @epochs
      expected_sats = ctx.multi_sats |> Enum.map(&elem(&1, 0)) |> Enum.sort()
      assert sol.used_sats == expected_sats
      assert sol.metadata.n_epochs == 3
      assert sol.metadata.n_observations == 24
      assert sol.metadata.converged
      assert sol.metadata.status == :state_tolerance

      for {clock, expected} <- Enum.zip(sol.epoch_clocks, @epoch_clocks_m) do
        assert abs(clock.rx_clock_m - expected) < 1.0e-4
        assert abs(clock.rx_clock_s - expected / @c) < 1.0e-13
      end

      for {sat, expected} <- true_ambiguities(ctx.multi_sats) do
        assert abs(sol.ambiguities_m[sat] - expected) < 1.0e-4
      end

      for residual <- sol.residuals_m do
        assert abs(residual.code_m) < 1.0e-4
        assert abs(residual.phase_m) < 1.0e-4
      end
    end

    test "can seed a multi-epoch arc from code-only SPP solutions", ctx do
      assert {:ok, sol} =
               PrecisePositioning.solve_float_epochs(ctx.sp3, ctx.epoch_observations,
                 spp_initial_guess: {3_513_400.0, 780_100.0, 5_249_000.0, 0.0}
               )

      assert position_error(sol.position, @truth) < 1.0e-3
      assert abs(hd(Enum.map(sol.epoch_clocks, & &1.rx_clock_m)) - hd(@epoch_clocks_m)) < 1.0e-4
    end

    test "multi-epoch input errors are tagged", ctx do
      one = [hd(ctx.epoch_observations)]

      assert {:error, {:too_few_epochs, 1, 2}} =
               PrecisePositioning.solve_float_epochs(ctx.sp3, one)

      duplicated = [hd(ctx.epoch_observations), hd(ctx.epoch_observations)]

      assert {:error, {:duplicate_epoch, @epoch}} =
               PrecisePositioning.solve_float_epochs(ctx.sp3, duplicated)

      assert {:error, {:invalid_option, :ambiguity_tolerance_m}} =
               PrecisePositioning.solve_float_epochs(ctx.sp3, ctx.epoch_observations,
                 initial_guess: {3_513_400.0, 780_100.0, 5_249_000.0, -20.0},
                 ambiguity_tolerance_m: -1.0
               )
    end
  end

  describe "solve_float/4 errors" do
    test "empty and too-few observation sets are tagged", ctx do
      assert {:error, :no_observations} = PrecisePositioning.solve_float(ctx.sp3, [], @epoch)

      few = Enum.take(ctx.observations, 3)

      assert {:error, {:too_few_satellites, 3, 4}} =
               PrecisePositioning.solve_float(ctx.sp3, few, @epoch)
    end

    test "duplicate and malformed observations are tagged", ctx do
      [first | rest] = ctx.observations
      first_sat = first.satellite_id

      assert {:error, {:duplicate_observation, ^first_sat}} =
               PrecisePositioning.solve_float(ctx.sp3, [first, first | rest], @epoch)

      assert {:error, {:invalid_observation, {"G01", :bad, 1.0}}} =
               PrecisePositioning.solve_float(ctx.sp3, [{"G01", :bad, 1.0}], @epoch)
    end

    test "bad initial guess and bad sigmas are tagged", ctx do
      assert {:error, :invalid_initial_guess} =
               PrecisePositioning.solve_float(ctx.sp3, ctx.observations, @epoch,
                 initial_guess: {:bad, 0.0, 0.0, 0.0}
               )

      assert {:error, {:invalid_sigma, :code_sigma_m}} =
               PrecisePositioning.solve_float(ctx.sp3, ctx.observations, @epoch,
                 initial_guess: {3_513_400.0, 780_100.0, 5_249_000.0, -20.0},
                 code_sigma_m: 0.0
               )

      assert {:error, {:invalid_option, :max_iterations}} =
               PrecisePositioning.solve_float(ctx.sp3, ctx.observations, @epoch,
                 initial_guess: {3_513_400.0, 780_100.0, 5_249_000.0, -20.0},
                 max_iterations: 0
               )

      assert {:error, {:invalid_option, :position_tolerance_m}} =
               PrecisePositioning.solve_float(ctx.sp3, ctx.observations, @epoch,
                 initial_guess: {3_513_400.0, 780_100.0, 5_249_000.0, -20.0},
                 position_tolerance_m: -1.0
               )
    end
  end

  defp synth_observations(sats) do
    sats
    |> Enum.with_index()
    |> Enum.map(fn {{sat, obs}, idx} ->
      code = obs.geometric_range_m - @c * obs.sat_clock_s + @clock_m
      ambiguity = ambiguity_m(idx)

      %{
        satellite_id: sat,
        code_m: code,
        phase_m: code + ambiguity
      }
    end)
  end

  defp synth_epoch_observations(multi_sats) do
    @epochs
    |> Enum.zip(@epoch_clocks_m)
    |> Enum.map(fn {epoch, clock_m} ->
      observations =
        Enum.map(multi_sats, fn {sat, predictions} ->
          {_epoch, obs} =
            Enum.find(predictions, fn {prediction_epoch, _obs} -> prediction_epoch == epoch end)

          code = obs.geometric_range_m - @c * obs.sat_clock_s + clock_m
          idx = Enum.find_index(multi_sats, fn {candidate, _} -> candidate == sat end)

          %{
            satellite_id: sat,
            code_m: code,
            phase_m: code + ambiguity_m(idx)
          }
        end)

      %{epoch: epoch, observations: observations}
    end)
  end

  defp true_ambiguities(sats) do
    sats
    |> Enum.with_index()
    |> Map.new(fn {{sat, _obs}, idx} -> {sat, ambiguity_m(idx)} end)
  end

  defp ambiguity_m(idx), do: 15_000.0 + idx * 17.25

  defp position_error(%{x_m: x, y_m: y, z_m: z}, {tx, ty, tz}) do
    :math.sqrt((x - tx) * (x - tx) + (y - ty) * (y - ty) + (z - tz) * (z - tz))
  end
end
