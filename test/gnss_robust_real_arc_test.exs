defmodule Orbis.GNSS.RobustRealArcTest do
  @moduledoc """
  Real-data guard for the opt-in `:robust` FDE solve path (spp-robustness defect
  3): on a real degraded GSDC Pixel-5 arc, `:robust` must never silently degrade
  a fix.

    * `:robust` with no noise model and no escape hatch REFUSES, before any
      solve, so it can never return a worse fix by accident.
    * `:robust` with a realistic noise model is a no-op-or-better: its 3D median
      error over the arc is <= the bare solve's 3D median (never worse).

  Raw phone observations are staged outside the repo (see spp-robustness-spec.md);
  the whole module is skipped when the staged arc is not present, so the
  in-repo battery is unaffected. The fast in-suite contract refusal is covered
  separately in `test/gnss_qc_test.exs`.
  """
  use ExUnit.Case, async: false

  alias Orbis.GNSS.Broadcast
  alias Orbis.GNSS.Positioning
  alias Orbis.GNSS.RINEX.Observations

  @work "/tmp/gsdc-work"
  @oracle Path.join(
            __DIR__,
            "fixtures/rtk/gsdc_2021_12_15_mtv1_pixel5_p222_demo5_rtklib_oracle.json"
          )
  @rover Path.join([@work, "train/2021-12-15-US-MTV-1/GooglePixel5/supplemental/gnss_rinex.21o"])
  @nav Path.join([@work, "cors/BRDC00WRD_R_20213490000_01D_MN.rnx"])

  @l1_codes %{"G" => ["C1C"], "R" => ["C1C"], "E" => ["C1C", "C1X"], "C" => ["C2I"]}
  @min_epochs 100
  @code_sigma_m 5.0
  @max_pdop 1000.0

  @staged? File.exists?(@rover) and File.exists?(@nav)

  # Large/licensed local corpus, excluded by default (see test_helper.exs); run
  # with: ORBIS_BUILD=1 mix test --include local_data test/gnss_robust_real_arc_test.exs
  @moduletag :local_data

  describe "robust no-silent-degrade on a real GSDC Pixel-5 arc" do
    if @staged? do
      setup do
        oracle = @oracle |> File.read!() |> Jason.decode!()
        rover = Observations.load!(@rover)
        nav = Broadcast.load!(@nav)
        b = oracle["truth"]["base_station"]["marker_ecef_m"]
        seed = {b["x"], b["y"], b["z"], 0.0}
        oracle_by_key = Map.new(oracle["per_epoch"], &{&1["time"], &1})
        {:ok, rover: rover, nav: nav, seed: seed, oracle_by_key: oracle_by_key}
      end

      test "robust with no noise model refuses (never a worse fix)", ctx do
        entry = first_usable(ctx)
        {:ok, obs} = Observations.pseudoranges(ctx.rover, entry.index, codes: @l1_codes)

        assert Positioning.solve(ctx.nav, obs, entry.epoch,
                 initial_guess: ctx.seed,
                 troposphere: true,
                 robust: true
               ) == {:error, {:robust_requires_noise_model, :no_weights}}
      end

      test "robust with a realistic noise model is no-op-or-better vs bare (>=100 epochs)", ctx do
        rows =
          ctx.rover
          |> Observations.epochs()
          |> Enum.flat_map(fn entry ->
            key = entry.epoch |> naive_datetime() |> epoch_key()

            case Map.get(ctx.oracle_by_key, key) do
              nil ->
                []

              oracle_epoch ->
                {:ok, obs} = Observations.pseudoranges(ctx.rover, entry.index, codes: @l1_codes)

                if length(obs) < 5 do
                  []
                else
                  weights =
                    Map.new(obs, fn {sat, _pr} -> {sat, 1.0 / (@code_sigma_m * @code_sigma_m)} end)

                  bare =
                    Positioning.solve(ctx.nav, obs, entry.epoch,
                      initial_guess: ctx.seed,
                      troposphere: true
                    )

                  robust =
                    Positioning.solve(ctx.nav, obs, entry.epoch,
                      initial_guess: ctx.seed,
                      troposphere: true,
                      robust: true,
                      weights: weights,
                      max_pdop: @max_pdop
                    )

                  truth = oracle_epoch["truth_ecef_m"]
                  t = {truth["x"], truth["y"], truth["z"]}

                  case {bare, robust} do
                    {{:ok, b}, {:ok, r}} -> [%{bare: err3d(b, t), robust: err3d(r, t)}]
                    _ -> []
                  end
                end
            end
          end)

        assert length(rows) >= @min_epochs,
               "underpowered: #{length(rows)} < #{@min_epochs} matched epochs"

        bare_med = median(Enum.map(rows, & &1.bare))
        robust_med = median(Enum.map(rows, & &1.robust))

        # The contract: robust-with-a-realistic-model never makes the arc worse.
        assert robust_med <= bare_med,
               "robust median #{robust_med} m exceeded bare median #{bare_med} m"
      end
    else
      @tag :skip
      test "staged GSDC arc not present at #{@rover}; skipped" do
        assert true
      end
    end
  end

  defp first_usable(ctx) do
    ctx.rover
    |> Observations.epochs()
    |> Enum.find(fn entry ->
      {:ok, obs} = Observations.pseudoranges(ctx.rover, entry.index, codes: @l1_codes)
      length(obs) >= 5
    end)
  end

  defp err3d(%{position: %{x_m: x, y_m: y, z_m: z}}, {tx, ty, tz}) do
    :math.sqrt((x - tx) * (x - tx) + (y - ty) * (y - ty) + (z - tz) * (z - tz))
  end

  defp median([]), do: nil

  defp median(values) do
    sorted = Enum.sort(values)
    n = length(sorted)
    mid = div(n, 2)

    if rem(n, 2) == 1,
      do: Enum.at(sorted, mid),
      else: (Enum.at(sorted, mid - 1) + Enum.at(sorted, mid)) / 2.0
  end

  defp naive_datetime({{year, month, day}, {hour, minute, second}}) do
    whole = trunc(second)
    micro = round((second - whole) * 1_000_000)
    NaiveDateTime.new!(Date.new!(year, month, day), Time.new!(hour, minute, whole, {micro, 6}))
  end

  defp epoch_key(%NaiveDateTime{} = ndt) do
    {micro, _p} = ndt.microsecond
    ms = round(micro / 1_000)
    {base, ms} = if ms == 1_000, do: {NaiveDateTime.add(ndt, 1, :second), 0}, else: {ndt, ms}
    base = %{base | microsecond: {0, 0}}
    NaiveDateTime.to_iso8601(base) <> "." <> String.pad_leading(Integer.to_string(ms), 3, "0")
  end
end
