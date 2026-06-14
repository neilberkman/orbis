defmodule HuberIrls202606 do
  @moduledoc false

  # GSDC truth metric for the opt-in crate-layer Huber/IRLS reweighting
  # (`:huber`): before/after single-frequency SPP on real degraded GSDC Pixel-5
  # phone observations, bare static-elevation-weighted crate solve vs the Huber-on
  # solve, on IDENTICAL inputs, measured against the vendored demo5/RTKLIB
  # position-domain oracles. See huber-irls-spec.md.
  #
  # The oracle JSONs carry RTKLIB output and GSDC truth per epoch but NO raw
  # observations; the raw phone L1 pseudoranges and broadcast NAV are staged
  # locally. orbis SPP is run per matched epoch as (A) bare and (B) Huber-on.
  #
  # Run from the orbis worktree root:
  #   ORBIS_BUILD=1 mix run test/fixtures/rtk/generators/huber_irls_2026_06.exs
  # Options: --work DIR (default /tmp/gsdc-work), --results FILE, --report FILE,
  #          --stride N.

  alias Orbis.GNSS.Broadcast
  alias Orbis.GNSS.Positioning
  alias Orbis.GNSS.RINEX.Observations

  @default_work "/tmp/gsdc-work"
  @default_results "/tmp/huber-irls-2026-06-results.json"

  @oracle_fixtures [
    "gsdc_2021_08_04_sjc1_pixel5_p222_demo5_rtklib_oracle.json",
    "gsdc_2021_08_24_svl1_pixel5_p222_demo5_rtklib_oracle.json",
    "gsdc_2021_12_15_mtv1_pixel5_p222_demo5_rtklib_oracle.json",
    "gsdc_2021_12_28_mtv1_pixel5_p222_demo5_rtklib_oracle.json"
  ]

  @l1_codes %{
    "G" => ["C1C"],
    "R" => ["C1C"],
    "E" => ["C1C", "C1X"],
    "C" => ["C2I"]
  }

  # Declared up-front gate parameters (see huber-irls-spec.md).
  @min_epochs 100
  @credibility_factor 2.0
  # Huber tuning, frozen in the spec.
  @huber_k 1.345
  # MAD scale floor (m): the realistic phone L1 code sigma, the load-bearing
  # scale source that sets where Huber engages on metre-class noise.
  @huber_sigma_m 5.0
  @huber_max_iter 5

  def main(args) do
    {opts, _argv, invalid} =
      OptionParser.parse(args,
        strict: [work: :string, results: :string, report: :string, stride: :integer],
        aliases: [w: :work, r: :results]
      )

    if invalid != [], do: raise(ArgumentError, "invalid arguments: #{inspect(invalid)}")

    generator_dir = __DIR__
    fixture_dir = Path.expand("..", generator_dir)
    work = Keyword.get(opts, :work, @default_work)
    stride = max(Keyword.get(opts, :stride, 1), 1)
    results_path = Keyword.get(opts, :results, @default_results)

    report_path =
      Keyword.get(
        opts,
        :report,
        Path.join(generator_dir, "huber-irls-measurement-2026-06.md")
      )

    arcs = Enum.map(@oracle_fixtures, &measure_arc(fixture_dir, &1, work, stride))
    pooled = pool(arcs)

    results = %{
      generated_at: DateTime.utc_now() |> DateTime.to_iso8601(),
      work_dir: work,
      stride: stride,
      gate_params: %{
        min_epochs: @min_epochs,
        credibility_factor: @credibility_factor,
        huber_k: @huber_k,
        huber_sigma_m: @huber_sigma_m,
        huber_max_iter: @huber_max_iter
      },
      arcs: arcs,
      pooled: pooled
    }

    File.write!(results_path, Jason.encode!(results, pretty: true))
    File.write!(report_path, report(results))
    IO.puts("wrote #{results_path}")
    IO.puts("wrote #{report_path}")
    print_summary(arcs, pooled)
  end

  defp measure_arc(fixture_dir, fixture, work, stride) do
    oracle_path = Path.join(fixture_dir, fixture)
    oracle = oracle_path |> File.read!() |> Jason.decode!()
    inputs = oracle["inputs"]
    drive = inputs["drive"]

    rover_path = Path.join([work, drive, "supplemental/gnss_rinex.21o"])
    nav_name = inputs["nav"] |> Path.basename() |> String.replace_suffix(".gz", "")
    nav_path = Path.join([work, "cors", nav_name])

    IO.puts("measuring #{oracle["reference"]["label"]}")

    cond do
      not File.exists?(rover_path) ->
        %{
          fixture: fixture,
          label: oracle["reference"]["label"],
          blocked: "missing rover RINEX #{rover_path}"
        }

      not File.exists?(nav_path) ->
        %{
          fixture: fixture,
          label: oracle["reference"]["label"],
          blocked: "missing broadcast NAV #{nav_path}"
        }

      true ->
        rover_obs = Observations.load!(rover_path)
        nav = Broadcast.load!(nav_path)
        oracle_by_key = Map.new(oracle["per_epoch"], &{&1["time"], &1})
        base_arp = base_arp(oracle)

        epochs =
          rover_obs
          |> matched_epochs(oracle_by_key)
          |> Enum.with_index()
          |> Enum.filter(fn {_e, i} -> rem(i, stride) == 0 end)
          |> Enum.map(&elem(&1, 0))

        rows = Enum.flat_map(epochs, &solve_set(rover_obs, nav, &1, base_arp))

        bare = summarize(Enum.map(rows, & &1.bare_3d), Enum.map(rows, & &1.bare_h))
        huber = summarize(Enum.map(rows, & &1.huber_3d), Enum.map(rows, & &1.huber_h))
        demo5 = demo5_summary(oracle["per_epoch"])

        %{
          fixture: fixture,
          label: oracle["reference"]["label"],
          matched_epochs: length(epochs),
          solved_epochs: length(rows),
          identical_off: Enum.all?(rows, & &1.identical_off?),
          # Fraction of solved epochs where the Huber fix differs from the bare
          # fix: direct proof the outer reweighting loop engaged (the crate-side
          # outer_iterations counter is not surfaced through the NIF, which keeps
          # the off-path solution tuple byte-identical).
          huber_changed_fraction: changed_fraction(rows),
          bare: bare,
          huber: huber,
          demo5: demo5,
          verdict: verdict(bare, huber, demo5, length(rows))
        }
    end
  end

  defp base_arp(oracle) do
    b = oracle["truth"]["base_station"]["marker_ecef_m"]
    {b["x"], b["y"], b["z"]}
  end

  defp matched_epochs(rover_obs, oracle_by_key) do
    rover_obs
    |> Observations.epochs()
    |> Enum.flat_map(fn entry ->
      key = entry.epoch |> naive_datetime() |> epoch_key()

      case Map.get(oracle_by_key, key) do
        nil -> []
        oracle_epoch -> [%{entry: entry, oracle_epoch: oracle_epoch}]
      end
    end)
  end

  defp solve_set(rover_obs, nav, %{entry: entry, oracle_epoch: oracle_epoch}, base_arp) do
    {:ok, observations} = Observations.pseudoranges(rover_obs, entry.index, codes: @l1_codes)

    if length(observations) < 5 do
      []
    else
      seed = with_clock(base_arp, 0.0)

      bare =
        Positioning.solve(nav, observations, entry.epoch, initial_guess: seed, troposphere: true)

      # Huber-on: identical inputs/seed/troposphere, opt-in crate IRLS only.
      huber =
        Positioning.solve(nav, observations, entry.epoch,
          initial_guess: seed,
          troposphere: true,
          huber: true,
          huber_k: @huber_k,
          huber_sigma: @huber_sigma_m,
          huber_max_iter: @huber_max_iter
        )

      case {bare, huber} do
        {{:ok, bare_sol}, {:ok, huber_sol}} ->
          truth = oracle_epoch["truth_ecef_m"]
          truth_tuple = {truth["x"], truth["y"], truth["z"]}

          [
            %{
              bare_3d: error_3d(bare_sol, truth_tuple),
              bare_h: error_h(bare_sol, truth_tuple),
              huber_3d: error_3d(huber_sol, truth_tuple),
              huber_h: error_h(huber_sol, truth_tuple),
              # Did the Huber fix move off the bare fix?
              huber_changed?: ecef(huber_sol) != ecef(bare_sol),
              # Additive-off cross-check: a SECOND bare solve with :huber false
              # must be byte-identical to the first bare solve.
              identical_off?: identical_off?(nav, observations, entry.epoch, seed, bare_sol)
            }
          ]

        _ ->
          []
      end
    end
  end

  # Re-solve with an explicit `huber: false` and assert byte-identical ECEF to
  # the bare solve, proving the default path is unchanged on real data.
  defp identical_off?(nav, observations, epoch, seed, bare_sol) do
    case Positioning.solve(nav, observations, epoch,
           initial_guess: seed,
           troposphere: true,
           huber: false
         ) do
      {:ok, off_sol} -> ecef(off_sol) == ecef(bare_sol)
      _ -> false
    end
  end

  defp error_3d(sol, truth), do: norm3(sub3(ecef(sol), truth))

  defp error_h(sol, {tx, ty, tz}) do
    {lat, lon, _h} = ecef_to_geodetic({tx, ty, tz})
    {e, n, _u} = ecef_delta_to_enu(sub3(ecef(sol), {tx, ty, tz}), lat, lon)
    :math.sqrt(e * e + n * n)
  end

  defp ecef(%{position: %{x_m: x, y_m: y, z_m: z}}), do: {x, y, z}

  defp summarize([], []), do: %{n: 0}

  defp summarize(threed, horiz) do
    %{
      n: length(threed),
      median_3d_m: percentile(threed, 50),
      p95_3d_m: percentile(threed, 95),
      median_h_m: percentile(horiz, 50),
      p95_h_m: percentile(horiz, 95)
    }
  end

  defp demo5_summary(per_epoch) do
    threed = Enum.map(per_epoch, & &1["error_3d_m"])
    horiz = Enum.map(per_epoch, & &1["horizontal_error_m"])

    %{
      n: length(threed),
      median_3d_m: percentile(threed, 50),
      p95_3d_m: percentile(threed, 95),
      median_h_m: percentile(horiz, 50),
      p95_h_m: percentile(horiz, 95)
    }
  end

  defp verdict(bare, huber, demo5, n) do
    powered? = n >= @min_epochs
    substrate_ok? = bare.median_3d_m <= @credibility_factor * demo5.median_3d_m

    improved? =
      powered? and huber.median_3d_m <= bare.median_3d_m and huber.p95_3d_m <= bare.p95_3d_m

    %{
      powered?: powered?,
      substrate_ok?: substrate_ok?,
      improved?: improved?,
      median_non_regress?: powered? and huber.median_3d_m <= bare.median_3d_m,
      delta_median_3d_m: bare.median_3d_m - huber.median_3d_m,
      delta_p95_3d_m: bare.p95_3d_m - huber.p95_3d_m,
      classification:
        cond do
          not powered? ->
            "underpowered (#{n} < #{@min_epochs} epochs)"

          improved? ->
            "improved: Huber non-regression on median and p95"

          huber.median_3d_m <= bare.median_3d_m ->
            "median-only: median non-regress, p95 regressed (null on strict bar)"

          true ->
            "null: Huber did not beat the bare elevation-weighted solve"
        end
    }
  end

  defp pool(arcs) do
    powered = Enum.filter(arcs, &(Map.get(&1, :verdict, %{})[:powered?] == true))

    %{
      powered_arc_count: length(powered),
      total_arc_count: length(arcs),
      all_off_identical?:
        powered != [] and Enum.all?(powered, &(Map.get(&1, :identical_off) == true)),
      all_powered_improved?:
        powered != [] and Enum.all?(powered, &(&1.verdict[:improved?] == true)),
      all_powered_median_non_regress?:
        powered != [] and Enum.all?(powered, &(&1.verdict[:median_non_regress?] == true))
    }
  end

  defp print_summary(arcs, pooled) do
    IO.puts("\n=== huber-irls GSDC truth metric ===")

    Enum.each(arcs, fn arc ->
      if Map.has_key?(arc, :blocked) do
        IO.puts("#{arc.label}: BLOCKED #{arc.blocked}")
      else
        v = arc.verdict

        IO.puts(
          "#{arc.label}: n=#{arc.solved_epochs} off-identical=#{arc.identical_off} huber-changed=#{fmt(arc.huber_changed_fraction)}\n" <>
            "  bare   med3D=#{fmt(arc.bare.median_3d_m)} p95=#{fmt(arc.bare.p95_3d_m)} medH=#{fmt(arc.bare.median_h_m)} p95H=#{fmt(arc.bare.p95_h_m)}\n" <>
            "  huber  med3D=#{fmt(arc.huber.median_3d_m)} p95=#{fmt(arc.huber.p95_3d_m)} medH=#{fmt(arc.huber.median_h_m)} p95H=#{fmt(arc.huber.p95_h_m)}\n" <>
            "  delta med3D=#{fmt(v.delta_median_3d_m)} p95=#{fmt(v.delta_p95_3d_m)} | demo5 med=#{fmt(arc.demo5.median_3d_m)} | #{v.classification}"
        )
      end
    end)

    IO.puts(
      "pooled: powered #{pooled.powered_arc_count}/#{pooled.total_arc_count}, off byte-identical? #{pooled.all_off_identical?}, all median non-regress? #{pooled.all_powered_median_non_regress?}, all (median AND p95) non-regress? #{pooled.all_powered_improved?}"
    )
  end

  defp report(results) do
    rows =
      results.arcs
      |> Enum.map_join("\n", fn arc ->
        if Map.has_key?(arc, :blocked) do
          "| #{arc.label} | BLOCKED | #{arc.blocked} | | | | | | |"
        else
          v = arc.verdict

          "| #{arc.label} | #{arc.solved_epochs} | #{fmt(arc.bare.median_3d_m)} | #{fmt(arc.bare.p95_3d_m)} | #{fmt(arc.huber.median_3d_m)} | #{fmt(arc.huber.p95_3d_m)} | #{fmt(v.delta_median_3d_m)} | #{fmt(v.delta_p95_3d_m)} | #{fmt(arc.demo5.median_3d_m)} | #{v.classification} |"
        end
      end)

    """
    # huber-irls GSDC truth metric (#{results.generated_at})

    Before/after single-frequency SPP on real GSDC Pixel-5 phone observations,
    bare static-elevation-weighted crate solve vs the opt-in crate-layer
    Huber/IRLS solve (`:huber`), on IDENTICAL inputs, vs the demo5/RTKLIB
    position-domain oracle. Truth is GSDC ground truth carried in the oracle. See
    huber-irls-spec.md (pre-registered).

    Gate params: min epochs #{results.gate_params.min_epochs}, credibility factor #{results.gate_params.credibility_factor}x demo5 median (absolute floor only), Huber k #{fmt(results.gate_params.huber_k)}, MAD scale floor #{fmt(results.gate_params.huber_sigma_m)} m, max outer #{results.gate_params.huber_max_iter}. Epoch stride #{results.stride}.

    | arc | n | bare med 3D | bare p95 3D | huber med 3D | huber p95 3D | delta med 3D | delta p95 3D | demo5 med 3D | classification |
    |---|---|---|---|---|---|---|---|---|---|
    #{rows}

    All values in metres. Delta is bare minus Huber (positive = Huber better).

    Pooled: powered arcs #{results.pooled.powered_arc_count}/#{results.pooled.total_arc_count}; Huber-off byte-identical to bare on every powered arc? #{results.pooled.all_off_identical?}; all powered arcs median non-regress (Huber <= bare on median)? #{results.pooled.all_powered_median_non_regress?}; all powered arcs median AND p95 non-regress? #{results.pooled.all_powered_improved?}.

    Strict bar (pre-registered): Huber 3D median <= bare AND p95 <= bare on EVERY
    powered arc, no slack. A non-positive delta is a null result, not massaged.
    The default path (`:huber` off) producing byte-identical numbers re-proves
    additive-off on real data.

    Reading: orbis runs an unaided single-frequency SPP per epoch; demo5 is a
    tuned multi-GNSS RTK reference and is the absolute context bar, not the
    comparand for the reweighting delta. The capability claim is the bare-vs-Huber
    delta on identical orbis SPP inputs.
    """
  end

  defp fmt(nil), do: "nil"
  defp fmt(x) when is_float(x), do: :erlang.float_to_binary(x, decimals: 3)
  defp fmt(x), do: to_string(x)

  defp changed_fraction([]), do: nil

  defp changed_fraction(rows) do
    changed = Enum.count(rows, & &1.huber_changed?)
    changed / length(rows)
  end

  defp percentile([], _p), do: nil

  defp percentile(values, p) do
    sorted = Enum.sort(values)
    n = length(sorted)
    rank = p / 100.0 * (n - 1)
    lo = trunc(rank)
    hi = min(lo + 1, n - 1)
    frac = rank - lo
    Enum.at(sorted, lo) * (1.0 - frac) + Enum.at(sorted, hi) * frac
  end

  defp sub3({ax, ay, az}, {bx, by, bz}), do: {ax - bx, ay - by, az - bz}
  defp norm3({x, y, z}), do: :math.sqrt(x * x + y * y + z * z)
  defp with_clock({x, y, z}, c), do: {x, y, z, c}

  @earth_a_m 6_378_137.0
  @earth_f 1.0 / 298.257_223_563

  defp ecef_to_geodetic({x, y, z}) do
    e2 = @earth_f * (2.0 - @earth_f)
    p = :math.sqrt(x * x + y * y)
    lon = :math.atan2(y, x)
    lat = :math.atan2(z, p * (1.0 - e2))
    lat = ecef_lat_iter(z, p, lat, e2, 0)
    {lat, lon, 0.0}
  end

  defp ecef_lat_iter(_z, _p, lat, _e2, 5), do: lat

  defp ecef_lat_iter(z, p, lat, e2, iter) do
    sin_lat = :math.sin(lat)
    n = @earth_a_m / :math.sqrt(1.0 - e2 * sin_lat * sin_lat)
    h = p / :math.cos(lat) - n
    new_lat = :math.atan2(z, p * (1.0 - e2 * n / (n + h)))
    ecef_lat_iter(z, p, new_lat, e2, iter + 1)
  end

  defp ecef_delta_to_enu({dx, dy, dz}, lat, lon) do
    sin_lat = :math.sin(lat)
    cos_lat = :math.cos(lat)
    sin_lon = :math.sin(lon)
    cos_lon = :math.cos(lon)
    e = -sin_lon * dx + cos_lon * dy
    n = -sin_lat * cos_lon * dx - sin_lat * sin_lon * dy + cos_lat * dz
    u = cos_lat * cos_lon * dx + cos_lat * sin_lon * dy + sin_lat * dz
    {e, n, u}
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

HuberIrls202606.main(System.argv())
