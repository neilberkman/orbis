defmodule Orbis.Passes do
  @moduledoc """
  Satellite pass prediction over a ground station.

  Finds time windows when a satellite is above the local horizon (or a
  minimum elevation threshold) as seen from a ground station. The algorithm:

  1. Propagate at coarse time steps, computing the elevation angle at each.
  2. When elevation crosses 0 deg (rise), bisect to refine the exact rise time.
  3. Track until elevation crosses 0 deg again (set), bisect to refine.
  4. Find the maximum elevation during the pass via golden-section search.
  5. Optionally discard passes whose max elevation is below `:min_elevation`.

  Elevation is computed through the full topocentric pipeline
  (TEME→GCRS→topocentric) via `Orbis.Coordinates.to_topocentric/3`.
  """

  alias Orbis.Elements

  @type ground_station :: %{
          latitude: float(),
          longitude: float(),
          altitude_m: float()
        }

  @type pass :: %{
          rise: DateTime.t(),
          set: DateTime.t(),
          max_elevation: float(),
          max_elevation_time: DateTime.t(),
          duration_seconds: float()
        }

  @doc """
  Predict visible passes of a satellite over a ground station.

  ## Parameters

    * `tle` — parsed `%Orbis.Elements{}` struct
    * `ground_station` — `%{latitude: deg, longitude: deg, altitude_m: m}`
    * `start_time` — `DateTime.t()` beginning of the search window
    * `end_time` — `DateTime.t()` end of the search window
    * `opts` — keyword list
      * `:min_elevation` — minimum peak elevation in degrees to keep a pass
        (default `0.0`)
      * `:step_seconds` — coarse propagation step in seconds (default `60`)

  ## Returns

  A list of pass maps sorted by rise time:

      [%{
        rise: DateTime.t(),
        set: DateTime.t(),
        max_elevation: float(),       # degrees
        max_elevation_time: DateTime.t(),
        duration_seconds: float()
      }]
  """
  @spec predict(Elements.t(), ground_station(), DateTime.t(), DateTime.t(), keyword()) :: [pass()]
  def predict(
        %Elements{} = tle,
        ground_station,
        %DateTime{} = start_time,
        %DateTime{} = end_time,
        opts \\ []
      ) do
    min_elevation = Keyword.get(opts, :min_elevation, 0.0)
    step = Keyword.get(opts, :step_seconds, 60)

    if step <= 0, do: raise(ArgumentError, "step_seconds must be positive, got #{step}")

    start_time
    |> coarse_scan(end_time, step, tle, ground_station)
    |> extract_passes(tle, ground_station)
    |> Enum.filter(fn pass -> pass.max_elevation >= min_elevation end)
  end

  # ---------------------------------------------------------------------------
  # Coarse scan — propagate at fixed intervals and record (time, elevation)
  # ---------------------------------------------------------------------------

  defp coarse_scan(start_time, end_time, step, tle, gs) do
    total_seconds = DateTime.diff(end_time, start_time, :second)
    num_steps = max(div(total_seconds, step), 0)

    Enum.map(0..num_steps, fn i ->
      dt = DateTime.add(start_time, i * step, :second)
      {dt, elevation_at(tle, dt, gs)}
    end)
  end

  # ---------------------------------------------------------------------------
  # Extract passes from the coarse elevation samples
  # ---------------------------------------------------------------------------

  defp extract_passes(samples, tle, gs) do
    # Determine initial state: if the first sample is above horizon,
    # we are already mid-pass — use the first sample time as rise.
    initial_rise =
      case samples do
        [{dt, el} | _] when el >= 0 -> dt
        _ -> nil
      end

    samples
    |> Enum.chunk_every(2, 1, :discard)
    |> find_crossings(tle, gs, initial_rise, [])
  end

  # Walk consecutive pairs looking for positive-going and negative-going zero crossings.
  defp find_crossings([], _tle, _gs, _rise, passes), do: Enum.reverse(passes)

  defp find_crossings([pair | rest], tle, gs, rise_time, passes) do
    [{dt_a, el_a}, {dt_b, el_b}] = pair

    cond do
      # Rising edge: elevation crosses from negative to non-negative
      is_nil(rise_time) and el_a < 0 and el_b >= 0 ->
        refined_rise = bisect_crossing(tle, gs, dt_a, dt_b)
        find_crossings(rest, tle, gs, refined_rise, passes)

      # Setting edge: elevation crosses from non-negative to negative
      not is_nil(rise_time) and el_a >= 0 and el_b < 0 ->
        refined_set = bisect_crossing(tle, gs, dt_a, dt_b)
        pass = build_pass(tle, gs, rise_time, refined_set)
        find_crossings(rest, tle, gs, nil, [pass | passes])

      true ->
        find_crossings(rest, tle, gs, rise_time, passes)
    end
  end

  # ---------------------------------------------------------------------------
  # Bisection — refine a zero crossing to ~0.1 s precision
  # ---------------------------------------------------------------------------

  @bisect_iterations 20

  defp bisect_crossing(tle, gs, dt_low, dt_high) do
    el_low = elevation_at(tle, dt_low, gs)

    Enum.reduce(1..@bisect_iterations, {dt_low, dt_high, el_low}, fn _i, {lo, hi, el_lo} ->
      mid = midpoint_datetime(lo, hi)
      el_mid = elevation_at(tle, mid, gs)

      if same_sign?(el_lo, el_mid) do
        {mid, hi, el_mid}
      else
        {lo, mid, el_lo}
      end
    end)
    |> then(fn {lo, hi, _} -> midpoint_datetime(lo, hi) end)
  end

  defp same_sign?(a, b), do: (a >= 0 and b >= 0) or (a < 0 and b < 0)

  defp midpoint_datetime(a, b) do
    diff_us = DateTime.diff(b, a, :microsecond)
    DateTime.add(a, div(diff_us, 2), :microsecond)
  end

  # ---------------------------------------------------------------------------
  # Build a pass struct — find max elevation via golden-section search
  # ---------------------------------------------------------------------------

  defp build_pass(tle, gs, rise, set) do
    {max_el, max_el_time} = find_max_elevation(tle, gs, rise, set)

    %Orbis.Pass{
      rise: rise,
      set: set,
      max_elevation: max_el,
      max_elevation_time: max_el_time
    }
  end

  # Golden-section search for the peak elevation within [rise, set].
  @golden_ratio (1 + :math.sqrt(5)) / 2
  @golden_resphi 2 - @golden_ratio
  @golden_iterations 30

  defp find_max_elevation(tle, gs, rise, set) do
    total_us = DateTime.diff(set, rise, :microsecond)

    {a, b} =
      Enum.reduce(1..@golden_iterations, {0, total_us}, fn _i, {a, b} ->
        x1 = round(a + @golden_resphi * (b - a))
        x2 = round(b - @golden_resphi * (b - a))

        dt1 = DateTime.add(rise, x1, :microsecond)
        dt2 = DateTime.add(rise, x2, :microsecond)

        el1 = elevation_at(tle, dt1, gs)
        el2 = elevation_at(tle, dt2, gs)

        if el1 > el2, do: {a, x2}, else: {x1, b}
      end)

    best_us = div(a + b, 2)
    best_dt = DateTime.add(rise, best_us, :microsecond)
    best_el = elevation_at(tle, best_dt, gs)

    {best_el, best_dt}
  end

  # Elevation via the full topocentric pipeline (TEME→GCRS→topocentric)
  defp elevation_at(tle, datetime, station) do
    case Orbis.SGP4.propagate(tle, datetime) do
      {:ok, teme_state} ->
        gcrs = Orbis.Coordinates.teme_to_gcrs(teme_state, datetime)
        look = Orbis.Coordinates.to_topocentric(gcrs, datetime, station)
        look.elevation

      {:error, _} ->
        -90.0
    end
  end
end
