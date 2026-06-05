defmodule Orbis.ReducedOrbit.Piecewise do
  @moduledoc """
  A long position track approximated by a sequence of contiguous, independently
  fitted `Orbis.ReducedOrbit` segments.

  A single `Orbis.ReducedOrbit` distills a whole track into one set of mean
  elements; extrapolated over a day it drifts (GPS ~thousands of km with the
  circular model, ~8 km eccentric). A `Piecewise` model instead splits the span
  `[t0, t1]` into contiguous segments of `:segment_s` seconds and fits each one
  with the **existing** `Orbis.ReducedOrbit.fit/2`. Every query then lands
  *inside* a fit window, so the error collapses to the in-window residual
  (sub-km to a few km) rather than the extrapolation error — at the cost of
  storing N small models. It is pure orchestration over the single-segment
  primitives; the orbit math, frames, and time scales are unchanged.

  This is the best accuracy-per-byte option for caching and transport. It is
  **not** orbit determination and **not** a substitute for SP3 or SGP4 — it is a
  compact approximation, honest about its residual through a source-backed
  `drift/3`.

  ## Models

  Each segment is one of the two `Orbis.ReducedOrbit` models, selected with the
  same `:model` option:

    * `:circular_secular` (the **default**) — a circular orbit;
    * `:eccentric_secular` — recovers the radial `a·e` signal.

  See `Orbis.ReducedOrbit` for the per-model details.

  ## Size / accuracy tradeoff

  Storage grows ~linearly with the number of segments: a span of `T` seconds
  split at `:segment_s` holds `ceil(T / segment_s)` models, each serialized via
  `Orbis.ReducedOrbit.to_map/1`. Shorter segments cost more bytes but keep every
  query closer to the centre of a fit window, shrinking the residual. The table
  below is the **max** position error over a full day (model-vs-source, measured
  through `drift/3` against the vendored MGEX fixtures, GPST), where both the
  single and piecewise models are fitted over that same full day and evaluated
  across it. (Fitting only part of the span and extrapolating is far worse for
  the single model — see `Orbis.ReducedOrbit`; here the single model gets its
  best case, the whole-span least-squares fit, and piecewise still wins.)

  ### `circular_secular`

  | Orbit class                  | single    | 2 h pw | 4 h pw | 6 h pw |
  |------------------------------|-----------|--------|--------|--------|
  | GPS, e ~ 0.024 (G21)         | ~1 437 km | ~331 km| ~653 km| ~830 km|
  | Galileo, e ~ 1e-4 (E01)      | ~7 km     | ~1 km  | ~3 km  | ~4 km  |
  | BeiDou MEO, e ~ 9e-4 (C21)   | ~54 km    | ~12 km | ~24 km | ~38 km |
  | BeiDou IGSO, e ~ 5e-3 (C08)  | ~533 km   | ~50 km | ~102 km| ~147 km|

  ### `eccentric_secular`

  | Orbit class                  | single  | 2 h pw  | 4 h pw  | 6 h pw  |
  |------------------------------|---------|---------|---------|---------|
  | GPS, e ~ 0.024 (G21)         | ~440 m  | ~90 m   | ~280 m  | ~310 m  |
  | Galileo, e ~ 1e-4 (E01)      | ~780 m  | ~90 m   | ~260 m  | ~380 m  |
  | BeiDou MEO, e ~ 9e-4 (C21)   | ~430 m  | ~120 m  | ~330 m  | ~420 m  |
  | BeiDou IGSO, e ~ 5e-3 (C08)  | ~870 m  | ~20 m   | ~70 m   | ~130 m  |

  For the eccentric model the single whole-day fit is already sub-km; piecewise
  still cuts the max error by roughly 3× to tens of × (most cells ~3-9×, the
  near-circular IGSO as much as ~40×). For the circular model the unmodelled
  `a·e` radial signal dominates and piecewise's benefit is largest in absolute
  terms (hundreds of km off GPS/IGSO). Note the *shorter* the segment the smaller
  the residual: the 2 h split beats the 4 h and 6 h splits in every cell, the
  monotonic accuracy-for-bytes tradeoff.

  These are characterisations, not guarantees; the exact numbers shift with the
  fit window, segment length, and drift cadence. Always measure a given fit with
  `drift/3` against the source.

  ## Segment selection

  Segments tile `[t0, t1]` with no gaps. A query epoch is resolved by finding the
  segment whose half-open interval `[seg_t0, seg_t1)` contains it; the final
  segment is treated as inclusive at the very end so the exact end-of-span epoch
  resolves to the last segment. An epoch exactly on an interior boundary resolves
  to the **later** segment (where it is the in-window start), which is
  deterministic. Selection is `O(segments)`; the segment count is modest and the
  ordered list is binary-searchable if it ever grows large. An epoch before `t0`
  or after `t1` returns `{:error, :out_of_range}`.

  ## Persistence

  `to_map/1` emits a stable, versioned container (string keys, JSON-safe) holding
  the per-segment maps via `Orbis.ReducedOrbit.to_map/1`; `from_map/1` validates
  the version and model and reconstructs, with the same tagged-error discipline
  as the single model.
  """
  alias Orbis.ReducedOrbit

  defstruct version: 1,
            model: "circular_secular",
            frame: "GCRS",
            time_scale: "UTC",
            window: nil,
            segment_s: nil,
            segments: []

  @type epoch :: ReducedOrbit.epoch()
  @type segment :: %{t0: NaiveDateTime.t(), t1: NaiveDateTime.t(), model: ReducedOrbit.t()}

  @type t :: %__MODULE__{
          version: pos_integer(),
          model: String.t(),
          frame: String.t(),
          time_scale: String.t(),
          window: {NaiveDateTime.t(), NaiveDateTime.t()},
          segment_s: number(),
          segments: [segment()]
        }

  # The crate requires four samples to fit one segment; mirror that count for the
  # whole-thing "nothing fit" surface.
  @min_samples 4

  @model_ids ~w(circular_secular eccentric_secular)

  # ------------------------------------------------------------------------
  # Fit
  # ------------------------------------------------------------------------

  @doc """
  Fit a piecewise model over a span, one contiguous `Orbis.ReducedOrbit` segment
  per `:segment_s` seconds.

  ## Sources

    * an `Orbis.SP3` handle — requires `:satellite_id` and `:window`; each
      segment is fitted with `Orbis.ReducedOrbit.fit/2` over its sub-window at
      `:cadence_s`;
    * a list of `{epoch, {x_m, y_m, z_m}}` ECEF samples — partitioned by segment
      interval, each sublist fitted directly.

  ## Options

    * `:window` — `{t0, t1}` epochs bounding the full span (`t1` strictly after
      `t0`, else `:invalid_window`)
    * `:segment_s` — positive segment length in seconds, e.g. `7200`
      (non-positive → `:invalid_segment`)
    * `:cadence_s` — positive SP3 sampling step in seconds
    * `:satellite_id` — e.g. `"G05"` (SP3 source)
    * `:model` — `:circular_secular` (**default**) or `:eccentric_secular`
    * `:time_scale` — for the sample-list source, the scale its epochs are in

  Segments are contiguous (`seg_t1` of one is `seg_t0` of the next); the final
  segment may be shorter. A `:segment_s` at least the full span yields a single
  segment equal to the whole window (piecewise with one segment ≡ single).

  Returns `{:ok, %Orbis.ReducedOrbit.Piecewise{}}` or a tagged error. The error
  set is exactly the single model's fit errors (`{:too_few_samples, got, req}`,
  `:invalid_window`, `:invalid_cadence`, `:satellite_id_required`,
  `{:unsupported_model, m}`, `{:unsupported_source_frame, f}`,
  `{:unsupported_time_scale, s}`, `:transform_unavailable`, …) plus
  `:invalid_segment`. A too-few-samples failure on a non-terminal segment is
  surfaced (a genuinely under-covered interior span is an error, not a silent
  hole); only the terminal short segment may be dropped. If nothing fits at all,
  `{:error, {:too_few_samples, 0, #{@min_samples}}}`.
  """
  @spec fit(Orbis.SP3.t() | [{epoch(), {number(), number(), number()}}], keyword()) ::
          {:ok, t()} | {:error, term()}
  def fit(source, opts \\ []) do
    with {:ok, {t0, t1}} <- validate_window(opts),
         {:ok, segment_s} <- validate_segment_s(Keyword.get(opts, :segment_s)),
         bounds = segment_bounds(t0, t1, segment_s),
         {:ok, segments} <- fit_segments(source, opts, bounds) do
      first = hd(segments)
      # If the trailing short segment was dropped (under-covered), the stored
      # window must shrink to the last kept segment's end so the advertised span
      # equals the actual coverage — no uncovered tail between the last segment
      # and the window end.
      coverage_end = List.last(segments).t1

      {:ok,
       %__MODULE__{
         version: 1,
         model: first.model.model,
         frame: first.model.frame,
         time_scale: first.model.time_scale,
         window: {t0, coverage_end},
         segment_s: segment_s,
         segments: segments
       }}
    end
  end

  # Mirror the single model's window discipline: {t0, t1} with t1 strictly after
  # t0, normalized to NaiveDateTime.
  defp validate_window(opts) do
    case Keyword.get(opts, :window) do
      {t0, t1} ->
        t0n = to_naive(t0)
        t1n = to_naive(t1)

        if NaiveDateTime.after?(t1n, t0n),
          do: {:ok, {t0n, t1n}},
          else: {:error, :invalid_window}

      _ ->
        {:error, :invalid_window}
    end
  end

  # The segment length drives a second-resolution tiling, so it must round to at
  # least one whole second; otherwise the step is 0 and the tiling cannot
  # advance. Return the rounded integer so the stored length is the one used.
  defp validate_segment_s(s) when is_number(s) and s > 0 do
    step = round(s)
    if step >= 1, do: {:ok, step}, else: {:error, :invalid_segment}
  end

  defp validate_segment_s(_), do: {:error, :invalid_segment}

  # Contiguous half-open intervals tiling [t0, t1]; the last may be shorter.
  defp segment_bounds(t0, t1, segment_s) do
    step = round(segment_s)
    do_bounds(t0, t1, step, [])
  end

  defp do_bounds(seg_t0, t1, step, acc) do
    if NaiveDateTime.before?(seg_t0, t1) do
      seg_t1 = NaiveDateTime.add(seg_t0, step, :second)
      seg_t1 = if NaiveDateTime.after?(seg_t1, t1), do: t1, else: seg_t1
      do_bounds(seg_t1, t1, step, [{seg_t0, seg_t1} | acc])
    else
      Enum.reverse(acc)
    end
  end

  defp fit_segments(source, opts, bounds) do
    last_index = length(bounds) - 1

    result =
      bounds
      |> Enum.with_index()
      |> Enum.reduce_while([], fn {{seg_t0, seg_t1}, index}, acc ->
        terminal? = index == last_index

        case fit_one(source, opts, seg_t0, seg_t1) do
          {:ok, model} ->
            {:cont, [%{t0: seg_t0, t1: seg_t1, model: model} | acc]}

          {:error, reason} ->
            cond do
              terminal? and too_few?(reason) ->
                # The trailing short segment can be undercovered; drop it.
                {:cont, acc}

              too_few?(reason) ->
                # An interior gap would break contiguity — surface it.
                {:halt, {:error, reason}}

              true ->
                {:halt, {:error, reason}}
            end
        end
      end)

    case result do
      {:error, _} = err -> err
      [] -> {:error, {:too_few_samples, 0, @min_samples}}
      segments -> {:ok, Enum.reverse(segments)}
    end
  end

  defp fit_one(%Orbis.SP3{} = sp3, opts, seg_t0, seg_t1) do
    forwarded = opts |> Keyword.delete(:segment_s) |> Keyword.put(:window, {seg_t0, seg_t1})
    ReducedOrbit.fit(sp3, forwarded)
  end

  defp fit_one(samples, opts, seg_t0, seg_t1) when is_list(samples) do
    # The list source ignores :window; partition by interval (half-open, the
    # terminal segment inclusive at the very end) and fit the sublist.
    forwarded = Keyword.delete(opts, :segment_s)

    sublist =
      Enum.filter(samples, fn {ep, _} ->
        e = to_naive(ep)

        not NaiveDateTime.before?(e, seg_t0) and
          (NaiveDateTime.before?(e, seg_t1) or NaiveDateTime.compare(e, seg_t1) == :eq)
      end)

    ReducedOrbit.fit(sublist, forwarded)
  end

  defp too_few?({:too_few_samples, _, _}), do: true
  # A zero-span tail window surfaces as :invalid_window from the inner fit.
  defp too_few?(:invalid_window), do: true
  defp too_few?(_), do: false

  # ------------------------------------------------------------------------
  # Evaluation
  # ------------------------------------------------------------------------

  @doc """
  Position of the piecewise model at `epoch`, ECEF (ITRF) meters by default.

  Selects the segment covering `epoch` and delegates to
  `Orbis.ReducedOrbit.position/3`. Pass `frame: :gcrs` for the inertial position.
  An epoch outside the full span returns `{:error, :out_of_range}`.
  """
  @spec position(t(), epoch(), keyword()) :: {:ok, ReducedOrbit.vec3()} | {:error, term()}
  def position(%__MODULE__{} = pw, epoch, opts \\ []) do
    with {:ok, seg} <- select_segment(pw, epoch) do
      ReducedOrbit.position(seg.model, epoch, opts)
    end
  end

  @doc """
  Position and velocity of the piecewise model at `epoch`.

  Selects the covering segment and delegates to
  `Orbis.ReducedOrbit.position_velocity/3`. An epoch outside the full span
  returns `{:error, :out_of_range}`.
  """
  @spec position_velocity(t(), epoch(), keyword()) :: {:ok, map()} | {:error, term()}
  def position_velocity(%__MODULE__{} = pw, epoch, opts \\ []) do
    with {:ok, seg} <- select_segment(pw, epoch) do
      ReducedOrbit.position_velocity(seg.model, epoch, opts)
    end
  end

  @doc """
  Select the segment whose coverage interval contains `epoch`.

  Returns `{:ok, segment}` or `{:error, :out_of_range}`. Interior boundaries
  resolve to the later segment; the exact end-of-span epoch resolves to the last
  segment.
  """
  @spec select_segment(t(), epoch()) :: {:ok, segment()} | {:error, term()}
  def select_segment(%__MODULE__{window: {t0, t1}, segments: segments}, epoch) do
    e = to_naive(epoch)

    if NaiveDateTime.before?(e, t0) or NaiveDateTime.after?(e, t1) do
      {:error, :out_of_range}
    else
      find_segment(segments, e)
    end
  end

  defp find_segment([], _e), do: {:error, :out_of_range}

  defp find_segment([seg], e) do
    # Terminal segment: inclusive upper bound so end-of-span resolves here.
    if not NaiveDateTime.before?(e, seg.t0) and not NaiveDateTime.after?(e, seg.t1),
      do: {:ok, seg},
      else: {:error, :out_of_range}
  end

  defp find_segment([seg | rest], e) do
    # Half-open [t0, t1): a boundary epoch falls through to the next segment.
    if not NaiveDateTime.before?(e, seg.t0) and NaiveDateTime.before?(e, seg.t1),
      do: {:ok, seg},
      else: find_segment(rest, e)
  end

  # ------------------------------------------------------------------------
  # Drift (source-backed, whole span)
  # ------------------------------------------------------------------------

  @doc """
  Evaluate the piecewise model error against the **source** ephemeris over the
  whole span.

  Samples the source across the span and compares each truth sample to the
  covering segment's ECEF position (the single-segment drift NIF is per-model, so
  the piecewise report is composed in Elixir from `position/3`). For an
  `Orbis.SP3` source it samples over `:window` (defaulting to the model's full
  span) at `:cadence_s` for `:satellite_id`; for a sample list it uses those
  directly. Returns

      {:ok, %{per_epoch: [%{epoch:, error_m:}], max_m:, rms_m:, threshold_horizon:,
              requested:, used:}}

  matching the single-segment `Orbis.ReducedOrbit.drift/3` report. Epochs outside
  the model's span are skipped (counted in `requested`, not `used`).
  `threshold_horizon` is the first epoch the ECEF error exceeds `:threshold_m`
  (or `nil`).
  """
  @spec drift(t(), Orbis.SP3.t() | [{epoch(), {number(), number(), number()}}], keyword()) ::
          {:ok, map()} | {:error, term()}
  def drift(%__MODULE__{} = pw, %Orbis.SP3{} = sp3, opts) do
    sat_id = Keyword.get(opts, :satellite_id)

    with :ok <- require_satellite_id(sat_id),
         :ok <- same_scale(pw, sp3),
         {:ok, cadence} <- valid_cadence(Keyword.get(opts, :cadence_s, 900)),
         {:ok, {t0, t1}} <- fetch_window(opts, pw.window),
         {:ok, samples, requested} <- sample_sp3(sp3, sat_id, t0, t1, cadence) do
      run_drift(pw, samples, requested, opts)
    end
  end

  def drift(%__MODULE__{} = pw, samples, opts) when is_list(samples) do
    run_drift(pw, samples, length(samples), opts)
  end

  defp run_drift(_pw, [], _requested, _opts), do: {:error, {:too_few_samples, 0, 1}}

  defp run_drift(pw, samples, requested, opts) do
    with {:ok, threshold_f} <- threshold_value(Keyword.get(opts, :threshold_m, :infinity)) do
      per_epoch =
        samples
        |> Enum.reduce([], fn {ep, {tx, ty, tz}}, acc ->
          epoch = to_naive(ep)

          case position(pw, epoch, frame: :ecef) do
            {:ok, %{x_m: x, y_m: y, z_m: z}} ->
              err = :math.sqrt((x - tx) ** 2 + (y - ty) ** 2 + (z - tz) ** 2)
              [%{epoch: epoch, error_m: err} | acc]

            # Epochs outside the model span (or an evaluation error) are skipped,
            # exactly as the single drift skips epochs off the product coverage.
            {:error, _} ->
              acc
          end
        end)
        |> Enum.reverse()

      finish_drift(per_epoch, threshold_f, requested)
    end
  end

  defp finish_drift([], _threshold_f, _requested), do: {:error, {:too_few_samples, 0, 1}}

  defp finish_drift(per_epoch, threshold_f, requested) do
    errors = Enum.map(per_epoch, & &1.error_m)
    max_m = Enum.max(errors)
    rms_m = :math.sqrt(Enum.reduce(errors, 0.0, fn e, acc -> acc + e * e end) / length(errors))

    horizon =
      Enum.find_value(per_epoch, fn %{epoch: ep, error_m: e} ->
        if e > threshold_f, do: ep
      end)

    {:ok,
     %{
       per_epoch: per_epoch,
       max_m: max_m,
       rms_m: rms_m,
       threshold_horizon: horizon,
       requested: requested,
       used: length(per_epoch)
     }}
  end

  # ------------------------------------------------------------------------
  # Persistence
  # ------------------------------------------------------------------------

  @doc """
  Serialize a piecewise model to a stable, versioned, JSON-safe map (string
  keys). Each segment's model is serialized via `Orbis.ReducedOrbit.to_map/1`.
  See `from_map/1` for the inverse.
  """
  @spec to_map(t()) :: map()
  def to_map(%__MODULE__{} = pw) do
    {t0, t1} = pw.window

    %{
      "version" => pw.version,
      "kind" => "piecewise",
      "model" => pw.model,
      "frame" => pw.frame,
      "time_scale" => pw.time_scale,
      "segment_s" => pw.segment_s,
      "window" => %{
        "start" => NaiveDateTime.to_iso8601(t0),
        "end" => NaiveDateTime.to_iso8601(t1)
      },
      "segments" =>
        Enum.map(pw.segments, fn seg ->
          %{
            "t0" => NaiveDateTime.to_iso8601(seg.t0),
            "t1" => NaiveDateTime.to_iso8601(seg.t1),
            "model" => ReducedOrbit.to_map(seg.model)
          }
        end)
    }
  end

  @doc """
  Reconstruct a piecewise model from a `to_map/1` map. Validates the version and
  model id.

  Returns `{:ok, %Orbis.ReducedOrbit.Piecewise{}}` or
  `{:error, {:unsupported_version, v}}` / `{:error, {:unsupported_model, m}}` /
  `{:error, :malformed_map}`. A segment whose inner model fails `from_map`, or
  whose model id differs from the container, makes the whole map malformed —
  never a raise, never a nil-filled struct.
  """
  @spec from_map(map()) :: {:ok, t()} | {:error, term()}
  def from_map(%{"version" => 1, "kind" => "piecewise", "model" => model} = map)
      when model in @model_ids do
    frame = Map.get(map, "frame", "GCRS")

    with %{"segments" => seg_maps, "window" => window_map, "segment_s" => segment_s} <- map,
         true <- is_list(seg_maps),
         # An empty segment list is a state fit/2 can never produce (it surfaces
         # {:too_few_samples, 0, _} when nothing fits); reject it as malformed.
         false <- seg_maps == [],
         true <- is_number(segment_s),
         {:ok, scale} <- valid_scale(Map.get(map, "time_scale", "UTC")),
         {:ok, {t0, t1}} <- window_from_map(window_map),
         # Each segment's inner model must agree with the container on model id,
         # frame, and time scale, or the persisted scale/frame contract is a lie
         # (position/3 would evaluate mixed-scale segments under a single
         # container scale). fit/2 only ever emits agreeing, gap-free segments.
         {:ok, segments} <- segments_from_maps(seg_maps, model, scale, frame),
         :ok <- validate_contiguity(segments, t0, t1) do
      {:ok,
       %__MODULE__{
         version: 1,
         model: model,
         frame: frame,
         time_scale: scale,
         window: {t0, t1},
         segment_s: segment_s,
         segments: segments
       }}
    else
      _ -> {:error, :malformed_map}
    end
  end

  def from_map(%{"version" => v, "kind" => "piecewise", "model" => model})
      when model in @model_ids, do: {:error, {:unsupported_version, v}}

  def from_map(%{"kind" => "piecewise", "model" => model}),
    do: {:error, {:unsupported_model, model}}

  def from_map(_), do: {:error, :malformed_map}

  defp segments_from_maps(seg_maps, model, scale, frame) do
    Enum.reduce_while(seg_maps, {:ok, []}, fn seg_map, {:ok, acc} ->
      case segment_from_map(seg_map, model, scale, frame) do
        {:ok, seg} -> {:cont, {:ok, [seg | acc]}}
        :error -> {:halt, :error}
      end
    end)
    |> case do
      {:ok, segs} -> {:ok, Enum.reverse(segs)}
      :error -> :error
    end
  end

  defp segment_from_map(
         %{"t0" => t0_iso, "t1" => t1_iso, "model" => model_map},
         container_model,
         container_scale,
         container_frame
       )
       when is_binary(t0_iso) and is_binary(t1_iso) do
    with {:ok, t0} <- NaiveDateTime.from_iso8601(t0_iso),
         {:ok, t1} <- NaiveDateTime.from_iso8601(t1_iso),
         {:ok,
          %ReducedOrbit{
            model: ^container_model,
            time_scale: ^container_scale,
            frame: ^container_frame
          } = model} <-
           ReducedOrbit.from_map(model_map) do
      {:ok, %{t0: t0, t1: t1, model: model}}
    else
      _ -> :error
    end
  end

  defp segment_from_map(_, _, _, _), do: :error

  # A genuine fit/2 product tiles [t0, t1] with gap-free, abutting segments whose
  # ends meet exactly and whose span equals the container window. Reject a
  # persisted map that violates this (an interior gap or a window that does not
  # match the segments).
  defp validate_contiguity([first | _] = segments, t0, t1) do
    last = List.last(segments)

    abutting? =
      segments
      |> Enum.chunk_every(2, 1, :discard)
      |> Enum.all?(fn [a, b] -> NaiveDateTime.compare(a.t1, b.t0) == :eq end)

    if abutting? and NaiveDateTime.compare(first.t0, t0) == :eq and
         NaiveDateTime.compare(last.t1, t1) == :eq do
      :ok
    else
      :error
    end
  end

  defp validate_contiguity(_, _, _), do: :error

  # ------------------------------------------------------------------------
  # Helpers (mirrored from the single model for parity)
  # ------------------------------------------------------------------------

  defp require_satellite_id(nil), do: {:error, :satellite_id_required}
  defp require_satellite_id(_), do: :ok

  defp valid_cadence(c) when is_number(c) and c > 0, do: {:ok, c}
  defp valid_cadence(_), do: {:error, :invalid_cadence}

  defp threshold_value(:infinity), do: {:ok, 1.0e308}
  defp threshold_value(t) when is_number(t) and t >= 0, do: {:ok, t / 1.0}
  defp threshold_value(_), do: {:error, :invalid_threshold}

  @time_scales ~w(UTC TAI TT TDB GPST GST BDT)
  defp valid_scale(s) when s in @time_scales, do: {:ok, s}
  defp valid_scale(_), do: :error

  defp same_scale(%{time_scale: model_scale}, %Orbis.SP3{time_scale: sp3_scale}) do
    if model_scale == sp3_scale,
      do: :ok,
      else: {:error, {:time_scale_mismatch, model_scale, sp3_scale}}
  end

  # drift defaults the window to the model's full span for convenience, but still
  # validates a supplied one with the single model's discipline.
  defp fetch_window(opts, default) do
    case Keyword.get(opts, :window, default) do
      {t0, t1} ->
        t0n = to_naive(t0)
        t1n = to_naive(t1)

        if NaiveDateTime.after?(t1n, t0n),
          do: {:ok, {t0n, t1n}},
          else: {:error, :invalid_window}

      _ ->
        {:error, :invalid_window}
    end
  end

  defp sample_sp3(sp3, sat_id, t0, t1, cadence_s) do
    steps = time_steps(t0, t1, cadence_s)

    samples =
      steps
      |> Enum.reduce([], fn ep, acc ->
        case Orbis.SP3.position(sp3, sat_id, ep) do
          {:ok, %{x_m: x, y_m: y, z_m: z}} -> [{ep, {x, y, z}} | acc]
          {:error, _} -> acc
        end
      end)
      |> Enum.reverse()

    {:ok, samples, length(steps)}
  end

  defp time_steps(t0, t1, cadence_s) do
    span = NaiveDateTime.diff(t1, t0)
    count = trunc(span / cadence_s)
    for k <- 0..count, do: NaiveDateTime.add(t0, round(k * cadence_s), :second)
  end

  defp window_from_map(%{"start" => s, "end" => e}) when is_binary(s) and is_binary(e) do
    with {:ok, t0} <- NaiveDateTime.from_iso8601(s), {:ok, t1} <- NaiveDateTime.from_iso8601(e) do
      {:ok, {t0, t1}}
    else
      _ -> :error
    end
  end

  defp window_from_map(_), do: :error

  defp to_naive(%NaiveDateTime{} = ndt), do: ndt

  defp to_naive({{y, mo, d}, {h, mi, s, us}}) do
    NaiveDateTime.new!(y, mo, d, h, mi, s, {us, 6})
  end

  defp to_naive({{y, mo, d}, {h, mi, s}}) when is_integer(s) do
    NaiveDateTime.new!(y, mo, d, h, mi, s, {0, 6})
  end

  defp to_naive({{y, mo, d}, {h, mi, s}}) when is_float(s) do
    sec = trunc(s)
    micro = round((s - sec) * 1_000_000)
    NaiveDateTime.new!(y, mo, d, h, mi, sec, {micro, 6})
  end
end
