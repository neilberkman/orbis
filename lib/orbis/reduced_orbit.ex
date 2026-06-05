defmodule Orbis.ReducedOrbit do
  @moduledoc """
  A compact, fitted mean-element approximation of a satellite's orbit.

  `Orbis.ReducedOrbit` distills a position track — from a precise SP3 product or
  a list of ECEF samples — into a handful of mean elements that reproduce the
  motion cheaply, for caching, transport, and quick visibility math. It is **not**
  orbit determination and **not** a substitute for SGP4 or precise ephemeris: it
  deliberately discards short-period structure and is honest about the error it
  leaves behind (`rms_m`/`max_m` on the fit, and a source-backed `drift/3`).

  ## Models

  Two models are available, chosen with the `:model` option on `fit/2`:

    * `:circular_secular` (the **default**) — a circular orbit; best for
      near-circular orbits (Galileo).
    * `:eccentric_secular` — adds eccentricity through a nonsingular `(h, k)`
      parameterization, recovering the radial `a·e` signal (~hundreds of km for
      GPS/BeiDou) while degrading smoothly to the circular model as `e -> 0`.

  ### `circular_secular`

  A circular orbit (eccentricity fixed at zero) whose orbital plane precesses at
  a constant nodal rate. At an offset `dt = t - t0` from the reference epoch the
  angles advance linearly,

      u(t)    = arg_lat0 + n * dt          # argument of latitude
      raan(t) = raan0    + raan_rate * dt
      e       = 0

  and the inertial (GCRS) position is the in-plane circle rotated by the node and
  inclination, `r = Rz(raan) * Rx(i) * a * [cos u, sin u, 0]`.

  The nodal rate `raan_rate` is **fitted**, but seeded from the J2 secular nodal
  regression (Vallado, *Fundamentals of Astrodynamics and Applications*):

      raan_rate_j2 = -1.5 * n * J2 * (Re / a)^2 * cos(i)

  Both the fitted value (`raan_rate_rad_s`) and the J2 seed (`raan_rate_j2_rad_s`)
  are kept; `raan_rate_mode` is `"fitted_j2_seeded"`. The model does not claim to
  be a pure J2 propagation.

  ### `eccentric_secular`

  Eight free elements: the four circular plane elements plus `h = e·sin ω`,
  `k = e·cos ω`, `L0` (mean argument of latitude at epoch), and `n`. Derived
  `e = sqrt(h² + k²)` and `ω = atan2(h, k)`. At an offset `dt` the model advances
  `λ = L0 + n·dt`, forms the mean anomaly `M = λ − ω`, solves Kepler's equation
  `E − e·sin E = M`, and places the satellite at radius `r = a(1 − e·cos E)` and
  argument of latitude `u = ω + ν`. The `(h, k)` form is nonsingular: at `e = 0`
  it reproduces `circular_secular` exactly with `arg_lat0 = L0`. The struct then
  carries `h`, `k`, `e`, and `arg_perigee_rad` (= ω).

  ## Frames

  Fitting and evaluation run internally in **GCRS**; positions are returned in
  **ECEF (ITRF) meters by default**, or GCRS via `frame: :gcrs`. ECEF velocity
  includes the Earth-rotation transport term. Sample/query epochs are interpreted
  consistently for the Earth-rotation conversion; the ECEF product (the primary
  output) is self-consistent across the fit, evaluation, and drift.

  ## Expected accuracy

  Representative drift (extrapolated model-vs-source position error): a fit to the
  first ~6 hours of an MGEX SP3 track, drifted over the rest of the day. Numbers
  are the **max** position error over the full day, measured against vendored MGEX
  products (GRG for GPS/Galileo, a trimmed GBM product for BeiDou):

  | Orbit class                   | `circular_secular` | `eccentric_secular` |
  |-------------------------------|--------------------|---------------------|
  | GPS, e ~ 0.024 (G21)          | ~8 100 km          | ~8 km               |
  | GPS, e ~ 0.020 (G02)          | ~9 400 km          | ~11 km              |
  | BeiDou IGSO, e ~ 5e-3 (C08)   | ~2 200 km          | ~12 km              |
  | BeiDou MEO, e ~ 9e-4 (C21)    | ~140 km            | ~5 km               |
  | Galileo, e ~ 1e-4 (E01)       | ~8 km              | ~8 km (≈ circular)  |

  Reading of the table:

    * For **eccentric** orbits the circular model is unusable — the unmodelled
      radial `a·e` signal compounds badly under extrapolation. `:eccentric_secular`
      is the **recommended** model and brings the error down by one-to-three orders
      of magnitude. Note this holds even for the *small* eccentricities of BeiDou
      MEO/IGSO (e ~ 1e-3 to 5e-3): a few-hundred-metre to ~200 km radial signal
      still wrecks the circular extrapolation, so eccentric wins there too.
    * For **near-circular** orbits (Galileo, e ~ 1e-4) both models are comparable:
      the eccentric model is essentially identical to the circular one (~8 km
      either way), so it does not regress and is a safe default when the orbit
      class is unknown.

  These figures are measured through `drift/3` (the public NIF-backed pipeline,
  GPST) against the vendored fixtures; the exact numbers shift slightly with the
  fit window and drift cadence.

  These are characterisations, not guarantees: always measure a given fit with
  `drift/3` against the source. LEO accuracy is characterised once the SGP4/TLE
  source lands. This is a compact approximation for caching / visibility, never a
  substitute for SP3 or SGP4.

  ## Persistence

  `to_map/1` emits a stable, versioned map (frame, model, units, epoch scale,
  elements, fit stats) that `from_map/1` reads back, for caching/transport.
  """
  # Serialization is the explicit, versioned `to_map/1`/`from_map/1` contract
  # (the `fit.window` tuple is not directly JSON-encodable), so no `@derive`.
  alias Orbis.NIF

  defstruct version: 1,
            model: "circular_secular",
            frame: "GCRS",
            raan_rate_mode: "fitted_j2_seeded",
            time_scale: "UTC",
            epoch: nil,
            a_m: nil,
            e: 0.0,
            i_rad: nil,
            raan_rad: nil,
            raan_rate_rad_s: nil,
            raan_rate_j2_rad_s: nil,
            arg_lat_rad: nil,
            mean_motion_rad_s: nil,
            h: nil,
            k: nil,
            arg_perigee_rad: nil,
            fit: nil

  @type epoch ::
          NaiveDateTime.t()
          | {{integer(), integer(), integer()}, {integer(), integer(), number()}}
  @type vec3 :: %{x_m: float(), y_m: float(), z_m: float()}

  @type t :: %__MODULE__{
          version: pos_integer(),
          model: String.t(),
          frame: String.t(),
          raan_rate_mode: String.t(),
          time_scale: String.t(),
          epoch: NaiveDateTime.t(),
          a_m: float(),
          e: float(),
          i_rad: float(),
          raan_rad: float(),
          raan_rate_rad_s: float(),
          raan_rate_j2_rad_s: float(),
          arg_lat_rad: float(),
          mean_motion_rad_s: float(),
          h: float() | nil,
          k: float() | nil,
          arg_perigee_rad: float() | nil,
          fit: map()
        }

  @default_cadence_s 900

  # ------------------------------------------------------------------------
  # Fit
  # ------------------------------------------------------------------------

  @doc """
  Fit a mean-element model to a source orbit.

  ## Sources

    * an `Orbis.SP3` handle — requires `:satellite_id` and `:window`; samples the
      product at `:cadence_s` over the window;
    * a list of `{epoch, {x_m, y_m, z_m}}` ECEF samples — the window is taken from
      the samples; `:frame` must be `:ecef` (the default).

  ## Options

    * `:model` — `:circular_secular` (**default**) or `:eccentric_secular`. The
      circular model fixes eccentricity at zero (best for near-circular orbits);
      the eccentric model recovers the `a·e` radial signal (recommended for GPS
      and other eccentric orbits). See the moduledoc accuracy table.
    * `:satellite_id` — e.g. `"G05"` (SP3 source)
    * `:window` — `{t0, t1}` epochs bounding the fit (SP3 source)
    * `:cadence_s` — positive sampling step in seconds (SP3 source, default `#{@default_cadence_s}`)
    * `:frame` — for the sample-list source, `:ecef` (default)
    * `:time_scale` — for the sample-list source, the scale its epochs are in
      (`"UTC"` default, e.g. `"GPST"`); SP3 sources use the product's own scale

  Epochs are interpreted in the model's time scale (recorded on the result). The
  reference epoch `t0` is the earliest sample, so the result is independent of the
  caller's sample order.

  Returns `{:ok, %Orbis.ReducedOrbit{}}` or a tagged error:
  `{:too_few_samples, got, required}`, `:invalid_window`, `:invalid_cadence`,
  `:satellite_id_required`, `:singular_plane_fit`, `:raan_ambiguous`,
  `{:unsupported_source_frame, frame}`, `{:unsupported_model, model}`,
  `:transform_unavailable`, `:fit_did_not_converge`.
  """
  @spec fit(Orbis.SP3.t() | [{epoch(), {number(), number(), number()}}], keyword()) ::
          {:ok, t()} | {:error, term()}
  def fit(source, opts \\ [])

  def fit(%Orbis.SP3{} = sp3, opts) do
    sat_id = Keyword.get(opts, :satellite_id)

    with :ok <- require_satellite_id(sat_id),
         {:ok, model} <- valid_model(Keyword.get(opts, :model, :circular_secular)),
         {:ok, cadence} <- valid_cadence(Keyword.get(opts, :cadence_s, @default_cadence_s)),
         {:ok, {t0, t1}} <- fetch_window(opts),
         {:ok, samples, requested} <- sample_sp3(sp3, sat_id, t0, t1, cadence) do
      # Epochs are interpreted in the product's own scale (typically GPST),
      # matching Orbis.SP3's contract, so the frame conversion is correct.
      meta = %{
        source: "sp3:#{sat_id}",
        window: {to_naive(t0), to_naive(t1)},
        cadence_s: cadence,
        scale: sp3.time_scale,
        model: model,
        requested: requested
      }

      run_fit(samples, meta)
    end
  end

  def fit(samples, opts) when is_list(samples) do
    with {:ok, model} <- valid_model(Keyword.get(opts, :model, :circular_secular)),
         :ecef <- Keyword.get(opts, :frame, :ecef) do
      meta = %{
        source: "samples",
        window: window_of(samples),
        cadence_s: nil,
        scale: Keyword.get(opts, :time_scale, "UTC"),
        model: model,
        requested: length(samples)
      }

      run_fit(samples, meta)
    else
      {:error, _} = err -> err
      other -> {:error, {:unsupported_source_frame, other}}
    end
  end

  defp run_fit([], _meta), do: {:error, {:too_few_samples, 0, 4}}

  defp run_fit(samples, meta) do
    tuples =
      Enum.map(samples, fn {ep, {x, y, z}} -> {datetime_tuple(ep), x / 1.0, y / 1.0, z / 1.0} end)

    model_str = Atom.to_string(meta.model)

    case safe_nif(fn -> NIF.reduced_orbit_fit(tuples, meta.scale, model_str) end) do
      {:ok, model_atom, epoch_tuple, elements, {rms, max, n_samples}} ->
        fit_stats = %{
          rms_m: rms,
          max_m: max,
          n_samples: n_samples,
          requested: meta.requested,
          window: meta.window,
          cadence_s: meta.cadence_s,
          source: meta.source
        }

        # The reference epoch is the fitter's t0 (earliest sample), returned from
        # Rust, so it is correct regardless of the caller's sample order.
        build_model(model_atom, to_naive(epoch_tuple), meta.scale, elements, fit_stats)

      {:error, reason} ->
        {:error, reason}

      {:nif_raised, _} ->
        {:error, :transform_unavailable}
    end
  end

  # Circular elements: eight floats, no eccentricity vector.
  defp build_model(
         :circular_secular,
         epoch,
         scale,
         [a_m, e, i_rad, raan, raan_rate, raan_rate_j2, arg_lat, n],
         fit_stats
       ) do
    {:ok,
     %__MODULE__{
       model: "circular_secular",
       epoch: epoch,
       time_scale: scale,
       a_m: a_m,
       e: e,
       i_rad: i_rad,
       raan_rad: raan,
       raan_rate_rad_s: raan_rate,
       raan_rate_j2_rad_s: raan_rate_j2,
       arg_lat_rad: arg_lat,
       mean_motion_rad_s: n,
       fit: fit_stats
     }}
  end

  # Eccentric elements: ten floats with h, k appended. `e` and arg_perigee are
  # derived for display; the load-bearing values are h and k.
  defp build_model(
         :eccentric_secular,
         epoch,
         scale,
         [a_m, e, i_rad, raan, raan_rate, raan_rate_j2, arg_lat, n, h, k],
         fit_stats
       ) do
    {:ok,
     %__MODULE__{
       model: "eccentric_secular",
       epoch: epoch,
       time_scale: scale,
       a_m: a_m,
       e: e,
       i_rad: i_rad,
       raan_rad: raan,
       raan_rate_rad_s: raan_rate,
       raan_rate_j2_rad_s: raan_rate_j2,
       arg_lat_rad: arg_lat,
       mean_motion_rad_s: n,
       h: h,
       k: k,
       arg_perigee_rad: arg_perigee(h, k),
       fit: fit_stats
     }}
  end

  # omega = atan2(h, k); undefined (pinned to 0) at e -> 0.
  defp arg_perigee(h, k) when h * h + k * k < 1.0e-24, do: 0.0
  defp arg_perigee(h, k), do: :math.atan2(h, k)

  # ------------------------------------------------------------------------
  # Evaluation
  # ------------------------------------------------------------------------

  @doc """
  Position of the model at `epoch`, ECEF (ITRF) meters by default.

  Pass `frame: :gcrs` for the inertial position. Returns
  `{:ok, %{x_m:, y_m:, z_m:}}` or `{:error, reason}`.
  """
  @spec position(t(), epoch(), keyword()) :: {:ok, vec3()} | {:error, term()}
  def position(%__MODULE__{} = model, epoch, opts \\ []) do
    with {:ok, frame} <- frame_string(Keyword.get(opts, :frame, :ecef)) do
      case safe_nif(fn ->
             NIF.reduced_orbit_position(
               datetime_tuple(model.epoch),
               model.time_scale,
               elements_tuple(model),
               datetime_tuple(epoch),
               frame
             )
           end) do
        {:nif_raised, _} -> {:error, :transform_unavailable}
        {x, y, z} -> {:ok, %{x_m: x, y_m: y, z_m: z}}
      end
    end
  end

  @doc """
  Position and velocity of the model at `epoch`.

  ECEF velocity includes the Earth-rotation transport term. Returns
  `{:ok, %{position: %{x_m:, y_m:, z_m:}, velocity: %{vx_m_s:, vy_m_s:, vz_m_s:}}}`.
  """
  @spec position_velocity(t(), epoch(), keyword()) :: {:ok, map()} | {:error, term()}
  def position_velocity(%__MODULE__{} = model, epoch, opts \\ []) do
    with {:ok, frame} <- frame_string(Keyword.get(opts, :frame, :ecef)) do
      case safe_nif(fn ->
             NIF.reduced_orbit_position_velocity(
               datetime_tuple(model.epoch),
               model.time_scale,
               elements_tuple(model),
               datetime_tuple(epoch),
               frame
             )
           end) do
        {:nif_raised, _} ->
          {:error, :transform_unavailable}

        {{x, y, z}, {vx, vy, vz}} ->
          {:ok,
           %{
             position: %{x_m: x, y_m: y, z_m: z},
             velocity: %{vx_m_s: vx, vy_m_s: vy, vz_m_s: vz}
           }}
      end
    end
  end

  # ------------------------------------------------------------------------
  # Drift (source-backed)
  # ------------------------------------------------------------------------

  @doc """
  Evaluate the model error against the **source** ephemeris over a horizon.

  This compares the model to fresh truth samples (not to itself): for an
  `Orbis.SP3` source it samples the product over `:window` at `:cadence_s` for
  `:satellite_id`; for a list of `{epoch, {x_m, y_m, z_m}}` samples it uses those
  directly. Returns

      {:ok, %{per_epoch: [%{epoch:, error_m:}], max_m:, rms_m:, threshold_horizon:}}

  where `threshold_horizon` is the first epoch the ECEF error exceeds
  `:threshold_m` (or `nil` if it never does / no threshold given).
  """
  @spec drift(t(), Orbis.SP3.t() | [{epoch(), {number(), number(), number()}}], keyword()) ::
          {:ok, map()} | {:error, term()}
  def drift(%__MODULE__{} = model, %Orbis.SP3{} = sp3, opts) do
    sat_id = Keyword.get(opts, :satellite_id)

    with :ok <- require_satellite_id(sat_id),
         :ok <- same_scale(model, sp3),
         {:ok, cadence} <- valid_cadence(Keyword.get(opts, :cadence_s, @default_cadence_s)),
         {:ok, {t0, t1}} <- fetch_window(opts),
         {:ok, samples, _requested} <- sample_sp3(sp3, sat_id, t0, t1, cadence) do
      run_drift(model, samples, opts)
    end
  end

  def drift(%__MODULE__{} = model, samples, opts) when is_list(samples) do
    run_drift(model, samples, opts)
  end

  defp run_drift(_model, [], _opts), do: {:error, {:too_few_samples, 0, 1}}

  defp run_drift(model, samples, opts) do
    with {:ok, threshold_f} <- threshold_value(Keyword.get(opts, :threshold_m, :infinity)) do
      epochs = Enum.map(samples, fn {ep, _} -> to_naive(ep) end)

      truth =
        Enum.map(samples, fn {ep, {x, y, z}} ->
          {datetime_tuple(ep), x / 1.0, y / 1.0, z / 1.0}
        end)

      run_drift_nif(model, truth, threshold_f, epochs)
    end
  end

  defp run_drift_nif(model, truth, threshold_f, epochs) do
    case safe_nif(fn ->
           NIF.reduced_orbit_drift(
             datetime_tuple(model.epoch),
             model.time_scale,
             elements_tuple(model),
             truth,
             threshold_f
           )
         end) do
      {:nif_raised, _} ->
        {:error, :transform_unavailable}

      {errors, max_m, rms_m, idx} ->
        per_epoch =
          epochs
          |> Enum.zip(errors)
          |> Enum.map(fn {ep, err} -> %{epoch: ep, error_m: err} end)

        horizon = if idx >= 0, do: Enum.at(epochs, idx)
        {:ok, %{per_epoch: per_epoch, max_m: max_m, rms_m: rms_m, threshold_horizon: horizon}}
    end
  end

  # ------------------------------------------------------------------------
  # Persistence
  # ------------------------------------------------------------------------

  @doc """
  Serialize a fitted model to a stable, versioned map (string keys) for caching
  or transport. See `from_map/1` for the inverse.
  """
  @spec to_map(t()) :: map()
  def to_map(%__MODULE__{} = m) do
    %{
      "version" => m.version,
      "model" => m.model,
      "frame" => m.frame,
      "time_scale" => m.time_scale,
      "epoch" => NaiveDateTime.to_iso8601(m.epoch),
      "elements" =>
        %{
          "a_m" => m.a_m,
          "e" => m.e,
          "i_rad" => m.i_rad,
          "raan_rad" => m.raan_rad,
          "raan_rate_rad_s" => m.raan_rate_rad_s,
          "raan_rate_j2_rad_s" => m.raan_rate_j2_rad_s,
          "raan_rate_mode" => m.raan_rate_mode,
          "arg_lat_rad" => m.arg_lat_rad,
          "mean_motion_rad_s" => m.mean_motion_rad_s
        }
        |> maybe_put_eccentric(m),
      "fit" => %{
        "rms_m" => m.fit.rms_m,
        "max_m" => m.fit.max_m,
        "n_samples" => m.fit.n_samples,
        "requested" => Map.get(m.fit, :requested),
        "cadence_s" => m.fit.cadence_s,
        "source" => m.fit.source,
        "window" => window_to_map(m.fit.window)
      },
      "units" => %{"length" => "m", "angle" => "rad", "rate" => "rad/s", "time" => "s"}
    }
  end

  # The eccentric model adds the eccentricity vector and derived argument of
  # perigee to the elements map; the circular model's map is unchanged.
  defp maybe_put_eccentric(elements, %__MODULE__{model: "eccentric_secular"} = m) do
    Map.merge(elements, %{
      "h" => m.h,
      "k" => m.k,
      "arg_perigee_rad" => m.arg_perigee_rad
    })
  end

  defp maybe_put_eccentric(elements, _m), do: elements

  @doc """
  Reconstruct a model from a `to_map/1` map. Validates the version and model id.

  Returns `{:ok, %Orbis.ReducedOrbit{}}` or `{:error, {:unsupported_version, v}}` /
  `{:error, {:unsupported_model, model}}` / `{:error, :malformed_map}`.
  """
  @spec from_map(map()) :: {:ok, t()} | {:error, term()}
  def from_map(%{"version" => 1, "model" => "circular_secular"} = map) do
    with %{"elements" => el, "fit" => fit, "epoch" => epoch_iso} <- map,
         {:ok, epoch} <- NaiveDateTime.from_iso8601(epoch_iso) do
      {:ok,
       %__MODULE__{
         version: 1,
         model: "circular_secular",
         frame: Map.get(map, "frame", "GCRS"),
         time_scale: Map.get(map, "time_scale", "UTC"),
         raan_rate_mode: Map.get(el, "raan_rate_mode", "fitted_j2_seeded"),
         epoch: epoch,
         a_m: el["a_m"],
         e: Map.get(el, "e", 0.0),
         i_rad: el["i_rad"],
         raan_rad: el["raan_rad"],
         raan_rate_rad_s: el["raan_rate_rad_s"],
         raan_rate_j2_rad_s: el["raan_rate_j2_rad_s"],
         arg_lat_rad: el["arg_lat_rad"],
         mean_motion_rad_s: el["mean_motion_rad_s"],
         fit: %{
           rms_m: fit["rms_m"],
           max_m: fit["max_m"],
           n_samples: fit["n_samples"],
           requested: fit["requested"],
           cadence_s: fit["cadence_s"],
           source: fit["source"],
           window: window_from_map(fit["window"])
         }
       }}
    else
      _ -> {:error, :malformed_map}
    end
  end

  def from_map(%{"version" => 1, "model" => "eccentric_secular"} = map) do
    with %{"elements" => el, "fit" => fit, "epoch" => epoch_iso} <- map,
         h when is_number(h) <- el["h"],
         k when is_number(k) <- el["k"],
         {:ok, epoch} <- NaiveDateTime.from_iso8601(epoch_iso) do
      {:ok,
       %__MODULE__{
         version: 1,
         model: "eccentric_secular",
         frame: Map.get(map, "frame", "GCRS"),
         time_scale: Map.get(map, "time_scale", "UTC"),
         raan_rate_mode: Map.get(el, "raan_rate_mode", "fitted_j2_seeded"),
         epoch: epoch,
         a_m: el["a_m"],
         e: el["e"],
         i_rad: el["i_rad"],
         raan_rad: el["raan_rad"],
         raan_rate_rad_s: el["raan_rate_rad_s"],
         raan_rate_j2_rad_s: el["raan_rate_j2_rad_s"],
         arg_lat_rad: el["arg_lat_rad"],
         mean_motion_rad_s: el["mean_motion_rad_s"],
         h: h,
         k: k,
         arg_perigee_rad: Map.get(el, "arg_perigee_rad") || arg_perigee(h, k),
         fit: %{
           rms_m: fit["rms_m"],
           max_m: fit["max_m"],
           n_samples: fit["n_samples"],
           requested: fit["requested"],
           cadence_s: fit["cadence_s"],
           source: fit["source"],
           window: window_from_map(fit["window"])
         }
       }}
    else
      _ -> {:error, :malformed_map}
    end
  end

  def from_map(%{"version" => v, "model" => "circular_secular"}),
    do: {:error, {:unsupported_version, v}}

  def from_map(%{"version" => v, "model" => "eccentric_secular"}),
    do: {:error, {:unsupported_version, v}}

  def from_map(%{"model" => model}), do: {:error, {:unsupported_model, model}}
  def from_map(_), do: {:error, :malformed_map}

  # ------------------------------------------------------------------------
  # Helpers
  # ------------------------------------------------------------------------

  defp require_satellite_id(nil), do: {:error, :satellite_id_required}
  defp require_satellite_id(_), do: :ok

  defp valid_model(:circular_secular), do: {:ok, :circular_secular}
  defp valid_model(:eccentric_secular), do: {:ok, :eccentric_secular}
  defp valid_model(other), do: {:error, {:unsupported_model, other}}

  # No threshold means no crossing is ever reported (an unreachable bound). A
  # finite threshold must be a non-negative number; NaN/negative are rejected.
  defp threshold_value(:infinity), do: {:ok, 1.0e308}
  defp threshold_value(t) when is_number(t) and t >= 0, do: {:ok, t / 1.0}
  defp threshold_value(_), do: {:error, :invalid_threshold}

  # The SP3 sampling step must be a positive, finite number of seconds.
  defp valid_cadence(c) when is_number(c) and c > 0, do: {:ok, c}
  defp valid_cadence(_), do: {:error, :invalid_cadence}

  # Drifting against an SP3 product only makes sense if the model's epochs are in
  # the same scale as the product, otherwise equal calendar values denote
  # different physical instants.
  defp same_scale(%{time_scale: model_scale}, %Orbis.SP3{time_scale: sp3_scale}) do
    if model_scale == sp3_scale,
      do: :ok,
      else: {:error, {:time_scale_mismatch, model_scale, sp3_scale}}
  end

  defp fetch_window(opts) do
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

  defp sample_sp3(sp3, sat_id, t0, t1, cadence_s) do
    steps = time_steps(t0, t1, cadence_s)

    samples =
      steps
      |> Enum.reduce([], fn ep, acc ->
        case Orbis.SP3.position(sp3, sat_id, ep) do
          {:ok, %{x_m: x, y_m: y, z_m: z}} -> [{ep, {x, y, z}} | acc]
          # An epoch outside the product's coverage is skipped; the requested
          # count is returned so the caller can see the horizon was clipped.
          {:error, _} -> acc
        end
      end)
      |> Enum.reverse()

    {:ok, samples, length(steps)}
  end

  defp time_steps(t0, t1, cadence_s) do
    t0n = to_naive(t0)
    t1n = to_naive(t1)
    span = NaiveDateTime.diff(t1n, t0n)
    count = trunc(span / cadence_s)
    for k <- 0..count, do: NaiveDateTime.add(t0n, round(k * cadence_s), :second)
  end

  defp window_of([]), do: nil

  defp window_of(samples) do
    epochs = Enum.map(samples, fn {ep, _} -> to_naive(ep) end)
    {Enum.min_by(epochs, &NaiveDateTime.to_erl/1), Enum.max_by(epochs, &NaiveDateTime.to_erl/1)}
  end

  # The flat element list the NIF consumes. Circular is eight floats (unchanged);
  # eccentric appends h, k (ten floats), which the NIF reads back by length.
  defp elements_tuple(%__MODULE__{model: "eccentric_secular"} = m) do
    [
      m.a_m,
      m.e,
      m.i_rad,
      m.raan_rad,
      m.raan_rate_rad_s,
      m.raan_rate_j2_rad_s,
      m.arg_lat_rad,
      m.mean_motion_rad_s,
      m.h,
      m.k
    ]
  end

  defp elements_tuple(%__MODULE__{} = m) do
    [
      m.a_m,
      m.e,
      m.i_rad,
      m.raan_rad,
      m.raan_rate_rad_s,
      m.raan_rate_j2_rad_s,
      m.arg_lat_rad,
      m.mean_motion_rad_s
    ]
  end

  defp frame_string(:ecef), do: {:ok, "ecef"}
  defp frame_string(:gcrs), do: {:ok, "gcrs"}
  defp frame_string(other), do: {:error, {:unsupported_frame, other}}

  defp to_naive(%NaiveDateTime{} = ndt), do: ndt

  defp to_naive({{y, mo, d}, {h, mi, s, us}}) do
    NaiveDateTime.new!(y, mo, d, h, mi, s, {us, 6})
  end

  defp to_naive({{y, mo, d}, {h, mi, s}}) do
    {sec, micro} = split_seconds(s)
    NaiveDateTime.new!(y, mo, d, h, mi, sec, {micro, 6})
  end

  defp datetime_tuple(%NaiveDateTime{} = ndt) do
    {{ndt.year, ndt.month, ndt.day}, {ndt.hour, ndt.minute, ndt.second, elem(ndt.microsecond, 0)}}
  end

  defp datetime_tuple({{_, _, _}, {_, _, _}} = tuple), do: datetime_tuple(to_naive(tuple))

  defp split_seconds(s) when is_integer(s), do: {s, 0}

  defp split_seconds(s) when is_float(s) do
    sec = trunc(s)
    {sec, round((s - sec) * 1_000_000)}
  end

  defp window_to_map(nil), do: nil

  defp window_to_map({t0, t1}),
    do: %{"start" => NaiveDateTime.to_iso8601(t0), "end" => NaiveDateTime.to_iso8601(t1)}

  defp window_from_map(nil), do: nil

  defp window_from_map(%{"start" => s, "end" => e}) do
    with {:ok, t0} <- NaiveDateTime.from_iso8601(s), {:ok, t1} <- NaiveDateTime.from_iso8601(e) do
      {t0, t1}
    else
      _ -> nil
    end
  end

  # Run a NIF call, mapping an unexpected native raise to a sentinel the callers
  # turn into `:transform_unavailable`.
  defp safe_nif(fun) do
    fun.()
  rescue
    e in ErlangError -> {:nif_raised, e.original}
  end
end
