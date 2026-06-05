defmodule Orbis.GnssGeometry do
  @moduledoc """
  Satellite-geometry and mission-planning layer above the GNSS observables:
  from a static receiver position and a precise (SP3) ephemeris, answer the
  three planning questions — which satellites are visible, how good is the
  geometry (dilution of precision), and when does each satellite rise and set.

  This module solves no positioning problem; it reads satellite states through
  `Orbis.GnssObservables` and applies standard textbook GNSS geometry.

  ## Visibility

  A satellite is *visible* when its topocentric elevation is at or above an
  elevation mask. Azimuth and elevation come from `Orbis.GnssObservables`, which
  rotates the receiver-to-satellite line of sight into the local east-north-up
  (ENU) frame at the receiver's geodetic latitude/longitude.

  ## Dilution of precision (DOP)

  Dilution of precision summarises how the receiver-to-satellite geometry maps
  range-measurement noise into solution uncertainty. From a design (geometry)
  matrix `G` whose rows are the line-of-sight unit vectors plus a receiver-clock
  column, and an optional diagonal weight matrix `W`, the cofactor matrix is

      Q = (G^T W G)^-1

  a 4x4 symmetric matrix ordered `[x, y, z, clock]`. The position block is in
  ECEF metres and the clock state is in the same length unit as the ranges.

  ### Sign and column convention

  Each row is `[-e_x, -e_y, -e_z, 1]`, where `e` is the **ECEF** receiver-to-
  satellite unit line of sight (the partial derivative of the predicted range
  with respect to the receiver position is `-e`; the clock column is `+1`). The
  geometry matrix is therefore built in ECEF, exactly as
  `Orbis.PointPositioning` builds it — the horizontal/vertical split is taken
  *after* inverting, by rotating the 3x3 position block into the local ENU frame
  at the receiver's geodetic latitude/longitude:

      R = [[-sin l,         cos l,        0   ],
           [-sin p cos l,  -sin p sin l,  cos p],
           [ cos p cos l,   cos p sin l,  sin p]]

  with `p` the geodetic latitude and `l` the longitude (radians); the rotated
  block is `Q_enu = R Q_pos R^T`. The DOP scalars are then

    * `pdop = sqrt(qE + qN + qU)` (the ENU position block),
    * `hdop = sqrt(qE + qN)`,
    * `vdop = sqrt(qU)`,
    * `tdop = sqrt(Q[3][3])` (the clock variance),
    * `gdop = sqrt(Q[0][0] + Q[1][1] + Q[2][2] + Q[3][3])` (the cofactor trace,
      which is rotation invariant, so it equals the ENU-frame trace).

  ### Weights

  The default is the unweighted geometric DOP (`W = I`), the standard textbook
  cofactor `(G^T G)^-1`. An elevation weighting (`weights: :elevation`, with
  `w = sin^2(elevation)`) is also available; it reproduces the weighting that a
  least-squares positioning solve applies, and is what lets the DOP here be
  cross-checked component-for-component against `Orbis.PointPositioning`'s
  reported DOP for the same satellite set and epoch.

  ### Limitation

  This is a single-receiver-clock (single-system) DOP. A mixed-constellation
  geometry with one receiver clock per system (extra clock columns) is out of
  scope; restrict the visible set to one system (e.g. `systems: ["G"]`) for a
  well-posed DOP.

  ## Passes

  A *pass* is a contiguous interval over which a satellite stays above the mask.
  Rise and set are detected by threshold-crossing on the sampled elevation, so
  they are resolved only to the sampling step `step_seconds`: a finer step gives
  finer rise/set epochs.
  """

  alias Orbis.GnssObservables
  alias Orbis.SP3

  @default_mask_deg 5.0
  @deg_to_rad :math.pi() / 180.0

  @type receiver ::
          {number(), number(), number()} | %{x_m: number(), y_m: number(), z_m: number()}

  @type visible_sat :: %{
          satellite_id: String.t(),
          elevation_deg: float(),
          azimuth_deg: float()
        }

  @type dop_result :: %{
          gdop: float(),
          pdop: float(),
          hdop: float(),
          vdop: float(),
          tdop: float(),
          n_satellites: non_neg_integer(),
          satellites: [String.t()]
        }

  # --- visibility -----------------------------------------------------------

  @doc """
  List the satellites visible from `receiver` at `epoch`, above the elevation
  mask, sorted by elevation descending.

  ## Options

    * `:elevation_mask_deg` - minimum elevation in degrees (default `5.0`); a
      satellite is included iff its elevation is at or above this value.
    * `:systems` - keep only these constellations, given as leading-letter
      strings (e.g. `["G"]` for GPS, `["G", "E"]` for GPS + Galileo). Default:
      keep all systems.

  Returns a list of `%{satellite_id, elevation_deg, azimuth_deg}` or
  `{:error, :invalid_receiver}` for a malformed receiver. Never raises.
  """
  @spec visible(SP3.t(), receiver(), NaiveDateTime.t(), keyword()) ::
          [visible_sat()] | {:error, :invalid_receiver}
  def visible(%SP3{} = sp3, receiver, %NaiveDateTime{} = epoch, opts \\ []) do
    with {:ok, rx} <- normalize_receiver(receiver) do
      mask = Keyword.get(opts, :elevation_mask_deg, @default_mask_deg)
      systems = Keyword.get(opts, :systems)

      sp3
      |> GnssObservables.predict_all(rx, epoch, light_time: false)
      |> Enum.flat_map(fn
        {sat_id, {:ok, obs}} ->
          if system_allowed?(sat_id, systems) and obs.elevation_deg >= mask do
            [
              %{
                satellite_id: sat_id,
                elevation_deg: obs.elevation_deg,
                azimuth_deg: obs.azimuth_deg
              }
            ]
          else
            []
          end

        {_sat_id, {:error, _}} ->
          []
      end)
      |> Enum.sort_by(& &1.elevation_deg, :desc)
    end
  end

  # --- dilution of precision ------------------------------------------------

  @doc """
  Dilution of precision for the visible satellites at `epoch`.

  Builds the geometry matrix from the visible satellites' ECEF line-of-sight
  unit vectors (rows `[-e_x, -e_y, -e_z, 1]`), forms `Q = (G^T W G)^-1`, rotates
  the position block into ENU, and returns all five DOP scalars plus the
  satellite count and ids.

  ## Options

  In addition to the `visible/4` options (`:elevation_mask_deg`, `:systems`):

    * `:weights` - `:unit` (default, `W = I`, the standard geometric DOP) or
      `:elevation` (`w = sin^2(elevation)`, the least-squares weighting).
    * `:light_time` - apply the light-time / Sagnac line-of-sight corrections
      when forming the geometry (default `false`, the planning value). Set
      `true` to match a converged positioning geometry exactly.
    * `:satellites` - an explicit list of satellite ids to use instead of the
      visibility scan (still subject to predicting successfully); useful to pin
      the geometry to a known set.

  Returns `%{gdop, pdop, hdop, vdop, tdop, n_satellites, satellites}` or a
  tagged error: `{:error, :invalid_receiver}`, `{:error, :too_few_satellites}`
  (fewer than four usable directions), or `{:error, :singular_geometry}`. Never
  raises.
  """
  @spec dop(SP3.t(), receiver(), NaiveDateTime.t(), keyword()) ::
          dop_result() | {:error, atom()}
  def dop(%SP3{} = sp3, receiver, %NaiveDateTime{} = epoch, opts \\ []) do
    with {:ok, rx} <- normalize_receiver(receiver) do
      weights = Keyword.get(opts, :weights, :unit)
      light_time? = Keyword.get(opts, :light_time, false)

      sats = dop_satellite_ids(sp3, rx, epoch, opts)

      rows_weights =
        Enum.flat_map(sats, fn sat_id ->
          case GnssObservables.predict(sp3, sat_id, rx, epoch,
                 light_time: light_time?,
                 sagnac: light_time?
               ) do
            {:ok, obs} ->
              {ex, ey, ez} = obs.los_unit
              w = weight_for(weights, obs.elevation_deg)
              [{sat_id, [-ex, -ey, -ez, 1.0], w}]

            {:error, _} ->
              []
          end
        end)

      compute_dop(rows_weights, rx)
    end
  end

  # Resolve the satellite set the DOP is built over: either an explicit
  # `:satellites` list, or the visibility scan under the mask/system filter.
  defp dop_satellite_ids(sp3, rx, epoch, opts) do
    case Keyword.get(opts, :satellites) do
      nil ->
        sp3
        |> visible(rx, epoch, opts)
        |> Enum.map(& &1.satellite_id)

      list when is_list(list) ->
        list
    end
  end

  defp compute_dop(rows_weights, rx) do
    n = length(rows_weights)

    if n < 4 do
      {:error, :too_few_satellites}
    else
      rows = Enum.map(rows_weights, fn {_id, row, _w} -> row end)
      ws = Enum.map(rows_weights, fn {_id, _row, w} -> w end)
      ids = Enum.map(rows_weights, fn {id, _row, _w} -> id end)

      a = normal_matrix(rows, ws)

      case inv4(a) do
        :singular ->
          {:error, :singular_geometry}

        {:ok, q} ->
          finalize_dop(q, rx, n, ids)
      end
    end
  end

  defp finalize_dop(q, rx, n, ids) do
    {lat_rad, lon_rad} = receiver_lat_lon_rad(rx)
    r = ecef_to_enu_rotation(lat_rad, lon_rad)
    enu = rotate_pos_block(q, r)

    qe = elem(elem(enu, 0), 0)
    qn = elem(elem(enu, 1), 1)
    qu = elem(elem(enu, 2), 2)
    qt = at(q, 3, 3)

    gdop_arg = at(q, 0, 0) + at(q, 1, 1) + at(q, 2, 2) + at(q, 3, 3)
    pdop_arg = qe + qn + qu
    hdop_arg = qe + qn
    vdop_arg = qu
    tdop_arg = qt

    if Enum.all?([gdop_arg, pdop_arg, hdop_arg, vdop_arg, tdop_arg], &finite_nonneg?/1) do
      %{
        gdop: :math.sqrt(gdop_arg),
        pdop: :math.sqrt(pdop_arg),
        hdop: :math.sqrt(hdop_arg),
        vdop: :math.sqrt(vdop_arg),
        tdop: :math.sqrt(tdop_arg),
        n_satellites: n,
        satellites: ids
      }
    else
      {:error, :singular_geometry}
    end
  end

  defp weight_for(:unit, _el_deg), do: 1.0

  defp weight_for(:elevation, el_deg) do
    s = :math.sin(el_deg * @deg_to_rad)
    s * s
  end

  # A variance argument is usable only when non-negative: a negative diagonal
  # entry signals a numerically singular geometry. (On the BEAM, float overflow
  # or division by zero raises rather than producing an infinity or NaN, so a
  # non-negativity test is sufficient here.)
  defp finite_nonneg?(arg), do: arg >= 0.0

  # --- time series ----------------------------------------------------------

  @doc """
  Per-epoch dilution of precision over a time window.

  Samples `{t0, t1}` (inclusive of `t0`, up to and including `t1`) every
  `step_seconds` and computes `dop/4` at each sample. Returns a list of
  `%{epoch, gdop, pdop, hdop, vdop, tdop, n_satellites, satellites}` for the
  epochs whose geometry yields a finite DOP; epochs with too few satellites or a
  singular geometry are skipped. An empty or inverted window returns `[]`.

  `opts` are the `dop/4` options. Errors from a malformed receiver propagate as
  `{:error, :invalid_receiver}`.
  """
  @spec dop_series(
          SP3.t(),
          receiver(),
          {NaiveDateTime.t(), NaiveDateTime.t()},
          pos_integer(),
          keyword()
        ) ::
          [map()] | {:error, :invalid_receiver}
  def dop_series(%SP3{} = sp3, receiver, {t0, t1}, step_seconds, opts \\ []) do
    with {:ok, rx} <- normalize_receiver(receiver) do
      {t0, t1}
      |> epochs(step_seconds)
      |> Enum.flat_map(fn epoch ->
        case dop(sp3, rx, epoch, opts) do
          {:error, _} -> []
          d -> [Map.put(d, :epoch, epoch)]
        end
      end)
    end
  end

  @doc """
  Per-epoch count of visible satellites over a time window.

  Samples `{t0, t1}` every `step_seconds` and returns a list of
  `%{epoch, n_visible}`. An empty or inverted window returns `[]`. `opts` are the
  `visible/4` options. A malformed receiver yields `{:error, :invalid_receiver}`.
  """
  @spec visibility_series(
          SP3.t(),
          receiver(),
          {NaiveDateTime.t(), NaiveDateTime.t()},
          pos_integer(),
          keyword()
        ) ::
          [%{epoch: NaiveDateTime.t(), n_visible: non_neg_integer()}]
          | {:error, :invalid_receiver}
  def visibility_series(%SP3{} = sp3, receiver, {t0, t1}, step_seconds, opts \\ []) do
    with {:ok, rx} <- normalize_receiver(receiver) do
      {t0, t1}
      |> epochs(step_seconds)
      |> Enum.map(fn epoch ->
        %{epoch: epoch, n_visible: length(visible(sp3, rx, epoch, opts))}
      end)
    end
  end

  # --- passes ---------------------------------------------------------------

  @doc """
  Rise / peak / set passes for each satellite over a time window.

  Samples `{t0, t1}` every `step_seconds` and, for each satellite, splits the
  samples into contiguous runs above the elevation mask. Each run is one pass:

      %{
        satellite_id: String.t(),
        rise_epoch: NaiveDateTime.t(),
        set_epoch: NaiveDateTime.t(),
        peak_elevation_deg: float(),
        peak_epoch: NaiveDateTime.t()
      }

  `rise_epoch` is the first sample above the mask and `set_epoch` the last; both
  are resolved only to the sampling step (the true crossing lies within one
  `step_seconds` of the reported epoch). A satellite already above the mask at
  `t0`, or still above it at `t1`, yields a pass clamped to the window. `opts`
  are the `visible/4` options. A malformed receiver yields
  `{:error, :invalid_receiver}`.
  """
  @spec passes(
          SP3.t(),
          receiver(),
          {NaiveDateTime.t(), NaiveDateTime.t()},
          pos_integer(),
          keyword()
        ) ::
          [map()] | {:error, :invalid_receiver}
  def passes(%SP3{} = sp3, receiver, {t0, t1}, step_seconds, opts \\ []) do
    with {:ok, rx} <- normalize_receiver(receiver) do
      mask = Keyword.get(opts, :elevation_mask_deg, @default_mask_deg)
      systems = Keyword.get(opts, :systems)
      epochs = epochs({t0, t1}, step_seconds)

      # satellite_id => [{epoch, elevation_deg}] in sampling order, above mask
      sp3
      |> SP3.satellite_ids()
      |> Enum.filter(&system_allowed?(&1, systems))
      |> Enum.flat_map(fn sat_id ->
        samples =
          Enum.flat_map(epochs, fn epoch ->
            case GnssObservables.predict(sp3, sat_id, rx, epoch, light_time: false) do
              {:ok, obs} when obs.elevation_deg >= mask -> [{epoch, obs.elevation_deg}]
              {:ok, _} -> [:below]
              {:error, _} -> [:below]
            end
          end)

        samples
        |> split_runs()
        |> Enum.map(&pass_from_run(sat_id, &1))
      end)
      |> Enum.sort_by(& &1.rise_epoch, NaiveDateTime)
    end
  end

  # Split a list of `{epoch, el}` | `:below` markers into the contiguous runs of
  # above-mask samples.
  defp split_runs(samples) do
    samples
    |> Enum.chunk_by(&(&1 == :below))
    |> Enum.reject(fn chunk -> match?([:below | _], chunk) end)
    |> Enum.reject(&(&1 == []))
  end

  defp pass_from_run(sat_id, run) do
    {rise_epoch, _} = hd(run)
    {set_epoch, _} = List.last(run)
    {peak_epoch, peak_el} = Enum.max_by(run, fn {_e, el} -> el end)

    %{
      satellite_id: sat_id,
      rise_epoch: rise_epoch,
      set_epoch: set_epoch,
      peak_elevation_deg: peak_el,
      peak_epoch: peak_epoch
    }
  end

  # --- epoch sampling -------------------------------------------------------

  defp epochs({t0, t1}, step_seconds) when is_integer(step_seconds) and step_seconds > 0 do
    if NaiveDateTime.after?(t0, t1) do
      []
    else
      build_epochs(t0, t1, step_seconds, [])
    end
  end

  defp build_epochs(t, t1, step_seconds, acc) do
    if NaiveDateTime.after?(t, t1) do
      Enum.reverse(acc)
    else
      build_epochs(NaiveDateTime.add(t, step_seconds, :second), t1, step_seconds, [t | acc])
    end
  end

  # --- linear algebra (4x4 cofactor inverse, ENU rotation) ------------------

  # Normal matrix A = G^T W G, accumulated left-to-right over the rows in input
  # order, the weight applied between the two design entries.
  defp normal_matrix(rows, weights) do
    pairs = Enum.zip(rows, weights)

    for i <- 0..3 do
      for j <- 0..3 do
        Enum.reduce(pairs, 0.0, fn {row, w}, acc ->
          acc + Enum.at(row, i) * w * Enum.at(row, j)
        end)
      end
      |> List.to_tuple()
    end
    |> List.to_tuple()
  end

  defp at(m, i, j), do: elem(elem(m, i), j)

  # Determinant of a 4x4 by first-row cofactor expansion, fixed term order.
  defp det4(a) do
    m01 = at(a, 2, 0) * at(a, 3, 1) - at(a, 2, 1) * at(a, 3, 0)
    m02 = at(a, 2, 0) * at(a, 3, 2) - at(a, 2, 2) * at(a, 3, 0)
    m03 = at(a, 2, 0) * at(a, 3, 3) - at(a, 2, 3) * at(a, 3, 0)
    m12 = at(a, 2, 1) * at(a, 3, 2) - at(a, 2, 2) * at(a, 3, 1)
    m13 = at(a, 2, 1) * at(a, 3, 3) - at(a, 2, 3) * at(a, 3, 1)
    m23 = at(a, 2, 2) * at(a, 3, 3) - at(a, 2, 3) * at(a, 3, 2)

    c0 = at(a, 1, 1) * m23 - at(a, 1, 2) * m13 + at(a, 1, 3) * m12
    c1 = at(a, 1, 0) * m23 - at(a, 1, 2) * m03 + at(a, 1, 3) * m02
    c2 = at(a, 1, 0) * m13 - at(a, 1, 1) * m03 + at(a, 1, 3) * m01
    c3 = at(a, 1, 0) * m12 - at(a, 1, 1) * m02 + at(a, 1, 2) * m01

    at(a, 0, 0) * c0 - at(a, 0, 1) * c1 + at(a, 0, 2) * c2 - at(a, 0, 3) * c3
  end

  # Determinant of the 3x3 obtained by deleting row `skip_r` and column
  # `skip_c`, surviving rows/cols in ascending order, fixed first-row expansion.
  defp minor3(a, skip_r, skip_c) do
    [r0, r1, r2] = Enum.reject(0..3, &(&1 == skip_r))
    [c0, c1, c2] = Enum.reject(0..3, &(&1 == skip_c))

    b00 = at(a, r0, c0)
    b01 = at(a, r0, c1)
    b02 = at(a, r0, c2)
    b10 = at(a, r1, c0)
    b11 = at(a, r1, c1)
    b12 = at(a, r1, c2)
    b20 = at(a, r2, c0)
    b21 = at(a, r2, c1)
    b22 = at(a, r2, c2)

    b00 * (b11 * b22 - b12 * b21) -
      b01 * (b10 * b22 - b12 * b20) +
      b02 * (b10 * b21 - b11 * b20)
  end

  @doc """
  Explicit 4x4 cofactor (adjugate / determinant) inverse.

  `a` is a 4x4 matrix as a tuple of four 4-tuples. Returns `{:ok, inverse}` (same
  shape) or `:singular` when the determinant is exactly zero. The `(i, j)`
  cofactor is `(-1)^(i+j)` times the `(i, j)` minor; the inverse is the transpose
  of the cofactor matrix over the determinant, so `inv[j][i] = cofactor(i, j) / det`.
  """
  @spec inv4(tuple()) :: {:ok, tuple()} | :singular
  def inv4(a) do
    det = det4(a)

    if det == 0.0 do
      :singular
    else
      # Build inv with inv[j][i] = sign(i,j) * minor3(a,i,j) / det.
      cells =
        for j <- 0..3 do
          for i <- 0..3 do
            sign = if rem(i + j, 2) == 0, do: 1.0, else: -1.0
            sign * minor3(a, i, j) / det
          end
          |> List.to_tuple()
        end
        |> List.to_tuple()

      {:ok, cells}
    end
  end

  # ECEF -> ENU rotation at geodetic latitude/longitude (radians).
  defp ecef_to_enu_rotation(lat_rad, lon_rad) do
    sphi = :math.sin(lat_rad)
    cphi = :math.cos(lat_rad)
    slam = :math.sin(lon_rad)
    clam = :math.cos(lon_rad)

    {
      {-slam, clam, 0.0},
      {-sphi * clam, -sphi * slam, cphi},
      {cphi * clam, cphi * slam, sphi}
    }
  end

  # Rotate the 3x3 position cofactor block: Q_enu = R Q_xyz R^T, formed as
  # (R Q) R^T with plain left-to-right inner sums in a fixed order.
  defp rotate_pos_block(q, r) do
    rq =
      for i <- 0..2 do
        for j <- 0..2 do
          Enum.reduce(0..2, 0.0, fn k, acc -> acc + at(r, i, k) * at(q, k, j) end)
        end
        |> List.to_tuple()
      end
      |> List.to_tuple()

    for i <- 0..2 do
      for j <- 0..2 do
        Enum.reduce(0..2, 0.0, fn k, acc -> acc + at(rq, i, k) * at(r, j, k) end)
      end
      |> List.to_tuple()
    end
    |> List.to_tuple()
  end

  # --- receiver geodetic ----------------------------------------------------

  # Receiver geodetic latitude/longitude in radians, from ECEF metres. Mirrors
  # the observables module: divide metres to km for the WGS84 inverse, which
  # returns degrees, then convert to radians.
  defp receiver_lat_lon_rad({x, y, z}) do
    geo = Orbis.Coordinates.to_geodetic({x / 1000.0, y / 1000.0, z / 1000.0})
    {geo.latitude * @deg_to_rad, geo.longitude * @deg_to_rad}
  end

  # --- helpers --------------------------------------------------------------

  defp system_allowed?(_sat_id, nil), do: true

  defp system_allowed?(<<letter::binary-size(1), _rest::binary>>, systems), do: letter in systems

  defp system_allowed?(_sat_id, _systems), do: false

  defp normalize_receiver({x, y, z}) when is_number(x) and is_number(y) and is_number(z),
    do: {:ok, {x * 1.0, y * 1.0, z * 1.0}}

  defp normalize_receiver(%{x_m: x, y_m: y, z_m: z})
       when is_number(x) and is_number(y) and is_number(z), do: {:ok, {x * 1.0, y * 1.0, z * 1.0}}

  defp normalize_receiver(_), do: {:error, :invalid_receiver}
end
