defmodule Orbis.GNSS.BroadcastComparison do
  @moduledoc """
  Broadcast-ephemeris accuracy: compare a broadcast navigation product against a
  precise SP3 product over a window.

  This is the standard broadcast-orbit / clock accuracy check (the orbit and clock
  pieces of the signal-in-space range error, SISRE). For each satellite at each
  epoch it differences the broadcast-evaluated ECEF position and clock against the
  precise SP3 values, decomposes the position error into radial / along-track /
  cross-track (RAC) components, and summarizes the differences as RMS and maximum
  statistics per satellite and overall.

  Both products are evaluated through `Orbis.GNSS.Ephemeris`, so the frame
  (ITRF/IGS ECEF, meters), the time scale (GPST), and the clock sign convention
  (positive = satellite clock ahead of system time) are exactly as documented
  there. Only epochs where **both** sources return a valid state contribute to the
  statistics; an epoch missing from either product is skipped, never extrapolated.

  ## RAC frame

  The radial/along-track/cross-track unit vectors are built from the **precise**
  state: radial along the position vector, cross-track along the orbital angular
  momentum `r x v`, and along-track completing the right-handed triad. The SP3
  position/clock interpolation does not expose a velocity, so the velocity is
  derived by a centered finite difference of the precise position
  (`(r(t+dt) - r(t-dt)) / 2dt`, falling back to a one-sided difference at a window
  edge). The position-difference vector `broadcast - precise` is projected onto
  this triad.

  ## Expected magnitudes

  GPS LNAV broadcast orbits differ from IGS precise orbits at roughly the 1-2 m
  RMS level (3D), dominated by the along-track and radial components; Galileo and
  BeiDou MEO are comparable. A result far outside this band (tens of meters or
  more) indicates a parse or evaluation defect rather than normal broadcast error.

  The broadcast models follow IS-GPS-200 (GPS LNAV), the Galileo OS-SIS-ICD, and
  the BeiDou BDS-SIS-ICD; the precise product is SP3-c / SP3-d (IGS).

  ## Example

      {:ok, broadcast} = Orbis.GNSS.Broadcast.load("BRDC.rnx")
      {:ok, sp3} = Orbis.GNSS.SP3.load("igs.sp3")

      report =
        Orbis.GNSS.BroadcastComparison.compare(broadcast, sp3, ["G01", "E11"], %{
          from: ~N[2020-06-25 02:00:00],
          to: ~N[2020-06-25 04:00:00],
          step_s: 300
        })

      report.overall.orbit_3d_rms_m   # 3D orbit RMS over all satellites, meters
      report.per_satellite["G01"].radial_rms_m
  """

  alias Orbis.GNSS.Broadcast
  alias Orbis.GNSS.Ephemeris
  alias Orbis.GNSS.SP3

  defmodule Stats do
    @moduledoc """
    Orbit and clock difference statistics for one satellite (or the overall set).

    All values are meters except `count` (the number of compared epochs).
    `orbit_3d_rms_m` / `orbit_3d_max_m` are the Euclidean position-difference
    magnitudes. `radial_*`, `along_*`, `cross_*` are the RMS and max of the signed
    RAC components of the position difference (`broadcast - precise`).
    `clock_rms_m` / `clock_max_m` are the satellite-clock differences scaled to
    meters by the speed of light; they are `nil` when neither product carried a
    clock estimate for any compared epoch.
    """
    @enforce_keys [
      :count,
      :orbit_3d_rms_m,
      :orbit_3d_max_m,
      :radial_rms_m,
      :radial_max_m,
      :along_rms_m,
      :along_max_m,
      :cross_rms_m,
      :cross_max_m,
      :clock_rms_m,
      :clock_max_m
    ]
    defstruct [
      :count,
      :orbit_3d_rms_m,
      :orbit_3d_max_m,
      :radial_rms_m,
      :radial_max_m,
      :along_rms_m,
      :along_max_m,
      :cross_rms_m,
      :cross_max_m,
      :clock_rms_m,
      :clock_max_m
    ]

    @type t :: %__MODULE__{
            count: non_neg_integer(),
            orbit_3d_rms_m: float() | nil,
            orbit_3d_max_m: float() | nil,
            radial_rms_m: float() | nil,
            radial_max_m: float() | nil,
            along_rms_m: float() | nil,
            along_max_m: float() | nil,
            cross_rms_m: float() | nil,
            cross_max_m: float() | nil,
            clock_rms_m: float() | nil,
            clock_max_m: float() | nil
          }
  end

  defmodule Report do
    @moduledoc """
    The result of a broadcast-vs-precise comparison.

    `per_satellite` maps each satellite id to its `Stats`; `overall` aggregates
    every compared epoch across all satellites. `missing` lists `{satellite_id,
    count}` pairs counting epochs that were skipped because one or both products
    had no valid state there.
    """
    @enforce_keys [:per_satellite, :overall, :missing]
    defstruct [:per_satellite, :overall, :missing]

    @type t :: %__MODULE__{
            per_satellite: %{String.t() => Stats.t()},
            overall: Stats.t(),
            missing: [{String.t(), non_neg_integer()}]
          }
  end

  # Speed of light (m/s), CODATA / IS-GPS-200 defining constant, to scale a
  # clock-offset difference (seconds) into a range error (meters).
  @speed_of_light_m_s 299_792_458.0

  @doc """
  Compare a broadcast product against a precise SP3 product over `window`.

  `broadcast` is an `Orbis.GNSS.Broadcast` handle and `precise` an
  `Orbis.GNSS.SP3` handle for the same day; `sat_ids` is a list of canonical
  RINEX tokens; `window` is the `Orbis.GNSS.Ephemeris` window map (`:from`,
  `:to`, `:step_s`). Returns an `Orbis.GNSS.BroadcastComparison.Report`.

  The window should sit within both the SP3 file span and the broadcast records'
  fit intervals; epochs missing from either product are counted in
  `report.missing` and excluded from the statistics.
  """
  @spec compare(Broadcast.t(), SP3.t(), [String.t()], Ephemeris.window()) :: Report.t()
  def compare(%Broadcast{} = broadcast, %SP3{} = precise, sat_ids, window)
      when is_list(sat_ids) do
    per_sat =
      Map.new(sat_ids, fn sat_id ->
        {sat_id, compare_satellite(broadcast, precise, sat_id, window)}
      end)

    overall = aggregate(Enum.flat_map(per_sat, fn {_id, {diffs, _missing}} -> diffs end))

    per_satellite = Map.new(per_sat, fn {id, {diffs, _missing}} -> {id, stats(diffs)} end)

    missing =
      per_sat
      |> Enum.map(fn {id, {_diffs, missing}} -> {id, missing} end)
      |> Enum.filter(fn {_id, n} -> n > 0 end)
      |> Enum.sort()

    %Report{per_satellite: per_satellite, overall: overall, missing: missing}
  end

  # Per-epoch differences for one satellite, plus the count of skipped epochs.
  defp compare_satellite(broadcast, precise, sat_id, window) do
    epochs = window_epochs(window)
    dt_s = window.step_s / 2.0

    {diffs, missing} =
      Enum.reduce(epochs, {[], 0}, fn epoch, {acc, miss} ->
        case diff_at(broadcast, precise, sat_id, epoch, dt_s) do
          {:ok, diff} -> {[diff | acc], miss}
          :skip -> {acc, miss + 1}
        end
      end)

    {Enum.reverse(diffs), missing}
  end

  # A single epoch's broadcast-minus-precise difference, decomposed into RAC and a
  # clock difference. Returns `:skip` when either product lacks a valid position
  # at the epoch (or the precise velocity finite difference is unavailable).
  defp diff_at(broadcast, precise, sat_id, epoch, dt_s) do
    with {:ok, b} <- position(broadcast, sat_id, epoch),
         {:ok, p} <- position(precise, sat_id, epoch),
         {:ok, vel} <- precise_velocity(precise, sat_id, epoch, dt_s) do
      d = sub(b.pos, p.pos)
      {radial, along, cross} = project_rac(d, p.pos, vel)

      clock_m =
        case {b.clock_s, p.clock_s} do
          {bc, pc} when is_float(bc) and is_float(pc) -> (bc - pc) * @speed_of_light_m_s
          _ -> nil
        end

      {:ok,
       %{
         orbit_3d: norm(d),
         radial: radial,
         along: along,
         cross: cross,
         clock_m: clock_m
       }}
    else
      _ -> :skip
    end
  end

  # Pull an ECEF position + clock from either product into a common shape.
  defp position(%Broadcast{} = src, sat_id, epoch) do
    case Broadcast.position(src, sat_id, epoch) do
      {:ok, %Broadcast.State{x_m: x, y_m: y, z_m: z, clock_s: clk}} ->
        {:ok, %{pos: {x, y, z}, clock_s: clk}}

      {:error, _} = err ->
        err
    end
  end

  defp position(%SP3{} = src, sat_id, epoch) do
    case SP3.position(src, sat_id, epoch) do
      {:ok, %SP3.State{x_m: x, y_m: y, z_m: z, clock_s: clk}} ->
        {:ok, %{pos: {x, y, z}, clock_s: clk}}

      {:error, _} = err ->
        err
    end
  end

  # Centered finite-difference velocity of the precise position, falling back to a
  # one-sided difference when a neighbor epoch lies outside the SP3 span. The SP3
  # interpolation API exposes position only, so velocity is derived here.
  defp precise_velocity(precise, sat_id, epoch, dt_s) do
    half = round(dt_s)
    plus = NaiveDateTime.add(epoch, half, :second)
    minus = NaiveDateTime.add(epoch, -half, :second)

    case {sp3_pos(precise, sat_id, plus), sp3_pos(precise, sat_id, minus)} do
      {{:ok, rp}, {:ok, rm}} ->
        {:ok, scale(sub(rp, rm), 1.0 / (2.0 * half))}

      {{:ok, rp}, _} ->
        with {:ok, r0} <- sp3_pos(precise, sat_id, epoch),
             do: {:ok, scale(sub(rp, r0), 1.0 / half)}

      {_, {:ok, rm}} ->
        with {:ok, r0} <- sp3_pos(precise, sat_id, epoch),
             do: {:ok, scale(sub(r0, rm), 1.0 / half)}

      _ ->
        :error
    end
  end

  defp sp3_pos(precise, sat_id, epoch) do
    case SP3.position(precise, sat_id, epoch) do
      {:ok, %SP3.State{x_m: x, y_m: y, z_m: z}} -> {:ok, {x, y, z}}
      {:error, _} -> :error
    end
  end

  # --- RAC projection ------------------------------------------------------

  # Project a difference vector onto the radial/along-track/cross-track triad of
  # the orbit defined by position `r` and velocity `v`. Radial is along `r`,
  # cross-track along the angular momentum `r x v`, along-track completes the
  # right-handed set (cross x radial).
  defp project_rac(d, r, v) do
    radial_hat = unit(r)
    cross_hat = unit(cross(r, v))
    along_hat = cross(cross_hat, radial_hat)

    {dot(d, radial_hat), dot(d, along_hat), dot(d, cross_hat)}
  end

  # --- statistics ----------------------------------------------------------

  defp stats(diffs), do: aggregate(diffs)

  defp aggregate([]) do
    %Stats{
      count: 0,
      orbit_3d_rms_m: nil,
      orbit_3d_max_m: nil,
      radial_rms_m: nil,
      radial_max_m: nil,
      along_rms_m: nil,
      along_max_m: nil,
      cross_rms_m: nil,
      cross_max_m: nil,
      clock_rms_m: nil,
      clock_max_m: nil
    }
  end

  defp aggregate(diffs) do
    clocks = diffs |> Enum.map(& &1.clock_m) |> Enum.reject(&is_nil/1)

    %Stats{
      count: length(diffs),
      orbit_3d_rms_m: rms(Enum.map(diffs, & &1.orbit_3d)),
      orbit_3d_max_m: max_abs(Enum.map(diffs, & &1.orbit_3d)),
      radial_rms_m: rms(Enum.map(diffs, & &1.radial)),
      radial_max_m: max_abs(Enum.map(diffs, & &1.radial)),
      along_rms_m: rms(Enum.map(diffs, & &1.along)),
      along_max_m: max_abs(Enum.map(diffs, & &1.along)),
      cross_rms_m: rms(Enum.map(diffs, & &1.cross)),
      cross_max_m: max_abs(Enum.map(diffs, & &1.cross)),
      clock_rms_m: rms_or_nil(clocks),
      clock_max_m: max_abs_or_nil(clocks)
    }
  end

  defp rms([]), do: nil

  defp rms(values) do
    sum_sq = Enum.reduce(values, 0.0, fn x, acc -> acc + x * x end)
    :math.sqrt(sum_sq / length(values))
  end

  defp rms_or_nil([]), do: nil
  defp rms_or_nil(values), do: rms(values)

  defp max_abs([]), do: nil
  defp max_abs(values), do: values |> Enum.map(&abs/1) |> Enum.max()

  defp max_abs_or_nil([]), do: nil
  defp max_abs_or_nil(values), do: max_abs(values)

  # --- vector helpers ------------------------------------------------------

  defp window_epochs(%{from: from, to: to, step_s: step_s}) do
    total_s = NaiveDateTime.diff(to, from, :second)
    n = div(total_s, step_s)
    for i <- 0..n//1, do: NaiveDateTime.add(from, i * step_s, :second)
  end

  defp sub({ax, ay, az}, {bx, by, bz}), do: {ax - bx, ay - by, az - bz}

  defp scale({x, y, z}, s), do: {x * s, y * s, z * s}

  defp dot({ax, ay, az}, {bx, by, bz}), do: ax * bx + ay * by + az * bz

  defp cross({ax, ay, az}, {bx, by, bz}) do
    {ay * bz - az * by, az * bx - ax * bz, ax * by - ay * bx}
  end

  defp norm(v), do: :math.sqrt(dot(v, v))

  defp unit(v) do
    n = norm(v)
    if n == 0.0, do: {0.0, 0.0, 0.0}, else: scale(v, 1.0 / n)
  end
end
