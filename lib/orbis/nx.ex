defmodule Orbis.Nx do
  @moduledoc """
  Batch/tensor analysis helpers for Orbis.

  This layer is for high-throughput workflows like visibility matrices,
  coverage grids, and Monte Carlo studies. It complements the exact scalar
  APIs in `Orbis`, rather than replacing them.
  """

  alias Orbis.Nx.{Coverage, Geometry, RF, Visibility}

  @doc """
  Compute topocentric look angles for many ITRS positions and stations.

  Expected shapes:
  - `sat_positions`: `[n, 3]` in ITRS km
  - `stations`: `[m, 3]` as `{lat_deg, lon_deg, alt_m}`

  Returns tensors shaped `[n, m]`.
  """
  defdelegate look_angles(sat_positions, stations, opts \\ []), to: Geometry

  @doc """
  Return a boolean visibility mask for `min_elevation` degrees.

  Result shape: `[n, m]`.
  """
  defdelegate visible_mask(sat_positions, stations, opts \\ []), to: Visibility

  @doc """
  Batch free-space path loss.

  `range_km` may be any broadcastable tensor.
  """
  defdelegate fspl(range_km, frequency_mhz), to: RF

  @doc """
  Batch link-margin calculation with broadcastable inputs.
  """
  defdelegate link_margin(params), to: RF

  @doc """
  Compute simple access counts over a time series.

  Expected shape for `elevation_series`: `[t, s, g]`.
  """
  defdelegate access_counts(elevation_series, opts \\ []), to: Coverage
end
