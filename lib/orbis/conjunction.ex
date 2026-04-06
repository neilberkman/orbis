defmodule Orbis.Conjunction do
  @moduledoc """
  Find close approaches between two satellites.

  Uses coarse-fine search: scans at a configurable step size, then
  refines with golden-section search within each candidate interval.

  ## Examples

      {:ok, el1} = Orbis.Format.TLE.parse(line1a, line2a)
      {:ok, el2} = Orbis.Format.TLE.parse(line1b, line2b)

      approaches = Orbis.Conjunction.find(el1, el2,
        start_min: 0.0,
        end_min: 1440.0,
        step_min: 1.0,
        threshold_km: 50.0
      )

      for {tca_min, distance_km} <- approaches do
        IO.puts("TCA: +\#{Float.round(tca_min / 60, 1)}h, miss: \#{Float.round(distance_km, 2)} km")
      end
  """

  alias Orbis.Elements

  @doc """
  Find all close approaches between two satellites within a time window.

  Times are in minutes from the first satellite's epoch.

  ## Options

    * `:start_min` - start of search window in minutes (default: 0.0)
    * `:end_min` - end of search window in minutes (required)
    * `:step_min` - coarse scan step size in minutes (default: 1.0)
    * `:threshold_km` - only report approaches closer than this (default: 50.0)

  ## Returns

  List of `{tca_min, distance_km}` tuples, sorted by time.
  """
  @spec find(Elements.t(), Elements.t(), keyword()) ::
          [{float(), float()}] | {:error, term()}
  def find(%Elements{} = el1, %Elements{} = el2, opts) do
    end_min = Keyword.fetch!(opts, :end_min)
    start_min = Keyword.get(opts, :start_min, 0.0)
    step_min = Keyword.get(opts, :step_min, 1.0)
    threshold_km = Keyword.get(opts, :threshold_km, 50.0)

    {l1a, l2a} = Orbis.Format.TLE.encode(el1)
    {l1b, l2b} = Orbis.Format.TLE.encode(el2)

    case Orbis.NIF.find_conjunctions(
           l1a,
           l2a,
           l1b,
           l2b,
           start_min,
           end_min,
           step_min,
           threshold_km
         ) do
      {:ok, results} -> results
      {:error, reason} -> {:error, reason}
    end
  end
end
