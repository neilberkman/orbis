defmodule Orbis.Constellation do
  @moduledoc """
  Manage and propagate satellite constellations.

  Load TLEs for an entire constellation and propagate all satellites
  to a given time. Useful for coverage analysis, constellation status
  monitoring, and visibility computations.

  ## Examples

      # Load from CelesTrak
      {:ok, constellation} = Orbis.Constellation.load("globalstar")
      constellation.count
      #=> 85

      # Propagate all satellites to now
      positions = Orbis.Constellation.propagate_all(constellation, DateTime.utc_now())
      Enum.each(positions, fn {norad_id, pos} ->
        IO.puts("\#{norad_id}: \#{inspect(pos)}")
      end)

      # Find visible satellites from a ground station
      visible = Orbis.Constellation.visible_from(constellation, station, datetime)
  """

  defstruct [:name, :satellites, :count]

  @type t :: %__MODULE__{
          name: String.t(),
          satellites: [Orbis.Elements.t()],
          count: non_neg_integer()
        }

  @doc """
  Load a constellation from CelesTrak by group name.

  ## Examples

      {:ok, c} = Orbis.Constellation.load("globalstar")
      c.count
      #=> 85
  """
  @spec load(String.t()) :: {:ok, t()} | {:error, term()}
  def load(group_name) do
    case Orbis.CelesTrak.fetch_group(group_name) do
      {:ok, tles} ->
        {:ok,
         %__MODULE__{
           name: group_name,
           satellites: tles,
           count: length(tles)
         }}

      {:error, reason} ->
        {:error, reason}
    end
  end

  @doc """
  Create a constellation from a list of TLEs.

  ## Examples

      constellation = Orbis.Constellation.from_tles("custom", tles)
  """
  @spec from_tles(String.t(), [Orbis.Elements.t()]) :: t()
  def from_tles(name, tles) do
    %__MODULE__{
      name: name,
      satellites: tles,
      count: length(tles)
    }
  end

  @doc """
  Propagate all satellites to a given time.

  Returns a list of `{catalog_number, {:ok, teme_state}}` or
  `{catalog_number, {:error, reason}}` tuples.

  ## Examples

      results = Orbis.Constellation.propagate_all(constellation, ~U[2024-07-04 00:00:00Z])
      for {id, {:ok, teme}} <- results do
        IO.puts("\#{id}: \#{inspect(teme.position)}")
      end
  """
  @spec propagate_all(t(), DateTime.t()) :: [
          {String.t(), {:ok, Orbis.TemeState.t()} | {:error, term()}}
        ]
  def propagate_all(%__MODULE__{satellites: sats}, datetime) do
    sats
    |> Task.async_stream(
      fn tle ->
        {tle.catalog_number, Orbis.SGP4.propagate(tle, datetime)}
      end,
      max_concurrency: System.schedulers_online() * 2,
      timeout: 5_000
    )
    |> Enum.flat_map(fn
      {:ok, result} -> [result]
      {:exit, _reason} -> []
    end)
  end

  @doc """
  Find satellites visible from a ground station at a given time.

  Returns satellites above `min_elevation` degrees, sorted by elevation
  (highest first).

  ## Options

    * `:min_elevation` - minimum elevation in degrees (default: 10.0)

  ## Examples

      visible = Orbis.Constellation.visible_from(constellation, station, datetime)
      for sat <- visible do
        IO.puts("\#{sat.catalog_number}: el=\#{sat.elevation}° range=\#{sat.range_km} km")
      end
  """
  @spec visible_from(t(), map(), DateTime.t(), keyword()) :: [map()]
  def visible_from(%__MODULE__{} = constellation, station, datetime, opts \\ []) do
    min_el = Keyword.get(opts, :min_elevation, 10.0)

    constellation
    |> propagate_all(datetime)
    |> Enum.filter(fn {_id, result} -> match?({:ok, _}, result) end)
    |> Enum.map(fn {id, {:ok, teme}} ->
      look =
        Orbis.Coordinates.to_topocentric(
          Orbis.Coordinates.teme_to_gcrs(teme, datetime),
          datetime,
          station
        )

      %{
        catalog_number: id,
        elevation: look.elevation,
        azimuth: look.azimuth,
        range_km: look.range_km,
        position: teme.position
      }
    end)
    |> Enum.filter(&(&1.elevation >= min_el))
    |> Enum.sort_by(& &1.elevation, :desc)
  end
end
