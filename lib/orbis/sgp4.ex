defmodule Orbis.SGP4 do
  @moduledoc """
  SGP4/SDP4 orbit propagation from Two-Line Element sets.
  """

  alias Orbis.Elements
  alias Orbis.TemeState

  @doc """
  Propagate a TLE to a specific datetime, returning a TEME state vector.

  Uses the sgp4 Rust crate in AFSPC compatibility mode. Elements are
  passed as individual fields, so this works for both TLE and OMM inputs.

  Returns `{:ok, %Orbis.TemeState{}}` with position in km and velocity in km/s,
  or `{:error, reason}`.
  """
  @spec propagate(Elements.t(), DateTime.t()) :: {:ok, TemeState.t()} | {:error, String.t()}
  def propagate(%Elements{} = tle, %DateTime{} = datetime) do
    datetime_tuple =
      {{datetime.year, datetime.month, datetime.day},
       {datetime.hour, datetime.minute, datetime.second, elem(datetime.microsecond, 0)}}

    with {:ok, elements_map} <- tle_to_elements_map(tle),
         {:ok, {position, velocity}} <-
           Orbis.NIF.propagate_with_elements(elements_map, datetime_tuple) do
      {:ok, %TemeState{position: position, velocity: velocity}}
    end
  rescue
    e -> {:error, "propagation failed: #{Exception.message(e)}"}
  end

  defp tle_to_elements_map(%Elements{epoch: nil}), do: {:error, "missing epoch"}
  defp tle_to_elements_map(%Elements{catalog_number: nil}), do: {:error, "missing catalog_number"}

  defp tle_to_elements_map(%Elements{} = tle) do
    year = tle.epoch.year
    jan1 = DateTime.new!(Date.new!(year, 1, 1), Time.new!(0, 0, 0, 0), "Etc/UTC")
    diff_us = DateTime.diff(tle.epoch, jan1, :microsecond)
    epochdays = 1.0 + diff_us / (86_400 * 1_000_000)
    epochyr = rem(year, 100)

    {:ok,
     %{
       catalog_number: String.trim(tle.catalog_number || "0"),
       bstar: tle.bstar || 0.0,
       mean_motion_dot: tle.mean_motion_dot || 0.0,
       mean_motion_double_dot: tle.mean_motion_double_dot || 0.0,
       eccentricity: tle.eccentricity || 0.0,
       arg_perigee_deg: tle.arg_perigee_deg || 0.0,
       inclination_deg: tle.inclination_deg || 0.0,
       mean_anomaly_deg: tle.mean_anomaly_deg || 0.0,
       mean_motion: tle.mean_motion || 0.0,
       raan_deg: tle.raan_deg || 0.0,
       epochyr: epochyr,
       epochdays: epochdays
     }}
  end
end
