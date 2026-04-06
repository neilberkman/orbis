defmodule Orbis.Format.OMM do
  @moduledoc """
  Parse and encode CCSDS Orbit Mean-Elements Messages (OMM).

  OMM is the modern standard format for orbital data, carrying the same
  elements as TLE but in structured JSON/XML. Used by CelesTrak and
  Space-Track.

  ## Examples

      {:ok, elements} = Orbis.Format.OMM.parse(json_map)
      json_map = Orbis.Format.OMM.encode(elements)
  """

  alias Orbis.Elements

  @doc """
  Parse an OMM JSON map into an `%Orbis.Elements{}` struct.

  Accepts maps with CelesTrak/Space-Track field names (e.g., `"NORAD_CAT_ID"`,
  `"INCLINATION"`, `"MEAN_MOTION"`). Handles both numeric and string values
  for numeric fields (Space-Track quirk).

  ## Examples

      iex> {:ok, el} = Orbis.Format.OMM.parse(%{
      ...>   "NORAD_CAT_ID" => 25544,
      ...>   "OBJECT_NAME" => "ISS (ZARYA)",
      ...>   "EPOCH" => "2024-01-01T00:00:00",
      ...>   "INCLINATION" => 51.6,
      ...>   "RA_OF_ASC_NODE" => 300.0,
      ...>   "ECCENTRICITY" => 0.0007,
      ...>   "ARG_OF_PERICENTER" => 90.0,
      ...>   "MEAN_ANOMALY" => 270.0,
      ...>   "MEAN_MOTION" => 15.5
      ...> })
      iex> el.catalog_number
      "25544"
      iex> el.object_name
      "ISS (ZARYA)"

  """
  @spec parse(map()) :: {:ok, Elements.t()} | {:error, String.t()}
  def parse(omm) when is_map(omm) do
    with {:ok, epoch} <- parse_epoch(omm["EPOCH"]),
         {:ok, ndot} <- to_float_field(omm, "MEAN_MOTION_DOT"),
         {:ok, nddot} <- to_float_field(omm, "MEAN_MOTION_DDOT"),
         {:ok, bstar} <- to_float_field(omm, "BSTAR") do
      {:ok,
       %Elements{
         object_name: omm["OBJECT_NAME"],
         catalog_number: to_string(omm["NORAD_CAT_ID"]),
         classification: omm["CLASSIFICATION_TYPE"] || "U",
         international_designator: omm["OBJECT_ID"] || "",
         epoch: epoch,
         mean_motion_dot: ndot,
         mean_motion_double_dot: nddot,
         bstar: bstar,
         ephemeris_type: omm["EPHEMERIS_TYPE"] || 0,
         elset_number: omm["ELEMENT_SET_NO"] || 999,
         inclination_deg: to_float(omm["INCLINATION"]),
         raan_deg: to_float(omm["RA_OF_ASC_NODE"]),
         eccentricity: to_float(omm["ECCENTRICITY"]),
         arg_perigee_deg: to_float(omm["ARG_OF_PERICENTER"]),
         mean_anomaly_deg: to_float(omm["MEAN_ANOMALY"]),
         mean_motion: to_float(omm["MEAN_MOTION"]),
         rev_number: omm["REV_AT_EPOCH"] || 0
       }}
    end
  rescue
    e -> {:error, "OMM parse error: #{Exception.message(e)}"}
  end

  @doc """
  Encode an `%Orbis.Elements{}` struct as an OMM JSON-compatible map.

  Returns a map with standard CCSDS OMM field names. Lossless for all
  f64 values.
  """
  @spec encode(Elements.t()) :: map()
  def encode(%Elements{} = el) do
    %{
      "OBJECT_NAME" => el.object_name,
      "OBJECT_ID" => el.international_designator,
      "NORAD_CAT_ID" => safe_int(el.catalog_number),
      "CLASSIFICATION_TYPE" => el.classification,
      "EPOCH" => DateTime.to_iso8601(el.epoch),
      "MEAN_MOTION_DOT" => el.mean_motion_dot,
      "MEAN_MOTION_DDOT" => el.mean_motion_double_dot,
      "BSTAR" => el.bstar,
      "EPHEMERIS_TYPE" => el.ephemeris_type,
      "ELEMENT_SET_NO" => el.elset_number,
      "INCLINATION" => el.inclination_deg,
      "RA_OF_ASC_NODE" => el.raan_deg,
      "ECCENTRICITY" => el.eccentricity,
      "ARG_OF_PERICENTER" => el.arg_perigee_deg,
      "MEAN_ANOMALY" => el.mean_anomaly_deg,
      "MEAN_MOTION" => el.mean_motion,
      "REV_AT_EPOCH" => el.rev_number
    }
  end

  # -- Private --

  defp parse_epoch(nil), do: {:error, "missing EPOCH"}

  defp parse_epoch(epoch_str) when is_binary(epoch_str) do
    # Try parsing as-is first (handles "...Z" and "...+05:00" correctly)
    case DateTime.from_iso8601(epoch_str) do
      {:ok, dt, _offset} ->
        # Shift to UTC if offset was present
        {:ok, DateTime.shift_zone!(dt, "Etc/UTC")}

      {:error, _} ->
        # No timezone info — treat as UTC (CelesTrak convention)
        case NaiveDateTime.from_iso8601(epoch_str) do
          {:ok, ndt} -> {:ok, DateTime.from_naive!(ndt, "Etc/UTC")}
          {:error, reason} -> {:error, "invalid epoch: #{reason}"}
        end
    end
  end

  defp to_float_field(omm, key) do
    case omm[key] do
      nil -> {:ok, 0.0}
      v when is_float(v) -> {:ok, v}
      v when is_integer(v) -> {:ok, v * 1.0}
      v when is_binary(v) -> {:ok, String.to_float(v)}
    end
  end

  defp to_float(v) when is_float(v), do: v
  defp to_float(v) when is_integer(v), do: v * 1.0
  defp to_float(v) when is_binary(v), do: String.to_float(v)
  defp to_float(nil), do: 0.0

  defp safe_int(s) when is_binary(s), do: s |> String.trim() |> String.to_integer()
  defp safe_int(n) when is_integer(n), do: n
end
