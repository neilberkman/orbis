defmodule Orbis.CCSDS.CDM do
  @moduledoc """
  Parse and encode CCSDS Conjunction Data Messages (CDM).

  Supports the KVN (Keyword=Value Notation) format per CCSDS 508.0-B-1.
  CDMs describe a predicted close approach between two space objects,
  including states, covariances, and collision probability.

  ## Examples

      {:ok, cdm} = Orbis.CCSDS.CDM.parse(kvn_string)
      cdm.tca                    # ~U[2010-03-13 22:37:52.618Z]
      cdm.miss_distance_m        # 715.0
      cdm.collision_probability  # 4.835e-05
      cdm.object1.state          # {r, v} in km and km/s
      cdm.object1.covariance     # 6x6 lower-triangle in m²
  """

  defmodule ObjectData do
    @moduledoc false

    @type t :: %__MODULE__{
            object_designator: String.t() | nil,
            catalog_name: String.t() | nil,
            object_name: String.t() | nil,
            international_designator: String.t() | nil,
            object_type: String.t() | nil,
            ref_frame: String.t() | nil,
            state: {{float(), float(), float()}, {float(), float(), float()}} | nil,
            covariance_rtn: list(float()) | nil
          }

    defstruct [
      :object_designator,
      :catalog_name,
      :object_name,
      :international_designator,
      :object_type,
      :ref_frame,
      :state,
      :covariance_rtn
    ]
  end

  @type t :: %__MODULE__{
          creation_date: DateTime.t() | nil,
          originator: String.t() | nil,
          message_id: String.t() | nil,
          tca: DateTime.t() | nil,
          miss_distance_m: float() | nil,
          relative_speed_m_s: float() | nil,
          collision_probability: float() | nil,
          collision_probability_method: String.t() | nil,
          hard_body_radius_m: float() | nil,
          object1: ObjectData.t() | nil,
          object2: ObjectData.t() | nil
        }

  defstruct [
    :creation_date,
    :originator,
    :message_id,
    :tca,
    :miss_distance_m,
    :relative_speed_m_s,
    :collision_probability,
    :collision_probability_method,
    :hard_body_radius_m,
    :object1,
    :object2
  ]

  @doc """
  Parse a CDM in KVN format.

  Returns `{:ok, %CDM{}}` or `{:error, reason}`.
  """
  @spec parse(String.t()) :: {:ok, t()} | {:error, String.t()}
  def parse(kvn_string) when is_binary(kvn_string) do
    lines =
      kvn_string
      |> String.split("\n")
      |> Enum.map(&String.trim/1)
      |> Enum.reject(&(&1 == "" or String.starts_with?(&1, "COMMENT")))

    kv = parse_kv_lines(lines)

    with {:ok, tca} <- parse_datetime(kv["TCA"]),
         {:ok, creation} <- parse_datetime(kv["CREATION_DATE"]) do
      # Parse HBR from COMMENT lines (NASA CARA convention)
      hbr = parse_hbr(kvn_string)

      # Split into object blocks
      {obj1_kv, obj2_kv} = split_object_blocks(lines)

      with {:ok, obj1} <- parse_object(obj1_kv),
           {:ok, obj2} <- parse_object(obj2_kv) do
        {:ok,
         %__MODULE__{
           creation_date: creation,
           originator: kv["ORIGINATOR"],
           message_id: kv["MESSAGE_ID"],
           tca: tca,
           miss_distance_m: parse_num(kv["MISS_DISTANCE"]),
           relative_speed_m_s: parse_num(kv["RELATIVE_SPEED"]),
           collision_probability: parse_num(kv["COLLISION_PROBABILITY"]),
           collision_probability_method: kv["COLLISION_PROBABILITY_METHOD"],
           hard_body_radius_m: hbr,
           object1: obj1,
           object2: obj2
         }}
      end
    end
  rescue
    e -> {:error, "CDM parse error: #{Exception.message(e)}"}
  end

  @doc """
  Convert a parsed CDM to inputs for `Orbis.Collision.probability/1`.

  Extracts states and covariances from the CDM, converting RTN covariance
  to position-only ECI covariance (simplified — uses diagonal approximation
  for the RTN-to-ECI transform).

  Returns the parameter map expected by `Orbis.Collision.probability/1`.
  """
  @spec to_collision_params(t()) :: map()
  def to_collision_params(%__MODULE__{} = cdm) do
    {r1, v1} = cdm.object1.state
    {r2, v2} = cdm.object2.state

    # Extract 3x3 position covariance from RTN (first 3x3 block)
    cov1_rtn = extract_pos_cov(cdm.object1.covariance_rtn)
    cov2_rtn = extract_pos_cov(cdm.object2.covariance_rtn)

    # Convert RTN covariance to ECI using the object's state
    cov1_eci = rtn_to_eci(cov1_rtn, r1, v1)
    cov2_eci = rtn_to_eci(cov2_rtn, r2, v2)

    # Convert m² to km²
    cov1_km2 = scale_cov(cov1_eci, 1.0e-6)
    cov2_km2 = scale_cov(cov2_eci, 1.0e-6)

    hbr_km = (cdm.hard_body_radius_m || 15.0) / 1000.0

    %{
      r1: r1,
      v1: v1,
      cov1: cov1_km2,
      r2: r2,
      v2: v2,
      cov2: cov2_km2,
      hard_body_radius_km: hbr_km
    }
  end

  # --- KVN Parsing ---

  defp parse_kv_lines(lines) do
    lines
    |> Enum.filter(&String.contains?(&1, "="))
    |> Enum.map(fn line ->
      case String.split(line, "=", parts: 2) do
        [k, v] -> {String.trim(k), String.trim(v) |> strip_units()}
        _ -> nil
      end
    end)
    |> Enum.reject(&is_nil/1)
    |> Map.new()
  end

  defp strip_units(val) do
    # Remove [m], [m/s], [km], [m**2] etc.
    Regex.replace(~r/\s*\[.*?\]\s*$/, val, "")
  end

  defp parse_datetime(nil), do: {:error, "missing datetime"}

  defp parse_datetime(str) do
    case DateTime.from_iso8601(str <> "Z") do
      {:ok, dt, _} ->
        {:ok, dt}

      _ ->
        case NaiveDateTime.from_iso8601(str) do
          {:ok, ndt} -> {:ok, DateTime.from_naive!(ndt, "Etc/UTC")}
          {:error, r} -> {:error, "bad datetime: #{r}"}
        end
    end
  end

  defp parse_num(nil), do: nil

  defp parse_num(str) do
    case Float.parse(str) do
      {f, _} -> f
      :error -> nil
    end
  end

  defp parse_hbr(text) do
    case Regex.run(~r/COMMENT\s+HBR\s*=\s*([\d.]+)/i, text) do
      [_, val] -> parse_num(val)
      _ -> nil
    end
  end

  defp split_object_blocks(lines) do
    # Find exact "OBJECT = OBJECT1" and "OBJECT = OBJECT2" markers
    # Must match key exactly to avoid false hits on OBJECT_DESIGNATOR etc.
    obj_marker_indices =
      lines
      |> Enum.with_index()
      |> Enum.filter(fn {line, _idx} ->
        case String.split(line, "=", parts: 2) do
          [k, _v] -> String.trim(k) == "OBJECT"
          _ -> false
        end
      end)
      |> Enum.map(fn {_line, idx} -> idx end)

    case obj_marker_indices do
      [i1, i2 | _] ->
        obj1_lines = Enum.slice(lines, i1, i2 - i1)
        obj2_lines = Enum.slice(lines, i2, length(lines) - i2)
        {parse_kv_lines(obj1_lines), parse_kv_lines(obj2_lines)}

      _ ->
        # Fallback: split in half
        mid = div(length(lines), 2)
        {parse_kv_lines(Enum.take(lines, mid)), parse_kv_lines(Enum.drop(lines, mid))}
    end
  end

  @state_keys ~w(X Y Z X_DOT Y_DOT Z_DOT)
  @cov_keys ~w(CR_R CT_R CT_T CN_R CN_T CN_N)

  defp parse_object(kv) do
    # Require state vector fields
    state_vals = Enum.map(@state_keys, &parse_num(kv[&1]))
    cov_vals = Enum.map(@cov_keys, &parse_num(kv[&1]))

    if Enum.any?(state_vals, &is_nil/1) do
      {:error, "incomplete state vector"}
    else
      [x, y, z, xd, yd, zd] = state_vals

      # Covariance defaults to zero (optional in some CDMs)
      [cr_r, ct_r, ct_t, cn_r, cn_t, cn_n] = Enum.map(cov_vals, &(&1 || 0.0))

      {:ok,
       %ObjectData{
         object_designator: kv["OBJECT_DESIGNATOR"],
         catalog_name: kv["CATALOG_NAME"],
         object_name: kv["OBJECT_NAME"],
         international_designator: kv["INTERNATIONAL_DESIGNATOR"],
         object_type: kv["OBJECT_TYPE"],
         ref_frame: kv["REF_FRAME"],
         state: {{x, y, z}, {xd, yd, zd}},
         covariance_rtn: [cr_r, ct_r, ct_t, cn_r, cn_t, cn_n]
       }}
    end
  end

  # --- Covariance transforms ---

  # Extract 3x3 position covariance from RTN lower-triangle
  defp extract_pos_cov([cr_r, ct_r, ct_t, cn_r, cn_t, cn_n | _]) do
    [
      [cr_r, ct_r, cn_r],
      [ct_r, ct_t, cn_t],
      [cn_r, cn_t, cn_n]
    ]
  end

  # Rotate RTN covariance to ECI: C_eci = R * C_rtn * R'
  defp rtn_to_eci(cov_rtn, {rx, ry, rz}, {vx, vy, vz}) do
    r = {rx, ry, rz}
    v = {vx, vy, vz}

    r_hat = normalize(r)
    h = cross(r, v)
    n_hat = normalize(h)
    t_hat = cross(n_hat, r_hat)

    # Rotation matrix: RTN → ECI (columns are R, T, N unit vectors)
    rot = [
      [elem(r_hat, 0), elem(t_hat, 0), elem(n_hat, 0)],
      [elem(r_hat, 1), elem(t_hat, 1), elem(n_hat, 1)],
      [elem(r_hat, 2), elem(t_hat, 2), elem(n_hat, 2)]
    ]

    rot_t = transpose(rot)
    mat_mul(mat_mul(rot, cov_rtn), rot_t)
  end

  defp mag({x, y, z}), do: :math.sqrt(x * x + y * y + z * z)
  defp normalize(v), do: {elem(v, 0) / mag(v), elem(v, 1) / mag(v), elem(v, 2) / mag(v)}

  defp cross({ax, ay, az}, {bx, by, bz}),
    do: {ay * bz - az * by, az * bx - ax * bz, ax * by - ay * bx}

  defp transpose(m), do: m |> Enum.zip() |> Enum.map(&Tuple.to_list/1)

  defp mat_mul(a, b) do
    bt = transpose(b)

    for row <- a,
        do:
          for(col <- bt, do: Enum.zip(row, col) |> Enum.map(fn {x, y} -> x * y end) |> Enum.sum())
  end

  defp scale_cov(m, s), do: for(row <- m, do: for(v <- row, do: v * s))
end
