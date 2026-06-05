defmodule Orbis.GNSS.Core.Types do
  @moduledoc false

  @type ecef :: {float(), float(), float()}
  @type ecef_input ::
          {number(), number(), number()} | %{x_m: number(), y_m: number(), z_m: number()}

  def normalize_ecef(value, error_tag \\ :invalid_receiver)

  def normalize_ecef({x, y, z}, _error_tag) when is_number(x) and is_number(y) and is_number(z),
    do: {:ok, {x * 1.0, y * 1.0, z * 1.0}}

  def normalize_ecef(%{x_m: x, y_m: y, z_m: z}, _error_tag)
      when is_number(x) and is_number(y) and is_number(z), do: {:ok, {x * 1.0, y * 1.0, z * 1.0}}

  def normalize_ecef(_value, error_tag), do: {:error, error_tag}

  def system_allowed?(_sat_id, nil), do: true
  def system_allowed?(<<letter::binary-size(1), _rest::binary>>, systems), do: letter in systems
  def system_allowed?(_sat_id, _systems), do: false
end
