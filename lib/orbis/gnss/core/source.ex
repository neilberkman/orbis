defmodule Orbis.GNSS.Core.Source do
  @moduledoc false

  alias Orbis.GNSS.SP3

  def same_time_scale(%{time_scale: model_scale}, %SP3{time_scale: sp3_scale}) do
    if model_scale == sp3_scale,
      do: :ok,
      else: {:error, {:time_scale_mismatch, model_scale, sp3_scale}}
  end
end
