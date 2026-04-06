defmodule Orbis.Pass do
  @moduledoc """
  A satellite pass over a ground station.
  """
  @derive Jason.Encoder
  @derive JSON.Encoder
  defstruct [:rise, :set, :max_elevation, :max_elevation_time]

  @type t :: %__MODULE__{
          rise: DateTime.t(),
          set: DateTime.t(),
          max_elevation: float(),
          max_elevation_time: DateTime.t()
        }
end
