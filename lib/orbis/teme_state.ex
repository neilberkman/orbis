defmodule Orbis.TemeState do
  @moduledoc """
  Represents the position and velocity of a satellite in the TEME coordinate system.

  TEME (True Equator Mean Equinox) is the reference frame used by SGP4 propagators.

  ## Units
  - Position: kilometers (km) from Earth's center  
  - Velocity: kilometers per second (km/s)
  """

  @type t :: %__MODULE__{
          position: {float(), float(), float()},
          velocity: {float(), float(), float()}
        }

  @derive Jason.Encoder
  @derive JSON.Encoder
  defstruct [
    :position,
    :velocity
  ]
end
