defmodule Orbis.Screening.Result do
  @moduledoc """
  Final result for a screened candidate pair.
  """

  @type t :: %__MODULE__{
          candidate: Orbis.Screening.Candidate.t(),
          collision: Orbis.Collision.Result.t() | nil,
          error: String.t() | nil
        }

  defstruct [
    :candidate,
    :collision,
    :error
  ]
end
