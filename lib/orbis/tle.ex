defmodule Orbis.TLE do
  @moduledoc """
  Backwards-compatible delegates for TLE/OMM parsing.

  Prefer `Orbis.Format.TLE` and `Orbis.Format.OMM` directly.
  """

  defdelegate parse(line1, line2), to: Orbis.Format.TLE
  defdelegate from_omm(omm), to: Orbis.Format.OMM, as: :parse
  defdelegate to_omm(elements), to: Orbis.Format.OMM, as: :encode
end
