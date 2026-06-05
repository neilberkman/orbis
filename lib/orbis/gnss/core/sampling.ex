defmodule Orbis.GNSS.Core.Sampling do
  @moduledoc false

  alias Orbis.GNSS.Core.Epoch
  alias Orbis.GNSS.SP3

  def sample_sp3(%SP3{} = sp3, sat_id, t0, t1, cadence_s) do
    steps = Epoch.steps(t0, t1, cadence_s)

    samples =
      steps
      |> Enum.reduce([], fn ep, acc ->
        case SP3.position(sp3, sat_id, ep) do
          {:ok, %{x_m: x, y_m: y, z_m: z}} -> [{ep, {x, y, z}} | acc]
          {:error, _} -> acc
        end
      end)
      |> Enum.reverse()

    {:ok, samples, length(steps)}
  end
end
