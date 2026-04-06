defmodule Orbis.ConstellationTest do
  use ExUnit.Case

  @fixtures_dir Path.join(__DIR__, "fixtures/celestrak")

  setup do
    body = File.read!(Path.join(@fixtures_dir, "stations.tle"))

    tles =
      body
      |> String.trim()
      |> String.split("\n")
      |> Enum.map(&String.trim/1)
      |> Enum.reject(&(&1 == ""))
      |> Enum.chunk_every(2, 1, :discard)
      |> Enum.reduce([], fn
        ["1 " <> _ = l1, "2 " <> _ = l2], acc ->
          case Orbis.Format.TLE.parse(l1, l2) do
            {:ok, tle} -> [tle | acc]
            _ -> acc
          end

        _, acc ->
          acc
      end)
      |> Enum.reverse()

    constellation = Orbis.Constellation.from_tles("stations", tles)
    %{constellation: constellation, epoch: hd(tles).epoch}
  end

  test "from_tles creates constellation", %{constellation: c} do
    assert c.name == "stations"
    assert c.count > 10
  end

  test "propagate_all propagates all satellites", %{constellation: c, epoch: dt} do
    results = Orbis.Constellation.propagate_all(c, dt)
    assert length(results) == c.count

    ok_count = Enum.count(results, fn {_, r} -> match?({:ok, _}, r) end)
    assert ok_count == c.count

    # All positions should be in valid orbit range
    for {_id, {:ok, teme}} <- results do
      {x, y, z} = teme.position
      radius = :math.sqrt(x * x + y * y + z * z)
      assert radius > 6300 and radius < 50_000
    end
  end
end
