defmodule Orbis.GNSS.Data.FtpTransportTest do
  # Not async: it stops/starts the global `:ftp` application.
  use ExUnit.Case, async: false

  alias Orbis.GNSS.Data.Download

  defp ftp_started? do
    Application.started_applications() |> Enum.any?(fn {app, _, _} -> app == :ftp end)
  end

  test "an FTP fetch starts the :ftp transport even if it was not pre-started" do
    # Simulate a host app (escript / bare script / release) that never started
    # the `:ftp` application: previously this crashed with `no process: :ftp_sup`.
    Application.stop(:ftp)
    refute ftp_started?(), "precondition: :ftp must be stopped"

    # Reserved `.invalid` TLD: the open fails fast with no real network, but the
    # transport must already have been started by the fetch path.
    result = Download.get("ftp://nonexistent.invalid/x.gz", :ftp, timeout_ms: 200)

    assert match?({:error, _}, result), "bogus host fetch should error, got #{inspect(result)}"
    assert ftp_started?(), ":ftp must be started by the fetch path"
  after
    # Leave the transport up for any later tests.
    Application.ensure_all_started(:ftp)
  end
end
