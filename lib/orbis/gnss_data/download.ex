defmodule Orbis.GnssData.Download do
  @moduledoc """
  Network transfer for GNSS products: HTTPS via `Req` and FTP via the Erlang
  stdlib `:ftp` application, with bounded retries, sensible timeouts, and a
  strict host allow-list.

  This module is only reached when `offline: true` is **not** set and the
  product is not already cached. It never accepts a caller-supplied URL: it is
  given a URL built by `Orbis.GnssData.Catalog`, and it re-checks the URL's host
  against `Catalog.allowed_hosts/0` before contacting it (defence in depth
  against SSRF). Cross-host HTTP redirects are **refused** rather than followed,
  so a compromised or on-path archive cannot redirect the fetch off the
  allow-list. Both transports also cap the number of compressed bytes pulled
  into memory (`:max_compressed_bytes`), since the remote file is untrusted (FTP
  is plaintext, anonymous). All failures are mapped to typed errors.
  """

  alias Orbis.GnssData.Catalog

  @default_timeout_ms 30_000
  @default_retries 3
  @default_backoff_ms 500

  # Generous cap on the compressed payload we will buffer. Real daily GNSS
  # products are a few MiB; this bounds memory against a hostile/oversized
  # response before the (output-side) decompression cap is even reached.
  @default_max_compressed_bytes 64 * 1024 * 1024

  @doc """
  Download the bytes at `url` using `protocol` (`:https` or `:ftp`).

  Options:

    * `:timeout_ms` — per-attempt timeout (default `30_000`)
    * `:retries` — number of attempts for transient errors (default `3`)
    * `:backoff_ms` — base backoff between retries, doubled each attempt
      (default `500`); set to `0` in tests to avoid sleeping
    * `:max_compressed_bytes` — cap on the compressed payload buffered into
      memory (default 64 MiB); exceeding it yields
      `{:error, {:download_size_exceeded, max, got}}`

  Returns `{:ok, compressed_bytes}` or a typed error. A `404`/file-not-found is
  treated as permanent and is **not** retried; transient network errors
  (including `408`/`429`) are.
  """
  @spec get(String.t(), :https | :ftp, keyword()) :: {:ok, binary()} | {:error, term()}
  def get(url, protocol, opts \\ []) when is_binary(url) and protocol in [:https, :ftp] do
    with :ok <- check_host(url) do
      retries = Keyword.get(opts, :retries, @default_retries)
      backoff = Keyword.get(opts, :backoff_ms, @default_backoff_ms)
      attempt(url, protocol, opts, retries, backoff)
    end
  end

  @doc """
  Whether `Req` is loaded and usable for HTTPS at runtime.
  """
  @spec req_available?() :: boolean()
  def req_available? do
    case Application.get_env(:orbis, :gnss_data_req_available) do
      nil -> Code.ensure_loaded?(Req) and function_exported?(Req, :get, 2)
      override -> override
    end
  end

  # --- retry loop ----------------------------------------------------------

  defp attempt(url, protocol, opts, retries_left, backoff) do
    case do_get(url, protocol, opts) do
      {:ok, _bytes} = ok ->
        ok

      {:error, reason} = err ->
        if retries_left > 1 and transient?(reason) do
          if backoff > 0, do: Process.sleep(backoff)
          attempt(url, protocol, opts, retries_left - 1, backoff * 2)
        else
          err
        end
    end
  end

  # 404 / missing files are permanent; everything else network-shaped retries.
  defp transient?({:file_not_found, _}), do: false
  # 408 (Request Timeout) and 429 (Too Many Requests) are the canonical
  # retry-me responses; the rest of 4xx is a permanent client error.
  defp transient?({:http_status, status}) when status in [408, 429], do: true
  defp transient?({:http_status, status}) when status in 400..499, do: false
  defp transient?({:http_status, _}), do: true
  # An oversized payload will be oversized again; do not retry it.
  defp transient?({:download_size_exceeded, _, _}), do: false
  defp transient?({:network, _}), do: true
  defp transient?({:ftp_error, _}), do: true
  defp transient?(_), do: false

  # --- HTTPS ---------------------------------------------------------------

  defp do_get(url, :https, opts) do
    if req_available?() do
      timeout = Keyword.get(opts, :timeout_ms, @default_timeout_ms)
      req_get(url, timeout, max_compressed_bytes(opts))
    else
      {:error, :req_not_available}
    end
  end

  defp do_get(url, :ftp, opts) do
    timeout = Keyword.get(opts, :timeout_ms, @default_timeout_ms)
    ftp_get(url, timeout, max_compressed_bytes(opts))
  end

  defp max_compressed_bytes(opts),
    do: Keyword.get(opts, :max_compressed_bytes, @default_max_compressed_bytes)

  # Redirects are disabled (`redirect: false`): a 3xx surfaces as a response we
  # reject, so a redirect can never carry the fetch to an off-allow-list host.
  # The body is streamed via `:into` and the transfer is halted the moment it
  # would exceed `max_bytes`, so an oversized response is never fully buffered.
  defp req_get(url, timeout, max_bytes) do
    collector = fn {:data, data}, {req, resp} ->
      acc = resp.private[:orbis_body] || []
      total = (resp.private[:orbis_total] || 0) + byte_size(data)

      if total > max_bytes do
        resp = Req.Response.put_private(resp, :orbis_over_limit, true)
        {:halt, {req, resp}}
      else
        resp =
          resp
          |> Req.Response.put_private(:orbis_body, [acc, data])
          |> Req.Response.put_private(:orbis_total, total)

        {:cont, {req, resp}}
      end
    end

    case Req.get(url,
           decode_body: false,
           redirect: false,
           into: collector,
           receive_timeout: timeout,
           connect_options: [timeout: timeout]
         ) do
      {:ok, %{status: 200} = resp} ->
        if resp.private[:orbis_over_limit] do
          {:error, {:download_size_exceeded, max_bytes, :over_limit}}
        else
          {:ok, IO.iodata_to_binary(resp.private[:orbis_body] || [])}
        end

      {:ok, %{status: 404}} ->
        {:error, {:file_not_found, url}}

      {:ok, %{status: status}} when status in 300..399 ->
        {:error, {:redirect_not_allowed, status}}

      {:ok, %{status: status}} ->
        {:error, {:http_status, status}}

      {:error, %{__exception__: true} = e} ->
        {:error, {:network, Exception.message(e)}}

      {:error, reason} ->
        {:error, {:network, reason}}
    end
  rescue
    e -> {:error, {:network, Exception.message(e)}}
  end

  # --- FTP -----------------------------------------------------------------

  defp ftp_get(url, timeout, max_bytes) do
    %URI{host: host, path: path} = URI.parse(url)

    cond do
      is_nil(host) -> {:error, {:ftp_error, :no_host}}
      is_nil(path) -> {:error, {:ftp_error, :no_path}}
      true -> ftp_fetch(host, path, timeout, max_bytes)
    end
  end

  defp ftp_fetch(host, path, timeout, max_bytes) do
    host_charlist = String.to_charlist(host)
    remote = String.to_charlist(path)

    case :ftp.open(host_charlist, mode: :passive, timeout: timeout) do
      {:ok, pid} ->
        try do
          ftp_login_and_recv(pid, remote, max_bytes)
        after
          :ftp.close(pid)
        end

      {:error, reason} ->
        {:error, {:ftp_error, reason}}
    end
  end

  # Pull the remote file in bounded chunks instead of `recv_bin/2`, which would
  # read the whole untrusted file into one binary with no cap.
  defp ftp_login_and_recv(pid, remote, max_bytes) do
    with :ok <- :ftp.user(pid, ~c"anonymous", ~c"orbis@example.com"),
         :ok <- :ftp.type(pid, :binary),
         :ok <- :ftp.recv_chunk_start(pid, remote) do
      ftp_recv_chunks(pid, max_bytes, [], 0)
    else
      {:error, :epath} -> {:error, {:file_not_found, to_string(remote)}}
      {:error, reason} -> {:error, {:ftp_error, reason}}
    end
  end

  defp ftp_recv_chunks(pid, max_bytes, acc, total) do
    case :ftp.recv_chunk(pid) do
      :ok ->
        {:ok, acc |> Enum.reverse() |> :erlang.iolist_to_binary()}

      {:ok, chunk} ->
        total = total + byte_size(chunk)

        if total > max_bytes do
          {:error, {:download_size_exceeded, max_bytes, total}}
        else
          ftp_recv_chunks(pid, max_bytes, [chunk | acc], total)
        end

      {:error, :epath} ->
        {:error, {:file_not_found, :ftp}}

      {:error, reason} ->
        {:error, {:ftp_error, reason}}
    end
  end

  # --- SSRF guard ----------------------------------------------------------

  defp check_host(url) do
    host = URI.parse(url).host

    if is_binary(host) and MapSet.member?(Catalog.allowed_hosts(), host) do
      :ok
    else
      {:error, {:unsupported_product, {:host_not_allowed, host}}}
    end
  end
end
