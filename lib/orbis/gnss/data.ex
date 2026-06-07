defmodule Orbis.GNSS.Data do
  @moduledoc """
  Optional fetch-and-cache layer for GNSS products (SP3, RINEX clock, broadcast
  navigation, IONEX).

  `Orbis.GNSS.Data` downloads, decompresses, checksums, and records provenance
  for the precise- and broadcast-product files that `Orbis.GNSS.SP3`,
  `Orbis.GNSS.Broadcast`, and `Orbis.GNSS.Positioning` consume — and then
  hands back a **local file path** (or a loaded handle). It is deliberately
  one-directional: the numerical layers never call into this module, so a solve
  never depends on network availability. You fetch once, then point the solver
  at the cached file.

  ## Quick start

      product = Orbis.GNSS.Data.mgex_sp3(:cod, ~D[2020-06-24])

      # Download (or reuse cache) and get a decompressed file path:
      {:ok, path} = Orbis.GNSS.Data.fetch(product)

      # Or fetch and load in one step:
      {:ok, sp3} = Orbis.GNSS.Data.sp3(product)
      {:ok, state} = Orbis.GNSS.SP3.position(sp3, "G01", ~N[2020-06-24 00:00:00])

  ## Catalog

  Products are identified by analysis center, content type, date, and sampling.

  Supported centers and what each publishes:

    * `:gfz` — GFZ operational rapid SP3/CLK over HTTPS (`isdc-data.gfz.de`)
    * `:cod`, `:grg`, `:wum` — CODE / CNES-CLS / Wuhan MGEX precise SP3/CLK over
      anonymous FTP (ESA GSSC mirror, `gssc.esa.int`); `:cod` also serves IONEX
    * `:igs` — the IGS merged broadcast navigation file (`:nav`) and the combined
      global ionosphere map (`:ionex`), over FTP

  Content types: `:sp3`, `:clk`, `:nav`, `:ionex`, `:obs` (station observation
  data, RINEX 3 / CRINEX). Precise products and IONEX
  follow the IGS long-name convention `AAAVPPPTTT_YYYYDDDHHMM_LEN_SMP_CNT.EXT`;
  broadcast navigation uses the no-sampling RINEX long-name
  `BRDC00IGS_R_YYYYDDDHHMM_01D_MN.rnx`. See `Orbis.GNSS.Data.Catalog`.

  ## The fetch pipeline

  `fetch/2` is cache-first:

    1. Resolve the canonical filename and cache path (pure, from the catalog).
    2. If the file is already cached, verify it: against the caller's `:sha256`
       when given, otherwise against the decompressed SHA-256 recorded in the
       file's provenance sidecar (every downloaded file has one). A verified hit
       returns with **no network**. A *corrupt* hit (checksum mismatch) or an
       *unverifiable* one (no sidecar — e.g. a hand-placed file) is, online,
       discarded and re-downloaded; offline, a corrupt hit is terminal and an
       unverifiable one is returned as the best available.
    3. Otherwise (and only when not `offline:`) download the `.gz` over HTTPS
       (`Req`, a required dependency) or FTP (`:ftp`) to memory, decompress with a gzip-bomb cap,
       verify any known checksum, and **atomically** commit the decompressed
       file into the cache (temp file + rename) together with its required
       `.provenance.json` sidecar (the commit fails if the sidecar cannot be
       written, so a cached file always carries its integrity hash).

  ## Offline mode

  Pass `offline: true` (or set `config :orbis, gnss_data_offline: true`) to
  forbid all network access: a verified cache hit is returned, a corrupt hit
  yields `{:error, {:checksum_mismatch, _, _}}`, and a miss returns
  `{:error, {:offline_miss, name}}`. This is how the test suite — and any user
  without connectivity — runs deterministically.

  ## Network tests

  Live-archive fetching is exercised by tests tagged `:network`, which are
  **excluded by default** (including in CI, which has no network); the rest of
  the suite is fully offline and deterministic. Run the live gate manually with
  `mix test --include network`.

  ## Options

    * `:offline` — when `true`, never touch the network (default from app config,
      else `false`)
    * `:cache_dir` — cache root (default `:filename.basedir(:user_cache,
      "orbis/gnss")`, overridable via `config :orbis, gnss_data_cache_dir:`)
    * `:sha256` — expected SHA-256 (hex) of the **decompressed** file; verified
      on both cache hits and fresh downloads
    * `:max_decompressed_bytes` — gzip-bomb cap (default 500 MiB)
    * `:timeout_ms` — per-attempt network timeout (default 30_000)
    * `:retries` — attempts for transient network errors (default 3)
    * `:backoff_ms` — base backoff between retries, doubled each attempt
      (default 500)
    * `:max_compressed_bytes` — cap on the compressed payload buffered into
      memory while downloading (default 64 MiB)

  ## Typed errors

  Every failure is a tagged tuple so callers can branch on it:

    * `{:error, {:offline_miss, name}}` — `offline: true` and not cached
    * `{:error, {:checksum_mismatch, expected, got}}` — digest verification failed
    * `{:error, {:unsupported_product, detail}}` — unknown center/content/sample,
      or a host outside the catalog
    * `{:error, :req_not_available}` — HTTPS downloads are disabled by config
    * `{:error, {:http_status, code}}` — non-2xx HTTP response
    * `{:error, {:redirect_not_allowed, code}}` — a 3xx redirect was refused
      (redirects are not followed, to keep the SSRF allow-list intact)
    * `{:error, {:file_not_found, url}}` — 404 / missing on the archive
    * `{:error, {:network, detail}}` — connection/timeout/DNS failure
    * `{:error, {:ftp_error, reason}}` — FTP-level failure
    * `{:error, {:download_size_exceeded, max, got}}` — compressed payload cap hit
    * `{:error, {:decompress_failed, reason}}` — corrupt gzip
    * `{:error, {:decompress_size_exceeded, max, got}}` — gzip-bomb cap hit
    * `{:error, {:cache_dir_not_writable, reason}}` — cannot create/write cache
    * `{:error, {:provenance_write_failed, reason}}` — the product downloaded but
      its required provenance sidecar could not be written (the product is rolled
      back so nothing unverifiable is left in the cache)
    * `{:error, {:unsafe_cache_name, name}}` — filename failed path-safety checks
    * `{:error, {:temp_file_error, reason}}` — temp write/rename failure
  """

  alias Orbis.GNSS.Broadcast
  alias Orbis.GNSS.Data.{Cache, Catalog, Download, Product}
  alias Orbis.GNSS.RINEX.Observations
  alias Orbis.GNSS.SP3

  @typedoc "A fetch error, always a tagged tuple. See the module docs."
  @type error :: {:error, term()}

  # --- product builders ----------------------------------------------------

  @doc """
  Build an MGEX SP3 (precise orbit) product for a center and date.

  Defaults to `05M` (5-minute) sampling; override with `sample:`.

  ## Examples

      iex> p = Orbis.GNSS.Data.mgex_sp3(:cod, ~D[2020-06-24])
      iex> p.center
      :cod
      iex> Orbis.GNSS.Data.Product.canonical_filename(p)
      {:ok, "COD0MGXFIN_20201760000_01D_05M_ORB.SP3"}
  """
  @spec mgex_sp3(atom(), Date.t(), keyword()) :: Product.t()
  def mgex_sp3(center, %Date{} = date, opts \\ []),
    do: build!(center, :sp3, date, Keyword.get(opts, :sample, "05M"))

  @doc """
  Build an MGEX clock (RINEX clock) product. Defaults to `30S` sampling.
  """
  @spec mgex_clk(atom(), Date.t(), keyword()) :: Product.t()
  def mgex_clk(center, %Date{} = date, opts \\ []),
    do: build!(center, :clk, date, Keyword.get(opts, :sample, "30S"))

  @doc """
  Build a broadcast-navigation (merged multi-GNSS RINEX NAV) product.

  Only `:igs` publishes this product (`BRDC00IGS_R_..._MN.rnx`). The RINEX
  navigation long-name carries no sampling field, so the `sample` argument is
  not part of the filename; it defaults to `01D` purely to satisfy validation.

  ## Examples

      iex> p = Orbis.GNSS.Data.mgex_nav(:igs, ~D[2020-06-25])
      iex> Orbis.GNSS.Data.Product.canonical_filename(p)
      {:ok, "BRDC00IGS_R_20201770000_01D_MN.rnx"}
  """
  @spec mgex_nav(atom(), Date.t(), keyword()) :: Product.t()
  def mgex_nav(center, %Date{} = date, opts \\ []),
    do: build!(center, :nav, date, Keyword.get(opts, :sample, "01D"))

  @doc """
  Build an IONEX (global ionosphere TEC map) product.

  Served by `:igs` (`IGS0OPSFIN`) and `:cod` (`COD0OPSFIN`). IONEX maps are
  sub-daily, so the sampling defaults to `01H`; pass `sample:` to override
  (e.g. `"02H"`).

  ## Examples

      iex> p = Orbis.GNSS.Data.mgex_ionex(:igs, ~D[2024-06-24])
      iex> Orbis.GNSS.Data.Product.canonical_filename(p)
      {:ok, "IGS0OPSFIN_20241760000_01D_01H_GIM.INX"}
  """
  @spec mgex_ionex(atom(), Date.t(), keyword()) :: Product.t()
  def mgex_ionex(center, %Date{} = date, opts \\ []),
    do: build!(center, :ionex, date, Keyword.get(opts, :sample, "01H"))

  @doc """
  Build a daily **station observation** product (RINEX 3 CRINEX, 30 s default).

  Station observation files are keyed by a 9-character site id (e.g.
  `"ESBC00DNK"`), not an analysis-center token, and resolve on the ESA GSSC
  anonymous archive's daily data tree. Override the sampling with `sample:`.

  ## Examples

      iex> p = Orbis.GNSS.Data.station_obs("ESBC00DNK", ~D[2020-06-25])
      iex> Orbis.GNSS.Data.Product.canonical_filename(p)
      {:ok, "ESBC00DNK_R_20201770000_01D_30S_MO.crx"}
  """
  @spec station_obs(String.t(), Date.t(), keyword()) :: Product.t()
  def station_obs(station, %Date{} = date, opts \\ []) when is_binary(station) do
    sample = Keyword.get(opts, :sample, "30S")

    case Product.new(:gssc, :obs, date, sample, station: station) do
      {:ok, p} -> p
      {:error, reason} -> raise ArgumentError, "invalid station OBS product: #{inspect(reason)}"
    end
  end

  @doc """
  Build a `Product` for any center/content/date/sample, returning a tuple.

  Use this instead of the bang builders when the inputs may be invalid.
  """
  @spec product(atom(), atom(), Date.t(), String.t()) ::
          {:ok, Product.t()} | {:error, {:unsupported_product, term()}}
  def product(center, content, %Date{} = date, sample),
    do: Product.new(center, content, date, sample)

  defp build!(center, content, date, sample) do
    case Product.new(center, content, date, sample) do
      {:ok, p} ->
        p

      {:error, reason} ->
        raise ArgumentError, "invalid GNSS product: #{inspect(reason)}"
    end
  end

  # --- fetch ---------------------------------------------------------------

  @doc """
  Fetch a product, returning the local path to its **decompressed** file.

  Cache-first: a verified cache hit returns immediately with no network. See the
  module docs for the full pipeline, options, and error taxonomy.

  Returns `{:ok, path}` or a typed `{:error, _}`.
  """
  @spec fetch(Product.t(), keyword()) :: {:ok, String.t()} | error()
  def fetch(%Product{} = product, opts \\ []) do
    cache_dir = cache_dir(opts)
    sha = Keyword.get(opts, :sha256)

    with {:ok, filename} <- Product.canonical_filename(product),
         {:ok, path} <- Cache.path_for(cache_dir, filename) do
      case Cache.classify(path, sha) do
        # Present and verified (against the caller's :sha256 if given, else the
        # provenance sidecar's stored hash).
        {:hit, ^path} ->
          {:ok, path}

        :absent ->
          fetch_miss(product, path, filename, sha, opts)

        # Present but unverifiable (no caller hash and no usable sidecar — a file
        # placed by hand). Online, fetch a verified, provenance-stamped copy;
        # offline, return it as the best available.
        :unverified ->
          if offline?(opts), do: {:ok, path}, else: download_and_cache(product, path, sha, opts)

        # Corrupt or stale. Offline it is terminal (nothing better to offer);
        # online we discard it and re-download — the atomic commit overwrites it,
        # and the fresh copy is verified before commit — so one bad file does not
        # permanently wedge the product.
        {:stale, mismatch} ->
          if offline?(opts),
            do: {:error, mismatch},
            else: download_and_cache(product, path, sha, opts)

        {:error, _} = err ->
          err
      end
    end
  end

  defp fetch_miss(product, path, filename, sha, opts) do
    if offline?(opts) do
      {:error, {:offline_miss, filename}}
    else
      download_and_cache(product, path, sha, opts)
    end
  end

  defp download_and_cache(product, path, sha, opts) do
    max_bytes = Keyword.get(opts, :max_decompressed_bytes, Cache.default_max_decompressed_bytes())

    with {:ok, url} <- Product.archive_url(product),
         {:ok, protocol} <- protocol_for(product),
         {:ok, compressed} <- Download.get(url, protocol, opts),
         {:ok, decompressed} <- Cache.gunzip(compressed, max_bytes),
         :ok <- verify(decompressed, sha),
         provenance = provenance(url, compressed, decompressed),
         {:ok, ^path} <- Cache.commit(path, decompressed, provenance) do
      {:ok, path}
    end
  end

  # Station observation products are not analysis-center products; they resolve
  # their protocol from the dedicated station archive path, not the @centers
  # token table.
  defp protocol_for(%Product{content: :obs}), do: {:ok, Catalog.station_obs_protocol()}
  defp protocol_for(%Product{center: center}), do: Catalog.protocol(center)

  defp verify(_decompressed, nil), do: :ok

  defp verify(decompressed, expected) do
    got = Cache.sha256(decompressed)

    if String.downcase(expected) == got do
      :ok
    else
      {:error, {:checksum_mismatch, expected, got}}
    end
  end

  defp provenance(url, compressed, decompressed) do
    %{
      "source_url" => url,
      "sha256_compressed" => Cache.sha256(compressed),
      "sha256_decompressed" => Cache.sha256(decompressed),
      "size_compressed" => byte_size(compressed),
      "size_decompressed" => byte_size(decompressed),
      "fetched_at" => DateTime.utc_now() |> DateTime.to_iso8601(),
      "fetcher" => "Orbis.GNSS.Data"
    }
  end

  # --- convenience loaders -------------------------------------------------

  @doc """
  Fetch an SP3 product and load it into an `Orbis.GNSS.SP3` handle.

  Equivalent to `fetch/2` followed by `Orbis.GNSS.SP3.load/1`. Returns
  `{:ok, %Orbis.GNSS.SP3{}}` or a typed error.
  """
  @spec sp3(Product.t(), keyword()) :: {:ok, SP3.t()} | error()
  def sp3(%Product{content: :sp3} = product, opts \\ []) do
    with {:ok, path} <- fetch(product, opts), do: SP3.load(path)
  end

  @doc """
  Fetch a broadcast-navigation product and load it into an
  `Orbis.GNSS.Broadcast` handle.

  Returns `{:ok, %Orbis.GNSS.Broadcast{}}` or a typed error.
  """
  @spec broadcast(Product.t(), keyword()) :: {:ok, Broadcast.t()} | error()
  def broadcast(%Product{content: :nav} = product, opts \\ []) do
    with {:ok, path} <- fetch(product, opts), do: Broadcast.load(path)
  end

  @doc """
  Fetch a station observation product and load it into an `Orbis.GNSS.RINEX.Observations`
  handle.

  `fetch/2` gunzips the `.gz`; the committed cache file is the (still Hatanaka)
  CRINEX text, which `Orbis.GNSS.RINEX.Observations.load/1` decodes before parsing. Returns
  `{:ok, %Orbis.GNSS.RINEX.Observations{}}` or a typed error.
  """
  @spec observations(Product.t(), keyword()) :: {:ok, Observations.t()} | error()
  def observations(%Product{content: :obs} = product, opts \\ []) do
    with {:ok, path} <- fetch(product, opts), do: Observations.load(path)
  end

  # --- option resolution ---------------------------------------------------

  defp cache_dir(opts) do
    Keyword.get(opts, :cache_dir) ||
      Application.get_env(:orbis, :gnss_data_cache_dir) ||
      Cache.default_dir()
  end

  defp offline?(opts) do
    case Keyword.fetch(opts, :offline) do
      {:ok, value} -> value
      :error -> Application.get_env(:orbis, :gnss_data_offline, false)
    end
  end
end
