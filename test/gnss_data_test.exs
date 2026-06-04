defmodule Orbis.GnssDataTest do
  # Not async: several tests exercise the app-config knobs (`:gnss_data_offline`,
  # `:gnss_data_req_available`) by setting global application env, which must not
  # race other tests.
  use ExUnit.Case, async: false

  alias Orbis.GnssData
  alias Orbis.GnssData.Cache

  @gz_fixture Path.join(__DIR__, "fixtures/gnss_data/GBM0MGXRAP_20201760000_01D_05M_ORB.SP3.gz")
  @nav_fixture Path.join(__DIR__, "fixtures/nav/ESBC00DNK_R_20201770000_01D_MN.rnx")

  setup do
    # A fresh, unique cache dir per test, removed afterwards. Never the user cache.
    cache_dir =
      Path.join(System.tmp_dir!(), "orbis_gnss_cache_#{System.unique_integer([:positive])}")

    on_exit(fn -> File.rm_rf(cache_dir) end)
    {:ok, cache_dir: cache_dir}
  end

  # The decompressed SP3 bytes the committed .gz fixture expands to.
  defp sp3_bytes, do: :zlib.gunzip(File.read!(@gz_fixture))

  # Seed `cache_dir` with `bytes` under the canonical name of `product`.
  defp seed(cache_dir, product, bytes) do
    {:ok, filename} = GnssData.Product.canonical_filename(product)
    {:ok, path} = Cache.path_for(cache_dir, filename)
    File.mkdir_p!(Path.dirname(path))
    File.write!(path, bytes)
    path
  end

  describe "gunzip/2 (decompression)" do
    test "expands the committed .gz fixture to the expected SP3 header" do
      compressed = File.read!(@gz_fixture)
      assert {:ok, decompressed} = Cache.gunzip(compressed)
      assert String.starts_with?(decompressed, "#cP2020")
      assert decompressed =~ "PG01  15000.000000"
    end

    test "rejects corrupt gzip data" do
      assert {:error, {:decompress_failed, _}} = Cache.gunzip(<<0, 1, 2, 3, 4, 5>>)
    end

    test "enforces the gzip-bomb size cap" do
      compressed = File.read!(@gz_fixture)
      # The fixture expands to ~944 bytes; cap at 100 to trip the guard.
      assert {:error, {:decompress_size_exceeded, 100, got}} = Cache.gunzip(compressed, 100)
      assert got > 100
    end

    test "aborts a real gzip bomb without materializing its full expansion" do
      # 256 MiB of zeros compresses to a tiny buffer but would dwarf a small cap.
      bomb = :zlib.gzip(:binary.copy(<<0>>, 256 * 1024 * 1024))
      cap = 1_000_000

      # Measure the peak memory the inflate uses: a streaming guard must stay
      # bounded near the cap, not balloon to the full 256 MiB expansion.
      before = :erlang.memory(:total)
      assert {:error, {:decompress_size_exceeded, ^cap, got}} = Cache.gunzip(bomb, cap)
      peak = :erlang.memory(:total) - before

      # We aborted as soon as we crossed the cap, so the reported size is just
      # over it (not the full 256 MiB), and we never allocated the remainder.
      assert got > cap
      assert got < 16 * 1024 * 1024
      assert peak < 64 * 1024 * 1024
    end
  end

  describe "fetch/2 offline cache hits" do
    test "returns a cached file without network", %{cache_dir: cache_dir} do
      product = GnssData.mgex_sp3(:cod, ~D[2020-06-24])
      seeded = seed(cache_dir, product, sp3_bytes())

      assert {:ok, ^seeded} = GnssData.fetch(product, offline: true, cache_dir: cache_dir)
    end

    test "verifies a known checksum on a cache hit", %{cache_dir: cache_dir} do
      product = GnssData.mgex_sp3(:cod, ~D[2020-06-24])
      bytes = sp3_bytes()
      seed(cache_dir, product, bytes)
      sha = Cache.sha256(bytes)

      assert {:ok, _path} =
               GnssData.fetch(product, offline: true, cache_dir: cache_dir, sha256: sha)
    end

    test "rejects a cache hit whose checksum does not match", %{cache_dir: cache_dir} do
      product = GnssData.mgex_sp3(:cod, ~D[2020-06-24])
      seed(cache_dir, product, sp3_bytes())
      wrong = String.duplicate("0", 64)

      assert {:error, {:checksum_mismatch, ^wrong, _got}} =
               GnssData.fetch(product, offline: true, cache_dir: cache_dir, sha256: wrong)
    end
  end

  describe "fetch/2 offline misses and errors" do
    test "returns an offline_miss when the file is absent", %{cache_dir: cache_dir} do
      product = GnssData.mgex_sp3(:cod, ~D[2020-06-24])

      assert {:error, {:offline_miss, "COD0MGXFIN_20201760000_01D_05M_ORB.SP3"}} =
               GnssData.fetch(product, offline: true, cache_dir: cache_dir)
    end

    test "an unsupported product is rejected before any network", %{cache_dir: cache_dir} do
      # Construct a Product struct directly to bypass the validating builder.
      bad = %GnssData.Product{center: :nope, content: :sp3, date: ~D[2020-06-24], sample: "05M"}

      assert {:error, {:unsupported_product, {:center, :nope}}} =
               GnssData.fetch(bad, offline: true, cache_dir: cache_dir)
    end

    test "respects app-config offline mode", %{cache_dir: cache_dir} do
      Application.put_env(:orbis, :gnss_data_offline, true)
      on_exit(fn -> Application.delete_env(:orbis, :gnss_data_offline) end)

      product = GnssData.mgex_sp3(:cod, ~D[2020-06-24])
      assert {:error, {:offline_miss, _}} = GnssData.fetch(product, cache_dir: cache_dir)
    end

    test "a corrupt cache hit is terminal offline but self-heals online", %{cache_dir: cache_dir} do
      product = GnssData.mgex_sp3(:gfz, ~D[2020-06-24], sample: "15M")
      seed(cache_dir, product, "corrupt bytes that do not match the digest")
      sha = Cache.sha256(sp3_bytes())

      # Offline: a mismatch is terminal (nothing better to offer).
      assert {:error, {:checksum_mismatch, ^sha, _got}} =
               GnssData.fetch(product, offline: true, cache_dir: cache_dir, sha256: sha)

      # Online: the poisoned entry is discarded and a fresh download is
      # attempted. We disable Req so the re-download stops at :req_not_available
      # rather than hitting the network — proving fetch did NOT terminate on the
      # stale checksum_mismatch.
      Application.put_env(:orbis, :gnss_data_req_available, false)
      on_exit(fn -> Application.delete_env(:orbis, :gnss_data_req_available) end)

      assert {:error, :req_not_available} =
               GnssData.fetch(product,
                 cache_dir: cache_dir,
                 sha256: sha,
                 retries: 1,
                 backoff_ms: 0
               )
    end
  end

  describe "convenience loaders (offline)" do
    test "sp3/2 returns a queryable SP3 handle", %{cache_dir: cache_dir} do
      product = GnssData.mgex_sp3(:cod, ~D[2020-06-24])
      seed(cache_dir, product, sp3_bytes())

      assert {:ok, %Orbis.SP3{} = sp3} =
               GnssData.sp3(product, offline: true, cache_dir: cache_dir)

      assert {:ok, state} = Orbis.SP3.position(sp3, "G01", ~N[2020-06-24 00:00:00])
      assert_in_delta state.x_m, 15_000_000.0, 1.0e-3
    end

    test "broadcast/2 returns a BroadcastEphemeris handle", %{cache_dir: cache_dir} do
      # :igs resolves to the real merged-nav name BRDC00IGS_R_20201770000_01D_MN.rnx.
      product = GnssData.mgex_nav(:igs, ~D[2020-06-25])

      assert {:ok, "BRDC00IGS_R_20201770000_01D_MN.rnx"} =
               GnssData.Product.canonical_filename(product)

      seed(cache_dir, product, File.read!(@nav_fixture))

      assert {:ok, %Orbis.BroadcastEphemeris{}} =
               GnssData.broadcast(product, offline: true, cache_dir: cache_dir)
    end

    test "an offline miss propagates through sp3/2", %{cache_dir: cache_dir} do
      product = GnssData.mgex_sp3(:cod, ~D[2020-06-24])

      assert {:error, {:offline_miss, _}} =
               GnssData.sp3(product, offline: true, cache_dir: cache_dir)
    end
  end

  describe "path-traversal safety" do
    test "rejects a cache filename containing a path separator", %{cache_dir: cache_dir} do
      assert {:error, {:unsafe_cache_name, _}} = Cache.path_for(cache_dir, "../escape")
      assert {:error, {:unsafe_cache_name, _}} = Cache.path_for(cache_dir, "a/b")
      assert {:error, {:unsafe_cache_name, _}} = Cache.path_for(cache_dir, "/etc/passwd")
    end

    test "accepts a valid canonical name" do
      assert {:ok, path} = Cache.path_for("/tmp/cache", "GBM0MGXRAP_20201760000_01D_05M_ORB.SP3")
      assert path == "/tmp/cache/GBM0MGXRAP_20201760000_01D_05M_ORB.SP3"
    end
  end

  describe "atomic commit + provenance" do
    test "commit writes the file and a provenance sidecar", %{cache_dir: cache_dir} do
      product = GnssData.mgex_sp3(:cod, ~D[2020-06-24])
      {:ok, filename} = GnssData.Product.canonical_filename(product)
      {:ok, path} = Cache.path_for(cache_dir, filename)

      bytes = sp3_bytes()

      provenance = %{
        "source_url" => "https://example/test",
        "sha256_decompressed" => Cache.sha256(bytes),
        "size_decompressed" => byte_size(bytes)
      }

      assert {:ok, ^path} = Cache.commit(path, bytes, provenance)
      assert File.read!(path) == bytes

      assert {:ok, decoded} = Cache.read_provenance(path)
      assert decoded["sha256_decompressed"] == Cache.sha256(bytes)
      # No leftover temp files in the cache directory.
      assert Enum.all?(File.ls!(cache_dir), &(not String.starts_with?(&1, ".tmp-")))
    end

    test "commit reports an unwritable cache directory" do
      # A path under a regular file (not a directory) cannot be created.
      blocker =
        Path.join(System.tmp_dir!(), "orbis_blocker_#{System.unique_integer([:positive])}")

      File.write!(blocker, "x")
      on_exit(fn -> File.rm(blocker) end)

      target = Path.join([blocker, "sub", "FILE"])
      assert {:error, {:cache_dir_not_writable, _}} = Cache.commit(target, "data", %{})
    end
  end

  describe "default cache integrity (no caller checksum)" do
    test "a default cache hit is verified against the provenance sidecar and self-heals online",
         %{cache_dir: cache_dir} do
      product = GnssData.mgex_sp3(:gfz, ~D[2020-06-24], sample: "15M")
      {:ok, filename} = GnssData.Product.canonical_filename(product)
      {:ok, path} = Cache.path_for(cache_dir, filename)
      good = sp3_bytes()

      # A committed file always carries its decompressed hash in the sidecar.
      assert {:ok, ^path} =
               Cache.commit(path, good, %{"sha256_decompressed" => Cache.sha256(good)})

      # A clean default (no :sha256) hit is returned with no network.
      assert {:ok, ^path} = GnssData.fetch(product, offline: true, cache_dir: cache_dir)

      # Corrupt the cached file but leave the sidecar: the stored hash no longer
      # matches the bytes, and the default hit must NOT trust it.
      File.write!(path, "corrupted after caching")

      # Offline: the sidecar mismatch is detected and is terminal.
      assert {:error, {:checksum_mismatch, _, _}} =
               GnssData.fetch(product, offline: true, cache_dir: cache_dir)

      # Online: the poisoned entry is discarded and a fresh download attempted.
      # Req disabled so it stops at :req_not_available — proving fetch did NOT
      # serve the corrupt cached bytes.
      Application.put_env(:orbis, :gnss_data_req_available, false)
      on_exit(fn -> Application.delete_env(:orbis, :gnss_data_req_available) end)

      assert {:error, :req_not_available} =
               GnssData.fetch(product, cache_dir: cache_dir, retries: 1, backoff_ms: 0)
    end

    test "an unverifiable cache hit (no sidecar) is served offline but refreshed online",
         %{cache_dir: cache_dir} do
      product = GnssData.mgex_sp3(:gfz, ~D[2020-06-24], sample: "15M")
      # seed/3 writes the product file but no provenance sidecar.
      path = seed(cache_dir, product, sp3_bytes())

      # Offline: nothing better to offer, so the unprovenanced file is returned.
      assert {:ok, ^path} = GnssData.fetch(product, offline: true, cache_dir: cache_dir)

      # Online: treated as a miss and re-downloaded rather than silently trusted
      # (Req disabled -> stops at :req_not_available).
      Application.put_env(:orbis, :gnss_data_req_available, false)
      on_exit(fn -> Application.delete_env(:orbis, :gnss_data_req_available) end)

      assert {:error, :req_not_available} =
               GnssData.fetch(product, cache_dir: cache_dir, retries: 1, backoff_ms: 0)
    end
  end

  describe "req availability" do
    test "fetch surfaces :req_not_available for an HTTPS center when Req is disabled",
         %{cache_dir: cache_dir} do
      Application.put_env(:orbis, :gnss_data_req_available, false)
      on_exit(fn -> Application.delete_env(:orbis, :gnss_data_req_available) end)

      product = GnssData.mgex_sp3(:gfz, ~D[2020-06-24], sample: "15M")

      # offline:false + cache miss + Req disabled -> :req_not_available.
      assert {:error, :req_not_available} =
               GnssData.fetch(product,
                 cache_dir: cache_dir,
                 retries: 1,
                 backoff_ms: 0
               )
    end
  end

  describe "real network fetch" do
    @tag :network
    test "downloads, decompresses, caches, and loads a real GFZ SP3", %{cache_dir: cache_dir} do
      # GFZ's operational rapid SP3 for 2020-06-24, served over HTTPS:
      # https://isdc-data.gfz.de/gnss/products/rapid/w2111/GFZ0OPSRAP_20201760000_01D_15M_ORB.SP3.gz
      product = GnssData.mgex_sp3(:gfz, ~D[2020-06-24], sample: "15M")

      assert {:ok, path} = GnssData.fetch(product, cache_dir: cache_dir)
      assert File.exists?(path)
      assert {:ok, _prov} = Cache.read_provenance(path)

      # A second fetch is served from the cache with no network.
      assert {:ok, ^path} = GnssData.fetch(product, offline: true, cache_dir: cache_dir)

      assert {:ok, %Orbis.SP3{} = sp3} = Orbis.SP3.load(path)
      assert {:ok, _state} = Orbis.SP3.position(sp3, "G01", ~N[2020-06-24 00:00:00])
    end

    @tag :network
    test "downloads a real MGEX SP3 over FTP (chunked recv)", %{cache_dir: cache_dir} do
      # CODE's MGEX final orbit for 2020-06-24, anonymous FTP on the ESA GSSC
      # mirror: ftp://gssc.esa.int/gnss/products/2111/COD0MGXFIN_20201760000_01D_05M_ORB.SP3.gz
      product = GnssData.mgex_sp3(:cod, ~D[2020-06-24])

      assert {:ok, path} = GnssData.fetch(product, cache_dir: cache_dir)
      assert File.exists?(path)
      assert {:ok, %Orbis.SP3{}} = Orbis.SP3.load(path)
    end

    @tag :network
    test "downloads the real merged broadcast navigation file over FTP", %{cache_dir: cache_dir} do
      # IGS merged multi-GNSS broadcast nav for 2020-06-25, anonymous FTP:
      # ftp://gssc.esa.int/gnss/data/daily/2020/177/BRDC00IGS_R_20201770000_01D_MN.rnx.gz
      product = GnssData.mgex_nav(:igs, ~D[2020-06-25])

      assert {:ok, path} = GnssData.fetch(product, cache_dir: cache_dir)
      assert File.exists?(path)
      assert {:ok, %Orbis.BroadcastEphemeris{}} = Orbis.BroadcastEphemeris.load(path)
    end
  end
end
