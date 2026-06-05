defmodule Orbis.GNSS.Data.Catalog do
  @moduledoc """
  Static, pure catalog of GNSS analysis centers and the rules that turn a
  product specification into a canonical filename and a full archive URL.

  Everything here is deterministic and network-free: GPS-week and day-of-year
  arithmetic, the IGS long-name filename convention, and the per-center archive
  layout. The fetch pipeline derives **all** network hosts and cache filenames
  from this module, which is what keeps the layer safe against SSRF (only known
  hosts are ever contacted) and path traversal (cache names come only from a
  validated canonical filename).

  ## Analysis centers

  Each center maps to a transfer protocol (`:https` or `:ftp`), an archive host,
  the product token it publishes, and the directory layout for each content
  type. Only centers and paths that resolve against a live, anonymous archive
  are listed; every entry below has been checked against its server.

  | Code   | Center                          | Protocol | Host                |
  |--------|---------------------------------|----------|---------------------|
  | `:gfz` | GFZ Potsdam (operational rapid) | HTTPS    | `isdc-data.gfz.de`  |
  | `:cod` | CODE / University of Bern (MGEX)| FTP      | `gssc.esa.int`      |
  | `:grg` | CNES/CLS (MGEX, final)          | FTP      | `gssc.esa.int`      |
  | `:wum` | Wuhan University (MGEX, final)  | FTP      | `gssc.esa.int`      |
  | `:igs` | IGS combined (broadcast / IONEX)| FTP      | `gssc.esa.int`      |

  The ESA Galileo Science Support Centre (`gssc.esa.int`) mirrors the IGS/MGEX
  archive over anonymous FTP and is the source for the MGEX precise products,
  the merged broadcast navigation file, and the global ionosphere maps. GFZ's
  own HTTPS server (`isdc-data.gfz.de`) carries its operational rapid SP3/CLK.

  ## Content types and what each center serves

    * `:sp3`, `:clk` — precise orbits and clocks. `:gfz` (operational rapid),
      `:cod`, `:grg`, `:wum` (MGEX). The GFZ token is `GFZ0OPSRAP`; the MGEX
      tokens are `COD0MGXFIN`, `GRG0MGXFIN`, `WUM0MGXFIN`.
    * `:nav` — the IGS merged multi-GNSS broadcast navigation file
      (`BRDC00IGS_R_..._MN.rnx`). Only `:igs` publishes it.
    * `:ionex` — the global ionosphere TEC map (`..._GIM.INX`). `:igs` serves the
      combined `IGS0OPSFIN` map; `:cod` serves `COD0OPSFIN`. IONEX cadence is
      sub-daily, so the default sampling is `01H`/`02H`, not `01D`.

  ## Filename conventions

  Precise products and IONEX follow the IGS long-name convention
  `AAAVPPPTTT_YYYYDDDHHMM_LEN_SMP_CNT.EXT` (e.g.
  `GFZ0OPSRAP_20201760000_01D_15M_ORB.SP3`). Broadcast navigation uses the
  RINEX long-name `SSSSMRCCC_R_YYYYDDDHHMM_LEN_CNT.fmt` with **no** sampling
  field and a lowercase extension (e.g. `BRDC00IGS_R_20201770000_01D_MN.rnx`).
  """

  @gps_epoch_jdn 2_444_245
  @seconds_per_day 86_400

  # Each center definition:
  #   :name      human-readable name
  #   :protocol  :https | :ftp
  #   :host      archive host (the SSRF allow-list source)
  #   :root      archive root URL (no trailing slash)
  #   :tokens    %{content => token_prefix} where the token already includes the
  #              solution code (e.g. "GFZ0OPSRAP", "COD0MGXFIN", "BRDC00IGS").
  #   :layouts   %{content => layout} describing the directory tree per content.
  @centers %{
    gfz: %{
      name: "GFZ (Deutsches GeoForschungsZentrum Potsdam) operational AC",
      protocol: :https,
      host: "isdc-data.gfz.de",
      root: "https://isdc-data.gfz.de/gnss/products",
      tokens: %{sp3: "GFZ0OPSRAP", clk: "GFZ0OPSRAP"},
      layouts: %{sp3: :gfz_class_week, clk: :gfz_class_week}
    },
    cod: %{
      name: "Center for Orbit Determination in Europe (CODE), University of Bern",
      protocol: :ftp,
      host: "gssc.esa.int",
      root: "ftp://gssc.esa.int/gnss",
      tokens: %{sp3: "COD0MGXFIN", clk: "COD0MGXFIN", ionex: "COD0OPSFIN"},
      layouts: %{sp3: :products_week, clk: :products_week, ionex: :ionex_year_doy}
    },
    grg: %{
      name: "CNES/CLS Analysis Center (MGEX)",
      protocol: :ftp,
      host: "gssc.esa.int",
      root: "ftp://gssc.esa.int/gnss",
      tokens: %{sp3: "GRG0MGXFIN", clk: "GRG0MGXFIN"},
      layouts: %{sp3: :products_week, clk: :products_week}
    },
    wum: %{
      name: "Wuhan University GNSS Analysis Center (MGEX)",
      protocol: :ftp,
      host: "gssc.esa.int",
      root: "ftp://gssc.esa.int/gnss",
      tokens: %{sp3: "WUM0MGXFIN", clk: "WUM0MGXFIN"},
      layouts: %{sp3: :products_week, clk: :products_week}
    },
    igs: %{
      name: "IGS Combined Analysis Center",
      protocol: :ftp,
      host: "gssc.esa.int",
      root: "ftp://gssc.esa.int/gnss",
      tokens: %{nav: "BRDC00IGS", ionex: "IGS0OPSFIN"},
      layouts: %{nav: :data_daily_year_doy, ionex: :ionex_year_doy}
    }
  }

  # Content type -> filename shape.
  #   :code   the 2- or 3-letter content code
  #   :ext    the file extension (case as published)
  #   :kind   :sampled (AAAVPPPTTT_DATE_LEN_SMP_CNT.EXT) or
  #           :nav (SSSSMRCCC_R_DATE_LEN_CNT.ext, no SMP field)
  @content %{
    sp3: %{code: "ORB", ext: "SP3", kind: :sampled},
    clk: %{code: "CLK", ext: "CLK", kind: :sampled},
    nav: %{code: "MN", ext: "rnx", kind: :nav},
    ionex: %{code: "GIM", ext: "INX", kind: :sampled},
    obs: %{code: "MO", ext: "crx", kind: :obs_station}
  }

  @doc """
  All supported analysis-center codes.
  """
  @spec centers() :: [atom()]
  def centers, do: Map.keys(@centers)

  @doc """
  All supported content-type codes.
  """
  @spec content_types() :: [atom()]
  def content_types, do: Map.keys(@content)

  @doc """
  Look up a center's static definition.

  Returns `{:ok, map}` or `{:error, {:unsupported_product, {:center, code}}}`.
  """
  @spec center(atom()) :: {:ok, map()} | {:error, {:unsupported_product, term()}}
  def center(code) when is_atom(code) do
    case Map.fetch(@centers, code) do
      {:ok, def} -> {:ok, def}
      :error -> {:error, {:unsupported_product, {:center, code}}}
    end
  end

  def center(code), do: {:error, {:unsupported_product, {:center, code}}}

  @doc """
  Human-readable center name, or `nil` if the code is unknown.
  """
  @spec center_name(atom()) :: String.t() | nil
  def center_name(code) do
    case center(code) do
      {:ok, def} -> def.name
      _ -> nil
    end
  end

  @doc """
  The content-type descriptor (`%{code:, ext:, kind:}`) for a content type.

  Returns `{:ok, map}` or `{:error, {:unsupported_product, {:content, type}}}`.
  """
  @spec content(atom()) :: {:ok, map()} | {:error, {:unsupported_product, term()}}
  def content(type) when is_atom(type) do
    case Map.fetch(@content, type) do
      {:ok, descriptor} -> {:ok, descriptor}
      :error -> {:error, {:unsupported_product, {:content, type}}}
    end
  end

  def content(type), do: {:error, {:unsupported_product, {:content, type}}}

  @doc """
  The GPS week number for a calendar date.

  GPS week 0 began on 1980-01-06. Uses exact integer day arithmetic, so it is
  leap-second-agnostic (week numbering is a calendar count, not a clock count).

  ## Examples

      iex> Orbis.GNSS.Data.Catalog.gps_week(~D[2020-06-24])
      2111
  """
  @spec gps_week(Date.t()) :: non_neg_integer()
  def gps_week(%Date{} = date) do
    div(days_since_gps_epoch(date), 7)
  end

  @doc """
  The GPS day-of-week for a calendar date (`0` = Sunday … `6` = Saturday).

  ## Examples

      iex> Orbis.GNSS.Data.Catalog.gps_day_of_week(~D[2020-06-24])
      3
  """
  @spec gps_day_of_week(Date.t()) :: 0..6
  def gps_day_of_week(%Date{} = date) do
    rem(days_since_gps_epoch(date), 7)
  end

  @doc """
  The day-of-year (`001`–`366`) for a calendar date.

  ## Examples

      iex> Orbis.GNSS.Data.Catalog.day_of_year(~D[2020-06-24])
      176
  """
  @spec day_of_year(Date.t()) :: 1..366
  def day_of_year(%Date{} = date) do
    jdn(date.year, date.month, date.day) - jdn(date.year, 1, 1) + 1
  end

  @doc """
  Build the canonical IGS long-name filename for a product.

  Precise products and IONEX use `AAAVPPPTTT_YYYYDDDHHMM_LEN_SMP_CNT.EXT`;
  broadcast navigation uses the no-sampling RINEX form
  `SSSSMRCCC_R_YYYYDDDHHMM_LEN_CNT.ext`. The center must actually publish the
  requested content type.

  Returns `{:ok, filename}` or an `{:error, {:unsupported_product, _}}` tuple.

  ## Examples

      iex> Orbis.GNSS.Data.Catalog.canonical_filename(:gfz, :sp3, ~D[2020-06-24], "15M")
      {:ok, "GFZ0OPSRAP_20201760000_01D_15M_ORB.SP3"}

      iex> Orbis.GNSS.Data.Catalog.canonical_filename(:igs, :nav, ~D[2020-06-25], "01D")
      {:ok, "BRDC00IGS_R_20201770000_01D_MN.rnx"}
  """
  @spec canonical_filename(atom(), atom(), Date.t(), String.t()) ::
          {:ok, String.t()} | {:error, {:unsupported_product, term()}}
  def canonical_filename(center, content, %Date{} = date, sample)
      when is_atom(center) and is_atom(content) and is_binary(sample) do
    with {:ok, cdef} <- center(center),
         {:ok, descriptor} <- content(content),
         {:ok, token} <- token_for(cdef, content),
         :ok <- validate_sample(sample) do
      {:ok, build_filename(descriptor, token, date, sample)}
    end
  end

  def canonical_filename(_center, _content, _date, _sample),
    do: {:error, {:unsupported_product, :bad_arguments}}

  @doc """
  Build the canonical IGS long-name filename for a daily station observation
  product (RINEX 3 CRINEX), e.g.
  `ESBC00DNK_R_20201770000_01D_30S_MO.crx`.

  Station observation files are keyed by a 9-character site id, not an
  analysis-center token, so they have their own builder.

  ## Examples

      iex> Orbis.GNSS.Data.Catalog.station_obs_filename("ESBC00DNK", ~D[2020-06-25], "30S")
      {:ok, "ESBC00DNK_R_20201770000_01D_30S_MO.crx"}
  """
  @spec station_obs_filename(String.t(), Date.t(), String.t()) ::
          {:ok, String.t()} | {:error, {:unsupported_product, term()}}
  def station_obs_filename(station, %Date{} = date, sample)
      when is_binary(station) and is_binary(sample) do
    with :ok <- validate_station(station),
         :ok <- validate_sample(sample),
         {:ok, descriptor} <- content(:obs) do
      {:ok, "#{station}_R_#{date_block(date)}_01D_#{sample}_#{descriptor.code}.#{descriptor.ext}"}
    end
  end

  def station_obs_filename(_station, _date, _sample),
    do: {:error, {:unsupported_product, :bad_arguments}}

  @doc """
  Build the full, compressed (`.gz`) archive URL for a daily station observation
  product on the ESA GSSC anonymous archive (the same daily data tree the
  broadcast navigation file uses).

  ## Examples

      iex> Orbis.GNSS.Data.Catalog.station_obs_url("ESBC00DNK", ~D[2020-06-25], "30S")
      {:ok, "ftp://gssc.esa.int/gnss/data/daily/2020/177/ESBC00DNK_R_20201770000_01D_30S_MO.crx.gz"}
  """
  @spec station_obs_url(String.t(), Date.t(), String.t()) ::
          {:ok, String.t()} | {:error, {:unsupported_product, term()}}
  def station_obs_url(station, %Date{} = date, sample) do
    with {:ok, filename} <- station_obs_filename(station, date, sample) do
      root = "ftp://gssc.esa.int/gnss"
      {:ok, "#{root}/#{dir_path(:data_daily_year_doy, date)}/#{filename}.gz"}
    end
  end

  @doc """
  Build the full, compressed (`.gz`) archive URL for a product.

  The directory follows the center/content layout; the filename is the canonical
  long-name plus a `.gz` suffix. The host is always one of the catalog hosts,
  never caller-supplied input.

  Returns `{:ok, url}` or an `{:error, {:unsupported_product, _}}` tuple.

  ## Examples

      iex> Orbis.GNSS.Data.Catalog.archive_url(:gfz, :sp3, ~D[2020-06-24], "15M")
      {:ok, "https://isdc-data.gfz.de/gnss/products/rapid/w2111/GFZ0OPSRAP_20201760000_01D_15M_ORB.SP3.gz"}

      iex> Orbis.GNSS.Data.Catalog.archive_url(:igs, :nav, ~D[2020-06-25], "01D")
      {:ok, "ftp://gssc.esa.int/gnss/data/daily/2020/177/BRDC00IGS_R_20201770000_01D_MN.rnx.gz"}
  """
  @spec archive_url(atom(), atom(), Date.t(), String.t()) ::
          {:ok, String.t()} | {:error, {:unsupported_product, term()}}
  def archive_url(center, content, %Date{} = date, sample) do
    with {:ok, cdef} <- center(center),
         {:ok, filename} <- canonical_filename(center, content, date, sample),
         {:ok, layout} <- layout_for(cdef, content) do
      {:ok, "#{cdef.root}/#{dir_path(layout, date)}/#{filename}.gz"}
    end
  end

  def archive_url(_center, _content, _date, _sample),
    do: {:error, {:unsupported_product, :bad_arguments}}

  @doc """
  The transfer protocol (`:https` or `:ftp`) for a center.
  """
  @spec protocol(atom()) :: {:ok, :https | :ftp} | {:error, {:unsupported_product, term()}}
  def protocol(center) do
    with {:ok, cdef} <- center(center), do: {:ok, cdef.protocol}
  end

  @doc """
  The set of hosts the layer is permitted to contact.

  Used by the download path as an allow-list so a malformed or unexpected URL
  can never cause a request to an off-catalog host.
  """
  @spec allowed_hosts() :: MapSet.t(String.t())
  def allowed_hosts do
    @centers |> Map.values() |> MapSet.new(& &1.host)
  end

  # --- filename construction -----------------------------------------------

  defp build_filename(%{kind: :sampled, code: code, ext: ext}, token, date, sample) do
    "#{token}_#{date_block(date)}_01D_#{sample}_#{code}.#{ext}"
  end

  defp build_filename(%{kind: :nav, code: code, ext: ext}, token, date, _sample) do
    "#{token}_R_#{date_block(date)}_01D_#{code}.#{ext}"
  end

  defp date_block(%Date{} = date), do: "#{date.year}#{pad3(day_of_year(date))}0000"

  defp token_for(%{tokens: tokens}, content) do
    case Map.fetch(tokens, content) do
      {:ok, token} -> {:ok, token}
      :error -> {:error, {:unsupported_product, {:content_not_served, content}}}
    end
  end

  defp layout_for(%{layouts: layouts}, content) do
    case Map.fetch(layouts, content) do
      {:ok, layout} -> {:ok, layout}
      :error -> {:error, {:unsupported_product, {:content_not_served, content}}}
    end
  end

  # --- directory layouts ---------------------------------------------------

  # GFZ HTTPS operational tree: <root>/rapid/w<gpsweek>/<file>.
  defp dir_path(:gfz_class_week, date), do: "rapid/w#{gps_week(date)}"

  # ESA GSSC flat MGEX products tree: <root>/products/<gpsweek>/<file>.
  defp dir_path(:products_week, date), do: "products/#{gps_week(date)}"

  # ESA GSSC daily broadcast tree: <root>/data/daily/<year>/<doy>/<file>.
  defp dir_path(:data_daily_year_doy, date),
    do: "data/daily/#{date.year}/#{pad3(day_of_year(date))}"

  # ESA GSSC IONEX tree: <root>/products/ionex/<year>/<doy>/<file>.
  defp dir_path(:ionex_year_doy, date),
    do: "products/ionex/#{date.year}/#{pad3(day_of_year(date))}"

  # --- validation ----------------------------------------------------------

  defp validate_sample(<<_::binary-size(3)>> = s) do
    if String.match?(s, ~r/\A[0-9]{2}[A-Z]\z/), do: :ok, else: bad_sample(s)
  end

  defp validate_sample(s), do: bad_sample(s)

  defp bad_sample(s), do: {:error, {:unsupported_product, {:sample, s}}}

  # A RINEX 3 site id is a 9-character SSSSMRCCC token (4-char monument, marker,
  # receiver, 3-char ISO country), upper-case alphanumeric. Validating it keeps
  # the cache filename path-safe and the archive URL on the known host.
  defp validate_station(<<_::binary-size(9)>> = s) do
    if String.match?(s, ~r/\A[A-Z0-9]{9}\z/), do: :ok, else: bad_station(s)
  end

  defp validate_station(s), do: bad_station(s)

  defp bad_station(s), do: {:error, {:unsupported_product, {:station, s}}}

  @doc """
  The transfer protocol for the daily station observation archive (`:ftp` on the
  ESA GSSC mirror).
  """
  @spec station_obs_protocol() :: :ftp
  def station_obs_protocol, do: :ftp

  # --- date arithmetic -----------------------------------------------------

  defp days_since_gps_epoch(%Date{} = date) do
    jdn(date.year, date.month, date.day) - @gps_epoch_jdn
  end

  defp pad3(n), do: n |> Integer.to_string() |> String.pad_leading(3, "0")

  defp jdn(year, month, day) do
    a = div(14 - month, 12)
    y = year + 4800 - a
    m = month + 12 * a - 3
    day + div(153 * m + 2, 5) + 365 * y + div(y, 4) - div(y, 100) + div(y, 400) - 32_045
  end

  @doc false
  def seconds_per_day, do: @seconds_per_day
end
