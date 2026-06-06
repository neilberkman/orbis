defmodule Orbis.GNSS.Constellation do
  @moduledoc """
  GNSS constellation identity catalogs and validation helpers.

  This module is a data/catalog layer: it builds normalized satellite identity
  records from public sources and compares those records with GNSS products. It
  does not alter positioning solves or infer application-specific health rules.

  GPS is supported first. CelesTrak `gps-ops` OMM/JSON is the base source for
  current NORAD catalog ids and PRN assignments; the PRN is parsed from
  `OBJECT_NAME` and rendered as the SP3/RINEX id (`"G13"`). NAVCEN's GPS
  constellation page can be parsed and merged as an optional status overlay for
  SVN and NANU usability details.

  ## Examples

      iex> omms = [
      ...>   %{"OBJECT_NAME" => "GPS BIIF-8  (PRN 03)", "NORAD_CAT_ID" => 40294}
      ...> ]
      iex> {:ok, [record]} = Orbis.GNSS.Constellation.from_celestrak_omm(omms)
      iex> {record.system, record.prn, record.norad_id, record.sp3_id}
      {:gps, 3, 40294, "G03"}

      iex> record = %Orbis.GNSS.Constellation.Record{
      ...>   system: :gps,
      ...>   prn: 3,
      ...>   svn: nil,
      ...>   norad_id: 40294,
      ...>   sp3_id: "G03",
      ...>   active?: true,
      ...>   usable?: true,
      ...>   source: %{}
      ...> }
      iex> Orbis.GNSS.Constellation.to_csv([record])
      "prn,norad_cat_id,active,sp3_id\\n3,40294,true,G03\\n"

  Live fetching is available through `fetch_gps/1`, but tests and production
  pipelines that need deterministic behavior should use `from_celestrak_omm/1`,
  `parse_navcen_html/1`, and `merge_navcen/2` with cached bytes.
  """

  alias Orbis.GNSS.SP3

  @celestrak_gps_group "gps-ops"
  @celestrak_gps_url "https://celestrak.org/NORAD/elements/gp.php?GROUP=gps-ops&FORMAT=json"
  @navcen_gps_url "https://www.navcen.uscg.gov/gps-constellation?order=field_gps_prn&sort=asc"

  defmodule Record do
    @moduledoc """
    A normalized GNSS satellite identity record.

    `active?` means the satellite is present in the base source used to build
    the catalog. `usable?` is an advisory health/status flag; with the current
    GPS implementation it is `true` unless a merged NAVCEN row has an active
    NANU that marks the PRN unusable/decommissioned. The `source` map preserves
    source-specific fields such as CelesTrak object names and NAVCEN NANU text.
    """

    @enforce_keys [:system, :prn, :norad_id, :sp3_id, :active?, :usable?, :source]
    defstruct [:system, :prn, :svn, :norad_id, :sp3_id, :active?, :usable?, :source]

    @type t :: %__MODULE__{
            system: :gps,
            prn: pos_integer(),
            svn: pos_integer() | nil,
            norad_id: pos_integer(),
            sp3_id: String.t(),
            active?: boolean(),
            usable?: boolean(),
            source: map()
          }
  end

  defmodule NavcenStatus do
    @moduledoc """
    A parsed row from NAVCEN's GPS constellation status table.
    """

    @enforce_keys [:system, :prn, :svn, :usable?, :active_nanu?, :source]
    defstruct [
      :system,
      :prn,
      :svn,
      :usable?,
      :active_nanu?,
      :nanu_type,
      :nanu_subject,
      :plane,
      :slot,
      :block_type,
      :clock,
      :source
    ]

    @type t :: %__MODULE__{
            system: :gps,
            prn: pos_integer(),
            svn: pos_integer() | nil,
            usable?: boolean(),
            active_nanu?: boolean(),
            nanu_type: String.t() | nil,
            nanu_subject: String.t() | nil,
            plane: String.t() | nil,
            slot: String.t() | nil,
            block_type: String.t() | nil,
            clock: String.t() | nil,
            source: map()
          }
  end

  defmodule Validation do
    @moduledoc """
    Validation report for a GNSS constellation catalog.
    """

    @enforce_keys [
      :missing_sp3_ids,
      :duplicate_prns,
      :duplicate_norad_ids,
      :inactive_unusable_prns,
      :extra_sp3_ids
    ]
    defstruct [
      :missing_sp3_ids,
      :duplicate_prns,
      :duplicate_norad_ids,
      :inactive_unusable_prns,
      :extra_sp3_ids
    ]

    @type t :: %__MODULE__{
            missing_sp3_ids: [String.t()],
            duplicate_prns: [pos_integer()],
            duplicate_norad_ids: [pos_integer()],
            inactive_unusable_prns: [pos_integer()],
            extra_sp3_ids: [String.t()]
          }
  end

  defmodule Diff do
    @moduledoc """
    Change report between two GNSS constellation catalog snapshots.
    """

    @enforce_keys [
      :added,
      :removed,
      :norad_reassigned,
      :sp3_id_changed,
      :svn_changed,
      :activity_changed,
      :usability_changed
    ]
    defstruct [
      :added,
      :removed,
      :norad_reassigned,
      :sp3_id_changed,
      :svn_changed,
      :activity_changed,
      :usability_changed
    ]

    @type record_change(value) :: %{
            required(:system) => :gps,
            required(:prn) => pos_integer(),
            required(:from) => value,
            required(:to) => value
          }

    @type t :: %__MODULE__{
            added: [Record.t()],
            removed: [Record.t()],
            norad_reassigned: [record_change(pos_integer())],
            sp3_id_changed: [record_change(String.t())],
            svn_changed: [record_change(pos_integer() | nil)],
            activity_changed: [record_change(boolean())],
            usability_changed: [record_change(boolean())]
          }
  end

  @type error ::
          {:error, {:unsupported_system, term()}}
          | {:error, {:bad_celestrak_record, term(), map()}}
          | {:error, {:bad_navcen_html, term()}}
          | {:error, {:navcen_fetch_failed, term()}}
          | {:error, :req_not_available}
          | {:error, term()}

  @doc """
  Fetch the current GPS catalog from public sources.

  CelesTrak `gps-ops` OMM/JSON is always fetched and used as the base identity
  source. By default, NAVCEN's GPS constellation table is also fetched and merged
  as a status/SVN overlay.

  Options:

    * `:include_navcen` - fetch and merge NAVCEN status rows (default: `true`)
    * `:navcen_html` - use this already-fetched NAVCEN HTML instead of network
      access; implies `include_navcen: true`
    * `:timeout_ms` - NAVCEN request timeout when fetching the overlay

  Returns `{:ok, [Record.t()]}` or a tagged `{:error, reason}`.
  """
  @spec fetch_gps(keyword()) :: {:ok, [Record.t()]} | error()
  def fetch_gps(opts \\ []) do
    with {:ok, omms} <- Orbis.CelesTrak.fetch_omm(@celestrak_gps_group),
         {:ok, records} <- from_celestrak_omm(omms) do
      maybe_merge_navcen(records, opts)
    end
  end

  @doc """
  Build GPS records from CelesTrak `gps-ops` OMM/JSON maps.

  CelesTrak does not publish SVN in this feed, so records built from this source
  alone have `svn: nil`.
  """
  @spec from_celestrak_omm([map()]) :: {:ok, [Record.t()]} | error()
  def from_celestrak_omm(omms) when is_list(omms) do
    omms
    |> Enum.reduce_while({:ok, []}, fn omm, {:ok, acc} ->
      case record_from_omm(omm) do
        {:ok, record} -> {:cont, {:ok, [record | acc]}}
        {:error, _} = err -> {:halt, err}
      end
    end)
    |> case do
      {:ok, records} -> {:ok, records |> Enum.reverse() |> Enum.sort_by(& &1.prn)}
      {:error, _} = err -> err
    end
  end

  def from_celestrak_omm(other),
    do: {:error, {:bad_celestrak_record, :not_a_list, %{value: other}}}

  @doc """
  Parse NAVCEN's GPS constellation status HTML.

  The parser targets the Drupal table field classes used by NAVCEN's public GPS
  constellation page. It returns status rows keyed by PRN/SVN; merge them into
  CelesTrak records with `merge_navcen/2`.
  """
  @spec parse_navcen_html(String.t()) :: {:ok, [NavcenStatus.t()]} | error()
  def parse_navcen_html(html) when is_binary(html) do
    statuses =
      html
      |> then(&Regex.scan(~r/<tr\b[^>]*>(.*?)<\/tr>/si, &1, capture: :all_but_first))
      |> Enum.map(fn [row] -> row end)
      |> Enum.filter(&String.contains?(&1, "views-field-field-gps-prn"))
      |> Enum.filter(&String.contains?(&1, "<td"))
      |> Enum.reduce_while([], fn row, acc ->
        case navcen_status_from_row(row) do
          {:ok, status} -> {:cont, [status | acc]}
          {:error, _} = err -> {:halt, err}
        end
      end)

    cond do
      match?({:error, _}, statuses) ->
        statuses

      statuses == [] ->
        {:error, {:bad_navcen_html, :no_gps_rows}}

      true ->
        {:ok, statuses |> Enum.reverse() |> Enum.sort_by(& &1.prn)}
    end
  end

  def parse_navcen_html(_), do: {:error, {:bad_navcen_html, :not_binary}}

  @doc """
  Merge NAVCEN status rows into normalized GPS records by PRN.

  NAVCEN does not publish NORAD catalog ids, so CelesTrak remains the identity
  base. When a PRN is present in both sources, this fills `svn`, updates
  `usable?`, and adds a `:navcen` entry to the record's source metadata, but
  only when the CelesTrak object-name block type and NAVCEN block type are
  compatible. During PRN transitions NAVCEN can still carry a NANU for an older
  vehicle on the same PRN; those rows are recorded under `:navcen_conflict`
  rather than merged into the normalized identity.
  """
  @spec merge_navcen([Record.t()], [NavcenStatus.t()]) :: [Record.t()]
  def merge_navcen(records, statuses) when is_list(records) and is_list(statuses) do
    by_prn = Map.new(statuses, &{&1.prn, &1})

    records
    |> Enum.map(fn record ->
      case Map.fetch(by_prn, record.prn) do
        {:ok, status} -> merge_status(record, status)
        :error -> record
      end
    end)
    |> Enum.sort_by(& &1.prn)
  end

  @doc """
  Export records as the compact mapping CSV:

      prn,norad_cat_id,active,sp3_id

  The `active` column is `true` only when `record.active?` and `record.usable?`
  are both true.

  ## Options

    * `:booleans` - how the `active` column is rendered: `:lower` (default —
      `true`/`false`, the conventional CSV form) or `:title` (`True`/`False`,
      for a consumer such as pandas that reads the column as Python booleans).
  """
  @spec to_csv([Record.t()], keyword()) :: String.t()
  def to_csv(records, opts \\ []) when is_list(records) do
    booleans = Keyword.get(opts, :booleans, :lower)

    body =
      records
      |> Enum.sort_by(& &1.prn)
      |> Enum.map(fn record ->
        active = format_bool(operational?(record), booleans)
        "#{record.prn},#{record.norad_id},#{active},#{record.sp3_id}\n"
      end)

    IO.iodata_to_binary(["prn,norad_cat_id,active,sp3_id\n", body])
  end

  defp format_bool(value, :title), do: if(value, do: "True", else: "False")
  defp format_bool(value, _lower), do: to_string(value)

  @doc """
  Compare two constellation catalog snapshots by `{system, prn}` identity.

  The result separates added/removed PRNs from identity/status changes on a PRN
  that exists in both snapshots.

  `diff/2` assumes each input has at most one record per `{system, prn}`. Run
  `validate/1` on both snapshots first when ingesting external or hand-edited
  catalogs, and treat duplicate findings as malformed input rather than a
  constellation change.

      iex> previous = [
      ...>   %Orbis.GNSS.Constellation.Record{
      ...>     system: :gps,
      ...>     prn: 3,
      ...>     svn: 69,
      ...>     norad_id: 40294,
      ...>     sp3_id: "G03",
      ...>     active?: true,
      ...>     usable?: true,
      ...>     source: %{}
      ...>   }
      ...> ]
      iex> current = [%{hd(previous) | norad_id: 99999, active?: false}]
      iex> diff = Orbis.GNSS.Constellation.diff(previous, current)
      iex> [norad] = diff.norad_reassigned
      iex> {norad.system, norad.prn, norad.from, norad.to}
      {:gps, 3, 40294, 99999}
      iex> [activity] = diff.activity_changed
      iex> {activity.system, activity.prn, activity.from, activity.to}
      {:gps, 3, true, false}
      iex> Orbis.GNSS.Constellation.changed?(diff)
      true
  """
  @spec diff([Record.t()], [Record.t()]) :: Diff.t()
  def diff(previous, current) when is_list(previous) and is_list(current) do
    previous_by_key = Map.new(previous, &{record_key(&1), &1})
    current_by_key = Map.new(current, &{record_key(&1), &1})

    previous_keys = previous_by_key |> Map.keys() |> MapSet.new()
    current_keys = current_by_key |> Map.keys() |> MapSet.new()

    common_keys =
      previous_keys
      |> MapSet.intersection(current_keys)
      |> MapSet.to_list()
      |> Enum.sort_by(&sort_key/1)

    %Diff{
      added:
        current
        |> Enum.reject(&Map.has_key?(previous_by_key, record_key(&1)))
        |> Enum.sort_by(&record_sort_key/1),
      removed:
        previous
        |> Enum.reject(&Map.has_key?(current_by_key, record_key(&1)))
        |> Enum.sort_by(&record_sort_key/1),
      norad_reassigned: changed_field(common_keys, previous_by_key, current_by_key, :norad_id),
      sp3_id_changed: changed_field(common_keys, previous_by_key, current_by_key, :sp3_id),
      svn_changed: changed_field(common_keys, previous_by_key, current_by_key, :svn),
      activity_changed: changed_field(common_keys, previous_by_key, current_by_key, :active?),
      usability_changed: changed_field(common_keys, previous_by_key, current_by_key, :usable?)
    }
  end

  def diff(_previous, _current) do
    raise ArgumentError, "Orbis.GNSS.Constellation.diff/2 expects two record lists"
  end

  @doc """
  Validate catalog identity without an SP3 product.

  Reports duplicate PRNs, duplicate NORAD ids, and PRNs that are inactive or
  unusable according to the normalized records.
  """
  @spec validate([Record.t()]) :: Validation.t()
  def validate(records) when is_list(records), do: validation(records, nil)

  @doc """
  Validate catalog identity against SP3 satellite ids.

  The second argument may be a loaded `%Orbis.GNSS.SP3{}` or a list of SP3/RINEX
  satellite tokens. `missing_sp3_ids` reports active+usable catalog GPS ids that
  are absent from the product; `extra_sp3_ids` reports GPS ids present in the SP3
  product but absent from the active+usable catalog.
  """
  @spec validate_sp3([Record.t()], SP3.t() | [String.t()]) :: Validation.t()
  def validate_sp3(records, %SP3{} = sp3), do: validate_sp3(records, SP3.satellite_ids(sp3))

  def validate_sp3(records, sp3_ids) when is_list(records) and is_list(sp3_ids) do
    validation(records, sp3_ids)
  end

  @doc """
  Returns `true` when a validation report has no findings.
  """
  @spec valid?(Validation.t()) :: boolean()
  def valid?(%Validation{} = report) do
    report.missing_sp3_ids == [] and report.duplicate_prns == [] and
      report.duplicate_norad_ids == [] and report.inactive_unusable_prns == [] and
      report.extra_sp3_ids == []
  end

  @doc """
  Returns `true` when a constellation diff has any findings.
  """
  @spec changed?(Diff.t()) :: boolean()
  def changed?(%Diff{} = diff) do
    diff.added != [] or diff.removed != [] or diff.norad_reassigned != [] or
      diff.sp3_id_changed != [] or diff.svn_changed != [] or diff.activity_changed != [] or
      diff.usability_changed != []
  end

  @doc """
  Validate against SP3 satellite ids and raise unless the catalog is clean.

  A build-time gate: returns `:ok` when `validate_sp3/2` reports no findings,
  otherwise raises `ArgumentError` describing them — for example an
  active-and-usable catalog PRN that is missing from a current SP3 product (a
  stale-active satellite). Intended for catalog-build / automation steps, not the
  positioning runtime.
  """
  @spec validate_sp3!([Record.t()], SP3.t() | [String.t()]) :: :ok
  def validate_sp3!(records, sp3_or_ids) do
    report = validate_sp3(records, sp3_or_ids)

    if valid?(report) do
      :ok
    else
      raise ArgumentError, "GNSS catalog failed SP3 validation — " <> describe_findings(report)
    end
  end

  defp describe_findings(%Validation{} = report) do
    [
      {"missing_sp3_ids", report.missing_sp3_ids},
      {"extra_sp3_ids", report.extra_sp3_ids},
      {"duplicate_prns", report.duplicate_prns},
      {"duplicate_norad_ids", report.duplicate_norad_ids},
      {"inactive_unusable_prns", report.inactive_unusable_prns}
    ]
    |> Enum.reject(fn {_label, values} -> values == [] end)
    |> Enum.map_join("; ", fn {label, values} -> "#{label}: #{inspect(values)}" end)
  end

  @doc """
  Render the canonical SP3/RINEX satellite token for a supported GNSS PRN.

      iex> Orbis.GNSS.Constellation.sp3_id(:gps, 7)
      "G07"
  """
  @spec sp3_id(:gps, pos_integer()) :: String.t()
  def sp3_id(:gps, prn) when is_integer(prn) and prn > 0, do: "G" <> pad_prn(prn)

  # --- source fetching -----------------------------------------------------

  defp maybe_merge_navcen(records, opts) do
    cond do
      html = Keyword.get(opts, :navcen_html) ->
        with {:ok, statuses} <- parse_navcen_html(html) do
          {:ok, merge_navcen(records, statuses)}
        end

      Keyword.get(opts, :include_navcen, true) ->
        with {:ok, statuses} <- fetch_navcen_statuses(opts) do
          {:ok, merge_navcen(records, statuses)}
        end

      true ->
        {:ok, records}
    end
  end

  defp fetch_navcen_statuses(opts) do
    if Code.ensure_loaded?(Req) and function_exported?(Req, :get, 2) do
      timeout = Keyword.get(opts, :timeout_ms, 30_000)

      case Req.get(@navcen_gps_url,
             receive_timeout: timeout,
             connect_options: [timeout: timeout]
           ) do
        {:ok, %{status: 200, body: body}} when is_binary(body) ->
          parse_navcen_html(body)

        {:ok, %{status: 200, body: body}} ->
          {:error, {:navcen_fetch_failed, {:unexpected_body, body}}}

        {:ok, %{status: status}} ->
          {:error, {:navcen_fetch_failed, {:http_status, status}}}

        {:error, reason} ->
          {:error, {:navcen_fetch_failed, reason}}
      end
    else
      {:error, :req_not_available}
    end
  end

  # --- CelesTrak parsing ---------------------------------------------------

  defp record_from_omm(%{} = omm) do
    with {:ok, prn} <- prn_from_object_name(omm["OBJECT_NAME"]),
         {:ok, norad_id} <- int_field(omm, "NORAD_CAT_ID") do
      {:ok,
       %Record{
         system: :gps,
         prn: prn,
         svn: nil,
         norad_id: norad_id,
         sp3_id: sp3_id(:gps, prn),
         active?: true,
         usable?: true,
         source: %{
           celestrak: %{
             group: @celestrak_gps_group,
             url: @celestrak_gps_url,
             object_name: omm["OBJECT_NAME"],
             object_id: omm["OBJECT_ID"],
             epoch: omm["EPOCH"],
             block_type: block_type_from_object_name(omm["OBJECT_NAME"])
           }
         }
       }}
    else
      {:error, reason} -> {:error, {:bad_celestrak_record, reason, omm}}
    end
  end

  defp record_from_omm(other), do: {:error, {:bad_celestrak_record, :not_a_map, %{value: other}}}

  defp prn_from_object_name(name) when is_binary(name) do
    case Regex.run(~r/\(PRN\s*0*([0-9]{1,3})\)/i, name) do
      [_, prn] -> parse_positive_int(prn, :prn)
      _ -> {:error, {:missing_prn, name}}
    end
  end

  defp prn_from_object_name(other), do: {:error, {:missing_prn, other}}

  defp block_type_from_object_name(name) when is_binary(name) do
    cond do
      Regex.match?(~r/\bBIIR-?M\b|\bBIIRM\b/i, name) -> "IIR-M"
      Regex.match?(~r/\bBIII\b/i, name) -> "III"
      Regex.match?(~r/\bBIIF\b/i, name) -> "IIF"
      Regex.match?(~r/\bBIIR\b/i, name) -> "IIR"
      true -> nil
    end
  end

  defp block_type_from_object_name(_), do: nil

  defp int_field(map, key) do
    case Map.fetch(map, key) do
      {:ok, value} -> parse_positive_int(value, key)
      :error -> {:error, {:missing_field, key}}
    end
  end

  # --- NAVCEN parsing ------------------------------------------------------

  defp navcen_status_from_row(row) do
    with {:ok, prn} <- navcen_int(row, "gps-prn"),
         {:ok, svn} <- navcen_optional_int(row, "gps-svn") do
      nanu_type = navcen_text(row, "nanu-type")
      active_nanu? = navcen_active?(row)
      usable? = not (active_nanu? and unusable_nanu_type?(nanu_type))
      subject = navcen_text(row, "nanu-subject")
      plane = navcen_text(row, "gps-con-plane")
      slot = navcen_text(row, "gps-con-slot")
      block_type = navcen_text(row, "gps-con-block-type")
      clock = navcen_text(row, "gps-con-clock")

      {:ok,
       %NavcenStatus{
         system: :gps,
         prn: prn,
         svn: svn,
         usable?: usable?,
         active_nanu?: active_nanu?,
         nanu_type: blank_to_nil(nanu_type),
         nanu_subject: blank_to_nil(subject),
         plane: blank_to_nil(plane),
         slot: blank_to_nil(slot),
         block_type: blank_to_nil(block_type),
         clock: blank_to_nil(clock),
         source: %{
           navcen: %{
             url: @navcen_gps_url,
             svn: svn,
             block_type: blank_to_nil(block_type),
             plane: blank_to_nil(plane),
             slot: blank_to_nil(slot),
             clock: blank_to_nil(clock),
             nanu_type: blank_to_nil(nanu_type),
             nanu_subject: blank_to_nil(subject),
             active_nanu?: active_nanu?
           }
         }
       }}
    else
      {:error, reason} -> {:error, {:bad_navcen_html, reason}}
    end
  end

  defp navcen_int(row, field) do
    row
    |> navcen_text(field)
    |> parse_positive_int(field)
  end

  defp navcen_optional_int(row, field) do
    case navcen_text(row, field) do
      "" -> {:ok, nil}
      text -> parse_positive_int(text, field)
    end
  end

  defp navcen_text(row, field) do
    pattern =
      ~r/<td\b[^>]*views-field-field-#{Regex.escape(field)}[^>]*>(.*?)<\/td>/si

    case Regex.run(pattern, row) do
      [_, text] -> clean_html(text)
      _ -> ""
    end
  end

  defp navcen_active?(row) do
    case Regex.run(~r/<td\b[^>]*nanu-active-check[^>]*>(.*?)<\/td>/si, row) do
      [_, text] -> clean_html(text) == "1"
      _ -> false
    end
  end

  defp unusable_nanu_type?(nil), do: false
  defp unusable_nanu_type?(""), do: false

  defp unusable_nanu_type?(type) do
    type = type |> String.trim() |> String.upcase()
    type in ~w(UNUSABLE DECOM FCSTDV FCSTMX FCSTEXTD)
  end

  defp clean_html(text) do
    text
    |> String.replace(~r/<[^>]+>/, "")
    |> html_unescape()
    |> String.replace(~r/\s+/u, " ")
    |> String.trim()
  end

  defp html_unescape(text) do
    text
    |> String.replace("&amp;", "&")
    |> String.replace("&lt;", "<")
    |> String.replace("&gt;", ">")
    |> String.replace("&quot;", "\"")
    |> String.replace("&#39;", "'")
    |> String.replace("&nbsp;", " ")
  end

  defp merge_status(%Record{} = record, %NavcenStatus{} = status) do
    if navcen_compatible?(record, status) do
      %{
        record
        | svn: status.svn,
          usable?: status.usable?,
          source: Map.merge(record.source, status.source)
      }
    else
      %{
        record
        | source:
            Map.put(
              record.source,
              :navcen_conflict,
              status.source.navcen
              |> Map.put(:svn, status.svn)
              |> Map.put(:block_type, status.block_type)
            )
      }
    end
  end

  defp navcen_compatible?(%Record{} = record, %NavcenStatus{} = status) do
    celestrak_block = get_in(record.source, [:celestrak, :block_type])
    navcen_block = normalize_block_type(status.block_type)

    is_nil(celestrak_block) or is_nil(navcen_block) or celestrak_block == navcen_block
  end

  defp normalize_block_type(nil), do: nil
  defp normalize_block_type(""), do: nil
  defp normalize_block_type(type), do: type |> String.trim() |> String.upcase()

  # --- validation ----------------------------------------------------------

  defp validation(records, nil) do
    %Validation{
      missing_sp3_ids: [],
      duplicate_prns: duplicates(records, & &1.prn),
      duplicate_norad_ids: duplicates(records, & &1.norad_id),
      inactive_unusable_prns: inactive_unusable_prns(records),
      extra_sp3_ids: []
    }
  end

  defp validation(records, sp3_ids) do
    base = validation(records, nil)

    catalog_ids =
      records
      |> Enum.filter(&operational?/1)
      |> MapSet.new(&String.upcase(&1.sp3_id))

    sp3_ids =
      sp3_ids
      |> Enum.map(&String.upcase/1)
      |> Enum.filter(&String.starts_with?(&1, "G"))
      |> MapSet.new()

    %{
      base
      | missing_sp3_ids:
          catalog_ids |> MapSet.difference(sp3_ids) |> MapSet.to_list() |> Enum.sort(),
        extra_sp3_ids:
          sp3_ids |> MapSet.difference(catalog_ids) |> MapSet.to_list() |> Enum.sort()
    }
  end

  defp duplicates(records, key_fun) do
    records
    |> Enum.group_by(key_fun)
    |> Enum.filter(fn {_key, rows} -> length(rows) > 1 end)
    |> Enum.map(fn {key, _rows} -> key end)
    |> Enum.sort()
  end

  defp inactive_unusable_prns(records) do
    records
    |> Enum.reject(&operational?/1)
    |> Enum.map(& &1.prn)
    |> Enum.uniq()
    |> Enum.sort()
  end

  defp operational?(%Record{} = record), do: record.active? and record.usable?

  # --- catalog diff --------------------------------------------------------

  defp changed_field(keys, previous_by_key, current_by_key, field) do
    Enum.flat_map(keys, fn key ->
      previous = Map.fetch!(previous_by_key, key)
      current = Map.fetch!(current_by_key, key)
      from = Map.fetch!(previous, field)
      to = Map.fetch!(current, field)

      if from == to do
        []
      else
        {system, prn} = key
        [%{system: system, prn: prn, from: from, to: to}]
      end
    end)
  end

  defp record_key(%Record{} = record), do: {record.system, record.prn}
  defp record_sort_key(%Record{} = record), do: sort_key(record_key(record))
  defp sort_key({system, prn}), do: {system, prn}

  # --- scalar helpers ------------------------------------------------------

  defp parse_positive_int(value, _field) when is_integer(value) and value > 0, do: {:ok, value}

  defp parse_positive_int(value, field) when is_binary(value) do
    value = String.trim(value)

    case Integer.parse(value) do
      {int, ""} when int > 0 -> {:ok, int}
      _ -> {:error, {:bad_integer, field, value}}
    end
  end

  defp parse_positive_int(value, field), do: {:error, {:bad_integer, field, value}}

  defp pad_prn(prn) when prn < 10, do: "0#{prn}"
  defp pad_prn(prn), do: Integer.to_string(prn)

  defp blank_to_nil(nil), do: nil
  defp blank_to_nil(""), do: nil
  defp blank_to_nil(value), do: value
end
