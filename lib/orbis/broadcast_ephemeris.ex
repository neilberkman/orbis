defmodule Orbis.BroadcastEphemeris do
  @moduledoc """
  A parsed RINEX broadcast-navigation product (GPS LNAV, Galileo I/NAV+F/NAV,
  BeiDou D1/D2).

  Holds the broadcast Keplerian elements and clock terms as a resource handle,
  the broadcast-ephemeris counterpart to `Orbis.SP3`. Pass a handle to
  `Orbis.PointPositioning.solve/4` to position from broadcast ephemeris instead
  of a precise SP3 product. The navigation file is parsed exactly once; the
  parsed product is held as a reference, not re-parsed per call.

  Parsing covers GPS, Galileo, and BeiDou records (including BeiDou
  geostationary satellites); other constellations in a mixed file are skipped.
  """

  alias Orbis.NIF

  @enforce_keys [:handle]
  defstruct [:handle]

  @type t :: %__MODULE__{handle: reference()}

  @doc """
  Parse a RINEX 3 navigation file from disk.

  Returns `{:ok, %Orbis.BroadcastEphemeris{}}` or `{:error, reason}`. The file
  is read and parsed once; the parsed product is held as a resource handle.
  """
  @spec load(String.t()) :: {:ok, t()} | {:error, term()}
  def load(path) when is_binary(path) do
    with {:ok, text} <- File.read(path) do
      parse(text)
    end
  end

  @doc """
  Like `load/1` but raises on failure.
  """
  @spec load!(String.t()) :: t()
  def load!(path) when is_binary(path) do
    case load(path) do
      {:ok, eph} ->
        eph

      {:error, reason} ->
        raise ArgumentError, "could not load RINEX NAV #{path}: #{inspect(reason)}"
    end
  end

  @doc """
  Parse an in-memory RINEX 3 navigation text buffer into a handle.
  """
  @spec parse(String.t()) :: {:ok, t()} | {:error, term()}
  def parse(text) when is_binary(text) do
    case NIF.broadcast_parse(text) do
      handle when is_reference(handle) -> {:ok, %__MODULE__{handle: handle}}
      {:error, _} = err -> err
    end
  rescue
    e in ErlangError -> {:error, e.original}
  end
end
