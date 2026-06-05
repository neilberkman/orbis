defmodule Orbis.GnssSignal.CA do
  @moduledoc """
  GPS L1 coarse/acquisition (C/A) code generation and correlation.

  The C/A code is the family of 1023-chip Gold codes broadcast by GPS
  satellites on the L1 carrier. Each satellite's code is the modulo-2 sum
  (XOR) of two 10-stage maximal-length linear feedback shift registers (LFSRs),
  G1 and G2, as defined in IS-GPS-200 (Section 3.3.2.3 and Figure 3-9):

    * G1 uses the polynomial `1 + x^3 + x^10` (feedback from stages 3 and 10).
    * G2 uses the polynomial `1 + x^2 + x^3 + x^6 + x^8 + x^9 + x^10`
      (feedback from stages 2, 3, 6, 8, 9, and 10).

  Both registers are initialized to the all-ones state and clocked at the
  1.023 Mcps chipping rate; the sequence repeats every 1023 chips (1 ms).

  The per-satellite code is `G1 XOR G2i`, where `G2i` is a delayed copy of the
  G2 sequence formed by the modulo-2 sum of two G2 register stages. The
  stage pair is fixed per PRN by the IS-GPS-200 Table 3-I code-phase-selection
  table (PRN 1 selects stages 2 and 6, PRN 2 selects 3 and 7, and so on);
  selecting a stage pair is equivalent to an integer cyclic delay of the G2
  m-sequence.

  ## Chip mapping (BPSK)

  Binary chips are mapped to bipolar `±1` values using the standard BPSK
  convention:

    * binary `0` maps to `+1`
    * binary `1` maps to `-1`

  With this mapping a code's circular autocorrelation peaks at `+1023` at zero
  lag, and the codes are nearly balanced (512 chips of `-1` and 511 of `+1`,
  so the chip sum is `-1`).

  Supported PRNs are the 32 GPS space-vehicle assignments, PRN 1 through 32.

  ## Examples

      iex> Orbis.GnssSignal.CA.code_length()
      1023

      iex> Orbis.GnssSignal.CA.chip_rate_hz()
      1_023_000

      iex> {:ok, chips} = Orbis.GnssSignal.CA.code(1)
      iex> length(chips)
      1023

      iex> Orbis.GnssSignal.CA.chip(1, 0)
      {:ok, -1}

      iex> Orbis.GnssSignal.CA.code(33)
      {:error, {:unsupported_prn, 33}}

  """

  @code_length 1023
  @chip_rate_hz 1_023_000

  # IS-GPS-200 Table 3-I code-phase-selection: PRN -> the pair of G2 register
  # stages (1-based) whose modulo-2 sum forms the delayed G2 output.
  @phase_select %{
    1 => {2, 6},
    2 => {3, 7},
    3 => {4, 8},
    4 => {5, 9},
    5 => {1, 9},
    6 => {2, 10},
    7 => {1, 8},
    8 => {2, 9},
    9 => {3, 10},
    10 => {2, 3},
    11 => {3, 4},
    12 => {5, 6},
    13 => {6, 7},
    14 => {7, 8},
    15 => {8, 9},
    16 => {9, 10},
    17 => {1, 4},
    18 => {2, 5},
    19 => {3, 6},
    20 => {4, 7},
    21 => {5, 8},
    22 => {6, 9},
    23 => {1, 3},
    24 => {4, 6},
    25 => {5, 7},
    26 => {6, 8},
    27 => {7, 9},
    28 => {8, 10},
    29 => {1, 6},
    30 => {2, 7},
    31 => {3, 8},
    32 => {4, 9}
  }

  # Precompute the bipolar (+/-1) code for every supported PRN at compile time
  # so runtime calls are pure lookups.
  @codes Map.new(@phase_select, fn {prn, taps} ->
           {prn, Orbis.GnssSignal.CA.Generator.bipolar_code(taps)}
         end)

  # Tuple form of each code for O(1) indexed chip access.
  @code_tuples Map.new(@codes, fn {prn, chips} -> {prn, List.to_tuple(chips)} end)

  @doc """
  Returns the number of chips in one C/A code period.

  ## Examples

      iex> Orbis.GnssSignal.CA.code_length()
      1023

  """
  @spec code_length() :: 1023
  def code_length, do: @code_length

  @doc """
  Returns the C/A chipping rate in hertz (1.023 Mcps).

  ## Examples

      iex> Orbis.GnssSignal.CA.chip_rate_hz()
      1_023_000

  """
  @spec chip_rate_hz() :: 1_023_000
  def chip_rate_hz, do: @chip_rate_hz

  @doc """
  Returns the 1023 bipolar (`±1`) C/A chips for a PRN.

  Binary `0` maps to `+1` and binary `1` maps to `-1`. Returns
  `{:error, {:unsupported_prn, prn}}` for any PRN outside 1..32.

  ## Examples

      iex> {:ok, chips} = Orbis.GnssSignal.CA.code(1)
      iex> length(chips)
      1023
      iex> Enum.take(chips, 4)
      [-1, -1, 1, 1]

  """
  @spec code(integer()) :: {:ok, [integer()]} | {:error, {:unsupported_prn, integer()}}
  def code(prn) do
    case Map.fetch(@codes, prn) do
      {:ok, chips} -> {:ok, chips}
      :error -> {:error, {:unsupported_prn, prn}}
    end
  end

  @doc """
  Returns the single bipolar (`±1`) chip at a 0-based, wrapping index.

  The index is taken modulo 1023, so negative indices and indices at or beyond
  the code length wrap around the period. Returns
  `{:error, {:unsupported_prn, prn}}` for any PRN outside 1..32.

  ## Examples

      iex> Orbis.GnssSignal.CA.chip(1, 0)
      {:ok, -1}

      iex> Orbis.GnssSignal.CA.chip(1, 1023)
      {:ok, -1}

  """
  @spec chip(integer(), integer()) ::
          {:ok, integer()} | {:error, {:unsupported_prn, integer()}}
  def chip(prn, index) when is_integer(index) do
    case Map.fetch(@code_tuples, prn) do
      {:ok, chips} ->
        {:ok, elem(chips, Integer.mod(index, @code_length))}

      :error ->
        {:error, {:unsupported_prn, prn}}
    end
  end

  @doc """
  Circular autocorrelation of a bipolar code over all 1023 lags.

  The result is the integer-valued sequence whose element at lag `k` is the
  sum of `code[i] * code[i + k mod n]`. For a C/A code the zero-lag value is
  1023 and every other value is one of `{-65, -1, 63}`.

  ## Examples

      iex> {:ok, code} = Orbis.GnssSignal.CA.code(1)
      iex> [peak | _] = Orbis.GnssSignal.CA.autocorrelation(code)
      iex> peak
      1023

  """
  @spec autocorrelation([integer()]) :: [integer()]
  def autocorrelation(code) when is_list(code), do: cross_correlation(code, code)

  @doc """
  Circular cross-correlation of two equal-length bipolar codes over all lags.

  The element at lag `k` is the sum of `code_a[i] * code_b[i + k mod n]`. For
  two distinct C/A codes every value is one of `{-65, -1, 63}`.
  """
  @spec cross_correlation([integer()], [integer()]) :: [integer()]
  def cross_correlation(code_a, code_b)
      when is_list(code_a) and is_list(code_b) and length(code_a) == length(code_b) do
    a = List.to_tuple(code_a)
    b = List.to_tuple(code_b)
    n = tuple_size(a)

    Enum.map(0..(n - 1), fn lag -> correlation_at(a, b, lag, n) end)
  end

  @doc """
  Single-lag circular correlation of two equal-length bipolar codes.

  Returns the integer sum of `code_a[i] * code_b[i + lag mod n]`.
  """
  @spec correlation_at([integer()], [integer()], integer()) :: integer()
  def correlation_at(code_a, code_b, lag)
      when is_list(code_a) and is_list(code_b) and length(code_a) == length(code_b) do
    a = List.to_tuple(code_a)
    b = List.to_tuple(code_b)
    correlation_at(a, b, lag, tuple_size(a))
  end

  defp correlation_at(a, b, lag, n) do
    Enum.reduce(0..(n - 1), 0, fn i, acc ->
      acc + elem(a, i) * elem(b, Integer.mod(i + lag, n))
    end)
  end
end
