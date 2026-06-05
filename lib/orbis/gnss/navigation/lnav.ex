defmodule Orbis.GNSS.Navigation.LNAV do
  @moduledoc """
  GPS L1 C/A LNAV navigation message synthesis and decoding (subframes 1-3).

  The legacy navigation (LNAV) message is the data stream modulated onto the
  GPS L1 C/A signal at 50 bits per second. Its structure is defined in
  IS-GPS-200 (Section 20.3): the message is organized into 1500-bit *frames*,
  each frame being five 300-bit *subframes*, and each subframe being ten 30-bit
  *words*. Every word carries 24 source data bits (most significant first)
  followed by 6 parity bits.

  This module covers the clock and ephemeris subframes:

    * Subframe 1 - SV clock correction and health (IS-GPS-200 Table 20-I).
    * Subframe 2 - first half of the ephemeris (IS-GPS-200 Table 20-II).
    * Subframe 3 - second half of the ephemeris (IS-GPS-200 Table 20-III).

  The first word of every subframe is the telemetry (TLM) word; the second is
  the hand-over word (HOW). Both are described in IS-GPS-200 Section 20.3.3.

  ## Words and bits

  A subframe is represented as a flat list of 300 bits (`0`/`1`), most
  significant bit first, with the ten words concatenated in transmission
  order. `word_length/0` is 30 and `subframe_length/0` is 300.

  ## Parameters and units

  `encode/2` and `decode/1` exchange an `Orbis.GNSS.Navigation.LNAV.Ephemeris` struct
  whose fields hold the natural engineering-unit values (the products of the
  transmitted integers and their IS-GPS-200 scale factors). See that struct's
  documentation for the per-field units. Angular ephemeris quantities use
  *semicircles* (and semicircles/second), the harmonic correction terms use
  radians, distances use meters, and clock/time quantities use seconds, exactly
  as tabulated in IS-GPS-200.

  ## Parity

  The 6 parity bits of each word are produced by the (32, 26) Hamming code of
  IS-GPS-200 Section 20.3.5.2 and Table 20-XIV, including the rule that the two
  trailing parity bits of the previous word (`D29*`, `D30*`) feed the current
  word and that `D30*` complements the 24 transmitted data bits. The last two
  data bits of the HOW and of word 10 are solved so that those words'
  `D29`/`D30` parity bits are zero, per IS-GPS-200 Section 20.3.3.2. At the
  start of each subframe the previous parity bits are seeded to zero, producing
  self-consistent stand-alone subframes.

  ## Examples

      iex> Orbis.GNSS.Navigation.LNAV.word_length()
      30

      iex> Orbis.GNSS.Navigation.LNAV.subframe_length()
      300

      iex> Orbis.GNSS.Navigation.LNAV.preamble()
      139

  """

  alias Orbis.GNSS.Navigation.LNAV.Ephemeris

  @word_length 30
  @subframe_length 300
  # IS-GPS-200 Section 20.3.3.1: TLM preamble is the bit pattern 1000 1011.
  @preamble 0b10001011

  @two_pow_4 16
  @two_pow_m5 :math.pow(2, -5)
  @two_pow_m19 :math.pow(2, -19)
  @two_pow_m29 :math.pow(2, -29)
  @two_pow_m31 :math.pow(2, -31)
  @two_pow_m33 :math.pow(2, -33)
  @two_pow_m43 :math.pow(2, -43)
  @two_pow_m55 :math.pow(2, -55)

  @doc """
  Bit length of a single LNAV word (IS-GPS-200 Section 20.3.2).

  ## Examples

      iex> Orbis.GNSS.Navigation.LNAV.word_length()
      30

  """
  @spec word_length() :: 30
  def word_length, do: @word_length

  @doc """
  Bit length of a single LNAV subframe (IS-GPS-200 Section 20.3.2).

  ## Examples

      iex> Orbis.GNSS.Navigation.LNAV.subframe_length()
      300

  """
  @spec subframe_length() :: 300
  def subframe_length, do: @subframe_length

  @doc """
  The 8-bit TLM preamble `1000 1011` as an integer (IS-GPS-200 Section 20.3.3.1).

  ## Examples

      iex> Orbis.GNSS.Navigation.LNAV.preamble()
      139

  """
  @spec preamble() :: 139
  def preamble, do: @preamble

  @doc """
  Extracts the 17-bit time-of-week count from a hand-over word.

  Accepts either a 30-bit HOW word or a full 300-bit subframe (whose word 2 is
  the HOW). The returned value is the truncated Z-count carried in the HOW
  (units of 6 seconds), per IS-GPS-200 Section 20.3.3.2.

  ## Examples

      iex> {:ok, sfs} = Orbis.GNSS.Navigation.LNAV.encode(Orbis.GNSS.Navigation.LNAV.Ephemeris.example(), tow: 12345)
      iex> Orbis.GNSS.Navigation.LNAV.tow(sfs[1])
      12345

  """
  @spec tow([0 | 1]) :: non_neg_integer()
  def tow(bits) when is_list(bits) do
    how = how_word(bits)
    bits_to_uint(Enum.slice(how, 0, 17))
  end

  @doc """
  Extracts the 3-bit subframe ID from a hand-over word.

  Accepts a 30-bit HOW word or a full 300-bit subframe. Returns the subframe
  identifier carried in HOW bits 20-22 (IS-GPS-200 Section 20.3.3.2).

  ## Examples

      iex> {:ok, sfs} = Orbis.GNSS.Navigation.LNAV.encode(Orbis.GNSS.Navigation.LNAV.Ephemeris.example(), tow: 0)
      iex> Orbis.GNSS.Navigation.LNAV.subframe_id(sfs[2])
      2

  """
  @spec subframe_id([0 | 1]) :: 1..5
  def subframe_id(bits) when is_list(bits) do
    how = how_word(bits)
    bits_to_uint(Enum.slice(how, 19, 3))
  end

  defp how_word(bits) do
    case length(bits) do
      @word_length -> bits
      @subframe_length -> Enum.slice(bits, @word_length, @word_length)
      _ -> raise ArgumentError, "expected a 30-bit word or 300-bit subframe"
    end
  end

  @doc """
  Computes the 6 parity bits of a word (IS-GPS-200 Table 20-XIV).

  `data24` is the list of 24 *source* data bits (most significant first, before
  the `D30*` complementation). `d29_prev` and `d30_prev` are the two trailing
  parity bits of the previous word. Returns `[D25, D26, D27, D28, D29, D30]`.
  """
  @spec parity([0 | 1], 0 | 1, 0 | 1) :: [0 | 1]
  def parity(data24, d29_prev, d30_prev) when is_list(data24) and length(data24) == 24 do
    d = List.to_tuple([nil | data24])
    # IS-GPS-200 Table 20-XIV. Index d[n] is the n-th source data bit.
    d25 =
      xor([
        d29_prev,
        e(d, 1),
        e(d, 2),
        e(d, 3),
        e(d, 5),
        e(d, 6),
        e(d, 10),
        e(d, 11),
        e(d, 12),
        e(d, 13),
        e(d, 14),
        e(d, 17),
        e(d, 18),
        e(d, 20),
        e(d, 23)
      ])

    d26 =
      xor([
        d30_prev,
        e(d, 2),
        e(d, 3),
        e(d, 4),
        e(d, 6),
        e(d, 7),
        e(d, 11),
        e(d, 12),
        e(d, 13),
        e(d, 14),
        e(d, 15),
        e(d, 18),
        e(d, 19),
        e(d, 21),
        e(d, 24)
      ])

    d27 =
      xor([
        d29_prev,
        e(d, 1),
        e(d, 3),
        e(d, 4),
        e(d, 5),
        e(d, 7),
        e(d, 8),
        e(d, 12),
        e(d, 13),
        e(d, 14),
        e(d, 15),
        e(d, 16),
        e(d, 19),
        e(d, 20),
        e(d, 22)
      ])

    d28 =
      xor([
        d30_prev,
        e(d, 2),
        e(d, 4),
        e(d, 5),
        e(d, 6),
        e(d, 8),
        e(d, 9),
        e(d, 13),
        e(d, 14),
        e(d, 15),
        e(d, 16),
        e(d, 17),
        e(d, 20),
        e(d, 21),
        e(d, 23)
      ])

    d29 =
      xor([
        d30_prev,
        e(d, 1),
        e(d, 3),
        e(d, 5),
        e(d, 6),
        e(d, 7),
        e(d, 9),
        e(d, 10),
        e(d, 14),
        e(d, 15),
        e(d, 16),
        e(d, 17),
        e(d, 18),
        e(d, 21),
        e(d, 22),
        e(d, 24)
      ])

    d30 =
      xor([
        d29_prev,
        e(d, 3),
        e(d, 5),
        e(d, 6),
        e(d, 8),
        e(d, 9),
        e(d, 10),
        e(d, 11),
        e(d, 13),
        e(d, 15),
        e(d, 19),
        e(d, 22),
        e(d, 23),
        e(d, 24)
      ])

    [d25, d26, d27, d28, d29, d30]
  end

  defp e(tuple, n), do: elem(tuple, n)
  defp xor(bits), do: Enum.reduce(bits, 0, &Bitwise.bxor/2)

  @doc """
  Verifies the parity of a single 30-bit word.

  `word30` is the 30-bit word as transmitted (data bits possibly complemented
  by `D30*`, followed by 6 received parity bits). `d29_prev`/`d30_prev` are the
  previous word's trailing parity bits. Returns `true` when the recomputed
  parity matches the received parity.
  """
  @spec parity_valid?([0 | 1], 0 | 1, 0 | 1) :: boolean()
  def parity_valid?(word30, d29_prev, d30_prev) when is_list(word30) and length(word30) == 30 do
    {transmitted_data, received_parity} = Enum.split(word30, 24)
    source = Enum.map(transmitted_data, &Bitwise.bxor(&1, d30_prev))
    parity(source, d29_prev, d30_prev) == received_parity
  end

  @doc """
  Encodes clock and ephemeris parameters into LNAV subframes 1-3.

  Returns `{:ok, %{1 => bits, 2 => bits, 3 => bits}}` where each value is a flat
  list of 300 bits (most significant first). Out-of-range parameters yield
  `{:error, {:out_of_range, field, value}}`; this function never raises on bad
  input.

  ## Options

    * `:tow` - the 17-bit time-of-week count placed in each HOW (0..131071,
      default 0).
    * `:alert` - HOW alert flag (`0`/`1`, default 0).
    * `:anti_spoof` - HOW anti-spoof flag (`0`/`1`, default 0).
    * `:integrity` - TLM integrity status flag (`0`/`1`, default 0).
    * `:tlm_message` - 14-bit TLM message field (default 0).

  """
  @spec encode(Ephemeris.t(), keyword()) ::
          {:ok, %{1 => [0 | 1], 2 => [0 | 1], 3 => [0 | 1]}} | {:error, term()}
  def encode(%Ephemeris{} = params, opts \\ []) do
    tow = Keyword.get(opts, :tow, 0)
    alert = Keyword.get(opts, :alert, 0)
    anti_spoof = Keyword.get(opts, :anti_spoof, 0)
    integrity = Keyword.get(opts, :integrity, 0)
    tlm_message = Keyword.get(opts, :tlm_message, 0)

    with :ok <- validate_uint(:tow, tow, 17),
         :ok <- validate_uint(:alert, alert, 1),
         :ok <- validate_uint(:anti_spoof, anti_spoof, 1),
         :ok <- validate_uint(:integrity, integrity, 1),
         :ok <- validate_uint(:tlm_message, tlm_message, 14),
         {:ok, w1} <- subframe1_words(params),
         {:ok, w2} <- subframe2_words(params),
         {:ok, w3} <- subframe3_words(params) do
      tlm = tlm_data(tlm_message, integrity)

      sf1 = assemble_subframe([tlm, how_data(tow, alert, anti_spoof, 1, :solve) | w1])
      sf2 = assemble_subframe([tlm, how_data(tow, alert, anti_spoof, 2, :solve) | w2])
      sf3 = assemble_subframe([tlm, how_data(tow, alert, anti_spoof, 3, :solve) | w3])

      {:ok, %{1 => sf1, 2 => sf2, 3 => sf3}}
    end
  end

  @doc """
  Decodes LNAV subframes 1-3 back into the engineering-unit parameter struct.

  Accepts `%{1 => bits, 2 => bits, 3 => bits}` of 300-bit subframes. Parity is
  verified on all 30 words first; a failure returns
  `{:error, {:parity_failed, subframe, word}}` (word is 1-based). On success
  returns `{:ok, %Orbis.GNSS.Navigation.LNAV.Ephemeris{}}`.
  """
  @spec decode(%{1 => [0 | 1], 2 => [0 | 1], 3 => [0 | 1]}) ::
          {:ok, Ephemeris.t()} | {:error, term()}
  def decode(%{1 => sf1, 2 => sf2, 3 => sf3}) do
    with :ok <- verify_subframe(sf1, 1),
         :ok <- verify_subframe(sf2, 2),
         :ok <- verify_subframe(sf3, 3) do
      w1 = source_words(sf1)
      w2 = source_words(sf2)
      w3 = source_words(sf3)

      params =
        %Ephemeris{}
        |> decode_subframe1(w1)
        |> decode_subframe2(w2)
        |> decode_subframe3(w3)

      {:ok, params}
    end
  end

  # --- Subframe 1 (clock/health), IS-GPS-200 Table 20-I -----------------------

  defp subframe1_words(p) do
    iodc = p.iodc || 0

    fields = [
      {:week_number, p.week_number, :uint, 10},
      {:l2_code, p.l2_code || 0, :uint, 2},
      {:ura_index, p.ura_index, :uint, 4},
      {:sv_health, p.sv_health, :uint, 6},
      {:iodc, iodc, :uint, 10},
      {:tgd, p.tgd, :sint, 8, @two_pow_m31},
      {:toc, p.toc, :uint, 16, @two_pow_4},
      {:af2, p.af2, :sint, 8, @two_pow_m55},
      {:af1, p.af1, :sint, 16, @two_pow_m43},
      {:af0, p.af0, :sint, 22, @two_pow_m31}
    ]

    with :ok <- validate_fields(fields) do
      iodc_msb = Bitwise.bsr(iodc, 8) |> Bitwise.band(0x3)
      iodc_lsb = Bitwise.band(iodc, 0xFF)

      l2_p_data_flag = p.l2_p_data_flag || 0

      # Word 3: WN(10) L2(2) URA(4) health(6) IODC-MSB(2)
      word3 =
        pack_uint(p.week_number, 10) ++
          pack_uint(p.l2_code || 0, 2) ++
          pack_uint(p.ura_index, 4) ++
          pack_uint(p.sv_health, 6) ++ pack_uint(iodc_msb, 2)

      # Word 4: L2-P data flag (bit 1) + 23 reserved bits.
      word4 = pack_uint(l2_p_data_flag, 1) ++ zeros(23)
      # Words 5, 6: reserved.
      word5 = zeros(24)
      word6 = zeros(24)
      # Word 7: 16 reserved bits then TGD(8).
      word7 = zeros(16) ++ pack_sint(p.tgd, 8, @two_pow_m31)
      # Word 8: IODC-LSB(8) toc(16).
      word8 = pack_uint(iodc_lsb, 8) ++ pack_uint_scaled(p.toc, 16, @two_pow_4)
      # Word 9: af2(8) af1(16).
      word9 = pack_sint(p.af2, 8, @two_pow_m55) ++ pack_sint(p.af1, 16, @two_pow_m43)
      # Word 10: af0(22) + 2 solved bits.
      word10 = pack_sint(p.af0, 22, @two_pow_m31)

      {:ok, [word3, word4, word5, word6, word7, word8, word9, {:solve, word10}]}
    end
  end

  defp decode_subframe1(p, w) do
    [word3, _w4, _w5, _w6, word7, word8, word9, word10] = w

    week_number = bits_to_uint(slice(word3, 1, 10))
    l2_code = bits_to_uint(slice(word3, 11, 2))
    ura_index = bits_to_uint(slice(word3, 13, 4))
    sv_health = bits_to_uint(slice(word3, 17, 6))
    iodc_msb = bits_to_uint(slice(word3, 23, 2))

    tgd = unpack_sint(slice(word7, 17, 8), @two_pow_m31)
    iodc_lsb = bits_to_uint(slice(word8, 1, 8))
    toc = unpack_uint_scaled(slice(word8, 9, 16), @two_pow_4)
    af2 = unpack_sint(slice(word9, 1, 8), @two_pow_m55)
    af1 = unpack_sint(slice(word9, 9, 16), @two_pow_m43)
    af0 = unpack_sint(slice(word10, 1, 22), @two_pow_m31)

    %{
      p
      | week_number: week_number,
        l2_code: l2_code,
        ura_index: ura_index,
        sv_health: sv_health,
        iodc: Bitwise.bor(Bitwise.bsl(iodc_msb, 8), iodc_lsb),
        tgd: tgd,
        toc: toc,
        af2: af2,
        af1: af1,
        af0: af0
    }
  end

  # --- Subframe 2 (ephemeris part 1), IS-GPS-200 Table 20-II ------------------

  defp subframe2_words(p) do
    fields = [
      {:iode, p.iode, :uint, 8},
      {:crs, p.crs, :sint, 16, @two_pow_m5},
      {:delta_n, p.delta_n, :sint, 16, @two_pow_m43},
      {:m0, p.m0, :sint, 32, @two_pow_m31},
      {:cuc, p.cuc, :sint, 16, @two_pow_m29},
      {:eccentricity, p.eccentricity, :uint, 32, @two_pow_m33},
      {:cus, p.cus, :sint, 16, @two_pow_m29},
      {:sqrt_a, p.sqrt_a, :uint, 32, @two_pow_m19},
      {:toe, p.toe, :uint, 16, @two_pow_4},
      {:fit_interval_flag, p.fit_interval_flag || 0, :uint, 1},
      {:aodo, p.aodo || 0, :uint, 5}
    ]

    with :ok <- validate_fields(fields) do
      m0 = pack_sint(p.m0, 32, @two_pow_m31)
      ecc = pack_uint_scaled(p.eccentricity, 32, @two_pow_m33)
      sqrt_a = pack_uint_scaled(p.sqrt_a, 32, @two_pow_m19)

      # Word 3: IODE(8) Crs(16).
      word3 = pack_uint(p.iode, 8) ++ pack_sint(p.crs, 16, @two_pow_m5)
      # Word 4: Delta-n(16) M0-MSB(8).
      word4 = pack_sint(p.delta_n, 16, @two_pow_m43) ++ Enum.slice(m0, 0, 8)
      # Word 5: M0-LSB(24).
      word5 = Enum.slice(m0, 8, 24)
      # Word 6: Cuc(16) e-MSB(8).
      word6 = pack_sint(p.cuc, 16, @two_pow_m29) ++ Enum.slice(ecc, 0, 8)
      # Word 7: e-LSB(24).
      word7 = Enum.slice(ecc, 8, 24)
      # Word 8: Cus(16) sqrtA-MSB(8).
      word8 = pack_sint(p.cus, 16, @two_pow_m29) ++ Enum.slice(sqrt_a, 0, 8)
      # Word 9: sqrtA-LSB(24).
      word9 = Enum.slice(sqrt_a, 8, 24)
      # Word 10: toe(16) fit(1) AODO(5) + 2 solved bits.
      word10 =
        pack_uint_scaled(p.toe, 16, @two_pow_4) ++
          pack_uint(p.fit_interval_flag || 0, 1) ++ pack_uint(p.aodo || 0, 5)

      {:ok, [word3, word4, word5, word6, word7, word8, word9, {:solve, word10}]}
    end
  end

  defp decode_subframe2(p, w) do
    [word3, word4, word5, word6, word7, word8, word9, word10] = w

    iode = bits_to_uint(slice(word3, 1, 8))
    crs = unpack_sint(slice(word3, 9, 16), @two_pow_m5)
    delta_n = unpack_sint(slice(word4, 1, 16), @two_pow_m43)
    m0_bits = slice(word4, 17, 8) ++ slice(word5, 1, 24)
    m0 = unpack_sint(m0_bits, @two_pow_m31)
    cuc = unpack_sint(slice(word6, 1, 16), @two_pow_m29)
    ecc_bits = slice(word6, 17, 8) ++ slice(word7, 1, 24)
    eccentricity = unpack_uint_scaled(ecc_bits, @two_pow_m33)
    cus = unpack_sint(slice(word8, 1, 16), @two_pow_m29)
    sqrt_a_bits = slice(word8, 17, 8) ++ slice(word9, 1, 24)
    sqrt_a = unpack_uint_scaled(sqrt_a_bits, @two_pow_m19)
    toe = unpack_uint_scaled(slice(word10, 1, 16), @two_pow_4)
    fit_interval_flag = bits_to_uint(slice(word10, 17, 1))
    aodo = bits_to_uint(slice(word10, 18, 5))

    %{
      p
      | iode: iode,
        crs: crs,
        delta_n: delta_n,
        m0: m0,
        cuc: cuc,
        eccentricity: eccentricity,
        cus: cus,
        sqrt_a: sqrt_a,
        toe: toe,
        fit_interval_flag: fit_interval_flag,
        aodo: aodo
    }
  end

  # --- Subframe 3 (ephemeris part 2), IS-GPS-200 Table 20-III -----------------

  defp subframe3_words(p) do
    fields = [
      {:cic, p.cic, :sint, 16, @two_pow_m29},
      {:omega0, p.omega0, :sint, 32, @two_pow_m31},
      {:cis, p.cis, :sint, 16, @two_pow_m29},
      {:i0, p.i0, :sint, 32, @two_pow_m31},
      {:crc, p.crc, :sint, 16, @two_pow_m5},
      {:omega, p.omega, :sint, 32, @two_pow_m31},
      {:omega_dot, p.omega_dot, :sint, 24, @two_pow_m43},
      {:iode, p.iode, :uint, 8},
      {:idot, p.idot, :sint, 14, @two_pow_m43}
    ]

    with :ok <- validate_fields(fields) do
      omega0 = pack_sint(p.omega0, 32, @two_pow_m31)
      i0 = pack_sint(p.i0, 32, @two_pow_m31)
      omega = pack_sint(p.omega, 32, @two_pow_m31)

      # Word 3: Cic(16) OMEGA0-MSB(8).
      word3 = pack_sint(p.cic, 16, @two_pow_m29) ++ Enum.slice(omega0, 0, 8)
      # Word 4: OMEGA0-LSB(24).
      word4 = Enum.slice(omega0, 8, 24)
      # Word 5: Cis(16) i0-MSB(8).
      word5 = pack_sint(p.cis, 16, @two_pow_m29) ++ Enum.slice(i0, 0, 8)
      # Word 6: i0-LSB(24).
      word6 = Enum.slice(i0, 8, 24)
      # Word 7: Crc(16) omega-MSB(8).
      word7 = pack_sint(p.crc, 16, @two_pow_m5) ++ Enum.slice(omega, 0, 8)
      # Word 8: omega-LSB(24).
      word8 = Enum.slice(omega, 8, 24)
      # Word 9: OMEGADOT(24).
      word9 = pack_sint(p.omega_dot, 24, @two_pow_m43)
      # Word 10: IODE(8) IDOT(14) + 2 solved bits.
      word10 = pack_uint(p.iode, 8) ++ pack_sint(p.idot, 14, @two_pow_m43)

      {:ok, [word3, word4, word5, word6, word7, word8, word9, {:solve, word10}]}
    end
  end

  defp decode_subframe3(p, w) do
    [word3, word4, word5, word6, word7, word8, word9, word10] = w

    cic = unpack_sint(slice(word3, 1, 16), @two_pow_m29)
    omega0_bits = slice(word3, 17, 8) ++ slice(word4, 1, 24)
    omega0 = unpack_sint(omega0_bits, @two_pow_m31)
    cis = unpack_sint(slice(word5, 1, 16), @two_pow_m29)
    i0_bits = slice(word5, 17, 8) ++ slice(word6, 1, 24)
    i0 = unpack_sint(i0_bits, @two_pow_m31)
    crc = unpack_sint(slice(word7, 1, 16), @two_pow_m5)
    omega_bits = slice(word7, 17, 8) ++ slice(word8, 1, 24)
    omega = unpack_sint(omega_bits, @two_pow_m31)
    omega_dot = unpack_sint(slice(word9, 1, 24), @two_pow_m43)
    idot = unpack_sint(slice(word10, 9, 14), @two_pow_m43)

    %{
      p
      | cic: cic,
        omega0: omega0,
        cis: cis,
        i0: i0,
        crc: crc,
        omega: omega,
        omega_dot: omega_dot,
        idot: idot
    }
  end

  # --- TLM / HOW --------------------------------------------------------------

  defp tlm_data(tlm_message, integrity) do
    # IS-GPS-200 Section 20.3.3.1: preamble(8) message(14) integrity(1) reserved(1).
    pack_uint(@preamble, 8) ++
      pack_uint(tlm_message, 14) ++ pack_uint(integrity, 1) ++ [0]
  end

  defp how_data(tow, alert, anti_spoof, sf_id, :solve) do
    # IS-GPS-200 Section 20.3.3.2: TOW(17) alert(1) A-S(1) SF-ID(3) + 2 solved.
    base =
      pack_uint(tow, 17) ++
        pack_uint(alert, 1) ++ pack_uint(anti_spoof, 1) ++ pack_uint(sf_id, 3)

    {:solve, base ++ [0, 0]}
  end

  # --- word/subframe assembly with parity -------------------------------------

  # Builds the 300-bit subframe from a list of eight 24-bit data words for
  # words 3..10 (the first two entries are the TLM and HOW data), chaining
  # parity through all ten words seeded with D29* = D30* = 0.
  defp assemble_subframe(data_words) do
    {words, _prev} =
      Enum.reduce(data_words, {[], {0, 0}}, fn entry, {acc, {d29_prev, d30_prev}} ->
        {data24, solve?} =
          case entry do
            {:solve, bits} -> {bits, true}
            bits -> {bits, false}
          end

        data24 = if solve?, do: solve_tbits(data24, d29_prev, d30_prev), else: data24

        source = pad24(data24)
        transmitted = Enum.map(source, &Bitwise.bxor(&1, d30_prev))
        par = parity(source, d29_prev, d30_prev)
        word = transmitted ++ par
        [d29, d30] = Enum.slice(par, 4, 2)
        {[word | acc], {d29, d30}}
      end)

    words |> Enum.reverse() |> List.flatten()
  end

  # Solve the two trailing data bits (positions 23, 24) so D29 = D30 = 0.
  defp solve_tbits(data24, d29_prev, d30_prev) do
    base = pad24(data24) |> List.replace_at(22, 0) |> List.replace_at(23, 0)
    [_, _, _, _, d29_k, d30_k] = parity(base, d29_prev, d30_prev)
    d24 = d29_k
    d23 = Bitwise.bxor(d30_k, d24)
    base |> List.replace_at(22, d23) |> List.replace_at(23, d24)
  end

  defp pad24(bits) do
    bits ++ zeros(24 - length(bits))
  end

  defp verify_subframe(bits, sf) when is_list(bits) and length(bits) == @subframe_length do
    words = Enum.chunk_every(bits, @word_length)

    {result, _prev} =
      Enum.reduce_while(Enum.with_index(words, 1), {:ok, {0, 0}}, fn {word, idx},
                                                                     {_acc, {d29p, d30p}} ->
        if parity_valid?(word, d29p, d30p) do
          [d29, d30] = Enum.slice(word, 28, 2)
          {:cont, {:ok, {d29, d30}}}
        else
          {:halt, {{:error, {:parity_failed, sf, idx}}, {d29p, d30p}}}
        end
      end)

    result
  end

  defp verify_subframe(_bits, sf), do: {:error, {:bad_subframe_length, sf}}

  # Returns words 3..10 as lists of 24 *source* data bits (D30* uncomplemented).
  defp source_words(bits) do
    words = Enum.chunk_every(bits, @word_length)

    {decoded, _prev} =
      Enum.map_reduce(words, {0, 0}, fn word, {_d29p, d30p} ->
        {data, par} = Enum.split(word, 24)
        source = Enum.map(data, &Bitwise.bxor(&1, d30p))
        [d29, d30] = Enum.slice(par, 4, 2)
        {source, {d29, d30}}
      end)

    # Drop TLM (word 1) and HOW (word 2); keep words 3..10.
    Enum.drop(decoded, 2)
  end

  # --- packing helpers --------------------------------------------------------

  defp pack_uint(value, bits) do
    for i <- (bits - 1)..0//-1, do: Bitwise.band(Bitwise.bsr(value, i), 1)
  end

  defp pack_uint_scaled(value, bits, scale) do
    pack_uint(round(value / scale), bits)
  end

  defp pack_sint(value, bits, scale) do
    int = round(value / scale)
    pack_twos_complement(int, bits)
  end

  defp pack_twos_complement(int, bits) do
    mask = Bitwise.bsl(1, bits) - 1
    pack_uint(Bitwise.band(int, mask), bits)
  end

  defp bits_to_uint(bits) do
    Enum.reduce(bits, 0, fn b, acc -> Bitwise.bor(Bitwise.bsl(acc, 1), b) end)
  end

  defp unpack_uint_scaled(bits, scale), do: bits_to_uint(bits) * scale

  defp unpack_sint(bits, scale), do: bits_to_sint(bits) * scale

  defp bits_to_sint(bits) do
    n = length(bits)
    raw = bits_to_uint(bits)

    if Bitwise.band(raw, Bitwise.bsl(1, n - 1)) == 0 do
      raw
    else
      raw - Bitwise.bsl(1, n)
    end
  end

  defp zeros(n), do: List.duplicate(0, n)
  defp slice(word, start_1based, len), do: Enum.slice(word, start_1based - 1, len)

  # --- range validation -------------------------------------------------------

  defp validate_fields(fields) do
    Enum.reduce_while(fields, :ok, fn field, :ok ->
      case validate_field(field) do
        :ok -> {:cont, :ok}
        err -> {:halt, err}
      end
    end)
  end

  defp validate_field({name, value, :uint, bits}) do
    validate_uint(name, value, bits)
  end

  defp validate_field({name, value, :uint, bits, scale}) do
    if is_number(value) and value >= 0 do
      validate_uint(name, round(value / scale), bits)
    else
      {:error, {:out_of_range, name, value}}
    end
  end

  defp validate_field({name, value, :sint, bits, scale}) do
    if is_number(value) do
      int = round(value / scale)
      limit = Bitwise.bsl(1, bits - 1)

      if int >= -limit and int <= limit - 1 do
        :ok
      else
        {:error, {:out_of_range, name, value}}
      end
    else
      {:error, {:out_of_range, name, value}}
    end
  end

  defp validate_uint(name, value, bits) do
    if is_integer(value) and value >= 0 and value < Bitwise.bsl(1, bits) do
      :ok
    else
      {:error, {:out_of_range, name, value}}
    end
  end
end
