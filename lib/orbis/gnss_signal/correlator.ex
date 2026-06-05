defmodule Orbis.GnssSignal.Correlator do
  @moduledoc """
  Baseband simulation, correlation, and acquisition for the GPS L1 C/A signal.

  This module builds the "generate a replica, correlate, and acquire a PRN"
  layer on top of `Orbis.GnssSignal.CA`. It works entirely in complex baseband
  (the carrier has already been removed down to a residual Doppler), using the
  standard time-domain signal-processing model described in the GNSS
  acquisition literature:

    * Kaplan & Hegarty, *Understanding GPS/GNSS: Principles and Applications*
      (3rd ed.), Ch. 5-8 (signal acquisition, coherent correlation, the
      sinc Doppler-mismatch loss, and the C/N0 to post-correlation SNR
      relation).
    * Misra & Enge, *Global Positioning System: Signals, Measurements, and
      Performance* (2nd ed.), Ch. 10-11.
    * Borre, Akos, Bertelsen, Rinder & Jensen, *A Software-Defined GPS and
      Galileo Receiver* (the 2D code-phase by Doppler search and the
      peak-to-noise acquisition metric).

  ## Sampled-code replica (zero-order-hold / nearest-chip)

  The 1023-chip C/A code is sampled to a sampling rate `fs` over an integration
  time `T` (so `N = round(fs * T)` samples). At sample `n` the code phase, in
  chips, advances at the code rate:

      code_rate = f_chip * (1 + fd_code / f_L1)
      sampled[n] = chip(prn, floor((code_phase + n * code_rate / fs) mod 1023))

  with `f_chip = 1_023_000` cps and `f_L1 = 1_575_420_000` Hz. `code_phase` is
  the initial offset in chips and `fd_code` is an optional code-rate Doppler
  scaling (defaults to `0.0`; code Doppler is negligible over one 1 ms period
  but the parameter is exposed). This is a clean nearest-chip (zero-order-hold)
  sampler: it picks the chip the sample instant falls within. `CA.chip/2`
  already wraps its index modulo 1023, so the sampler wraps for free.

  ## Coherent correlation

  For a complex baseband record `x[n] = xI[n] + j*xQ[n]` (a real-valued test
  signal has `xQ = 0`), correlation against a local carrier wipe-off
  `exp(-j 2*pi*f_d*n/fs)` times the real bipolar code `c[n]` is the coherent
  sum over the integration window:

      S = sum_{n=0}^{N-1} x[n] * c[n] * exp(-j 2*pi*f_d*n/fs)
      I = sum ( xI[n]*c[n]*cos(th_n) + xQ[n]*c[n]*sin(th_n) )
      Q = sum ( xQ[n]*c[n]*cos(th_n) - xI[n]*c[n]*sin(th_n) ),  th_n = 2*pi*f_d*n/fs
      power = I^2 + Q^2

  The amplitude is recovered as `sqrt(I^2 + Q^2)`; for a real input the carrier
  wipe-off still spreads energy across both I and Q so the magnitude is the
  meaningful quantity.

  ## Acquisition metric

  `acquire/3` performs a 2D search over code-phase bins (at sample resolution)
  and Doppler bins, computing the correlator `power` on each cell. The detection
  metric is the standard peak-to-noise ratio used in software receivers:

      metric = peak_power / mean(off_peak_powers)

  where the off-peak set excludes the peak cell and an exclusion zone of one
  code-phase bin on either side of the peak (at the peak's Doppler bin), so the
  main correlation lobe is not mistaken for the noise floor. For a clean signal
  this metric is large (order of the number of samples); for noise or the wrong
  PRN it is close to one. The alternative peak-to-second-peak ratio is also a
  standard choice; this module exposes peak-to-mean-off-peak as `metric`.

  ## Coherent integration loss and post-correlation SNR

  A residual frequency error `f` over a coherent integration time `T` attenuates
  the correlation by the sinc-squared Doppler-mismatch loss (Kaplan & Hegarty,
  the residual-carrier loss whose discrete form
  `|sin(N*pi*df*Ts) / (N*sin(pi*df*Ts))|` tends to `sinc(pi*f*T)` in the
  continuous limit):

      coherent_loss(f, T) = sinc^2(pi*f*T) = ( sin(pi*f*T) / (pi*f*T) )^2

  which is `1` (0 dB) at `f = 0` and has its first null at `f = 1/T`. The
  correlation *amplitude* scales as `|sinc(pi*f*T)| = sqrt(coherent_loss)`.

  For the post-correlation signal-to-noise ratio this module exposes only the
  relation it can state cleanly from a standard reference: coherent integration
  over `T` seconds gives the predetection SNR

      snr_post_db(cn0_dbhz, T) = cn0_dbhz + 10*log10(T)

  i.e. the carrier-to-noise-density ratio C/N0 (dB-Hz) plus the processing gain
  `10*log10(T)` corresponding to a `1/T` effective noise bandwidth (Kaplan &
  Hegarty; Misra & Enge). The exact noise-bandwidth convention (`1/T` versus
  `1/(2T)`) differs between texts; this module uses the `1/T` form. No second,
  competing gain formula is exposed.
  """

  alias Orbis.GnssSignal.CA

  @code_length 1023
  @chip_rate_hz 1_023_000
  @l1_hz 1_575_420_000

  @two_pi 2.0 * :math.pi()

  @default_doppler_min_hz -2500.0
  @default_doppler_max_hz 2500.0
  @default_doppler_step_hz 500.0
  @default_sample_rate_hz 2.046e6

  @doc """
  Builds a sampled `±1` C/A code replica.

  Options:

    * `:sample_rate_hz` - sampling rate in Hz (default `2.046e6`, 2 samples/chip).
    * `:integration_time_s` - integration time in seconds (default one code
      period, `1023 / 1_023_000` = 1 ms). Determines `N` together with the
      sample rate, unless `:num_samples` is given.
    * `:num_samples` - explicit sample count `N` (overrides `:integration_time_s`).
    * `:code_phase_chips` - initial code phase offset in chips (default `0.0`).
    * `:code_doppler_hz` - code-rate Doppler scaling in Hz (default `0.0`).

  Returns `{:ok, samples}` with `samples` a list of `±1` integers of length `N`,
  or propagates `{:error, {:unsupported_prn, prn}}` from `Orbis.GnssSignal.CA`.

  ## Examples

      iex> {:ok, s} = Orbis.GnssSignal.Correlator.replica(1, num_samples: 4, sample_rate_hz: 1.023e6)
      iex> s
      [-1, -1, 1, 1]

  """
  @spec replica(integer(), keyword()) ::
          {:ok, [integer()]} | {:error, {:unsupported_prn, integer()}}
  def replica(prn, opts \\ []) do
    fs = Keyword.get(opts, :sample_rate_hz, @default_sample_rate_hz)
    code_phase = Keyword.get(opts, :code_phase_chips, 0.0) * 1.0
    code_doppler = Keyword.get(opts, :code_doppler_hz, 0.0) * 1.0
    n = num_samples(opts, fs)

    with {:ok, code} <- CA.code(prn) do
      {:ok, sample_code(code, n, fs, code_phase, code_doppler)}
    end
  end

  @doc """
  Coherently correlates a complex baseband record against a PRN replica.

  `iq` is a list of samples, each either a `{i, q}` tuple or a bare real number
  (interpreted as `i` with `q = 0`). Options:

    * `:sample_rate_hz` - default `2.046e6`.
    * `:doppler_hz` - residual carrier Doppler to wipe off (default `0.0`).
    * `:code_phase_chips` - replica code phase offset in chips (default `0.0`).
    * `:code_doppler_hz` - replica code-rate Doppler (default `0.0`).

  Returns `{:ok, %{i: i, q: q, power: i*i + q*q}}`, or propagates the PRN error.
  The replica is generated at `length(iq)` samples to match the record.
  """
  @spec correlate(list(), integer(), keyword()) ::
          {:ok, %{i: float(), q: float(), power: float()}}
          | {:error, {:unsupported_prn, integer()} | :empty_samples}
  def correlate(iq, prn, opts \\ [])
  def correlate([], _prn, _opts), do: {:error, :empty_samples}

  def correlate(iq, prn, opts) when is_list(iq) do
    fs = Keyword.get(opts, :sample_rate_hz, @default_sample_rate_hz)
    doppler = Keyword.get(opts, :doppler_hz, 0.0) * 1.0
    code_phase = Keyword.get(opts, :code_phase_chips, 0.0) * 1.0
    code_doppler = Keyword.get(opts, :code_doppler_hz, 0.0) * 1.0
    n = length(iq)

    with {:ok, code} <- CA.code(prn) do
      sampled = sample_code(code, n, fs, code_phase, code_doppler)
      {i, q} = correlate_against(iq, sampled, fs, doppler)
      {:ok, %{i: i, q: q, power: i * i + q * q}}
    end
  end

  @doc """
  Low-level coherent correlation of a baseband record against an explicit
  sampled `±1` code.

  `iq` is a list of `{i, q}` tuples or bare reals; `code` is the sampled `±1`
  vector (same length); `fs` is the sample rate and `doppler_hz` the residual
  carrier to wipe off. Returns `{i, q}`, the real and imaginary parts of the
  coherent sum.
  """
  @spec correlate_against(list(), [integer()], number(), number()) :: {float(), float()}
  def correlate_against(iq, code, fs, doppler_hz) when is_list(iq) and is_list(code) do
    w = @two_pi * doppler_hz / fs

    {i, q, _n} =
      Enum.zip(iq, code)
      |> Enum.reduce({0.0, 0.0, 0}, fn {sample, c}, {acc_i, acc_q, n} ->
        {xi, xq} = to_iq(sample)
        theta = w * n
        cos = :math.cos(theta)
        sin = :math.sin(theta)
        cc = c * 1.0
        # x * exp(-j theta): real = xi*cos + xq*sin, imag = xq*cos - xi*sin
        di = (xi * cos + xq * sin) * cc
        dq = (xq * cos - xi * sin) * cc
        {acc_i + di, acc_q + dq, n + 1}
      end)

    {i, q}
  end

  @doc """
  Acquires a PRN by a 2D search over code phase and Doppler.

  `samples` is the complex baseband record (list of `{i, q}` tuples or bare
  reals). Options:

    * `:sample_rate_hz` - default `2.046e6`.
    * `:doppler_min_hz` / `:doppler_max_hz` / `:doppler_step_hz` - the Doppler
      search grid (defaults `-2500`, `2500`, `500` Hz).

  The code-phase axis is searched at sample resolution over one code period.
  Returns

      {:ok, %{
        code_phase_chips: float,
        doppler_hz: float,
        peak_metric: float,
        metric: float,
        peak_power: float,
        grid: %{
          doppler_hz: [float],
          code_phase_bins: integer,
          doppler_step_hz: float,
          samples_per_chip: float
        }
      }}

  `metric` (and its alias `peak_metric`) is the peak-to-mean-off-peak power
  ratio described in the module docs. Errors:

    * `{:error, :empty_samples}` for an empty record,
    * `{:error, :too_short}` if the record is shorter than one code period,
    * `{:error, {:unsupported_prn, prn}}` propagated from the code generator.
  """
  @spec acquire(list(), integer(), keyword()) ::
          {:ok, map()}
          | {:error, :empty_samples | :too_short | {:unsupported_prn, integer()}}
  def acquire(samples, prn, opts \\ [])
  def acquire([], _prn, _opts), do: {:error, :empty_samples}

  def acquire(samples, prn, opts) when is_list(samples) do
    fs = Keyword.get(opts, :sample_rate_hz, @default_sample_rate_hz)
    samples_per_chip = fs / @chip_rate_hz
    samples_per_code = round(samples_per_chip * @code_length)
    n = length(samples)

    if n < samples_per_code do
      {:error, :too_short}
    else
      with {:ok, code} <- CA.code(prn) do
        do_acquire(samples, code, prn, fs, samples_per_chip, samples_per_code, opts)
      end
    end
  end

  defp do_acquire(samples, code, _prn, fs, samples_per_chip, samples_per_code, opts) do
    dmin = Keyword.get(opts, :doppler_min_hz, @default_doppler_min_hz) * 1.0
    dmax = Keyword.get(opts, :doppler_max_hz, @default_doppler_max_hz) * 1.0
    dstep = Keyword.get(opts, :doppler_step_hz, @default_doppler_step_hz) * 1.0

    doppler_bins = doppler_grid(dmin, dmax, dstep)

    # Use exactly one code period of samples and the matching base replica.
    record = Enum.take(samples, samples_per_code)
    base_code = sample_code(code, samples_per_code, fs, 0.0, 0.0)
    base_tuple = List.to_tuple(base_code)

    iq_tuple =
      record
      |> Enum.map(&to_iq/1)
      |> List.to_tuple()

    # For each Doppler bin: wipe off the carrier once, then slide the base code
    # over every sample offset (circular) to get the full code-phase axis.
    grid =
      for d <- doppler_bins do
        wiped = carrier_wipeoff(iq_tuple, fs, d, samples_per_code)
        powers = code_phase_powers(wiped, base_tuple, samples_per_code)
        {d, powers}
      end

    # Find the global peak across the grid.
    {peak_power, peak_doppler, peak_offset} =
      Enum.reduce(grid, {-1.0, 0.0, 0}, fn {d, powers}, acc ->
        Enum.reduce(0..(samples_per_code - 1), acc, fn off, {bp, _bd, _bo} = a ->
          p = elem(powers, off)
          if p > bp, do: {p, d, off}, else: a
        end)
      end)

    metric = peak_to_mean_off_peak(grid, peak_power, peak_doppler, peak_offset, samples_per_code)
    code_phase_chips = peak_offset / samples_per_chip

    {:ok,
     %{
       code_phase_chips: code_phase_chips,
       doppler_hz: peak_doppler,
       peak_metric: metric,
       metric: metric,
       peak_power: peak_power,
       grid: %{
         doppler_hz: doppler_bins,
         code_phase_bins: samples_per_code,
         doppler_step_hz: dstep,
         samples_per_chip: samples_per_chip
       }
     }}
  end

  # Carrier wipe-off: returns a tuple of {i, q} after multiplying by exp(-j w n).
  defp carrier_wipeoff(iq_tuple, fs, doppler_hz, n) do
    w = @two_pi * doppler_hz / fs

    0..(n - 1)
    |> Enum.map(fn k ->
      {xi, xq} = elem(iq_tuple, k)
      theta = w * k
      cos = :math.cos(theta)
      sin = :math.sin(theta)
      {xi * cos + xq * sin, xq * cos - xi * sin}
    end)
    |> List.to_tuple()
  end

  # For a carrier-wiped record (tuple of {i,q}), compute the correlator power at
  # every circular code-phase offset of the base code. Returns a tuple of powers
  # indexed by sample offset.
  defp code_phase_powers(wiped, base_tuple, n) do
    0..(n - 1)
    |> Enum.map(fn offset ->
      {i, q} =
        Enum.reduce(0..(n - 1), {0.0, 0.0}, fn k, {ai, aq} ->
          {xi, xq} = elem(wiped, k)
          c = elem(base_tuple, Integer.mod(k + offset, n)) * 1.0
          {ai + xi * c, aq + xq * c}
        end)

      i * i + q * q
    end)
    |> List.to_tuple()
  end

  defp peak_to_mean_off_peak(grid, peak_power, peak_doppler, peak_offset, n) do
    {sum, count} =
      Enum.reduce(grid, {0.0, 0}, fn {d, powers}, {s, c} ->
        Enum.reduce(0..(n - 1), {s, c}, fn off, {s2, c2} ->
          if d == peak_doppler and abs_circular_diff(off, peak_offset, n) <= 1 do
            {s2, c2}
          else
            {s2 + elem(powers, off), c2 + 1}
          end
        end)
      end)

    cond do
      count == 0 -> 0.0
      sum <= 0.0 and peak_power > 0.0 -> 1.0e12
      sum <= 0.0 -> 0.0
      true -> peak_power / (sum / count)
    end
  end

  defp doppler_grid(dmin, dmax, dstep) do
    count = round((dmax - dmin) / dstep)

    0..count
    |> Enum.map(fn k -> dmin + k * dstep end)
    |> Enum.filter(fn d -> d <= dmax + 1.0e-9 end)
  end

  defp abs_circular_diff(a, b, n) do
    d = Integer.mod(abs(a - b), n)
    min(d, n - d)
  end

  @doc """
  Coherent integration loss from a residual frequency error.

  Returns the linear loss `sinc^2(pi*f*T) = (sin(pi*f*T)/(pi*f*T))^2` in
  `[0, 1]`: `1.0` at `f = 0`, with its first null at `f = 1/T`. The correlation
  amplitude scales as `sqrt` of this value.

  ## Examples

      iex> Orbis.GnssSignal.Correlator.coherent_loss(0.0, 1.0e-3)
      1.0

  """
  @spec coherent_loss(number(), number()) :: float()
  def coherent_loss(freq_error_hz, integration_time_s) do
    x = :math.pi() * freq_error_hz * integration_time_s

    if x == 0.0 do
      1.0
    else
      s = :math.sin(x) / x
      s * s
    end
  end

  @doc """
  Coherent integration loss in decibels, `10*log10(coherent_loss(f, T))`.

  Returns `:neg_infinity` at an exact null (loss of zero).
  """
  @spec coherent_loss_db(number(), number()) :: float() | :neg_infinity
  def coherent_loss_db(freq_error_hz, integration_time_s) do
    loss = coherent_loss(freq_error_hz, integration_time_s)

    if loss <= 0.0 do
      :neg_infinity
    else
      10.0 * :math.log10(loss)
    end
  end

  @doc """
  Post-correlation (predetection) SNR in dB from C/N0 and integration time.

  Uses the standard relation `snr_post_db = cn0_dbhz + 10*log10(T)` (Kaplan &
  Hegarty; Misra & Enge), corresponding to a `1/T` effective noise bandwidth.

  ## Examples

      iex> Float.round(Orbis.GnssSignal.Correlator.snr_post_db(40.0, 1.0e-3), 6)
      10.0

  """
  @spec snr_post_db(number(), number()) :: float()
  def snr_post_db(cn0_dbhz, integration_time_s) when integration_time_s > 0 do
    cn0_dbhz + 10.0 * :math.log10(integration_time_s)
  end

  # --- internal helpers ---

  defp num_samples(opts, fs) do
    case Keyword.get(opts, :num_samples) do
      nil ->
        t = Keyword.get(opts, :integration_time_s, @code_length / @chip_rate_hz)
        round(fs * t)

      n when is_integer(n) and n > 0 ->
        n
    end
  end

  defp sample_code(code, n, fs, code_phase, code_doppler) do
    code_tuple = List.to_tuple(code)
    len = tuple_size(code_tuple)
    code_rate = @chip_rate_hz * (1.0 + code_doppler / @l1_hz)
    per_sample = code_rate / fs

    Enum.map(0..(n - 1), fn k ->
      pos = code_phase + k * per_sample
      idx = Integer.mod(floor(pos), len)
      elem(code_tuple, idx)
    end)
  end

  defp to_iq({i, q}), do: {i * 1.0, q * 1.0}
  defp to_iq(i) when is_number(i), do: {i * 1.0, 0.0}
end
