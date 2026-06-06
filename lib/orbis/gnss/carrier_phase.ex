defmodule Orbis.GNSS.CarrierPhase do
  @moduledoc """
  Dual-frequency carrier-phase linear combinations and the precise-positioning
  prep tooling built on them: geometry-free and wide-lane phase, the
  narrow-lane code, Melbourne-Wubbena, cycle-slip detection, and Hatch
  carrier-smoothed code.

  This is pure Elixir over the existing observation primitives. The phase and
  code observations come from `Orbis.GNSS.RINEX.Observations` (`values/3`,
  `phases/3`). `phases/3` includes `frequency_hz`, `wavelength_m`, and
  `value_m` whenever the band is known; use
  `Orbis.GNSS.RINEX.Observations.band_frequency_hz/2` (or `/3` for GLONASS
  FDMA channels) when building arcs manually. A clean synthetic arc for testing
  can be built from
  `Orbis.GNSS.Observables.predict/5` over an `Orbis.GNSS.SP3` product.

  ## Notation

  For a satellite, band `i` has carrier frequency `f_i` (Hz), wavelength
  `lambda_i = c / f_i`, carrier-phase observation `phi_i` (cycles) so the phase
  in metres is `L_i = lambda_i * phi_i`, and code pseudorange `P_i` (metres),
  with `c = 299_792_458 m/s`.

  ## Combinations

    * **Geometry-free phase** `L_GF = L1 - L2` (metres). The geometry and the
      receiver/satellite clocks are common to both bands and cancel, leaving the
      (dispersive) ionospheric delay plus a constant ambiguity. On a continuous
      arc `L_GF` varies only with the slow ionosphere, so a cycle slip on either
      band makes it jump by (close to) an integer number of band wavelengths.

    * **Wide-lane wavelength** `lambda_WL = c / (f1 - f2)` (about 0.8619 m for
      GPS L1/L2). The long wavelength is why the wide-lane combination is robust.

    * **Wide-lane phase** `phi_WL = phi1 - phi2` (cycles), `L_WL = lambda_WL *
      phi_WL` (metres).

    * **Narrow-lane code** `P_NL = (f1*P1 + f2*P2) / (f1 + f2)` (metres). The
      frequency-weighted code combination paired with the wide-lane phase in
      Melbourne-Wubbena.

    * **Melbourne-Wubbena** `MW = L_WL - P_NL = lambda_WL*(phi1 - phi2) -
      (f1*P1 + f2*P2)/(f1 + f2)` (metres). Sign convention: **wide-lane phase
      minus narrow-lane code**. MW is geometry-, clock-, and ionosphere-free, so
      on a continuous arc it equals the constant wide-lane ambiguity times
      `lambda_WL` plus noise (it is roughly constant); a cycle slip shifts it by
      an integer number of wide-lane cycles, i.e. by `k * lambda_WL` metres.

    * **Code-minus-carrier** `CMC_i = P_i - L_i` (metres). This is the standard
      single-band diagnostic for code multipath/noise plus roughly twice the
      first-order ionosphere, with the carrier ambiguity left as a constant
      offset.

    * **Ionosphere-free Hatch smoothing** combines code and phase into `P_IF`
      and `L_IF` first, then applies the carrier-smoothed-code recursion to the
      ionosphere-free pair. It is a dual-frequency smoothing primitive: it still
      resets on detected cycle slips, but avoids the ionospheric divergence of a
      single-frequency Hatch filter.

  ## Arc shape

  `detect_cycle_slips/2` and `smooth_code/2` take an `arc`: a time-ordered list
  of per-epoch maps for **one** satellite. Each epoch map is a convenient subset
  of

      %{
        epoch: term(),        # opaque, passed through to the output
        phi1: float | nil,    # band-1 carrier phase, cycles
        phi2: float | nil,    # band-2 carrier phase, cycles
        p1:   float | nil,    # band-1 code, metres
        p2:   float | nil,    # band-2 code, metres
        lli1: integer | nil,  # band-1 LLI (bit 0 = loss of lock)
        lli2: integer | nil,  # band-2 LLI
        f1:   float | nil,    # band-1 carrier frequency, Hz (nil => skip)
        f2:   float | nil     # band-2 carrier frequency, Hz
      }

  An epoch with an unknown band frequency is **skipped and reported**, never
  raised: such epochs come back with `skipped: true` and no slip flags. GLONASS
  G1/G2 are usable when the RINEX observation header carries the satellite's
  FDMA frequency-channel number; `Observations.phases/3` applies that map
  automatically.

  ## Non-goals

  Out of scope here (documented so the boundary is explicit): integer ambiguity
  resolution (LAMBDA), double differencing / RTK, triple-frequency combinations,
  and any change to the underlying native primitives.
  """

  import Bitwise, only: [band: 2]

  alias Orbis.GNSS.IonosphereFree

  @c 299_792_458.0

  # Below this absolute frequency separation (Hz) the wide-lane / narrow-lane
  # denominators are treated as degenerate.
  @freq_epsilon 1.0

  # Default cycle-slip thresholds (see detect_cycle_slips/2 docs).
  @default_gf_threshold_m 0.05
  @default_mw_threshold_cycles 4.0

  # Default Hatch smoothing window cap (epochs).
  @default_hatch_window_cap 100

  @type epoch_map :: %{optional(atom()) => term()}
  @type slip_reason :: :lli | :geometry_free | :melbourne_wubbena
  @type slip_result :: %{
          epoch: term(),
          slip: boolean(),
          reasons: [slip_reason()],
          gf: float() | nil,
          mw: float() | nil,
          skipped: boolean()
        }
  @type smooth_result :: %{
          epoch: term(),
          p_smooth: float() | nil,
          window: non_neg_integer(),
          reset: boolean()
        }
  @type iono_free_smooth_result :: %{
          epoch: term(),
          p_smooth: float() | nil,
          p_if: float() | nil,
          l_if: float() | nil,
          window: non_neg_integer(),
          reset: boolean()
        }

  @doc """
  Carrier phase in metres, `L = c / f * phi`.

  `phi_cyc` is carrier phase in cycles and `f_hz` is the carrier frequency in
  hertz. Returns `{:ok, l_m}` or `{:error, :invalid_frequency}` when the
  frequency is not positive. Never raises.
  """
  @spec phase_meters(number(), number()) :: {:ok, float()} | {:error, :invalid_frequency}
  def phase_meters(phi_cyc, f_hz) when is_number(phi_cyc) and is_number(f_hz) do
    if f_hz > 0.0 do
      {:ok, @c / f_hz * phi_cyc}
    else
      {:error, :invalid_frequency}
    end
  end

  def phase_meters(_phi_cyc, _f_hz), do: {:error, :invalid_frequency}

  @doc """
  Geometry-free phase combination `L_GF = l1_m - l2_m` (metres).

  Both inputs are carrier phase already expressed in metres (`L_i = lambda_i *
  phi_i`). The result cancels geometry and clocks, leaving the ionosphere plus a
  constant ambiguity.
  """
  @spec geometry_free(number(), number()) :: float()
  def geometry_free(l1_m, l2_m) when is_number(l1_m) and is_number(l2_m) do
    l1_m - l2_m
  end

  @doc """
  Wide-lane wavelength `lambda_WL = c / (f1 - f2)` (metres).

  Returns `{:ok, lambda_WL}` or `{:error, :equal_frequencies}` when the two
  band frequencies are equal (within a small epsilon), which would divide by
  zero. For GPS L1/L2 this is about `0.8619 m`.
  """
  @spec wide_lane_wavelength(number(), number()) ::
          {:ok, float()} | {:error, :equal_frequencies}
  def wide_lane_wavelength(f1, f2) when is_number(f1) and is_number(f2) do
    if abs(f1 - f2) < @freq_epsilon do
      {:error, :equal_frequencies}
    else
      {:ok, @c / (f1 - f2)}
    end
  end

  def wide_lane_wavelength(_f1, _f2), do: {:error, :equal_frequencies}

  @doc """
  Narrow-lane code `P_NL = (f1*p1 + f2*p2) / (f1 + f2)` (metres).

  The frequency-weighted code combination. Returns `{:ok, P_NL}` or
  `{:error, :equal_frequencies}` when `f1 + f2` is degenerate (zero within
  epsilon).
  """
  @spec narrow_lane_code(number(), number(), number(), number()) ::
          {:ok, float()} | {:error, :equal_frequencies}
  def narrow_lane_code(p1_m, p2_m, f1, f2)
      when is_number(p1_m) and is_number(p2_m) and is_number(f1) and is_number(f2) do
    if abs(f1 + f2) < @freq_epsilon do
      {:error, :equal_frequencies}
    else
      {:ok, (f1 * p1_m + f2 * p2_m) / (f1 + f2)}
    end
  end

  def narrow_lane_code(_p1, _p2, _f1, _f2), do: {:error, :equal_frequencies}

  @doc """
  Melbourne-Wubbena combination (metres).

      MW = L_WL - P_NL
         = lambda_WL*(phi1 - phi2) - (f1*P1 + f2*P2)/(f1 + f2)

  Sign convention: wide-lane phase **minus** narrow-lane code. Inputs are the
  band carrier phases in cycles (`phi1`, `phi2`), the band codes in metres
  (`p1_m`, `p2_m`), and the band frequencies in Hz. Geometry-, clock-, and
  ionosphere-free.

  Returns `{:ok, MW}` or `{:error, :equal_frequencies}` (when `f1 == f2`, so the
  wide-lane wavelength is undefined, or `f1 + f2` is degenerate).
  """
  @spec melbourne_wubbena(number(), number(), number(), number(), number(), number()) ::
          {:ok, float()} | {:error, :equal_frequencies}
  def melbourne_wubbena(phi1_cyc, phi2_cyc, p1_m, p2_m, f1, f2)
      when is_number(phi1_cyc) and is_number(phi2_cyc) and is_number(p1_m) and is_number(p2_m) and
             is_number(f1) and is_number(f2) do
    with {:ok, lambda_wl} <- wide_lane_wavelength(f1, f2),
         {:ok, p_nl} <- narrow_lane_code(p1_m, p2_m, f1, f2) do
      l_wl = lambda_wl * (phi1_cyc - phi2_cyc)
      {:ok, l_wl - p_nl}
    end
  end

  def melbourne_wubbena(_phi1, _phi2, _p1, _p2, _f1, _f2), do: {:error, :equal_frequencies}

  @doc """
  Code-minus-carrier diagnostic `CMC = P - L` (metres).

  `p_m` is code pseudorange in metres, `phi_cyc` is carrier phase in cycles, and
  `f_hz` is the carrier frequency. The carrier phase is converted to metres with
  `phase_meters/2`, then subtracted from the code. Returns `{:ok, cmc_m}` or
  `{:error, :invalid_frequency}`. Never raises.
  """
  @spec code_minus_carrier(number(), number(), number()) ::
          {:ok, float()} | {:error, :invalid_frequency}
  def code_minus_carrier(p_m, phi_cyc, f_hz)
      when is_number(p_m) and is_number(phi_cyc) and is_number(f_hz) do
    case phase_meters(phi_cyc, f_hz) do
      {:ok, l_m} -> {:ok, p_m - l_m}
      {:error, _} = err -> err
    end
  end

  def code_minus_carrier(_p_m, _phi_cyc, _f_hz), do: {:error, :invalid_frequency}

  @doc """
  Detect cycle slips on a single-satellite arc.

  `arc` is a time-ordered list of per-epoch maps (see the module docs for the
  shape). One output map is produced per input epoch, in order:

      %{epoch: term(), slip: boolean(), reasons: [reason],
        gf: float | nil, mw: float | nil, skipped: boolean()}

  with `reason in [:lli, :geometry_free, :melbourne_wubbena]`. A slip at epoch
  `k` is flagged when **any** of:

    * an LLI loss-of-lock bit is set — bit 0 of `lli1` or `lli2`
      (`:lli`, parameter-free);
    * the epoch-to-epoch geometry-free step exceeds `:gf_threshold_m`
      (`:geometry_free`);
    * the epoch-to-epoch Melbourne-Wubbena step, expressed in wide-lane cycles,
      exceeds `:mw_threshold_cycles` (`:melbourne_wubbena`).

  The GF and MW steps require a usable predecessor; the first epoch (and any
  epoch whose predecessor lacked the needed observations) cannot flag `:gf` /
  `:mw` but can still flag `:lli`.

  An epoch with an unknown band frequency (`f1` or `f2` `nil`) is reported with
  `skipped: true`, `slip: false`, `reasons: []`, `gf: nil`,
  `mw: nil`, and never breaks the arc for its neighbours' LLI check.

  ## Options / default thresholds

    * `:gf_threshold_m` (default `#{@default_gf_threshold_m}` m): flags
      `:geometry_free` when `abs(L_GF(k) - L_GF(k-1)) > gf_threshold_m`. At a
      typical 30 s cadence the epoch-to-epoch ionospheric change is sub-cm,
      whereas a single L1 cycle slip jumps `L_GF` by ~0.19 m, so 5 cm separates
      the two cleanly.
    * `:mw_threshold_cycles` (default `#{@default_mw_threshold_cycles}` wide-lane
      cycles): flags `:melbourne_wubbena` when
      `abs(MW(k) - MW(k-1)) / lambda_WL > mw_threshold_cycles`. MW noise is
      dominated by code noise (~0.3-0.5 wide-lane cycles 1-sigma on raw code), so
      4 wide-lane cycles is a robust gate.

  Both thresholds must be non-negative numbers (a negative threshold would flag
  every epoch); any other value raises `ArgumentError`.
  """
  @spec detect_cycle_slips([epoch_map()], keyword()) :: [slip_result()]
  def detect_cycle_slips(arc, opts \\ []) when is_list(arc) do
    {gf_threshold_m, mw_threshold_cycles} = validate_slip_opts!(opts)

    {results, _prev} =
      Enum.map_reduce(arc, nil, fn ep, prev ->
        result = classify_epoch(ep, prev, gf_threshold_m, mw_threshold_cycles)
        # Carry forward only usable (non-skipped) epochs as the predecessor for
        # GF/MW continuity; a skipped epoch leaves the previous good one in place.
        next_prev = if result.skipped, do: prev, else: %{gf: result.gf, mw: result.mw}
        {result, next_prev}
      end)

    results
  end

  @doc """
  Single-frequency Hatch carrier-smoothed code on band 1.

  `arc` is the same per-epoch list as `detect_cycle_slips/2`; this filter uses
  `p1`, `phi1`, `lli1`, and `f1` (for `lambda_1 = c / f1`). The recursion is

      P_smooth(1) = P(1),                         N(1) = 1
      N(k)        = min(N(k-1) + 1, N_cap)
      P_smooth(k) = P(k)/N(k)
                    + ((N(k)-1)/N(k)) * (P_smooth(k-1) + (L1(k) - L1(k-1)))

  with `L1 = lambda_1 * phi1`. The window `N` grows to a cap (`:hatch_window_cap`,
  default `#{@default_hatch_window_cap}`), which must be a positive integer
  (`ArgumentError` otherwise). The slip-threshold options are validated as in
  `detect_cycle_slips/2`.

  **Reset rule.** On a detected cycle slip / LLI loss-of-lock at epoch `k`, or
  when the previous usable sample is missing (start of arc, a gap, or missing
  code/phase), the filter restarts: `N(k) = 1`, `P_smooth(k) = P(k)`, and no
  smoothing carries across the discontinuity. Cycle slips are detected with the
  same logic as `detect_cycle_slips/2`; pass `:gf_threshold_m` /
  `:mw_threshold_cycles` to tune it.

  Returns one map per epoch:

      %{epoch: term(), p_smooth: float | nil, window: non_neg_integer,
        reset: boolean()}

  `p_smooth` is `nil` (with `window: 0`) when the epoch lacks the band-1 code or
  phase, or is skipped (unknown frequency).

  ## Divergence note

  Single-frequency Hatch smoothing accumulates **ionospheric divergence**: the
  ionospheric delay has opposite sign on code and carrier, so as the window `N`
  grows the smoothed code is pulled away from the true range by roughly twice the
  ionospheric change accumulated over the window. Capping `N` bounds this bias;
  the dual-frequency (divergence-free) variant is out of scope here.
  """
  @spec smooth_code([epoch_map()], keyword()) :: [smooth_result()]
  def smooth_code(arc, opts \\ []) when is_list(arc) do
    cap = validate_cap!(opts)
    slips = detect_cycle_slips(arc, opts)

    {results, _state} =
      arc
      |> Enum.zip(slips)
      |> Enum.map_reduce(nil, fn {ep, slip}, state ->
        hatch_step(ep, slip, state, cap)
      end)

    results
  end

  @doc """
  Dual-frequency ionosphere-free Hatch carrier-smoothed code.

  This uses the same arc shape and cycle-slip detector as `smooth_code/2`, but it
  forms the ionosphere-free code and carrier phase at each epoch first:

      P_IF = gamma * P1 - (gamma - 1) * P2
      L_IF = gamma * L1 - (gamma - 1) * L2

  The Hatch recursion then smooths `P_IF` with epoch-to-epoch changes in `L_IF`.
  This avoids the ionospheric divergence of the single-frequency Hatch filter,
  while still resetting on LLI/cycle-slip detections and data gaps.

  Returns one map per epoch:

      %{epoch: term(), p_smooth: float | nil, p_if: float | nil,
        l_if: float | nil, window: non_neg_integer, reset: boolean()}

  `p_smooth`, `p_if`, and `l_if` are `nil` (with `window: 0`) when an epoch is
  missing a required dual-frequency code/phase/frequency value.
  """
  @spec smooth_iono_free_code([epoch_map()], keyword()) :: [iono_free_smooth_result()]
  def smooth_iono_free_code(arc, opts \\ []) when is_list(arc) do
    cap = validate_cap!(opts)
    slips = detect_cycle_slips(arc, opts)

    {results, _state} =
      arc
      |> Enum.zip(slips)
      |> Enum.map_reduce(nil, fn {ep, slip}, state ->
        iono_free_hatch_step(ep, slip, state, cap)
      end)

    results
  end

  # --- option validation --------------------------------------------------

  # Slip thresholds gate an absolute step, so a negative value would flag every
  # epoch (abs(...) is always >= 0 > negative); reject non-numbers and negatives.
  defp validate_slip_opts!(opts) do
    gf = Keyword.get(opts, :gf_threshold_m, @default_gf_threshold_m)
    mw = Keyword.get(opts, :mw_threshold_cycles, @default_mw_threshold_cycles)
    {non_negative!(:gf_threshold_m, gf), non_negative!(:mw_threshold_cycles, mw)}
  end

  defp non_negative!(_name, value) when is_number(value) and value >= 0.0, do: value

  defp non_negative!(name, value) do
    raise ArgumentError, "#{inspect(name)} must be a non-negative number, got: #{inspect(value)}"
  end

  # The Hatch window divides by N, which starts at 1 and is capped at this value,
  # so a cap below 1 would divide by zero; require a positive integer.
  defp validate_cap!(opts) do
    cap = Keyword.get(opts, :hatch_window_cap, @default_hatch_window_cap)

    if is_integer(cap) and cap >= 1 do
      cap
    else
      raise ArgumentError, ":hatch_window_cap must be a positive integer, got: #{inspect(cap)}"
    end
  end

  # --- cycle-slip classification ------------------------------------------

  defp classify_epoch(ep, prev, gf_threshold_m, mw_threshold_cycles) do
    f1 = Map.get(ep, :f1)
    f2 = Map.get(ep, :f2)
    epoch = Map.get(ep, :epoch)

    # An unknown OR non-numeric frequency (e.g. a GLONASS satellite whose band
    # frequency is nil, or a malformed entry) is skipped, matching the docs and
    # smooth_code/2 — and the good predecessor is preserved for the continuity
    # check rather than being reset by an unusable epoch.
    if not is_number(f1) or not is_number(f2) do
      %{epoch: epoch, slip: false, reasons: [], gf: nil, mw: nil, skipped: true}
    else
      phi1 = Map.get(ep, :phi1)
      phi2 = Map.get(ep, :phi2)
      p1 = Map.get(ep, :p1)
      p2 = Map.get(ep, :p2)

      gf = current_gf(phi1, phi2, f1, f2)
      mw = current_mw(phi1, phi2, p1, p2, f1, f2)

      lli_reason = if loss_of_lock?(ep), do: [:lli], else: []
      gf_reason = gf_reason(gf, prev, gf_threshold_m)
      mw_reason = mw_reason(mw, prev, f1, f2, mw_threshold_cycles)

      reasons = lli_reason ++ gf_reason ++ mw_reason

      %{
        epoch: epoch,
        slip: reasons != [],
        reasons: reasons,
        gf: gf,
        mw: mw,
        skipped: false
      }
    end
  end

  defp current_gf(phi1, phi2, f1, f2)
       when is_number(phi1) and is_number(phi2) and is_number(f1) and is_number(f2) do
    with {:ok, l1} <- phase_meters(phi1, f1),
         {:ok, l2} <- phase_meters(phi2, f2) do
      geometry_free(l1, l2)
    else
      {:error, _} -> nil
    end
  end

  defp current_gf(_phi1, _phi2, _f1, _f2), do: nil

  defp current_mw(phi1, phi2, p1, p2, f1, f2)
       when is_number(phi1) and is_number(phi2) and is_number(p1) and is_number(p2) do
    case melbourne_wubbena(phi1, phi2, p1, p2, f1, f2) do
      {:ok, mw} -> mw
      {:error, _} -> nil
    end
  end

  defp current_mw(_phi1, _phi2, _p1, _p2, _f1, _f2), do: nil

  defp gf_reason(gf, prev, threshold_m) when is_number(gf) and is_map(prev) do
    case prev do
      %{gf: prev_gf} when is_number(prev_gf) ->
        if abs(gf - prev_gf) > threshold_m, do: [:geometry_free], else: []

      _ ->
        []
    end
  end

  defp gf_reason(_gf, _prev, _threshold_m), do: []

  defp mw_reason(mw, prev, f1, f2, threshold_cycles) when is_number(mw) and is_map(prev) do
    case {prev, wide_lane_wavelength(f1, f2)} do
      {%{mw: prev_mw}, {:ok, lambda_wl}} when is_number(prev_mw) ->
        step_cycles = abs(mw - prev_mw) / abs(lambda_wl)
        if step_cycles > threshold_cycles, do: [:melbourne_wubbena], else: []

      _ ->
        []
    end
  end

  defp mw_reason(_mw, _prev, _f1, _f2, _threshold_cycles), do: []

  defp loss_of_lock?(ep) do
    lli_set?(Map.get(ep, :lli1)) or lli_set?(Map.get(ep, :lli2))
  end

  defp lli_set?(lli) when is_integer(lli), do: band(lli, 1) == 1
  defp lli_set?(_lli), do: false

  # --- Hatch smoothing -----------------------------------------------------

  defp hatch_step(ep, slip, state, cap) do
    epoch = Map.get(ep, :epoch)
    f1 = Map.get(ep, :f1)
    p1 = Map.get(ep, :p1)
    phi1 = Map.get(ep, :phi1)

    if slip.skipped or not is_number(f1) or not is_number(p1) or not is_number(phi1) do
      # No usable band-1 sample: emit nil and drop the running state so the
      # next usable epoch restarts cleanly.
      {%{epoch: epoch, p_smooth: nil, window: 0, reset: false}, nil}
    else
      case phase_meters(phi1, f1) do
        {:ok, l1} -> do_hatch(epoch, p1, l1, slip.slip, state, cap)
        {:error, _} -> {%{epoch: epoch, p_smooth: nil, window: 0, reset: false}, nil}
      end
    end
  end

  defp do_hatch(epoch, p1, l1, slip?, state, cap) do
    if slip? or is_nil(state) do
      # Reset: start a new smoothing window at the raw code.
      result = %{epoch: epoch, p_smooth: p1, window: 1, reset: state != nil and slip?}
      {result, %{p_smooth: p1, l1: l1, window: 1}}
    else
      %{p_smooth: prev_smooth, l1: prev_l1, window: prev_window} = state
      n = min(prev_window + 1, cap)
      p_smooth = p1 / n + (n - 1) / n * (prev_smooth + (l1 - prev_l1))
      result = %{epoch: epoch, p_smooth: p_smooth, window: n, reset: false}
      {result, %{p_smooth: p_smooth, l1: l1, window: n}}
    end
  end

  # --- ionosphere-free Hatch smoothing ------------------------------------

  defp iono_free_hatch_step(ep, slip, state, cap) do
    epoch = Map.get(ep, :epoch)

    case current_iono_free_code_phase(ep) do
      {:ok, p_if, l_if} ->
        do_iono_free_hatch(epoch, p_if, l_if, slip.slip, state, cap)

      {:error, _} ->
        result = %{epoch: epoch, p_smooth: nil, p_if: nil, l_if: nil, window: 0, reset: false}
        {result, nil}
    end
  end

  defp current_iono_free_code_phase(ep) do
    f1 = Map.get(ep, :f1)
    f2 = Map.get(ep, :f2)
    p1 = Map.get(ep, :p1)
    p2 = Map.get(ep, :p2)
    phi1 = Map.get(ep, :phi1)
    phi2 = Map.get(ep, :phi2)

    with true <- Enum.all?([f1, f2, p1, p2, phi1, phi2], &is_number/1),
         {:ok, p_if} <- IonosphereFree.iono_free(p1, p2, f1, f2),
         {:ok, l_if} <- IonosphereFree.iono_free_phase_cycles(phi1, phi2, f1, f2) do
      {:ok, p_if, l_if}
    else
      _ -> {:error, :missing_observation}
    end
  end

  defp do_iono_free_hatch(epoch, p_if, l_if, slip?, state, cap) do
    if slip? or is_nil(state) do
      result = %{
        epoch: epoch,
        p_smooth: p_if,
        p_if: p_if,
        l_if: l_if,
        window: 1,
        reset: state != nil and slip?
      }

      {result, %{p_smooth: p_if, l_if: l_if, window: 1}}
    else
      %{p_smooth: prev_smooth, l_if: prev_l_if, window: prev_window} = state
      n = min(prev_window + 1, cap)
      p_smooth = p_if / n + (n - 1) / n * (prev_smooth + (l_if - prev_l_if))

      result = %{
        epoch: epoch,
        p_smooth: p_smooth,
        p_if: p_if,
        l_if: l_if,
        window: n,
        reset: false
      }

      {result, %{p_smooth: p_smooth, l_if: l_if, window: n}}
    end
  end
end
