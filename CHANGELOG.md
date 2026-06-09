# Changelog

All notable changes to this project are documented here. The format is based on
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project
adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.14.0] - 2026-06-08

### Added

- `Orbis.GNSS.SP3.to_iodata/2` serializes an `%Orbis.GNSS.SP3{}` product back to
  standard SP3-c / SP3-d text — the inverse of the reader, so a read → `merge/2`
  → write pipeline emits a single standard SP3 file any reader consumes. Pure and
  deterministic; header fields are derived from the product; a satellite absent
  at an epoch is written as the SP3 missing-orbit sentinel (so a quarantined
  merge cell re-reads as missing, never a fabricated position). Round-trips to
  SP3 format precision (mm / sub-ns) for position-only and position+velocity,
  multi-constellation products.
- `Orbis.GNSS.Data.write_sp3/3` writes a product to disk with the fetch layer's
  atomic-commit discipline (same-directory temp file + `File.rename/2`), with an
  optional `gzip: true` for the gzipped-archive shape. Unblocks persisting a
  merged product, which was otherwise only an in-memory handle.

## [0.13.0] - 2026-06-08

### Added

- `Orbis.GNSS.Data.ops_ultra_sp3/3` and `ops_ultra_clk/3` add the ultra-rapid
  precise-product tier to the offline-safe catalog/fetch layer. The catalog now
  derives anonymous GSSC archive names and URLs for `IGS0OPSULT`, `COD0OPSULT`,
  `ESA0OPSULT`, `GFZ0OPSULT`, and `GRG0OPSULT` SP3 products (plus `GRG0OPSULT`
  clocks), including sub-daily issue times, `02D` spans, per-center sampling,
  and latest-available issue fallback before a target epoch.
- `Orbis.GNSS.Data.fetch_merged_sp3/3` fetches the same SP3 product from several
  centers in precedence order, tolerates not-yet-published or missing centers,
  and returns one merged `Orbis.GNSS.SP3` plus provenance and merge-audit
  metadata. One available center is returned as a flagged single-source result;
  zero available centers returns `{:error, {:no_products, reasons}}`; centers
  that cannot be combined (mismatched time scale / coordinate-system frame)
  return `{:error, {:incompatible_sources, %{centers:, reason:}}}` rather than
  leaking a raw merge error.

## [0.12.0] - 2026-06-08

### Added

- `Orbis.GNSS.SP3.merge/2` merges several SP3 products from different analysis
  centers into one consistent precise-ephemeris dataset. Coverage is the union
  across satellite×epoch (a satellite present in any input is present in the
  output, filling a single center's dropouts); overlapping records are resolved
  by robust consensus — the largest subset of centers agreeing within tolerance
  is combined (`:mean`, `:median`, or `:precedence`), disagreeing centers are
  recorded as outliers, and a cell with no agreeing subset is quarantined rather
  than averaged. Pure and deterministic; returns the merged product plus an audit
  report (`:quarantined`, `:single_source`, `:position_outliers`).
- `Orbis.GNSS.SP3.clock_reference_offset/3` and
  `Orbis.GNSS.SP3.align_clock_reference/3` expose the clock-datum primitive:
  precise clock products from different centers are referenced to different
  station/ensemble clocks, so their raw clocks differ by a per-epoch common
  offset. The first estimates that offset (robust median over common satellites);
  the second returns a copy of a product with its clocks shifted onto a
  reference's datum so the two are directly comparable. Positions need no such
  treatment.
- `Orbis.GNSS.BroadcastComparison` now reports `clock_datum_removed_rms_m` /
  `clock_datum_removed_max_m` alongside the raw clock statistics: the per-epoch
  common reference-clock offset (median over satellites) is removed to give the
  actual signal-in-space clock error, several times smaller than the raw value.
- `Orbis.GNSS.Ephemeris.sample/3` samples a precise (`Orbis.GNSS.SP3`) or
  broadcast (`Orbis.GNSS.Broadcast`) ephemeris over an epoch window into a
  unified per-satellite, per-epoch table of ECEF position and clock bias — the
  same call shape for either source, with out-of-coverage cells reported as an
  explicit `:no_ephemeris` gap rather than extrapolated.
- `Orbis.GNSS.Broadcast.position/3` evaluates a single satellite's broadcast
  ECEF position and clock at an epoch (IS-GPS-200 LNAV, Galileo OS-SIS-ICD,
  BeiDou BDS-SIS-ICD).
- `Orbis.GNSS.BroadcastComparison.compare/4` (and the `mix gnss.broadcast_diff`
  task, with a `--system` selector) computes per-satellite broadcast-vs-precise
  orbit and clock differences (3D plus radial/along/cross RMS and max) over a
  window — the standard broadcast ephemeris accuracy check. Validated over a full
  UTC day against the IGS combined broadcast (`BRDC00IGS`) and CODE MGEX final
  precise orbits (`COD0MGXFIN`): GPS LNAV ~1.4 m, Galileo I/NAV ~0.9 m, BeiDou
  ~2.5 m orbit RMS.
- `Orbis.GNSS.RTK.solve_float_baseline_epochs/3` and fixed RTK solvers now
  accept `code_smoothing: true` to apply per-receiver/per-ambiguity-arc Hatch
  carrier smoothing to code observations before forming double differences.
  The real Wettzell RTK gate verifies the smoothing reduces code residual RMS
  while still refusing unsafe integer fixes.
- `Orbis.GNSS.RTK.solve_fixed_baseline_epochs/3` now supports opt-in partial
  ambiguity resolution with `partial_ambiguity_resolution: true`. When the full
  ambiguity set fails the ratio test, Orbis tries confidence-ranked subsets and
  re-solves with the accepted subset fixed while rejected ambiguities remain
  float-estimated. The real Wettzell RTK gate now verifies a safe four-ambiguity
  partial fix improves the L1 baseline while the unsafe full-set fix remains
  rejected.

### Changed

- GNSS integer ambiguity fixing now uses a complete bounded integer
  least-squares scan over the caller's `integer_search_radius_cycles`, scored by
  the exact ambiguity covariance inverse. Fixed-solution metadata reports
  `integer_method: :bounded_ils` (or
  `:widelane_narrowlane_bounded_ils`) for this path.
- The default integer candidate cap for precise positioning and RTK fixed
  solvers is now `200_000`, enough for the default radius-1 search with up to 11
  ambiguities.
- `Orbis.GNSS.RTK.solve_float_baseline_epochs/3` and fixed RTK solvers now use
  non-reference satellites on the epochs where they are available instead of
  dropping a satellite from the entire arc when it is absent from one epoch. The
  reference satellite is still required across the arc.

### Fixed

- GNSS integer ambiguity fixing no longer treats a missing runner-up lattice
  candidate as infinite ratio confidence; one-candidate searches now return
  `integer_status: :not_fixed`.
- `Orbis.GNSS.SP3.position/3` (and everything built on it, including
  `Orbis.GNSS.Observables` and the ephemeris sampler) now refuses an epoch
  beyond the product's node coverage with an `epoch out of range` error instead
  of silently extrapolating the interpolation spline to a non-physical position.
  Queries within one sampling step of the ends still interpolate; in-coverage
  results are bit-for-bit unchanged.

## [0.11.0] - 2026-06-08

### Added

- `Orbis.GNSS.PrecisePositioning.solve_fixed_epochs/3` now reports
  `metadata.ambiguity_search` diagnostics (satellite order, float ambiguities,
  ambiguity covariance, and inverse covariance in cycles) so callers can audit
  the LAMBDA integer decision against the same lattice metric.
- `Orbis.GNSS.PrecisePositioning` now accepts `elevation_weighting: true` on
  float, multi-epoch, and fixed solves, scaling code and phase row sigmas by
  `1 / sin(elevation)` for a simple real-data stochastic model that down-weights
  low-elevation observations.
- `Orbis.GNSS.RTK.double_differences/3` for deterministic base/rover
  code-and-carrier double differences, the RTK measurement primitive that
  cancels receiver clocks and common short-baseline satellite errors before
  baseline estimation.
- `Orbis.GNSS.RTK.solve_float_baseline_epochs/3` for static float RTK baseline
  estimation from supplied satellite ECEF positions and multi-epoch
  code/carrier double differences, holding one float ambiguity per
  non-reference double-difference arc. The float solution now exposes the
  double-difference ambiguity covariance and inverse covariance in metres.
- `Orbis.GNSS.RTK.solve_fixed_baseline_epochs/3` for LAMBDA-fixed RTK baseline
  estimation. It starts from the float RTK baseline, fixes double-difference
  carrier ambiguities with the same correlated covariance used by the float
  solve, and re-solves the baseline with those integers held fixed.
- `Orbis.GNSS.RTK.solve_fixed_baseline_epochs/3` now accepts
  `ambiguity_offset_m`, so fixed RTK ambiguities can be modeled as
  `offset + integer * wavelength`. This is the hook needed for
  wide-lane-fixed / narrow-lane dual-frequency RTK workflows.
- `Orbis.GNSS.RTK.solve_widelane_fixed_baseline_epochs/3` for dual-frequency
  RTK fixing. It estimates Melbourne-Wubbena wide-lane double-difference
  integers, converts the arc to ionosphere-free narrow-lane measurements, then
  runs the existing correlated LAMBDA baseline solve with the wide-lane offsets
  held fixed.
- `Orbis.GNSS.RTK.solve_float_baseline_epochs/3` and
  `solve_fixed_baseline_epochs/3` now understand carrier-phase arc identities:
  map observations may carry `:ambiguity_id`, and LLI loss-of-lock can be
  handled with `on_cycle_slip: :error | :drop_satellite | :split_arc`. Split
  arcs reset the affected double-difference ambiguity while residuals keep the
  physical satellite id.
- `Orbis.GNSS.RTK.solve_float_baseline_epochs/3` and
  `solve_fixed_baseline_epochs/3` now accept `elevation_weighting: true`, which
  scales each undifferenced measurement sigma by
  `1 / max(sin(elevation), 0.05)` before propagating the correlated
  double-difference covariance.
- `Orbis.GNSS.PrecisePositioning.solve_widelane_fixed_epochs/3` now supports
  `on_cycle_slip: :split_arc`, which resets a satellite's carrier ambiguity at
  detected cycle slips and keeps any post-slip fragments long enough for
  wide-lane fixing. Split fragments are reported in
  `metadata.split_cycle_slip_arcs` and use suffixed ambiguity ids such as
  `"G21#2"` in `used_sats` and the ambiguity maps.

### Changed

- `Orbis.GNSS.PrecisePositioning.solve_fixed_epochs/3` now uses an
  LDL-consistent forward recursion for the decorrelated LAMBDA sphere search.
  This fixes the zero-candidate search miss on noisy real arcs without an
  original-space substitute path: those arcs now return a `FixedSolution` with
  `metadata.integer_status == :not_fixed` when candidates exist but fail the
  ratio test.
- `Orbis.GNSS.RTK.solve_float_baseline_epochs/3` now propagates the
  non-diagonal double-difference measurement covariance into the normal
  equations and ambiguity covariance instead of treating DD rows that share a
  reference satellite as independent.
- `Orbis.GNSS.RTK.solve_float_baseline_epochs/3` now chooses the
  highest-average-elevation common satellite as the default reference, with a
  deterministic satellite-id tie-break. `double_differences/3` still defaults to
  the lexicographically first common satellite because it has no geometry.

## [0.10.0] - 2026-06-07

### Added

- `Orbis.GNSS.IonosphereFree.iono_free_phase/4` and
  `iono_free_phase_cycles/4` for PPP/RTK-facing first-order ionosphere-free
  carrier-phase combinations, plus `Orbis.GNSS.CarrierPhase.phase_meters/2`,
  `code_minus_carrier/3`, and `smooth_iono_free_code/2` for code-carrier
  diagnostics and dual-frequency divergence-free Hatch smoothing.
- `Orbis.GNSS.PrecisePositioning.solve_float/4`, a first float-ambiguity
  carrier-phase estimator for one SP3-backed epoch from ionosphere-free code and
  phase observations. It estimates receiver ECEF position, clock, and one float
  ambiguity per satellite, exposing residuals and metadata for later PPP/RTK
  layers.
- `Orbis.GNSS.PrecisePositioning.solve_float_epochs/3`, a static multi-epoch
  float carrier-phase estimator that holds one ambiguity per satellite across an
  arc while estimating one receiver clock per epoch. This is the bridge from
  single-epoch float positioning toward PPP/RTK ambiguity fixing.
- `Orbis.GNSS.PrecisePositioning.solve_fixed_epochs/3`, an integer-fixed
  multi-epoch carrier-phase estimator. It starts from the float arc, builds the
  ambiguity covariance from the float normal matrix, runs LAMBDA integer
  decorrelation plus a covariance-weighted integer sphere search on explicit
  caller-supplied wavelengths, then re-solves receiver position and epoch clocks
  with the selected ambiguities held fixed. The fixed solution reports the
  integer method, ratio-test status, weighted scores, and evaluated candidate
  count.
- `Orbis.GNSS.PrecisePositioning.solve_widelane_fixed_epochs/3`, a
  dual-frequency convenience layer that fixes Melbourne-Wubbena wide-lane
  integers first, then uses LAMBDA on the remaining narrow-lane integer while
  returning both ambiguity sets.
- `Orbis.GNSS.PrecisePositioning` can now apply an opt-in a-priori
  Saastamoinen/Niell tropospheric slant delay to ionosphere-free code and phase
  observations (`troposphere: true` with surface meteorology options), including
  the float, multi-epoch, and fixed-ambiguity solve paths.
- `Orbis.GNSS.PrecisePositioning.solve_float_epochs/3` and
  `solve_fixed_epochs/3` can now estimate one residual zenith troposphere delay
  over a static arc (`estimate_ztd: true`, with `troposphere: true`), reporting
  `ztd_residual_m` and `metadata.ztd_estimated`.
- `Orbis.GNSS.PrecisePositioning.solve_widelane_fixed_epochs/3` accepts
  `on_cycle_slip: :drop_satellite` to remove slipped satellite arcs before the
  wide-lane / narrow-lane solve. The default remains `:error`; dropped satellites
  are reported in `metadata.dropped_cycle_slip_sats`.

### Changed

- `Req` is now a required dependency. Network-backed features (`CelesTrak`,
  `Orbis.GNSS.Data`, NAVCEN constellation status) are first-class Orbis
  capabilities, and making the HTTP client required keeps consumer compiles
  warning-free.
- The LAMBDA integer search now shrinks its live search bound to the current
  second-best candidate, so `solve_fixed_epochs/3` keeps the same integer
  decision and ratio-test semantics while visiting far fewer complete
  candidates.
- `Orbis.GNSS.PrecisePositioning.solve_fixed_epochs/3` now reports an empty
  LAMBDA sphere-search result as `{:error, {:no_integer_candidates, count}}`
  instead of conflating it with the `:too_many_integer_candidates` cap.

## [0.9.2] - 2026-06-06

### Added

- `Orbis.GNSS.Constellation.diff/2` and `changed?/1` for deterministic
  snapshot-to-snapshot catalog comparisons keyed by `{system, prn}`. The diff
  reports added/removed PRNs plus NORAD, SP3 id, SVN, activity, and usability
  changes in structured lists.
- GLONASS FDMA carrier-phase wavelengths. `Orbis.GNSS.RINEX.Observations`
  exposes the parsed `GLONASS SLOT / FRQ #` channel map and `phases/3` now
  derives carrier frequency, G1/G2 wavelengths, and metre phases for GLONASS
  satellites with a channel entry, so `Orbis.GNSS.CarrierPhase` can process
  real GLONASS phase arcs instead of skipping them.
- `Orbis.GNSS.ReducedOrbit` and `Orbis.GNSS.ReducedOrbit.Piecewise` can now fit
  and drift against `%Orbis.Elements{}` TLE/OMM sources by sampling SGP4 over the
  requested window (TEME → GCRS → ECEF, UTC scale). This closes the LEO reduced
  orbit source path without changing the Rust reduced-orbit numerics.

## [0.9.1] - 2026-06-05

### Added

- Rustler precompiled-NIF packaging support. Release tags now build GitHub
  Release archives for common Linux/macOS/Windows targets, and the Hex package
  will include `checksum-*.exs` so supported users do not need a local Rust
  toolchain. If no checksum file is present, Orbis source-builds instead of
  trying to download missing assets; `ORBIS_BUILD=1` remains the explicit
  source-build escape hatch.
- **`Orbis.GNSS.CarrierPhase`** — dual-frequency carrier-phase combinations and
  the quality tooling on them: geometry-free (`L1 - L2`), wide-lane wavelength,
  narrow-lane code, Melbourne-Wübbena, arc-wise cycle-slip detection (LLI bit,
  geometry-free step, and Melbourne-Wübbena step, with documented thresholds),
  and the single-frequency Hatch carrier-smoothed code (with slip/LLI reset).
  GPS/Galileo/BeiDou; GLONASS satellites are skipped (FDMA wavelengths not yet
  derived). Builds on the newly exposed phase observations; no crate change.
- `Orbis.GNSS.RINEX.Observations.values/3` and `phases/3` — expose the raw RINEX
  observations for an epoch (pseudorange, carrier phase, Doppler, signal strength
  with their LLI/SSI), and a carrier-phase convenience that adds the wavelength
  and the phase in metres for GPS/Galileo/BeiDou bands (`band_frequency_hz/2` is
  public; GLONASS FDMA wavelengths are not yet derived). `values/3` takes a
  `:codes` per-system filter so only the requested systems/codes cross the NIF
  boundary. This unlocks carrier-phase combinations without a parser change.
- `Orbis.GNSS.Constellation.validate_sp3!/2` — a build-time validation gate that
  returns `:ok` or raises `ArgumentError` describing the findings (e.g. a
  stale-active PRN that is active and usable in the catalog but missing from a
  current SP3 product). Intended for catalog-build automation, not the runtime.
- Python/georinex/scipy oracle gates for the recent Orbis-only GNSS layer:
  raw RINEX `values/3` / `phases/3`, `CarrierPhase` combinations/slip/Hatch
  smoothing, `IonosphereFree` coefficients and combinations, `GNSS.QC`
  weighting/chi-square thresholds, `GNSS.Observables.predict/5`, C/A
  code/correlation/acquisition, LNAV parity/subframe synthesis,
  visibility/DOP, velocity, DGNSS, `SolutionReport`, and `ReducedOrbit` /
  `ReducedOrbit.Piecewise` fit/evaluation/drift against Astropy/scipy.

### Changed

- `Orbis.GNSS.Constellation.to_csv/2` gains a `:booleans` option: `:lower`
  (default, conventional `true`/`false`) or `:title` (`True`/`False`, for a
  pandas-style consumer that reads the `active` column as Python booleans).
- `Orbis.GNSS.QC.chi2_inv/2` now inverts the regularized-gamma chi-square CDF
  and is checked against `scipy.stats.chi2.ppf`, replacing the older
  Wilson-Hilferty approximation.

## [0.9.0] - 2026-06-05

A large GNSS expansion — signal generation, measurement modelling, velocity,
quality control, and differential positioning — alongside a consolidation of
the whole GNSS surface under the `Orbis.GNSS.*` namespace.

### Added

- **`Orbis.GNSS.Signal.CA`** — GPS L1 C/A Gold-code generation, chip indexing,
  and auto/cross-correlation (IS-GPS-200 G1/G2 generators and per-PRN taps).
- **`Orbis.GNSS.Signal.Correlator`** — C/A code+carrier replica, coherent
  correlation, a 2-D code-phase/Doppler acquisition search, and the
  coherent-integration (sinc²) loss model.
- **`Orbis.GNSS.Navigation.LNAV`** — GPS LNAV subframe synthesis and decoding:
  TLM/HOW, time-of-week, subframe parity (IS-GPS-200 Table 20-XIV), and
  ephemeris bit-packing.
- **`Orbis.GNSS.Observables`** — predicted geometric range, range-rate, Doppler,
  satellite clock, elevation, and azimuth from a receiver position and an SP3
  ephemeris, with light-time (transmit-time) and Sagnac corrections.
- **`Orbis.GNSS.Geometry`** — satellite visibility above an elevation mask,
  dilution of precision (GDOP/PDOP/HDOP/VDOP/TDOP), DOP/visibility time series,
  and rise/set passes.
- **`Orbis.GNSS.Velocity`** — receiver velocity and clock drift from Doppler or
  pseudorange-rate measurements by least squares over the line-of-sight geometry.
- **`Orbis.GNSS.QC`** — measurement quality control: residual-based RAIM fault
  detection, leave-one-out fault detection and exclusion (FDE), and
  elevation/C-N₀ measurement weighting.
- **`Orbis.GNSS.IonosphereFree`** — the dual-frequency ionosphere-free
  pseudorange combination, with standard per-system frequency pairs
  (GPS L1/L2, Galileo E1/E5a, BeiDou B1I/B3I).
- **`Orbis.GNSS.DGNSS`** — code-differential positioning: base-station
  pseudorange corrections and corrected rover solves that cancel the errors
  common to both receivers (satellite clock, ephemeris, short-baseline
  atmosphere).
- **`Orbis.GNSS.SolutionReport`** — a per-satellite and summary diagnostic over
  a position solution: elevation/azimuth, post-fit and RAIM-normalized
  residuals, DOP, residual RMS, and the integrity verdict.
- **`Orbis.GNSS.ReducedOrbit.Piecewise`** — a piecewise (segmented)
  reduced-orbit model that tiles a span into contiguous fitted segments for
  tighter caching/transport accuracy than a single mean-element fit.

### Changed

- **Breaking:** GNSS modules now live under the `Orbis.GNSS.*` namespace. The
  old top-level GNSS names (`Orbis.SP3`, `Orbis.PointPositioning`,
  `Orbis.GnssData`, etc.) were removed instead of retained as aliases, matching
  the library's current single-client / pre-broad-adoption status. Examples:
  `Orbis.GNSS.SP3`, `Orbis.GNSS.Positioning`, `Orbis.GNSS.Data`,
  `Orbis.GNSS.RINEX.Observations`, `Orbis.GNSS.ReducedOrbit`,
  `Orbis.GNSS.Signal.CA`, and `Orbis.GNSS.Navigation.LNAV`.
- Internal GNSS implementation helpers were consolidated under
  `Orbis.GNSS.Core` for shared constants, ECEF input normalization,
  epoch/window handling, validation, source sampling, and versioned-map guards.
- Hardened public-API input validation across the GNSS modules: malformed
  receiver/base positions, out-of-range RAIM options, sub-second piecewise
  segment lengths, out-of-range LNAV flags, and duplicate observations now
  return tagged errors (or raise a clear `ArgumentError` for invalid options)
  instead of crashing, looping, or silently truncating.

## [0.8.0] - 2026-06-05

Observation parsing and a compact orbit model. Orbis can now read a station's
RINEX observation file end-to-end into pseudoranges, and distill a position
track into a tiny, transportable mean-element model.

### Added

- **`Orbis.GNSS.RINEX.Observations`** — RINEX 3 observation parsing with Hatanaka (CRINEX 1.0
  and 3.0) decoding. Decodes `.crx`/`.rnx`, exposes the header (incl. the
  surveyed `APPROX POSITION`), observation codes, and epochs, and extracts
  single-frequency pseudoranges (`pseudoranges/3`) in the
  `[{satellite_id, range_m}]` shape `Orbis.GNSS.Positioning.solve/4` consumes —
  closing the loop from a station's observation file to a recovered position.
  `Orbis.GNSS.Data` gains a station observation product fetch and an
  `observations/2` loader. CRINEX decoding is verified byte-for-byte against
  `crx2rnx`; an end-to-end test recovers a surveyed station position to metre
  level from real GPS observations.

- **`Orbis.GNSS.ReducedOrbit`** — a compact, fitted mean-element approximation of an
  orbit for caching, transport, and quick visibility math (not orbit
  determination). Fits from an `Orbis.GNSS.SP3` track or a list of ECEF samples;
  evaluates position/velocity (ECEF by default, GCRS on request); reports a
  source-backed `drift/3` against the source ephemeris; and serialises to a
  stable, versioned map (`to_map/1`/`from_map/1`). Two models: `:circular_secular`
  (default) and `:eccentric_secular` (nonsingular `h = e·sin ω`, `k = e·cos ω`),
  the latter recovering the radial `a·e` signal that the circular model discards —
  cutting full-day extrapolation error by one-to-three orders of magnitude for
  GPS and BeiDou while matching the circular model on near-circular Galileo.

## [0.7.0] - 2026-06-04

GNSS positioning. Orbis can now recover a receiver position from pseudoranges
against precise or broadcast ephemeris, with the supporting ephemeris,
correction, time, and data-fetch layers.

### Added

- **`Orbis.GNSS.Positioning`** — single-point positioning (SPP). Solves a
  receiver position, clock, and geometry diagnostics from one epoch of
  pseudoranges against either an `Orbis.GNSS.SP3` precise product or an
  `Orbis.GNSS.Broadcast` handle. Multi-constellation
  (GPS / Galileo / BeiDou / GLONASS) solves carry one receiver clock per system;
  the solution reports position, geodetic position, per-system clocks, DOP,
  residuals, used/rejected satellites, and solver metadata.
- **`Orbis.GNSS.SP3`** — SP3-c/SP3-d precise orbit/clock loading and arbitrary-epoch
  satellite position/clock interpolation, plus `satellite_ids/1` to read the
  product's declared satellite set.
- **`Orbis.GNSS.Constellation`** — a GPS constellation catalog built from
  CelesTrak `gps-ops` OMM identity and an optional NAVCEN status/SVN overlay
  (PRN ↔ SVN ↔ NORAD ↔ SP3 id, active/usable flags). Merges sources only when
  the block type matches, recording PRN-transition disagreements as conflicts
  rather than corrupting identity; exports the compact mapping CSV and validates
  a catalog (duplicate PRNs/NORAD ids, inactive/unusable PRNs, and missing/extra
  satellites against a loaded `Orbis.GNSS.SP3` product).
- **`Orbis.GNSS.Broadcast`** — RINEX 3.x and 4.xx navigation parsing and
  broadcast orbit/clock evaluation: GPS LNAV, Galileo I/NAV and F/NAV, BeiDou
  D1/D2 (including geostationary satellites), and GLONASS (PZ-90.11 state-vector
  propagation by Runge–Kutta integration).
- **`Orbis.GNSS.Ionosphere`** (broadcast Klobuchar, frequency-aware across L1/E1/B1I)
  and **`Orbis.GNSS.Troposphere`** (Saastamoinen zenith delay + Niell mapping)
  correction models.
- **`Orbis.GNSS.Data`** — an optional product fetch/cache layer: a catalog over
  public archives, HTTPS (`Req`) and FTP downloads, an atomic on-disk cache with
  SHA-256 integrity and provenance sidecars, a gzip-bomb guard, and an offline
  mode. Includes convenience loaders that return `Orbis.GNSS.SP3` /
  `Orbis.GNSS.Broadcast` handles. `Req` is an optional dependency.
- **`Orbis.GNSS.Time`** — GNSS epoch/seconds-of-week and day-of-year helpers.

### Notes

- The GNSS numerical core lives in the Rust `astrodynamics` / `astrodynamics-gnss`
  crate layer. Its libm-bound components (orbit and clock evaluation, ionosphere,
  troposphere, dilution of precision) are held to bit-exact (0 ULP) parity
  against pinned Python references; broadcast orbits are additionally validated
  against precise SP3 products. The least-squares solver's final position is a
  sub-micron solver-agreement result, not a 0-ULP claim.

---

Releases before 0.7.0 predate this changelog.
