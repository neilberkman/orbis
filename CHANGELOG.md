# Changelog

All notable changes to this project are documented here. The format is based on
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project
adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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
