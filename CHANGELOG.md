# Changelog

All notable changes to this project are documented here. The format is based on
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project
adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

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
