# Orbis Architecture

Satellite toolkit for Elixir — orbit propagation, coordinate transformations, orbit determination, conjunction assessment, and ground station operations. Rust NIF backend with bit-exact Skyfield parity on coordinate transforms.

## Module Structure

```
Orbis                         — public API entry point
Orbis.TLE                     — TLE struct and parser
Orbis.SGP4                    — orbit propagation (TLE → TEME state vectors)
Orbis.Coordinates             — frame transformations (TEME → GCRS → ITRS → geodetic → topocentric)
Orbis.Passes                  — ground station visibility / access windows
Orbis.Doppler                 — range rate and Doppler shift
Orbis.Eclipse                 — sunlit/penumbra/umbra status
Orbis.Atmosphere              — NRLMSISE-00 atmospheric density model
Orbis.Ephemeris               — JPL SPK/BSP reader (Sun, Moon, planets)
Orbis.Angles                  — Sun/Moon angles, phase angle, Earth angular radius
Orbis.NIF                     — private Rust NIF bindings (not part of public API)
```

## Rust NIF Architecture

The NIF is implemented in Rust via Rustler (`native/orbis_nif/`):

```
src/
  lib.rs                      — NIF entry points
  propagation.rs              — SGP4 propagation (wraps sgp4 crate)
  coordinates.rs              — TEME→GCRS→ITRS pipeline
  nutation.rs                 — IAU2000A nutation (678 lunisolar + 687 planetary terms)
  precession.rs               — IAU2006 precession + frame bias
  matrix.rs                   — matrix operations (Kahan-compensated triple product)
  time_scales.rs              — UTC→TAI→TT→TDB→UT1 conversions
  iod.rs                      — Gibbs and Herrick-Gibbs orbit determination
  gauss.rs                    — Gauss angles-only IOD
  lambert.rs                  — Lambert solver (Battin method)
  conjunction.rs              — Closest approach finder
  atmosphere.rs               — NRLMSISE-00 atmospheric density
  doppler.rs                  — Doppler shift computation
  ephemeris.rs                — JPL SPK/BSP ephemeris reader
  iau2000a_data.rs            — nutation coefficient tables (IERS/ERFA source)
  iers_data.rs                — embedded IERS delta-T / UT1-UTC table
```

### FMA Discipline

Bit-exact Skyfield parity requires specific floating-point behavior:

- **Most code**: normal Rust arithmetic. Rust does NOT fuse multiply-add by default.
- **mat3_vec3_mul only**: uses explicit `f64::mul_add()` for fused multiply-add, matching numpy's vectorized behavior.

## Key Principles

1. **Skyfield parity**: Coordinate transforms (TEME→GCRS→ITRS) produce 0 ULP output vs Skyfield.

2. **Vallado parity**: IOD methods (Gibbs, Herrick-Gibbs) match Vallado's Python at 0 ULP. Iterative methods (Gauss, Lambert) match at 1e-12 relative tolerance.

3. **NIF is private**: Users interact with `Orbis` and submodules. `Orbis.NIF` is `@moduledoc false`.

4. **Units**: positions in km, velocities in km/s, angles in radians internally, degrees at the public API boundary.

## Build

```
mix deps.get
mix test                                                    # default tests
mix test --include skyfield_parity --exclude spk_file       # parity tests
```

## Test Tags

- `:skyfield_parity` — bit-exact 0 ULP coordinate transform verification
- `:spk_file` — requires JPL DE421 BSP file at `/tmp/de421.bsp`
