# Accuracy & Validation

Orbis validates against established reference implementations at every
layer. Here is exactly what is tested, to what precision, and against
what oracle.

## Coordinate Transforms: 0 ULP vs Skyfield

The TEME→GCRS→ITRS pipeline produces **IEEE 754 bit-identical output**
(0 ULP — zero Units in the Last Place) to Python
[Skyfield](https://rhodesmill.org/skyfield/) on all tested platforms.

| Transform | Reference | Precision | Verified in CI |
|-----------|-----------|-----------|----------------|
| TEME→GCRS | Skyfield 1.49 | 0 ULP | Yes |
| GCRS→ITRS | Skyfield 1.49 | 0 ULP | Yes |
| Geodetic | Skyfield 1.49 | 0 ULP | Yes |
| Topocentric | Skyfield 1.49 | 0 ULP | Yes |

The transform includes IAU2000A nutation (1365 terms), IAU2006
precession, frame bias, and precise time scale conversions
(UTC→TAI→TT→TDB→UT1). The FMA discipline in `mat3_vec3_mul` matches
numpy's vectorized behavior.

**Test tag:** `:skyfield_parity` — run with `mix test --include skyfield_parity`

## SGP4 Propagation: < 1 mm vs Skyfield

SGP4 propagation uses the published
[`sgp4`](https://crates.io/crates/sgp4) Rust crate in AFSPC
compatibility mode. This is a clean-room Rust implementation that
differs from Skyfield's bundled C extension by sub-nanometer amounts
(different floating-point expression evaluation order and FMA
contractions).

The oracle test verifies position distance is < 1 mm for the ISS at
274 minutes from epoch.

## Orbit Determination: Vallado Reference

IOD methods are validated against David Vallado's
[valladopy](https://github.com/CelesTrak/fundamentals-of-astrodynamics)
Python implementation using the textbook examples.

| Method | Reference | Precision |
|--------|-----------|-----------|
| Gibbs (Algorithm 54) | Vallado Example 7-3 | 0 ULP |
| Herrick-Gibbs (Algorithm 55) | Vallado Example 7-4 | 0 ULP |
| Gauss angles-only (Algorithm 52) | Vallado Example 7-2 | 1e-12 relative |
| Lambert/Battin (Algorithm 61) | Vallado test suite | 1e-12 relative |

The deterministic methods (Gibbs, Herrick-Gibbs) match at the bit level.
The iterative methods (Gauss, Lambert) converge to the same result
within 1e-12 relative tolerance, which is the tolerance used by
Vallado's own test suite.

## Conjunction Assessment

Validated against the Iridium 33 / Cosmos 2251 collision of 2009-02-10:

- **Time of closest approach:** within 1 minute of the known collision time
- **Miss distance:** ~1.9 km (consistent with SGP4/TLE accuracy limits)

## RF Primitives

FSPL uses the standard inverse square law formula:

    FSPL = 32.45 + 20·log₁₀(f_MHz) + 20·log₁₀(d_km)

Verified from first principles (`20·log₁₀(4πd/λ)`) and independently
computed in Python. All other RF functions (EIRP, C/N₀, link margin)
are textbook dB arithmetic with exact analytical solutions.

## JPL Ephemeris

SPK/BSP reader tested at 0 ULP vs Skyfield for Mars, Venus, Sun, and
Moon positions from Earth.

**Test tag:** `:spk_file` — requires `/tmp/de421.bsp`

## What Is NOT Validated

- **Atmospheric density (NRLMSISE-00):** Implemented from the public-domain
  C translation but not cross-validated against a reference implementation.
  Results are physically reasonable but not precision-verified.

- **Pass prediction timing:** Uses the full topocentric pipeline for
  elevation, but rise/set times depend on the scan resolution
  (`step_seconds` option) and bisection precision.
