# Accuracy & Validation

Orbis validates against established reference implementations at every
layer. Here is exactly what is tested, to what precision, and against
what oracle.

## Coordinate Transforms vs Skyfield

The rotation stages of the TEME→GCRS→ITRS pipeline produce **IEEE 754
bit-identical output** (0 ULP — zero Units in the Last Place) to Python
[Skyfield](https://rhodesmill.org/skyfield/) on all tested platforms.
The derived geodetic and topocentric outputs are validated against
Skyfield to tight numerical tolerances, but are **not** bit-exact
today.

| Transform   | Reference     | Precision                              | Verified in CI |
| ----------- | ------------- | -------------------------------------- | -------------- |
| TEME→GCRS   | Skyfield 1.49 | 0 ULP (hex-float refs, pos+vel)        | Yes            |
| GCRS→ITRS   | Skyfield 1.49 | 0 ULP (hex-float refs, isolated input) | Yes            |
| Geodetic    | Skyfield 1.49 | tolerance: 1e-8° lat/lon, 1e-9 km alt  | Yes            |
| Topocentric | Skyfield 1.49 | tolerance: 1e-6° az/el, 1 mm range     | Yes            |

The 0-ULP rows assert IEEE-754 ULP distance against hex-captured
Skyfield reference values. The geodetic and topocentric rows assert
`assert_in_delta` against decimal Skyfield references at the tolerances
shown (`oracle_test.exs`); they are tolerance-tested, not bit-exact.

The transform includes IAU2000A nutation (1365 terms), IAU2006
precession, frame bias, and precise time scale conversions
(UTC→TAI→TT→TDB→UT1). The FMA discipline in `mat3_vec3_mul` matches
numpy's vectorized behavior.

**Test tag:** `:skyfield_parity` — run with `mix test --include skyfield_parity`.
All four rows above are tagged `:skyfield_parity`, and CI runs
`mix test --include skyfield_parity --exclude spk_file`
(see `.github/workflows/ci.yml`), so all four are exercised on every
push and pull request.

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

| Method                           | Reference           | Precision                      |
| -------------------------------- | ------------------- | ------------------------------ |
| Gibbs (Algorithm 54)             | Vallado Example 7-3 | velocity ULP=0\*, angles ULP≤2 |
| Herrick-Gibbs (Algorithm 55)     | Vallado Example 7-4 | velocity ULP=0\*, angles ULP≤2 |
| Gauss angles-only (Algorithm 52) | Vallado Example 7-2 | 1e-12 relative                 |
| Lambert/Battin (Algorithm 61)    | Vallado test suite  | 1e-12 relative                 |

The deterministic methods (Gibbs, Herrick-Gibbs) assert ULP distance on
their outputs (`iod_test.exs`): the velocity components assert `max_ulp = 0`
and the inter-vector angles `theta12`/`theta23` assert `max_ulp = 2`
(near-0-ULP, not strict 0 ULP).

\* The velocity reference values are **decimal literals** transcribed from
valladopy output (e.g. `5.5311472050176125`), not captured `float.hex()`
values. A decimal round-trip backs the ULP=0 velocity assertions, which is
weaker than the hex-float capture discipline used for the coordinate
transforms; treat the velocity ULP=0 as near-bit-exact rather than a
hex-float oracle.

The iterative methods (Gauss, Lambert) converge to the same result
within 1e-12 relative tolerance, which is the tolerance used by
Vallado's own test suite.

## Conjunction Assessment

Validated against the Iridium 33 / Cosmos 2251 collision of 2009-02-10:

- **Time of closest approach:** within ~1 hour of the known collision time (`conjunction_test.exs` asserts `tca_hours` within 1.0 of 22.1)
- **Miss distance:** under 10 km (event-level sanity check, `min_dist < 10.0` km; consistent with SGP4/TLE accuracy limits — not a precise-value oracle)

## RF Primitives

FSPL uses the standard inverse square law formula:

    FSPL = 32.45 + 20·log₁₀(f_MHz) + 20·log₁₀(d_km)

These are tolerance-only checks with **no committed oracle or golden
vector**. `fspl/2` asserts `assert_in_delta 0.01` dB against hand-computed
numbers written in test comments (`rf_test.exs`); `eirp/2` asserts exact
equality on trivial dB arithmetic; the remaining RF functions (C/N₀, link
margin, wavelength, dish gain) assert only sign and monotonicity. All RF
math is textbook dB arithmetic, but it is not cross-validated against an
external reference implementation.

## JPL Ephemeris

The SPK/BSP reader is tested at **0 ULP vs Skyfield for Mars and Venus
only** (Earth-relative positions, 6 components, `max_ulp = 0` against
Skyfield reference values in `oracle_test.exs`; the references are stored
as decimal literals). Sun, Moon, and Earth positions are covered only by
distance-magnitude **sanity bounds** (within ~20,000 km of the expected
lunar distance / 0.02 AU of 1 AU in `ephemeris_test.exs`); they are
**not** 0-ULP verified.

**Test tag:** `:spk_file` — requires `/tmp/de421.bsp`.

**Not verified in CI.** CI runs `mix test --include skyfield_parity
--exclude spk_file`, so the Mars/Venus 0-ULP test is **excluded** from CI.
The `de421.bsp` fixture it depends on is not committed to the repository,
so the only way to exercise this row is locally with the file present:
`mix test --include spk_file` after placing `de421.bsp` at `/tmp/de421.bsp`.

## What Is NOT Validated

- **Atmospheric density (NRLMSISE-00):** Implemented from the public-domain
  C translation but not cross-validated against a reference implementation.
  Results are physically reasonable but not precision-verified.

- **Pass prediction timing:** Uses the full topocentric pipeline for
  elevation, but rise/set times depend on the scan resolution
  (`step_seconds` option) and bisection precision.
