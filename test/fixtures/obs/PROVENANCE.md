# Observation fixtures — provenance

Small, trimmed RINEX 3 observation fixtures for the CRINEX round-trip and the
end-to-end single-point-positioning test. Committed offline; no test fetches
the network by default.

## ESBC00DNK_R_20201770000_01D_30S_MO_trim.crx / .rnx

- **Station / day:** ESBC00DNK (Esbjerg, Denmark), 2020 day-of-year 177
  (2020-06-25), 30 s sampling, MIXED observation, RINEX 3.05.
- **Upstream source:** the full daily file
  `CRNX/V3/ESBC00DNK_R_20201770000_01D_30S_MO.crx.gz` from the
  `nav-solutions/data` redistribution
  (`https://github.com/nav-solutions/data`), which mirrors the IGS/MGEX
  archive. The same station and day as the committed broadcast navigation
  fixture `test/fixtures/nav/ESBC00DNK_R_20201770000_01D_MN.rnx`.
- **Trim:** decoded with `crx2rnx`, kept the verbatim header plus the first two
  epochs (epoch boundary trim — a CRINEX body cannot be cut mid-stream because
  the difference engines are stateful), then re-compressed with `rnx2crx` so the
  committed `.crx` is self-consistent and re-initializes cleanly at epoch 1. The
  `.rnx` is the `crx2rnx` decode of the committed `.crx`.
- **Reference decoder:** RNXCMP `crx2rnx` / `rnx2crx` version 4.1.0 (the
  `hatanaka` Python package's bundled RNXCMP binaries).
- **Trim commands (equivalent):**

  ```
  gunzip -k ESBC00DNK_R_20201770000_01D_30S_MO.crx.gz
  crx2rnx - < ESBC00DNK_R_20201770000_01D_30S_MO.crx > full.rnx
  head -n 143 full.rnx > trim.rnx        # header + first 2 epochs
  rnx2crx - < trim.rnx > ESBC00DNK_R_20201770000_01D_30S_MO_trim.crx
  crx2rnx - < ESBC..._trim.crx > ESBC00DNK_R_20201770000_01D_30S_MO_trim.rnx
  ```

- **sha256:**
  - `.crx`: `73b2294711f317c20a043c290c5d590917b037caac8feb21d74dc7600f55f5c2`
  - `.rnx`: `c1fbc120be90d7498b3ff138f28ceb7865ce050ddd9cd2f55ce413e553f2e7e0`
- **Header APPROX POSITION XYZ (ECEF m):** `3582105.2910  532589.7313  5232754.8054`
  — the surveyed receiver position the SPP test recovers to metre level.

## algo0010_2015001_v1_trim.crx / .rnx

- **Station / day:** ALGO (Algonquin Park, Canada), 2015 day-of-year 001
  (2015-01-01), 30 s sampling, MIXED (GPS+GLONASS), RINEX 2.11 / CRINEX 1.0.
- **Purpose:** exercises the CRINEX 1.0 (RINEX 2) decode path — the
  12-satellite epoch-line wrap (20 satellites in epoch 1) and the
  five-observations-per-line wrap (8 observation types).
- **Upstream source:** the full daily file
  `gnss/data/daily/2015/001/algo0010.15d.Z` from the ESA GSSC archive
  (`ftp://gssc.esa.int`).
- **Trim:** decompressed, decoded, kept the verbatim header plus the first two
  epochs, then re-compressed with `rnx2crx` so the committed `.crx`
  re-initializes cleanly at epoch 1. The `.rnx` is the `crx2rnx` decode of the
  committed `.crx`.
- **Reference decoder:** RNXCMP `crx2rnx` / `rnx2crx` version 4.1.0 (the
  `hatanaka` Python package's bundled RNXCMP binaries).
- **sha256:**
  - `.crx`: `acc0d16347d28fb5911798f792046b1d32b8177a73b8b8fb4e521fa1fcf0af38`
  - `.rnx`: `f2eae58b37fa267b6f64549de8eb1504473057b46fba1d039b9d2f063b536f22`

The crate carries identical copies under
`crates/astrodynamics-gnss/tests/fixtures/obs/` for the crate's own
CRINEX round-trip and RINEX observation parser tests.
