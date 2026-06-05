# SP3 test-fixture provenance

## `GRG0MGXFIN_20201760000_01D_15M_ORB.SP3`

IGS MGEX final combined precise orbit + clock product (CNES/CLS/GRGS), 2020
day-of-year 176, 15-minute grid, GPS time, carrying GPS / GLONASS / Galileo (no
BeiDou). Redistributed public IGS product. Used by the SP3 interpolation and
reduced-orbit tests.

## `GBM_BDS_C21_C08_trim.sp3`

Derived from the GFZ rapid MGEX product
`GBM0MGXRAP_20201770000_01D_05M_ORB.SP3` (2020 day-of-year 177, 5-minute grid,
GPS time) by keeping the verbatim header and only the position records for two
BeiDou satellites — **C21** (MEO, e ≈ 9e-4) and **C08** (IGSO, e ≈ 5e-3) — across
all 288 epochs. Other satellites' records were dropped (the header still lists
the full constellation, which the SP3 reader tolerates); no values were altered.

- size 72293 bytes, sha256
  `f77d83a0da91e7112c2890ba7aae29326b8c621cfee58ac18e4243d86e40238b`.
- Source product: `ftp://ftp.gfz-potsdam.de/pub/GNSS/products/mgex/2111/`.
- Purpose: the real BeiDou drift gate for `Orbis.GNSS.ReducedOrbit`'s
  `:eccentric_secular` model (the GRG product carries no BeiDou). GEO satellites
  (C01–C05) are excluded — near-equatorial and not orbital-element-friendly. An
  identical copy lives in the `astrodynamics-gnss` crate fixtures for the same
  gate at the Rust layer.

## `degenerate_coincident_5sat.sp3`

Hand-authored rank-deficient fixture (five GPS satellites at one ECEF point) for
the graceful-degeneracy path; not a redistributed product.
