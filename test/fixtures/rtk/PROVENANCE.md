# RTKLIB Wettzell Oracle Fixture

`wtzr_wtzz_rtklib_oracle.json` records RTKLIB `rnx2rtkp` output for the
vendored WTZR/WTZZ 2020-06-25 120-epoch short-baseline arc.

- **Reference implementation:** RTKLIB `rnx2rtkp` version 2.4.2 from the local
  RTKLIB checkout at commit `71db0ff`.
- **Input observations:**
  - `test/fixtures/obs/WTZR00DEU_R_20201770000_01D_30S_MO_120epoch.rnx`
  - `test/fixtures/obs/WTZZ00DEU_R_20201770000_01D_30S_MO_120epoch.rnx`
- **Reference orbit/clock product for the committed RTKLIB solution:** broadcast
  navigation from `test/fixtures/nav/ESBC00DNK_R_20201770000_01D_MN.rnx`.
- **Navigation product used by RTKLIB:** `test/fixtures/nav/ESBC00DNK_R_20201770000_01D_MN.rnx`
- **Primary reference command:** `rnx2rtkp -p 3 -f 1 -h -m 10 ...`
  (kinematic/RTK, L1, fix-and-hold, 10 degree mask). The generated `.pos`
  source was `_planning/rtklib_wettzell/l1_brdc.pos`.

The JSON fixture includes the full per-epoch L1+broadcast fix-and-hold solution
and compact summaries for L1 instantaneous, L1 float, and L1/L2 broadcast
variants.
Truth is the antenna-reference-point baseline derived from the marker ECEF
coordinates and RINEX antenna-height fields already used by
`test/gnss_rtk_real_arc_test.exs`.

Important: passing an SP3 file to `rnx2rtkp` does not by itself select precise
ephemeris; RTKLIB defaults to broadcast (`pos1-sateph = brdc`). The similarly
named `_planning/rtklib_wettzell/modes/l1_sp3.pos` and `l12_sp3.pos` outputs are
numerically identical to the broadcast outputs except for the extra input-file
header line. An explicit `pos1-sateph = precise` run against the trimmed GBM SP3
fixture produced no valid RTK solutions, so this fixture is a broadcast-mode
RTKLIB reference rather than an SP3 parity oracle.
