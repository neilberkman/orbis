# RTK oracle generators

Reproducible RTKLIB references for the RTK filter kernel parity gates. The
oracle JSONs in the parent directory are generated from the vendored RINEX
observations and broadcast nav by RTKLIB's `rnx2rtkp`, then converted to the
shared per-epoch JSON shape by `pos_to_oracle.py`.

## Provenance

- **Oracle binary:** RTKLIB `rnx2rtkp` **v2.4.2-p13**, commit `71db0ff`, built in
  this workspace at `_tools/RTKLIB/app/rnx2rtkp/gcc/rnx2rtkp`.
- **Inputs** (all vendored under `test/fixtures/`):
  - rover obs: `obs/WTZZ00DEU_R_20201770000_01D_30S_MO_120epoch.rnx`
  - base obs: `obs/WTZR00DEU_R_20201770000_01D_30S_MO_120epoch.rnx`
  - nav: `nav/ESBC00DNK_R_20201770000_01D_MN.rnx` (mixed broadcast)
- **Base reference position** (ARP): `49.144200524 12.878913935 666.0917`, the
  same fixed base used by the canonical `wtzr_wtzz_rtklib_oracle.json`.
- **Truth:** the static antenna baseline of the co-located WTZR/WTZZ pair
  (receivers do not move) — copied verbatim into every oracle's `truth` block.

Each config below changes exactly **one** variable from the canonical
`l1_brdc_fix_and_hold` reference, so a parity failure isolates to that variable.

## Oracles

| JSON                                            | config                             | changed variable                                      | result                                                                      |
| ----------------------------------------------- | ---------------------------------- | ----------------------------------------------------- | --------------------------------------------------------------------------- |
| `wtzr_wtzz_kinematic_gps_rtklib_oracle.json`    | `track_a_kinematic_gps_l1.conf`    | `posmode` static→**kinematic** (GPS L1)               | 119/120 fixed, ~7mm converged; epoch-0 kinematic cold-start transient (~1m) |
| `wtzr_wtzz_multignss_static_rtklib_oracle.json` | `track_b_static_multignss_l1.conf` | `navsys` GPS→**GPS+GLO+GAL+BDS** (GLONASS float-only) | 120/120 fixed, ~2mm, 9–11 sats                                              |

**GLONASS is float-only** (`pos2-gloarmode=off`): FDMA inter-channel biases break
the clean double-difference integer assumption, so GLONASS contributes to the
float solution but not to ambiguity resolution. FDMA AR is a non-goal until a
gate proves the win.

## Reproduce

```sh
RNX=../../../../_tools/RTKLIB/app/rnx2rtkp/gcc/rnx2rtkp
OBS=../../obs; NAV=../../nav/ESBC00DNK_R_20201770000_01D_MN.rnx
$RNX -k track_a_kinematic_gps_l1.conf -o track_a.pos \
  $OBS/WTZZ00DEU_R_20201770000_01D_30S_MO_120epoch.rnx \
  $OBS/WTZR00DEU_R_20201770000_01D_30S_MO_120epoch.rnx $NAV
python3 pos_to_oracle.py track_a.pos track_a_kinematic_gps_l1.conf \
  kinematic_gps_l1_fix_and_hold "..." ../wtzr_wtzz_kinematic_gps_rtklib_oracle.json
```

## Next: kinematic moving-rover arcs (queued, not yet vendored)

The real Track A gate is genuinely-moving rovers (GSDC drive logs: moving phones,
published ground-truth CSVs). Those arcs use the **demo5** RTKLIB fork as the
oracle binary — the de-facto reference for low-cost/phone kinematic RTK that
vanilla 2.4.2 handles poorly — vendored with its own version/commit provenance.
Truth swaps from a surveyed baseline to a per-epoch trajectory; the JSON shape is
unchanged.
