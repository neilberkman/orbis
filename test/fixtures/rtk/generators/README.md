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
  - canonical/Track A nav: `nav/ESBC00DNK_R_20201770000_01D_MN.rnx`
    (mixed broadcast, filtered to GPS/Galileo/BeiDou)
  - Track B nav: `nav/BRDC00WRD_R_20201770000_01D_GREC.rnx`
    (BKG combined broadcast, filtered to GPS/GLONASS/Galileo/BeiDou)
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
| `wtzr_wtzz_multignss_static_rtklib_oracle.json` | `track_b_static_multignss_l1.conf` | `navsys` GPS→**GPS+GLO+GAL+BDS** (GLONASS float-only) | 120/120 fixed, ~1.8mm, 14–17 sats                                           |

**GLONASS is float-only** (`pos2-gloarmode=off`): FDMA inter-channel biases break
the clean double-difference integer assumption, so GLONASS contributes to the
float solution but not to ambiguity resolution. FDMA AR is a non-goal until a
gate proves the win.

## Reproduce

```sh
RNX=../../../../../../_tools/RTKLIB/app/rnx2rtkp/gcc/rnx2rtkp
OBS=../../obs; NAV=../../nav/ESBC00DNK_R_20201770000_01D_MN.rnx
$RNX -k track_a_kinematic_gps_l1.conf -o track_a.pos \
  $OBS/WTZZ00DEU_R_20201770000_01D_30S_MO_120epoch.rnx \
  $OBS/WTZR00DEU_R_20201770000_01D_30S_MO_120epoch.rnx $NAV
python3 pos_to_oracle.py track_a.pos track_a_kinematic_gps_l1.conf \
  kinematic_gps_l1_fix_and_hold "..." ../wtzr_wtzz_kinematic_gps_rtklib_oracle.json

NAV=../../nav/BRDC00WRD_R_20201770000_01D_GREC.rnx
$RNX -k track_b_static_multignss_l1.conf -o track_b_static_multignss_l1.pos \
  $OBS/WTZZ00DEU_R_20201770000_01D_30S_MO_120epoch.rnx \
  $OBS/WTZR00DEU_R_20201770000_01D_30S_MO_120epoch.rnx $NAV
python3 pos_to_oracle.py track_b_static_multignss_l1.pos track_b_static_multignss_l1.conf \
  static_multignss_grec_l1_fix_and_hold "..." ../wtzr_wtzz_multignss_static_rtklib_oracle.json
```

Track B was regenerated with `rnx2rtkp -y 2` as a scratch audit while keeping
the committed JSON shape unchanged. The RTKLIB status trace reports 14-17 total
satellites over the 120 epochs, including 5-6 GLONASS satellites per epoch
(GPS 5-6, Galileo 3-5, BeiDou 0 on this RTKLIB 2.4.2 L1 arc).

## Next: kinematic moving-rover arcs (queued, not yet vendored)

The real Track A gate is genuinely-moving rovers (GSDC drive logs: moving phones,
published ground-truth CSVs). Those arcs use the **demo5** RTKLIB fork as the
oracle binary — the de-facto reference for low-cost/phone kinematic RTK that
vanilla 2.4.2 handles poorly — vendored with its own version/commit provenance.
Truth swaps from a surveyed baseline to a per-epoch trajectory; the JSON shape is
unchanged.
