# spp-robustness pre-registration spec

Pre-registered BEFORE measuring. Capability: single-frequency (L1 C/A) SPP
robustness on degraded receiver data (poor and variable C/N0, multipath
outliers, weak or near-degenerate geometry).

## Scope (house rules)

Elixir reference and measurement ONLY. The Rust kernel, the NIF, and the
astrodynamics crate are out of scope and are not modified. Every new option is
additive and defaults to current behaviour (`nil`/`off`), so the existing
bit-exact solve path is unchanged when no robustness option is passed.

## Gap under test

The robustness machinery already exists in `Orbis.GNSS.QC` (RAIM chi-square
fault detection, leave-one-out FDE in `QC.fde/4`, elevation and C/N0 variance
and inverse-variance weights). It is NOT wired into the `Orbis.GNSS.Positioning.solve/4`
path that callers use, and `solve/4` has no integrity check at all: a multipath
outlier that survives the crate elevation mask biases the position and is
returned as `{:ok, solution}` with no flag. There is also no Elixir guard that
turns a numerically near-degenerate but full-rank geometry (huge PDOP) into a
tagged refusal instead of a silent bad fix.

The classical solve-domain robustness (Huber IRLS reweighting of the solve
residuals, C/N0 into the WLS normal equations) requires the crate to accept
per-iteration weights and is therefore OUT OF SCOPE. It is declared here and not
faked. The implementable Elixir-side win is FDE in the solve path plus a
degenerate-geometry guard.

## Implementation (additive, default-off)

In `Orbis.GNSS.Positioning.solve/4`:

  * `:robust` (default `false`). When set, the solve routes through the existing
    `QC.fde/4` leave-one-out loop instead of a bare crate call, returns the
    cleaned `Solution`, and records the excluded satellites in
    `solution.metadata.fde`. `:p_fa` and `:weights` are forwarded to the RAIM
    test inside FDE. `false` reproduces the current single crate solve
    bit-for-bit.
  * `:max_pdop` (default `nil`). After the solve, if `solution.dop` is `nil`
    (rank-deficient) or `solution.dop.pdop > max_pdop`, return
    `{:error, {:degenerate_geometry, pdop}}` (`pdop` is `nil` for the
    rank-deficient case) instead of `{:ok, bad_fix}`. `nil` disables the guard.

No positioning math is added; FDE and the chi-square test already exist and are
tested on synthetic SP3 data.

## Oracles

1. In-repo, deterministic SP3-synthesized labelled fault (exclusion
   correctness). Clean pseudoranges are synthesized from
   `test/fixtures/sp3/GRG0MGXFIN_20201760000_01D_15M_ORB.SP3` via
   `Observables.predict_all`, exactly as `test/gnss_qc_test.exs` does today, and
   a known bias is injected on one chosen satellite. The biased satellite is the
   ground-truth fault the robust solve must isolate.

2. Position-domain GSDC Pixel-5 demo5/RTKLIB oracles (before/after on real
   degraded data). The four vendored
   `test/fixtures/rtk/gsdc_*_pixel5_p222_demo5_rtklib_oracle.json` files carry
   RTKLIB output position, GSDC truth position, and per-epoch 3D/horizontal error
   vs truth. They are the BAR (median 3D ~3.6-4.0 m, p95 ~7.9 m on mtv1), not
   the input. Raw phone L1 observations come from the staged
   `/tmp/gsdc-work/<drive>/supplemental/gnss_rinex.21o` plus the staged
   broadcast NAV, fed to `Positioning.solve` per epoch and matched to the oracle
   truth by GPST time. If the staged raw observations are not reachable the
   real-data gate is reported BLOCKED; gate 1 is fully in-repo and unaffected.

## Gates (declared up front, never loosened)

GATE 1 (FDE exclusion correctness, deterministic SP3 oracle). Inject a known
bias on one satellite, sweep `{50, 100, 200, 500}` m, in an otherwise clean
GPS-only set (~7 sats at the fixture epoch). PASS bar, declared:

  * the robust solve excludes EXACTLY the biased satellite at the bias levels
    where the fault localizes (`excluded == [{biased, :raim_excluded}]`);
  * the recovered 3D position error is within `clean_3d_error + 0.5 m`;
  * zero false exclusions on the clean set (`:robust` returns
    `excluded == []`, position bit-identical to the bare solve).

GATE 2 (before/after on real GSDC degraded data, position-domain oracle). For
each vendored arc, run orbis SPP per matched epoch as (A) bare crate solve and
(B) robust solve (`:robust true`). Match to oracle truth by GPST. Metrics per
arc and pooled: 3D median, 3D p95, horizontal median and p95 error vs truth.

  * Population floor: >= 100 matched epochs per arc (arcs carry ~1465).
  * Credibility floor: orbis bare 3D median must be within 2x the demo5 oracle
    3D median, else the arc is reported as an unusable substrate, not massaged.
  * Robustness claim (declared): robust 3D median <= bare 3D median AND robust
    3D p95 <= bare 3D p95 on every usable arc, with the gain expected in the
    p95 tail (multipath/outlier removal). Equal-or-worse is a documented NULL
    result (the crate elevation weight already handles most gross multipath),
    reported, not hidden.

GATE 3 (degenerate-geometry guard). On a full-rank but high-PDOP synthetic
geometry, `:max_pdop` below the realized PDOP returns
`{:error, {:degenerate_geometry, pdop}}`; a generous threshold returns
`{:ok, solution}` identical to the unguarded solve. `:max_pdop nil` (default)
never refuses.

## Never weaken

The elevation mask, SIGMA0, and every existing assertion are untouched. An
out-of-tolerance or no-improvement outcome is a fail/null result reported with
oracle, metric, sample size, and tolerance. No tolerance gate is loosened to
manufacture a pass.
