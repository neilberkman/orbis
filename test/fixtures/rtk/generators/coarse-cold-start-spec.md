# Coarse cold-start convergence basin: pre-registered spec

Pre-registered before measuring. Capability: SPP convergence from a coarse or
degraded position prior, with no hardcoded seed.

## Question

A low-cost tracker often starts from only a rough position prior: region level,
last known hundreds of km away, or none. Characterize `Orbis.GNSS.Positioning`
convergence from a coarse or degraded initial prior, and if the convergence
basin is too small, widen it with an additive, default-off robustification that
needs no hardcoded answer.

## Oracle

ESBC00DNK (Esbjerg, Denmark) real IGS station, single epoch at 2020-06-25
00:00 GPST, GPS L1 C1C only (12 satellites).

- Observations: `test/fixtures/obs/ESBC00DNK_R_20201770000_01D_30S_MO_trim.crx`
- Broadcast nav: `test/fixtures/nav/ESBC00DNK_R_20201770000_01D_MN.rnx`
- Truth: the RINEX header APPROX POSITION XYZ via
  `Orbis.GNSS.RINEX.Observations.approx_position`, ECEF
  `(3582105.291, 532589.731, 5232754.805)`, radius 6.364e6 m. This is the same
  truth the committed SPP test (`rinex_obs_spp_test.exs`) trusts.

This is a static, single-frequency ground receiver, the closest vendored
analogue to a low-cost tracker epoch. Single-frequency broadcast SPP here is
metre class (about 2 to 4 m), which is the right accuracy floor for this
capability.

The GSDC Pixel-5 demo5/RTKLIB JSON oracles are RTK carrier-phase rover
references; the raw phone pseudoranges needed to drive `Positioning.solve`
directly are not vendored, so they are not used as the cold-start driver here. A
convergence-basin test needs only the per-epoch truth ECEF plus driveable
pseudoranges, and the vendored ESBC00DNK station fixture is fully self contained.

## Metric

3D ECEF position error vs the APPROX POSITION truth, in metres. Per prior we also
record `metadata.converged`, `metadata.iterations`, `metadata.status`, post-fit
residual RMS, and `n_used`.

## Sample size

n = 1 epoch for the headline characterization. This is a characterization, not a
gated multi-epoch claim; a single epoch is acknowledged as underpowered for any
basin-width number. The ESBC00DNK 120-epoch obs fixture exists to extend n later.

## Characterization (before)

For each prior in {truth+45km, surface 100/500/1000/3000 km tangential, regional
`{4.5e6, 0.5e6, 4.5e6}`, antipodal `-truth`, earth-center `{0,0,0,0}` = the
module default}, record the metric tuple above. Expected from the design pass:
near-surface seeds converge to a few metres; antipodal returns
`{:too_few_satellites, 0, 4}`; earth-center returns `converged=true` at about
6.36e6 m error with `iterations=0` (a silent false positive).

## Robustification (additive, default off)

New option `:coarse_search` on `Positioning.solve/4`, default `nil` = exact
current single-solve behavior, zero behavior change. When set (`true`, an integer
seed count, or `[seeds: N]`):

1. Generate N deterministic near-surface seeds on a golden-spiral lattice at mean
   Earth radius (6.371e6 m), clock 0. No hardcoded answer; the lattice is fixed
   and independent of truth.
2. Run the existing single solve once per seed (and once more from the caller's
   `:initial_guess` if one was given).
3. Score candidates with a redundancy-gated post-fit scorer: keep only
   candidates with `metadata.converged` AND `n_used >= 5` (redundancy required),
   then pick minimum residual RMS, tie-broken by GDOP.

The redundancy gate is load bearing: an exactly determined 4-satellite fit has
near-zero post-fit RMS regardless of correctness (zero degrees of freedom), so a
naive min-RMS scorer can pick a degenerate 4-sat candidate. Gating on n_used
fixes this. `residuals_m`, `used_sats`, `dop`, and `metadata.converged` are all
already returned by `Solution`, so the scorer is pure Elixir with no kernel
change.

The Rust SPP iteration core (NIF `spp_solve` / `spp_solve_broadcast`) is out of
scope per house rules and is not modified. The elevation-mask/weight freeze that
starves a far-off single seed lives in the kernel; the Elixir coarse search works
around it by supplying near-surface seeds, it does not change kernel behavior.

## Pass bar (after, with :coarse_search on)

From a no-prior cold start (golden-spiral surface seeds, no hardcoded answer),
the redundancy-gated scorer must return a fix with:

- err_3d <= 5.0 m (declared tolerance; the single-frequency SPP floor measured at
  about 2 to 4 m, so 5 m gives margin), AND
- `metadata.converged = true`, AND
- `n_used >= 5`.

The seed count N is a gate parameter pinned by the measurement (the design pass
saw 6 seeds already meet 5 m on this epoch; the committed measurement pins the
chosen default and reports the count sweep).

## Invariant (never weakened)

With `:coarse_search` unset, `solve/4` is byte-for-byte identical to the current
single solve. Proven by keeping every existing assertion in
`point_positioning_test.exs` and `rinex_obs_spp_test.exs` green, unchanged. Those
assertions are never loosened. Out-of-tolerance or underpowered is a fail, not a
pass.
