# spp-robustness Gate 2 measurement (2026-06-13T15:46:43.289233Z)

Before/after single-frequency SPP on real GSDC Pixel-5 phone observations,
bare crate solve vs opt-in robust (FDE) solve, vs the demo5/RTKLIB
position-domain oracle. Truth is GSDC ground truth carried in the oracle.

Gate params: min epochs 100, credibility factor 2.0x demo5 median (absolute-accuracy floor only), max_pdop 1.0e3.

| arc | n | bare med 3D | bare p95 3D | robust-unit med 3D | robust-unit p95 3D | robust-wtd med 3D | robust-wtd p95 3D | demo5 med 3D | wtd sats excl |
|---|---|---|---|---|---|---|---|---|---|
| gsdc_2021_08_04_sjc1_pixel5_p222_grec_l1_demo5 | 1453 | 10.213 | 37.535 | 13.054 | 64.905 | 9.680 | 42.778 | 4.522 | 1204 |
| gsdc_svl1_pixel5_p222_grec_l1_demo5 | 3136 | 9.241 | 26.877 | 10.596 | 39.410 | 8.837 | 24.673 | 3.977 | 1171 |
| gsdc_2021_12_15_mtv1_pixel5_p222_grec_l1_demo5 | 1465 | 8.062 | 29.437 | 9.327 | 37.561 | 7.729 | 28.366 | 3.653 | 367 |
| gsdc_2021_12_28_mtv1_pixel5_p222_grec_l1_demo5 | 1610 | 10.599 | 27.830 | 12.185 | 48.816 | 10.145 | 26.768 | 3.974 | 553 |

All values in metres. robust-unit is RAIM/FDE with unit weights (sigma=1 m
assumed); robust-wtd is RAIM/FDE with a realistic uniform phone code sigma of
5.000 m.

Pooled: powered arcs 4/4; substrate-floor-passing 0/4; all powered arcs improved (robust-wtd <= bare on median and p95)? false.

Reading: orbis runs an unaided single-frequency SPP per epoch; demo5 is a
tuned multi-GNSS RTK reference and is the absolute bar (the 2x floor), not the
comparand for the robustness delta. The robustness claim is the bare-vs-robust
delta on identical orbis SPP inputs. The unit-weight FDE over-excludes on real
phone noise (it treats several-metre code noise as faults under a 1 m sigma
assumption) and degrades the fix: reported as a null/negative result, not
hidden. A non-positive delta is a null result, not massaged.

## Findings

Gate 1 (FDE exclusion correctness, in-repo SP3 oracle, `test/gnss_qc_test.exs`):
PASS. A labelled bias sweep of 50, 100, 200, 500 m on one satellite in an
otherwise clean GPS-only set is each excluded as exactly the biased satellite
(`excluded == [{biased, :raim_excluded}]`), the recovered 3D error returns to
within `clean_error + 0.5 m`, and the clean set produces zero false exclusions
with a position bit-identical to the bare solve.

Gate 3 (degenerate-geometry guard): PASS. `:max_pdop` below the realized PDOP
returns `{:error, {:degenerate_geometry, pdop}}`; a generous threshold returns
the same solution as the unguarded solve; the default (`nil`) never refuses; the
guard composes with `:robust`.

Gate 2 (before/after on real GSDC phone data): the decisive finding is that
robustness depends entirely on the RAIM noise model.

  * Unit-weight FDE (the naive default, sigma = 1 m) is WORSE than the bare solve
    on all four arcs, on both median and p95, excluding roughly 4 to 6
    satellites per epoch. Under a 1 m sigma assumption the chi-square test reads
    ordinary several-metre phone code noise as a fault and strips good
    satellites, degrading geometry. This is the headline negative result.
  * Weighted FDE with a realistic uniform phone code sigma of 5 m improves the
    median on all four arcs and improves p95 on three of four (svl1, both mtv1),
    excluding only the genuine outliers (about 0.2 to 0.8 satellites per epoch).
    On sjc1 the median improves (10.21 to 9.68 m) but p95 regresses
    (37.5 to 42.8 m), so the strict "median AND p95 down on every arc" gate is
    NOT met (`all powered improved? false`); 3 of 4 arcs meet it.

The 2x credibility floor against demo5 fails on all arcs and is expected to:
orbis here is unaided single-frequency SPP, demo5 is tuned multi-GNSS RTK. The
floor gates an absolute-accuracy claim, which is not made; the in-scope claim is
the bare-vs-robust delta on identical orbis inputs.

Verdict: the robustness machinery (FDE in the solve path, degenerate guard) is
correctly wired and demonstrably isolates injected faults (Gate 1) and rejects
weak geometry (Gate 3). On real degraded phone data the gain is conditional: it
requires passing a realistic detection noise model through `:weights`. With the
default unit weights FDE hurts; this is the reportable result and is the reason
`:robust` is opt-in and pairs with `:weights`. The solve-domain Huber/IRLS
reweighting that would handle the noise floor without exclusion lives in the
crate and is out of scope for this Elixir-only task.
