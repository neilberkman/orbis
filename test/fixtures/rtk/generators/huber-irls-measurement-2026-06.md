# huber-irls GSDC truth metric (2026-06-14T05:24:39.872725Z)

Before/after single-frequency SPP on real GSDC Pixel-5 phone observations,
bare static-elevation-weighted crate solve vs the opt-in crate-layer
Huber/IRLS solve (`:huber`), on IDENTICAL inputs, vs the demo5/RTKLIB
position-domain oracle. Truth is GSDC ground truth carried in the oracle. See
huber-irls-spec.md (pre-registered).

Gate params: min epochs 100, credibility factor 2.0x demo5 median (absolute floor only), Huber k 1.345, MAD scale floor 5.000 m, max outer 5. Epoch stride 3.

| arc | n | bare med 3D | bare p95 3D | huber med 3D | huber p95 3D | delta med 3D | delta p95 3D | demo5 med 3D | classification |
|---|---|---|---|---|---|---|---|---|---|
| gsdc_2021_08_04_sjc1_pixel5_p222_grec_l1_demo5 | 485 | 9.815 | 34.532 | 9.719 | 34.591 | 0.096 | -0.059 | 4.522 | median-only: median non-regress, p95 regressed (null on strict bar) |
| gsdc_svl1_pixel5_p222_grec_l1_demo5 | 1046 | 9.267 | 27.134 | 8.742 | 23.595 | 0.524 | 3.539 | 3.977 | improved: Huber non-regression on median and p95 |
| gsdc_2021_12_15_mtv1_pixel5_p222_grec_l1_demo5 | 489 | 8.046 | 28.170 | 7.347 | 25.668 | 0.698 | 2.502 | 3.653 | improved: Huber non-regression on median and p95 |
| gsdc_2021_12_28_mtv1_pixel5_p222_grec_l1_demo5 | 537 | 11.063 | 27.426 | 10.646 | 24.430 | 0.418 | 2.996 | 3.974 | improved: Huber non-regression on median and p95 |

All values in metres. Delta is bare minus Huber (positive = Huber better).

Pooled: powered arcs 4/4; Huber-off byte-identical to bare on every powered arc? true; all powered arcs median non-regress (Huber <= bare on median)? true; all powered arcs median AND p95 non-regress? false.

Strict bar (pre-registered): Huber 3D median <= bare AND p95 <= bare on EVERY
powered arc, no slack. A non-positive delta is a null result, not massaged.
The default path (`:huber` off) producing byte-identical numbers re-proves
additive-off on real data.

Reading: orbis runs an unaided single-frequency SPP per epoch; demo5 is a
tuned multi-GNSS RTK reference and is the absolute context bar, not the
comparand for the reweighting delta. The capability claim is the bare-vs-Huber
delta on identical orbis SPP inputs.
