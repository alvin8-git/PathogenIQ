# T10 Benchmark — first run (2026-06-21)

Grading vs raw Kraken2 on the ZymoBIOMICS D6300 standard
(`V350234554_L01_UDB-32.Zymo_Std.fq.gz`, ~22M reads, single community).
Configs scored with `scripts/06_benchmark.py` against the known Zymo composition.

## Standard-8 DB (full RefSeq, 893 species called) — the realistic case

| config | PR-AUC | precision | recall | TP | FP | FN |
|---|---|---|---|---|---|---|
| kraken-raw | 0.658 | **0.009** | 0.80 | 8 | 885 | 2 |
| kraken+grading | 0.658 | **0.051** | 0.80 | 8 | 148 | 2 |
| kraken+grading+bg | 0.480 | 0.041 | **0.60** | 6 | 142 | 4 |

## tier1 DB (custom, 74 species called)

| config | PR-AUC | precision | recall | TP | FP | FN |
|---|---|---|---|---|---|---|
| kraken-raw | 1.000 | 0.108 | 1.00 | 8 | 66 | 0 |
| kraken+grading | 1.000 | 0.111 | 1.00 | 8 | 64 | 0 |
| kraken+grading+bg | 0.750 | 0.171 | **0.75** | 6 | 29 | 2 |

## Findings

1. **The grading read-floor works, and the effect is large on a realistic DB.**
   On Standard-8, grading lifts precision **0.009 → 0.051 (5.7x)** while keeping
   recall at 0.80 — it removes ~737 of 885 false positives (the low-read spurious
   tail) without dropping a true organism. This is the wedge demonstrating value.

2. **The current Salter kitome background is net-harmful here.** Adding it drops
   recall (0.80 → 0.60 on Std-8; 1.00 → 0.75 on tier1) by suppressing **real
   E. coli / Klebsiella** — organisms that are both genuine Zymo members and
   present in the kitome background — for only a marginal FP reduction. This
   empirically confirms two documented issues: the Enterobacteriaceae confound in
   the Salter-derived background, and the crude single-NTC dispersion prior
   (Plan-4), which makes the background distribution so wide that even abundant
   real organisms fall "within background."

3. **Metric caveats (single community).** PR-AUC is insensitive when FPs are
   bottom-ranked by read count (grading's value shows in operating-point
   precision, not AUC). precision@95%-recall is undefined here because recall caps
   at 0.80 (2 Zymo organisms undetected). One community cannot anchor a PR curve —
   a real validation needs >=3-5 labeled communities (outside-voice c).

## Conclusions / next

- **Ship the grading layer** (read floor / CI / tier) — it measurably improves
  precision at preserved recall.
- **Do NOT auto-apply the current kitome background** to samples likely to contain
  Enterobacteriaceae — it costs real detections. The shipped Tier-2 default needs:
  (a) Plan-3D dedup (E. coli/Shigella), (b) non-spiked blank datasets, and
  (c) a better dispersion model (Plan-4) so high-abundance real organisms are not
  suppressed.
- Add CAMI communities + a held-out split for a defensible PR-AUC.
