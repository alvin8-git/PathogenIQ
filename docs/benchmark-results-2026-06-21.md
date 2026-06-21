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

## Update — Plan-3D cross-mapping dedup (same day)

Added `crossmap.py` and wired it into grading + the benchmark. On the Standard-8
Zymo report it flagged all 4 Shigella species (flexneri, sonnei, dysenteriae,
boydii) as E. coli cross-mappers and dropped them (Grade X):

| config (Std-8) | precision | recall | TP | FP |
|---|---|---|---|---|
| kraken+grading (pre-dedup) | 0.051 | 0.80 | 8 | 148 |
| kraken+grading (with dedup) | 0.053 | 0.80 | 8 | 144 |

Key contrast with the NTC background: the dedup is **precision-positive with zero
recall cost** (it removes only phantom relatives of a dominant real organism),
whereas the kitome background dropped *real* E. coli/Klebsiella. The dedup is the
right tool for this artifact class; the background still needs Plan-4 + non-spiked
blanks before it's a net positive.

## Update — held-out across Zymo D6300 runs (`scripts/08_heldout_pr_auc.py`)

Three sequencing runs of the same standard (UDB-32/-40/-48), classified with the
Standard-8 DB. Held-out split: the read floor is calibrated (on F1) on r1, then
reported on the held-out r2 + r3.

| | raw precision | graded precision | graded recall |
|---|---|---|---|
| fixed grading (all runs out-of-sample) | 0.009 | 0.051 | ~1.0 |
| **calibrated floor=500 (held-out r2+r3)** | 0.009 | **0.875** | **0.700** |

A read floor learned on one run generalizes to the held-out two: precision
0.875, recall 0.700 (F1 0.778) — a ~90x precision lift over raw, at a recall cost
(3 of 10 Zymo organisms fall below the floor). Graded PR-AUC is *lower* (0.578 vs
0.651 raw) because aggressive filtering trades ranking-recall for operating-point
precision; both are reported.

## Update — first CROSS-community held-out (Zymo mock -> CAMI mouse gut)

Added a genuinely different community: CAMI II mouse-gut sample_0 (real CAMISIM
data, per-sample gold standard = 64 species present), classified with Standard-8.

Fixed grading, per community:

| community | raw precision | graded precision | graded recall |
|---|---|---|---|
| zymo_r1/2/3 | ~0.009 | ~0.05 | 0.80 |
| cami_mousegut | 0.006 | 0.015 | 0.53 |

**Cross-community held-out** — read floor calibrated on the 3 Zymo runs, then
applied to the held-out **mouse-gut** community:

| held-out community | raw precision | graded precision | graded recall |
|---|---|---|---|
| cami_mousegut (floor=500 from Zymo) | 0.006 | **0.463** | 0.484 |

A threshold learned on a mock-standard community **generalizes to a completely
different gut microbiome**: precision 0.006 -> 0.463 (77x), recall 0.484
(F1 ~0.47). This is the cross-community generalization the held-out PR-AUC was
for. Caveats: mouse-gut recall is bounded by Standard-8 DB coverage of gut
anaerobes + taxid-version matching between CAMISIM and the DB; still only 2
community *types* (marine/strain/HMP would strengthen it).

**Scope caveat:** the Zymo runs alone are one community (D6300) across runs — they
demonstrate run-to-run generalization; the mouse-gut adds the first true
cross-community test.
True multi-community PR-AUC needs CAMI reads (portal-gated at
data.cami-challenge.org/participate, multi-GB). The harness and the CAMI
gold-standard truth parser (`parse_cami_profile`, validated on the mouse-gut
profile -> 549 species) are built and tested, so CAMI communities slot straight
in once their reads are classified to a .kreport.
