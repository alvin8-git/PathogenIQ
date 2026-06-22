# Dispersion-prior validation (Plan-4 follow-up) — 2026-06-22

**Question (todo.md):** with one NTC the per-taxon NB dispersion `r` isn't
estimable, so `background.py` uses a shared prior (`_DEFAULT_DISPERSION = 2.0`).
How sensitive is the background test's false-positive rate to that prior, and
what default is defensible? If FPR is highly sensitive, the single-NTC (Tier-2)
"controlled-α" claim is hollow.

**Method (`scripts/10_validate_dispersion.py`):** leave-one-out across the 15
spike-free Salter kitome blanks (cached GCF-keyed counts from script 05). Build
the background from 14 blanks, test every taxon in the held-out blank. The
held-out run is *itself* a blank, so any taxon called "above background"
(p < α=0.01) is a **false positive** — there is no real signal in a blank (the
spike taxon is excluded). Sweep `r`; also report the limit of detection (min
reads for an NTC-absent taxon to clear the test at 1M depth) as the sensitivity
cost.

## Result

| dispersion r | LOO FPR | FPR ≤ α | LoD (reads) |
|---|---|---|---|
| 0.1 | 0.302 | no | 10 |
| 0.5 | 0.326 | no | 6 |
| 2.0 (current) | 0.488 | no | 4 |
| 10 | 0.581 | no | 4 |
| 1000 (≈Poisson) | 0.651 | no | 4 |

(`r=0.001` → FPR 1.0 is NB numerical degeneracy, LoD=1 read — not a usable regime.)

## Findings

1. **No dispersion controls FPR.** Across six orders of magnitude, LOO FPR moves
   only 0.30 → 0.65 and never approaches α=0.01 — it sits **30–65×** the nominal
   rate. The single-NTC NB test's controlled-α guarantee is **hollow** on real
   kitome.

2. **Dispersion is not the lever — NTC coverage is.** FPR has an irreducible
   floor ≈0.30 because **10 of 18 kitome taxa occur in only one blank**. A taxon
   absent from the training pool defaults to the pseudocount floor (0.5 RPM) and
   reads as real at *any* `r` — the dispersion prior governs variance around a
   *known* background rate, not unseen taxa. Sweeping `r` cannot fix a coverage
   gap.

3. **The `r` tradeoff is real but secondary.** Smaller `r` (wider null) lowers
   FPR toward the floor and raises LoD (4 → 10 reads); it also *widens* the null
   for taxa that ARE in the background, increasing suppression of genuine
   organisms — the over-suppression the T10 benchmark already flagged (real
   Enterobacteriaceae dropped). So lowering `r` cannot rescue α-control and would
   worsen recall.

## Decision

- **Keep `_DEFAULT_DISPERSION = 2.0`.** Lowering it can't deliver α-control (the
  limit is coverage) and would deepen the documented real-organism
  over-suppression. The validation's value is the knowledge, not a knob change.
- **Keep Tier-2 capped at Grade B.** This result is the empirical justification:
  a single pooled NTC cannot support a controlled-α "Grade A" detection.
- **The real fix is coverage:** a batch-matched NTC (Tier 1), or many more
  pooled blanks so run-unique contaminants stop hitting the pseudocount floor.
  Per-taxon dispersion only becomes estimable — and worth modeling — with ≥2–3
  batch-matched NTCs.

Reproduce: `python scripts/10_validate_dispersion.py`
