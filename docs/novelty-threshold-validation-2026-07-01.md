# R1 novelty threshold — validation on real air data (2026-07-01)

Calibrates the Kraken2 unclassified-fraction trigger (`novelty.py`,
`_DEFAULT_FLAG_THRESHOLD = 0.5`) on the PRJNA1228129 aircraft dataset.

## Method

One uniform Kraken2 (Standard 8 GB DB) pass over a fixed 500 000-read subsample of
each sample. The deployed screen runs on **non-host** reads and flags when
`unclassified_fraction >= threshold`. Rather than re-run host removal, the
non-host metric is recovered exactly as `U / (total − Homo sapiens)` from the
report (air samples carry ≤10 human reads / 500 k, so raw ≈ non-host anyway).

Groups: **spike** = Zymo D6300 (every organism in-DB → low-novelty positive
control, must NOT flag); **filter** = 5 environmental aircraft filters; **ntc** =
6 aircraft negative controls. Driver: `/data/alvin/tmp/novelty_validation/`.

## Result — bare unclassified fraction (the deployed metric)

| group  | n | min | median | max |
|--------|---|-----|--------|-----|
| spike  | 2 | 0.295 | 0.298 | 0.298 |
| filter | 5 | 0.319 | 0.421 | 0.478 |
| ntc    | 6 | 0.219 | 0.324 | 0.472 |

Two hard findings:

1. **0.5 fires on nothing.** The most-unclassified sample in the whole set is a
   *filter* at 0.478 — below 0.5. As calibrated, the trigger is inert on this air
   data.
2. **In-DB reads have a ~0.30 unclassified floor.** The Zymo spike — 100 % catalogued
   organisms — still sits at 0.30 (shotgun reads vs an 8 GB DB: intergenic / strain-
   variant / low-complexity sequence never classifies). So the usable dynamic range
   above the floor is small, and no air sample crosses 0.5.

## The metric does not separate "assemble-worthy" from "kitome"

The operational question is not spike-vs-filter, it is *does this sample hold novel
content worth the assembly arm?* — i.e. filter (real air) vs ntc (nothing to
recover). It does not:

- NTC `SRR32514323` scores **0.472 unclassified — higher than 4 of 5 real filters.**
  Low-biomass negative controls alias as "novel" because their sparse, contaminated
  read pool is just as hard to classify as environmental air.

We tested the rank-resolution refinement the code comment proposed (novelty =
reads that hit a genus but resolve to no species) and two combined variants. **All
overlap:**

| candidate metric | filter_min | ntc_max | verdict |
|------------------|-----------|---------|---------|
| unclassified fraction        | 0.319 | 0.472 | OVERLAP |
| genus-only (stuck at genus)  | 0.069 | 0.193 | OVERLAP |
| unclassified + genus-only    | 0.388 | 0.665 | OVERLAP |
| species-resolved fraction    | 0.079 | 0.576 | OVERLAP |

`genus-only` is worse than useless: the *spike* (all known) has the **highest**
value (0.30) — Zymo strains hit near-relatives and stick at genus — so it is a
confounded "known-organism" signal, not a novelty signal.

## Verdict & decision

- **Keep `_DEFAULT_FLAG_THRESHOLD = 0.5`.** It sits just above the observed NTC
  ceiling (0.472), so it is correctly calibrated to *never false-trigger on the air
  kitome*. The cost — it also never fires on the catalogued environmental filters —
  is acceptable: those organisms (*Sphingomonas*, *Methylobacterium*, …) **are** in
  the DB (consistent with the earlier 6/10 read-based concordance), so there is
  genuinely little read-level novelty to discover in this dataset.
- **Do NOT lower it.** A threshold ~0.31 (to catch filters) would flag 4/6 NTCs —
  strictly worse than inert.
- **The bare read-fraction gate cannot be both sensitive and specific on low-biomass
  air**, and neither can its rank-resolution refinement. Reliable novelty here is a
  *post-assembly, per-MAG* judgement (completeness + GTDB placement — R3–R5, already
  built), not a pre-assembly read-fraction pre-filter. Since `--assemble` already
  runs regardless of the trigger, the novelty gate is advisory, not gating — which
  is the correct role for it at this reliability.
- **Next refinement if the trigger must gate:** precede it with a biomass floor
  (NTCs are low-biomass by definition) so kitome can't alias as novelty.

Raw table: `/data/alvin/tmp/novelty_validation/novelty_validation.tsv`.
