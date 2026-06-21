# Breadth-of-coverage exploration (2026-06-21)

Explored the air-surveillance design's "coverage breadth over read depth" idea.
Built `pathogeniq/coverage.py` (breadth, depth, Lander-Waterman expected breadth,
depth-normalized `breadth_ratio = observed/expected`) from the minimap2 PAF the
alignment stage already produces.

## What breadth catches — and what it doesn't

**Synthetic (validated):** same read count, opposite verdict.
- 50 reads piled at one 150 bp locus → `breadth_ratio ≈ 0.02` (clumped artifact)
- 50 reads spread across the genome → `breadth_ratio ≈ 1.0` (real)

**Real kitome PAFs:** every taxon `breadth_ratio ≈ 1.0` even at 1–161 reads — those
reads are genuine genomic DNA *spread* across the contaminant genome. So breadth
does NOT flag low-abundance spread reads; it complements the NTC background.

**Cross-mapping test (the case I predicted breadth would solve — IT DOES NOT):**
Zymo reads aligned to E. coli vs Shigella refs:

| ref | reads | depth | breadth | breadth_ratio |
|---|---|---|---|---|
| E. coli (real) | 1.96M | 21x | 0.916 | 0.916 |
| Shigella (shadow) | 1.38M | 14x | 0.743 | 0.743 |

The Shigella shadow is **74% covered, not clumped** — because the genomes are >95%
identical, E. coli reads map across almost the whole Shigella genome. Breadth is
genome-wide for both, so `breadth_ratio` cannot cleanly separate them. (At 14–21x
depth the expected breadth saturates to 1.0, so the ratio ≈ breadth and the shadow
still scores high.)

## Conclusion — breadth is a real but narrow tool

- **Use it for:** clumped artifacts — PCR/optical duplicates, low-complexity
  pileups, and *partial/distant* cross-mapping where only a conserved region
  matches (the air design's "100 reads at one B. anthracis locus" example). Best
  as a secondary verification gate on flagged high-threat taxa.
- **Do NOT use it for:** near-identical-genome cross-mapping (E. coli/Shigella).
  Those shadows are broadly covered; the real discriminators are competitive read
  assignment (EM / Kraken LCA — already in the pipeline) or per-read alignment
  identity, not breadth. Note our existing read-count dedup works on the *Kraken*
  path precisely because Kraken's LCA assigns shared reads to one taxon; on a raw
  *alignment* path neither read count nor breadth separates them (E. coli only 1.4x
  Shigella's reads here).

## If/when integrated
- A `breadth_ratio` field on GradingInput, low ratio → demote/drop, as a
  verification gate after alignment (align.py already has the PAF; it currently
  discards the positions coverage.py needs).
- A relative breadth/depth ranking *within* a cross-mapping group (E. coli 0.92 >
  Shigella 0.74) is a marginal refinement over the read-count dedup, not a
  replacement.
- Net: worth adding as a clumped-artifact verifier; it is NOT the E. coli/Shigella
  fix (that stays with EM/LCA + the crossmap dedup).
