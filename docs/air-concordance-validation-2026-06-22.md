# Aircraft-air concordance validation (vs Jeilu et al. 2025)

Validate PathogenIQ on the **6 environmental aircraft-filter** samples of
PRJNA1228129 (`databases/aircraft/air-aircraft-filter/`) by **concordance** with
the taxa the paper documents — not a scored precision/recall benchmark.

## Why this is *concordance*, not a scored benchmark

The paper's calls are themselves predictions (MetaPhlAn4 / GTDB-Tk), not gold
truth — two tools agreeing can both be wrong (esp. on contaminants like
*C. acnes*). And false positives can't be scored: a taxon we call that the paper
didn't may be a real low-abundance hit *or* an FP. So this measures **agreement
with the published analysis**, which is still a meaningful real-data check (it is
how the paper itself cross-validated MetaPhlAn vs GTDB-Tk). The spike-in runs
(`air-aircraft-spike`) remain the only *scored* (known-truth) detection test.

## Paper's documented taxa — aircraft-filter ENVIRONMENTAL samples

From Fig 3D / Fig 4 / Table 1 (direct-extraction aircraft filters):

| Taxon | Paper | In PathogenIQ tier-1 clinical DB? | Recoverable by |
|-------|-------|-----------------------------------|----------------|
| *Sphingomonas hankookensis* | genus ~34.7% (dominant) | **No** (environmental) | assembly + GTDB-Tk only |
| *Escherichia coli* | genus ~22.6% | **Yes** | read-based |
| *Pseudomonas aeruginosa* | present | **Yes** | read-based |
| *Cutibacterium acnes* | common across samples | maybe (contaminant prior) | read-based (likely flagged contaminant) |
| *Staphylococcus epidermidis / hominis* | present | partial (*S. aureus* yes) | read-based (overlap only) |
| *Methylobacterium radiotolerans* | present | **No** (environmental) | assembly + GTDB-Tk only |
| *Geobacillus*, *Paenibacillus humicus*, *Microsaccharimonas* | enrichment-dominant | **No** | assembly only |

**Key point:** PathogenIQ's tier-1 DB is ~110 *clinical pathogens*, so the
environmental dominants (*Sphingomonas*, *Methylobacterium*) are **not targets** —
read-based profiling will "miss" them by design. They are exactly what the
**assembly/MAG arm (`--assemble`)** exists to recover, and the reason that arm
matters for air. So the two checks are complementary:

- **Read-based concordance** = recall of the *in-DB* overlap (E. coli, P.
  aeruginosa, staph) under the air NTC background.
- **Assembly concordance** = recover the *environmental* dominants de novo and
  GTDB-Tk-classify them, vs the paper's 9 MAGs.

## Method

1. Build an **air-specific NTC background** from the 6 `air-aircraft-ntc` controls
   (`_classify_taxon_counts` → pool) so the air kitome is subtracted (clinical
   reagent background ≠ air background).
2. **Read-based:** `pathogeniq run --background air_ntc.tsv` on each filter; record
   detected taxa + grades; check recall of the in-DB paper taxa.
3. **Assembly:** `pathogeniq run --assemble …` → MAGs; GTDB-Tk lineage vs the
   paper's environmental species.

## Tooling status (2026-06-22)

- Read-based path: **available** (fastp, bwa-mem2, minimap2, samtools, sourmash).
- `abricate` **missing** → AMR/virulence skip (irrelevant to taxon concordance).
- Assembly arm: **megahit, metabat2, checkm, gtdbtk all missing.** The assembly
  concordance is **blocked** until they (and the CheckM + ~100 GB GTDB-Tk DBs) are
  installed. megahit + metabat2 alone would recover bins but without taxonomy the
  species concordance can't be done.

## Results

_(read-based results appended after the run; assembly results pending tool install)_
