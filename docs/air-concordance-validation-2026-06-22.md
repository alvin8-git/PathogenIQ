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

## Paper supplements — checked 2026-06-22

The article ships **exactly one** supplement: *Additional file 1 — Supplementary
Table 1 (metadata)* = `MOESM1_ESM.xlsx`, a sample metadata + qPCR table (sample
name, type, processing, DNA conc, log-normalized 16S copies). **There is no
downloadable per-sample species table** (no MOESM2+, no data-repo link). So
species concordance is limited to figure-level data. The xlsx is still useful:
confirms AF1–6 = the aircraft filters, and its copy numbers can sanity-check the
spike-in absolute quantification.

**Most granular paper species set = Fig 3G top-20 heatmap** (across all samples):
Staphylococcus epidermidis, Cutibacterium acnes, GGB2722_SGB3663, S. hominis,
Sphingomonas hankookensis, Paenibacillus humicus, **Pseudomonas aeruginosa**,
**Escherichia coli**, Methylobacterium radiotolerans, S. capitis, Streptococcus
mitis, S. australis, **Bacillus cereus**, **Bacillus subtilis**, S. warneri,
Priestia megaterium, **Clostridium perfringens**, Limosilactobacillus fermentum,
Brevibacterium parabrevis. (**bold** = clinically-relevant + plausibly in our
tier-1 DB → the read-based concordance targets; the rest are commensal/
environmental, recoverable only via the triggered assembly arm.)

## What the "spike" is, and how the paper handled controls (from the PDF)

- **Spike = wet-lab ZymoBIOMICS D6300, not in-silico.** Per Methods ("Spike-in
  controls"): a 2 cm × 2 cm sterile face-mask section was inoculated with 100 µL
  of the ZymoBIOMICS Microbial Community Standard (Zymo D6300) and run through the
  same extraction + sequencing. So our Zymo-truth scoring of SpikeC1/SpikeC2 is
  correct; the two differ in input DNA (1.2 vs 10 ng/µL) — a sensitivity check.
  *(Distinct from the "synthetic spike-in DNA standard, 10⁹ copies/µL" added to
  every sample — that is a qPCR normalization/recovery spike, not the mock.)*
- **The paper did NOT subtract NTC/control species.** Its only "contaminant
  removal" is read-level: Trimmomatic + KneadData (human + **PhiX**) + over-
  represented-sequence removal. The **NTC is a qPCR-only** monitoring control.
  The control/unworn masks are sequenced and reported as their *own category*
  (Fig 3), never used to subtract taxa from samples. The paper instead argues
  biologically (e.g. *Sphingomonas hankookensis* is "valid biological signal
  rather than contamination ... a bona fide component").
- **Our "air NTC" is built from those control masks — which are environmental
  exposure controls, not no-template blanks.** The 6 pooled controls (C1–3
  unworn-mask environmental + CE enrichment) legitimately carry airborne + skin/
  handling microbiota, so using them as a kitome background subtracts *real air
  signal*, not just reagent noise. That is the root reason the air NTC erased
  Zymo *E. coli*/*P. aeruginosa*. The flag-not-subtract fix (2026-06-23) mitigates
  it for dual-use pathogens; the cleaner fix is a true reagent/no-template NTC.

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

**Update (2026-07-02):** the assembly arm is no longer blocked. `scripts/15_setup_mag_env.sh`
installs metabat2 + CheckM and GTDB-Tk (release r232) in isolated conda envs with
PATH wrappers into the core env; the CheckM and GTDB-Tk DBs are built on /data. The
first `--assemble` run on the 5 aircraft filters (scripts/16) is underway — MAG
taxonomy/concordance results are not yet finalized.

## Results — read-based, `--specimen air` (scored 2026-06-23)

Run: `air_all.sh` (headless, setsid). 5 filters + 2 Zymo spike-ins, air NTC
background, no PDF. AMR/VFDB columns **empty in every report** — `abricate` is
not installed (fields are `[]`, no crash; the screen is non-blocking by design),
so AMR/virulence concordance is out of scope for this run.

### Spike-ins — scored validation (Zymo D6300, known truth)

Truth = 8 bacteria @ ~12% + 2 yeast @ ~2% (10 organisms). Both replicates
(SRR32514283, SRR32514284) agree. Detected **6/10**, all Grade B (correctly
capped at B by the Tier-2 pooled background):

| Zymo member | In DB? | Detected | Note |
|-------------|--------|----------|------|
| Salmonella enterica | Y | ✅ ~22% | |
| Lactobacillus fermentum | Y | ✅ ~25% | |
| Enterococcus faecalis | Y | ✅ ~7% | |
| Staphylococcus aureus | Y | ✅ ~5% | pathogen, **contam=False** ✓ |
| Listeria monocytogenes | Y | ✅ ~2% | |
| Cryptococcus neoformans | Y | ✅ ~0.03% | yeast, recovered at <0.1% |
| Escherichia coli | Y | ❌ | **subtracted by air NTC** (E.coli 346 788 rpm in background) |
| Pseudomonas aeruginosa | Y | ❌ | **subtracted by air NTC** (140 308 rpm in background) |
| Saccharomyces cerevisiae | Y | ❌ | **subtracted by air NTC** (35 717 rpm in background) |
| Bacillus subtilis | **N** | ❌ | not in DB; surrogate *B. anthracis* also in background |

**The 4 misses are not detection failures — they are background over-subtraction.**
Raw read-based detection recovers ≥9/10 (all but the *B. subtilis* DB gap); the
air NTC background then zeroes out 3 genuine pathogens (*E. coli*, *P. aeruginosa*,
*S. cerevisiae*) because **those exact taxa dominate the air "control" runs** used
to build `air_ntc.tsv`.

### The air NTC background is contaminated (root cause)

`air_ntc.tsv` (6 aircraft-control runs) subtracts, by descending rpm:
*Shigella sonnei* 419 713 · *E. coli* 346 788 · *S. epidermidis* 309 000 ·
*P. aeruginosa* 140 308 · *C. acnes* 56 869 · *S. cerevisiae* 35 717 ·
*B. anthracis* 34 857 · *S. flexneri* 19 090.

This is the **same failure mode we rejected for the HUNT blanks** (E. coli–inflated):
the controls are not clean blanks — they carry kitome + carryover at pathogen-level
abundance, so the background subtracts real signal. *E. coli / Shigella / Pseudomonas*
are dual-use (kitome **and** clinical target); subtracting them at species level is
exactly the mistake the AIR contaminant prior deliberately avoids (flag, don't demote).

**Fix options (not yet done):** (a) vet the 6 control runs — confirm they are true
sampling blanks, not environmental samples; (b) for dual-use taxa, **flag-not-subtract**
in the air background (mirror the contaminant-prior policy) so a real *E. coli*/
*P. aeruginosa* hit is demoted to a capped grade rather than erased.

### Contaminant-prior check — PASS

Across all 7 reports, **no real pathogen was demoted** by the air prior:
*S. aureus, Salmonella, Listeria, Enterococcus, S. maltophilia, S. flexneri,
C. perfringens* all `contaminant_risk=False`. Only *C. acnes* (filter SRR32514297)
was flagged — correct. Confirms `test_air_does_not_demote_real_pathogens`.

### Filters — concordance vs paper (Fig 3G **bold** in-DB targets)

| Paper bold target | Recovered? | Where |
|-------------------|-----------|-------|
| *E. coli* (genus ~22.6%) | ✅ as ***Shigella flexneri*** | 4/5 filters (E.coli↔Shigella cross-map; genomically near-identical) |
| *Clostridium perfringens* | ✅ | SRR32514308 (Grade B, 4.5%) |
| *C. acnes* (common) | ✅ flagged contaminant | SRR32514297 |
| *Pseudomonas aeruginosa* | ❌ | **subtracted by air NTC** (in background at 140 308 rpm) |
| *Bacillus cereus / subtilis* | ❌ | DB gap + *B. anthracis* surrogate in background |

Also recovered, not in paper bold set: *Stenotrophomonas maltophilia* (3/5 filters,
Grade B) — plausible real environmental Gram-negative, unscorable (agreement-not-truth).

**Net:** read-based detection itself is healthy (6 Zymo + cross-mapped E. coli + C.
perfringens, contaminant prior behaving). The **single dominant problem is the air
NTC background**: built from contaminated controls, it erases *E. coli*, *P.
aeruginosa*, and *S. cerevisiae* from both the scored spike and the filter
concordance. Rebuild/flag-not-subtract before trusting air negatives. Assembly arm
(environmental dominants *Sphingomonas*/*Methylobacterium*) still blocked on tool
install.
