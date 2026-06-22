# PathogenIQ TODO — Next Session

## Critical Bug: Pipeline Produces Broken Output  — RESOLVED (verified 2026-06-21)

**Status:** FIXED in commits 8c77b68 (resolve names from GCF accessions) and
4906f09 (genome path resolution). `sketch.py::_build_md5_map` now maps md5 →
(name, taxon_id, genome) from the `.sig` files + `name_map.json`, and the DB is
built (`databases/tier1/` with sigs/ + name_map.json). The original symptom
below is historical. T1 of Plan 4 then added taxon_id threading on top.

**Original report (historical):**

**Problem:** Running the full pipeline on Zymo data produces reports with:
- Empty organism names
- NaN abundances and CI bounds
- INT64_MIN read counts (`-9223372036854775808`)
- All Grade X

**Root Cause (in `pathogeniq/sketch.py`):**
When `sourmash search` queries an SBT (Sequence Bloom Tree) index, the CSV output has:
- `name` column: **empty** (signatures were built without `--name`)
- `filename` column: the SBT zip file path, not individual genome paths
- `md5` column: valid md5 hash of the matched signature

This causes `SketchHit.name=""` and `SketchHit.genome_path=SBT.zip`, so:
1. `align.py` passes empty names to the report
2. `align.py` runs `minimap2` against a `.zip` file → empty PAFs
3. EM receives zero-alignment matrix → NaN/underflow values

**Fix Required:**
Build an md5-to-genome mapping from the signature files in `databases/tier1/sigs/`.
Each `.sig` file's md5 (from `sourmash sig describe`) maps to:
- Genome name = sig filename stem (e.g. `GCF_000005845.2_ASM584v2_genomic`)
- Genome path = `databases/tier1/genomes/{stem}.fna`

In `sketch.py`, after parsing CSV rows, look up each row's `md5` in the mapping
to get the correct `name` and `genome_path`.

## Already Completed This Session

- [x] README.md updated with Plan 2 features (AMR, PDF, contaminants, pipeline diagram)
- [x] README.md committed and pushed to GitHub
- [x] Version.md and Documentation.md created in previous session
- [x] Zymo pipeline output inspected — bug identified

## Next Steps (Priority Order)

1. **Fix `pathogeniq/sketch.py`** — md5-to-genome mapping for SBT search results
2. **Re-run Zymo pipeline** to verify fix
3. **Run unit tests** (`pytest tests/ -v`) to ensure no regressions
4. **Commit and push** the fix

---

## Plan 3 — Precision Upgrades (from pipeline review)

Three improvements identified from comparing PathogenIQ against an alternative DADA2-based
pipeline design. DADA2 itself was rejected (amplicon-only tool, incompatible with shotgun WGS).
The following are the transferable improvements:

### 3A — Targeted Read Extraction (Stage 3.5)

**What:** After sourmash shortlists candidate organisms, extract only reads that could plausibly
map to those taxa before running minimap2 alignment. Tool: `seqtk subseq` or `bbmap`.

**Why:** Currently PathogenIQ passes the entire non-human FASTQ to minimap2 against all candidate
genomes. Extracting targeted reads first:
- Reduces minimap2 wall time proportionally to the fraction of off-target reads
- Lowers false-positive alignments from reads that happen to seed to wrong genomes
- Enables per-organism parallelism in Stage 4

**Implementation sketch:**
```
sourmash gather → candidate list → seqtk subseq (reads with k-mers in candidates) → minimap2
```
New module: `pathogeniq/extract.py`, function `extract_candidate_reads(cfg, nonhuman_fastq, hits) -> Path`

**Where it fits:** Insert between `run_sketch_screen()` and `run_targeted_alignment()` in `cli.py`.
Design spec: new Stage 3.5.

---

### 3B — AMRFinderPlus (replace abricate)

**What:** Replace `abricate` with NCBI AMRFinderPlus as the primary AMR detection tool.

**Why:**
- abricate wraps CARD/ResFinder but is community-maintained and lags updates
- AMRFinderPlus is NCBI-maintained, updated quarterly, covers CARD + ResFinder + stress/virulence genes
- AMRFinderPlus outputs organism-level context when given a taxon flag (e.g. `--organism Staphylococcus_aureus`), reducing false positives
- Required for regulatory submissions (FDA ARGOS, NCBI SRA pathogen detection pipeline uses it)

**Implementation:** New `pathogeniq/amr_finder.py` module alongside existing `amr.py`.
CLI flag: `--amr-tool [abricate|amrfinderplus]` (default: amrfinderplus once installed).
Fallback gracefully to abricate if amrfinderplus not on PATH.

**Install:** `conda install -c bioconda ncbi-amrfinderplus && amrfinder --update`

---

### 3C — MLST / Strain Typing (new Stage 5 sub-step)

**What:** Add Multi-Locus Sequence Typing (MLST) for key clinical pathogens to provide
strain-level resolution beyond species ID.

**Why:** Species-level calls are insufficient for outbreak response and infection control.
Strain typing answers clinically critical questions:
- Is this *E. coli* ST131 (high-risk, fluoroquinolone-resistant clone)?
- Is this *S. aureus* ST8 (USA300 community MRSA)?
- Is this *K. pneumoniae* ST258 (carbapenem-resistant pandemic clone)?

**Tool:** `mlst` (Torsten Seemann) — fastest, uses PubMLST schemes, supports 100+ species.
Runs on assembled contigs or directly on reads (via `--minid`/`--mincov` thresholds).

**Implementation:**
- New `pathogeniq/mlst.py`, function `run_mlst(cfg, contigs_or_reads, organism_names) -> list[MLSTResult]`
- `MLSTResult` dataclass: `organism`, `scheme`, `sequence_type`, `allele_profile: dict[str, str]`
- Only runs for organisms with a PubMLST scheme (skip fungi, viruses, novel organisms)
- Output: MLST column in TSV report + "Sequence Type: ST131" in PDF report

**Priority organisms to cover:** *E. coli*, *S. aureus*, *K. pneumoniae*, *P. aeruginosa*,
*S. pneumoniae*, *N. meningitidis*, *L. monocytogenes*, *E. faecalis/faecium*

**Install:** `conda install -c bioconda mlst`

---

### 3D — False Positive Deduplication (Shigella/E.coli cross-mapping) — RESOLVED 2026-06-21

**Done:** `pathogeniq/crossmap.py` (`find_crossmappers` / `deduplicate_closely_related`)
flags the minor member of a known cross-mapping group when a relative outnumbers
it >=10x, and grade() drops it as X. Wired into build_entries + the benchmark.
Benchmark (Zymo/Standard-8): removed all 4 Shigella spp. as E. coli cross-mappers,
precision up, zero recall cost. Deviation from the spec below: drops (Grade X)
rather than demoting via contaminant_risk — the benchmark showed removal is needed.

Original spec (historical):


**What:** Post-EM deduplication step to collapse phylogenetically redundant detections.

**Why:** Zymo validation shows *Shigella sonnei* and *Shigella flexneri* appearing as false positives
alongside *E. coli* — caused by 95% genome identity. EM distributes reads but doesn't eliminate the calls.

**Implementation:** Add a `deduplicate_closely_related(entries, ani_threshold=0.97)` function in
`pathogeniq/report.py`. When two organisms exceed ANI threshold and one has >10× the read count,
flag the minor one as `contaminant_risk=True` with note "Likely cross-mapping from {major_organism}".

**Scope:** Only needed for a handful of known problematic pairs — hardcode a `KNOWN_CROSSMAPPERS`
dict initially rather than computing ANI at runtime.

---

## Plan 4 — NTC-Aware Grading (from eng review 2026-06-21)

See design doc `~/.gstack/projects/alvin8-git-PathogenIQ/alvin-master-design-20260621-100903.md`
for the full plan and GSTACK REVIEW REPORT.

**Method code: COMPLETE** (T1-T7 + T3, committed). taxon_id join key → NB
background model → tier-capped grading → consolidated renderers → Tier-2
default mechanism + curation recipe. 111 unit tests pass.

### 4-BLOCKER — Tier-2 default needs MORE blank datasets (do NOT ship one blank)

**What:** Build `pathogeniq/data/background_default.tsv` from a pool of several
real no-template / water-blank datasets, not a single control. Keep the shipped
file an empty placeholder until this is done.

**Why (learned 2026-06-21):** the obvious move — pool ENA project ERP006808 — is
WRONG twice over:
- ERP006808 is the Salter *S. bongori* serial-dilution experiment, not blanks.
  Pooling all 33 runs made a "background" dominated by the *Salmonella* spike
  (43% RPM), which would suppress real Enterobacteriaceae. Only `ERR588954`
  ("Water") is a true negative control.
- That one water blank has just 157 classified reads → *Pseudomonas* at 98%
  (real kitome but over-aggressive from one blank) and single-read detections of
  serious pathogens (*Rickettsia*, *Ehrlichia*) that became spurious background
  floors. The `min_reads` floor (now in build_background, default 2) removes the
  singletons, but one thin blank still can't anchor a credible default.

**Status (v1 SHIPPED):** `background_default.tsv` is now populated (14 taxa)
from 15 spike-free (<=10% Salmonella) runs of Salter ERP006808, selected by
`scripts/05_select_kitome_controls.py`. Tier 2 auto-activates. Caveat header
in-file: the Enterobacteriaceae are genuine kitome (verified stable across spike
thresholds — NOT spike cross-mapping), but the Shigella entries are E. coli/
Shigella reference redundancy (Plan-3D dedup), and it's a foreign 2014 prior.

**Still to improve (not blocking):** (a) pool non-spiked blank datasets (other
kitome studies) for a cleaner, less Enterobacteriaceae-heavy prior;
(b) Plan-3D dedup to collapse the E. coli/Shigella artifact; (c) periodically
refresh from the lab's own NTCs. The cleanest results still need a per-batch
NTC (Tier 1) — the shipped prior is a documented fallback, capped at Grade B.

### 4-FollowUp — Validate the single-NTC NB dispersion prior — RESOLVED 2026-06-22

**Done:** `scripts/10_validate_dispersion.py` +
`docs/dispersion-validation-2026-06-22.md`. Leave-one-out FPR across the 15
spike-free Salter blanks, swept over dispersion.

**Finding:** the single-NTC controlled-α guarantee IS hollow — but not because
the prior is mis-set. LOO FPR is 30-65× α (0.01) at *every* dispersion
(0.30→0.65 across 6 orders of magnitude); it never approaches α. The limit is
**NTC coverage, not the prior**: ~10/18 kitome taxa occur in only one blank, so
they hit the pseudocount floor in leave-one-out and read as real at any r.
Sweeping dispersion cannot fix a coverage gap.

**Decision:** keep `_DEFAULT_DISPERSION = 2.0` (lowering it can't deliver
α-control and worsens the documented real-organism over-suppression); keep
Tier-2 capped at Grade B (this is its empirical justification). Real fix =
batch-matched NTC (Tier 1) or many more pooled blanks; per-taxon dispersion only
becomes worth modeling with ≥2–3 batch-matched NTCs.

### T10 — Benchmark (DONE — validated across 9 held-out communities)

**Built + tested:** `pathogeniq/benchmark.py` (Kraken2 parser, grading adapter
with a Wilson-CI bootstrap stand-in, CAMI-profile + reads_mapping truth loaders,
PR scoring) + `scripts/06_benchmark.py` (single community) + `08_heldout_pr_auc.py`
(multi-community held-out, per-row rank) + `07_build_kraken_db.py` +
`09_reads_mapping_truth.py` (reads_mapping -> truth at a chosen rank).

**Result (`docs/benchmark-results-2026-06-21.md`):** floor=500 calibrated ONCE
on 3 Zymo runs (F1) then applied to all held-out communities — none seen in
calibration, spanning two taxonomic ranks:
- **CAMI species panel** (gut/marine/strain): mean precision 0.011 -> 0.494 (45x)
  at recall 0.684. Strain-madness weakest (0.326) — near-identical-strain
  cross-mapping, the Plan-3D case.
- **HMP genus panel** (skin/airway/oral): mean precision 0.008 -> 0.352 (44x) at
  recall 0.952. Scored at genus rank because the 97%-OTU truth resolves only to
  genus; kept a separate panel (mixing ranks in one mean is apples-to-oranges).

The grading wedge is community- and rank-agnostic, not Zymo-overfit.

**Metric finding (baked into the harness):** grading's value is FP removal, which
moves operating-point precision; PR-AUC stays flat when those FPs are already
bottom-ranked by read count. Both deltas reported.

**Known net-negative:** the Salter kitome NB background over-suppresses real
Enterobacteriaceae (drops recall) — do NOT auto-apply; needs Plan-4 + non-spiked
blanks (below). Cross-map dedup IS a clean precision win (zero recall cost).

**Still open (not blockers, enhancements):**
- CZID-export parser for a third external-baseline config.
- taxid<->GCF map so the GCF-keyed NB background can join Kraken taxids (or
  rebuild the kitome background taxid-keyed through Kraken2).
- Plan-4: sweep the dispersion prior, measure FPR at fixed α on the held-out
  split, report sensitivity — the path to making the background a net positive.

---

## Future Exploration

### Air-sample / bioaerosol surveillance — TO EXPLORE LATER

**What:** Adapt the pipeline to bioaerosol/air-sample surveillance. Air is even
lower-biomass + more kitome-contaminated than clinical specimens, so the NTC work
is the backbone. Source: Gemini design (https://gemini.google.com/share/1aee75f44cf8),
content captured below.

**Design (from the Gemini conversation):**
- Pipeline: QC + low-complexity filter -> host depletion -> Kraken2+Bracken ->
  statistical NTC subtraction (drop taxa with Z < 10) -> high-threat flag ->
  minimap2 verification + coverage check -> alert.
- NTC thresholding: RPM-normalize (or spike-in scale); three options — fold-change
  (sample > K x NTC, K>=5-10), mean+3sd, or Z = (rpm_sample - mu_ntc)/sigma_ntc
  with Z>10 to alert. Zero-deviation fix: pseudocount + a hard floor (>=3 unique
  non-overlapping reads) for high-consequence pathogens.
- **Coverage breadth over read depth:** 100 reads at one genomic position = PCR/
  contamination artifact; 100 reads across 80% of the genome = real threat.

**Already built here (the NTC design is largely DONE):**
- RPM normalization, pseudocount, hard read floor (min_reads), `--ntc` flag,
  build-an-NTC-dictionary-and-subtract, tiered statistical subtraction — all in
  `background.py` + `cli.py`. Our NB tail test is a more principled version of
  their Z-score; the Plan-4 dispersion work refines exactly their mu/sigma model.
- Kraken2 adapter + grading + PR benchmark (`benchmark.py`, `06/08` scripts).
- Cross-mapping dedup (`crossmap.py`).

**NEW ideas worth pursuing (gaps the air design highlights):**
1. **Coverage breadth** — EXPLORED 2026-06-21 (`pathogeniq/coverage.py`,
   `docs/coverage-exploration-2026-06-21.md`). Built + tested. Finding: breadth_ratio
   (Lander-Waterman normalized) catches CLUMPED artifacts (PCR dupes, low-complexity,
   partial/distant cross-mapping — the air "one B. anthracis locus" case) but does
   NOT solve near-identical cross-mapping (E. coli/Shigella both ~75-92% covered,
   genomes >95% identical). To integrate: a verification gate on flagged taxa
   (align.py already has the PAF; it discards the positions). NOT the E. coli/Shigella
   fix — that stays with EM/LCA + the crossmap dedup.
2. **Kaiju** (protein-level, six-frame) as a complementary classifier for divergent/
   novel viruses that nucleotide k-mers miss — relevant for air viral surveillance.
3. **Air-specific contaminant priors** — Cutibacterium, Ralstonia, Pseudomonas kit
   contaminants (extend `contaminants.py` / a non-spiked air-NTC background).
4. **Bracken** abundance re-estimation after Kraken2; **spike-in (ERCC/Zymo) scaling**
   for cross-run RPM comparability.

Not scoped yet — parked for a future session.

---

## Plan 6 — Incorporations from Jeilu et al. 2025 (aircraft air metagenomics)

Source: Jeilu et al. 2025, *Microbiome* 13:249 (doi:10.1186/s40168-025-02276-7),
PDF in `databases/CAMI/`. Their pipeline: KneadData QC → MetaPhlAn4 + HUMAnN3 +
StrainPhlAn (read-based) **and** MEGAHIT → MetaBAT2 → CheckM → GTDB-Tk → Prokka →
Abricate/CARD (assembly-based), with BHI **enrichment** culture, **qPCR + synthetic
spike-in** recovery normalization, and diversity/ordination stats.

**Key contrast (validates our design):** their shotgun profiling does NOT
statistically subtract a kitome background — they report raw MetaPhlAn profiles
and their top taxa (C. acnes, S. epidermidis) are exactly the skin/kit
contaminants our NTC background + contaminant priors + grading are built to flag.
Our NTC-aware grading, bootstrap CIs, and cross-map dedup are rigor their analysis
lacks. We are weaker on assembly, functional profiling, strain ID, and absolute
quantification — the gaps below.

| # | Incorporation | Value | Relevance | Effort |
|---|---------------|-------|-----------|--------|
| 1 | ~~**PhiX removal in QC**~~ — **DONE 2026-06-22**: `run_phix_removal` in `host_remove.py` (minimap2 vs bundled `data/phix.fasta` NC_001422.1, sr/map-ont), wired after host removal in the sample + NTC paths; non-blocking if ref missing. | High | All samples | Low |
| 2 | **Spike-in recovery normalization / absolute quantification** — synthetic DNA spike at known copies → recovery efficiency → absolute load (not just RPM). Enables cross-run/cross-site comparability + pathogen load. (Supersedes the parked "ERCC/Zymo scaling" note.) | High | Clinical load + air surveillance | Med |
| 3 | **Optional assembly + MAG stage** — MEGAHIT → MetaBAT2 → CheckM (≥50%/≤10%) → GTDB-Tk → Prokka. Recovers novel/divergent/uncultured organisms read-based profiling misses, and links ARGs to a genome. Air has many unknown environmental taxa. | High | Air surveillance (novel orgs) | High |
| 4 | **Enrichment (BHI culture) as a documented low-biomass wet-lab option** — boosts DNA yield, MAG + ARG recovery; cost = culture bias toward Firmicutes. Not pipeline code — a sample-prep SOP for Plan 5 air NTC/spike runs. | Med | Air (very low biomass) | Low (doc) |
| 5 | **Strain-level resolution** — StrainPhlAn-style reference strain tracking; complements the parked MLST (Plan 3C). Outbreak / source attribution. | Med | Clinical + surveillance | Med |
| 6 | ~~**Virulence-factor screen**~~ — **DONE 2026-06-22**: `run_virulence_screen` (ABRicate vs VFDB, captures PRODUCT as the factor) runs alongside the CARD AMR screen; rendered as a separate "Virulence Factors (VFDB)" section in JSON/PDF/HTML. | Med | Pathogen interpretation | Low |
| 7 | **Multi-sample diversity + ordination** — alpha (Shannon), beta (Bray-Curtis/UniFrac), PCoA, PERMANOVA across samples/sites/time for a surveillance dashboard. | Med | Air surveillance (fleets/time) | Med |
| 8 | **GTDB taxonomy option** — GTDB-Tk gives consistent genome-based names for environmental taxa NCBI lacks; offer as an alt to NCBI taxids for air. | Low | Air | Med |
| 9 | **Overrepresented-sequence removal** (FastQC flag) in QC — minor extra contaminant trim they apply. | Low | All | Low |

**Suggested order:** #1 (PhiX, quick correctness) → #6 (VFDB, near-free) →
#2 (spike-in absolute quant) → #3 (assembly/MAG, the big surveillance unlock,
pairs with the PRJNA1228129 air data) → the rest as surveillance scope grows.

---

## Wet-lab experimental design — air metagenome runs (Plan 5 prerequisite)

Air is ultra-low-biomass and kitome-dominated, so the **controls are not optional
extras — they are what makes a PathogenIQ result interpretable.** Each control
below unlocks a specific pipeline feature. Numbers are starting points; the
true-positive thresholds and spike dose MUST be titrated on the lab's own
reagents/sampler (a real workflow needs the calibration knob, not just defaults).

### Controls matrix — what to run and why

| Control | What it is | How many | Unlocks in PathogenIQ |
|---------|------------|----------|-----------------------|
| **Batch-matched NTC** (extraction blank) | Buffer/water through the *entire* extraction + library + sequencing protocol, same reagent lot, same run | **≥1 per extraction batch/run — mandatory** | Tier 1 → the *only* path to **Grade A**. Feeds `--ntc` for this batch's NB background |
| **Field / collection blank** | A sterile filter/cassette opened at the sampling site but not actively sampling, then processed identically | **≥1 per sampling session** | Separates *sampling/transport* contamination from *reagent* kitome; without it, field contaminants masquerade as real air taxa |
| **No-template library control** | Water into library prep only (no extraction) | 1 per library batch | Isolates library/index contamination + index hopping from extraction kitome |
| **Pooled-background NTCs** | The accumulating set used to build `background_default.tsv` | **≥8–12 across reagent lots**, then add each new batch NTC over time | The Tier-2 default. Dispersion study: *coverage* is the limit → more NTCs, spanning lots, is what lowers the false-positive floor |
| **Internal spike-in** (absolute quant) | Whole cells absent from air/clinical — **ZymoBIOMICS Spike-in Control I/II** (*Imtechella halochordis* + *Allobacillus halotolerans*), added **before extraction** | **Every sample** (it is the per-sample ruler) | Absolute load (copies/m³ air) + extraction-recovery normalization (Plan 6 #2). Add the spike genomes to the Tier-1 DB |
| **Positive process control** | **ZymoBIOMICS D6300** mock community spiked onto a blank matrix (à la Jeilu et al.) | **1–3 per batch** | End-to-end sensitivity + abundance-accuracy check; the integration-test analogue for the air domain |
| **Technical replicates** | Same sample split/extracted ≥2× | **2–3 per sample** at low biomass | Distinguishes stochastic single-read contamination from reproducible real signal (grading already de-weights singletons via `min_reads`) |

### Spike-in dose (the part people get wrong)

- **Target the spike at ~0.5–5% of total reads.** Too little → noisy denominator;
  too much → it swamps the tiny real air signal and wastes depth. Low-biomass air
  sits at the dangerous end — start LOW and titrate.
- **Quantify the exact input by ddPCR/qPCR** (don't trust the nominal tube conc).
  Absolute load = `reads_target / reads_spike × spike_copies_added / air_volume`.
- **Add before lysis/extraction** so the spike experiences the same extraction
  bias as the sample (Gram± lyse differently) — DNA-only spikes miss this.

### Sequencing / sampling notes

- **qPCR (16S) pre-check** each sample before sequencing — confirm any microbial
  DNA is present and gate samples below the limit of detection.
- **Low-input library prep** + extra PCR cycles; **1–2% PhiX** sequencing spike
  (also satisfies Plan 6 #1 PhiX removal). Budget depth for the spike + kitome
  fraction so the real signal still gets enough reads.
- **Sampling captures biomass:** HEPA/cabin filters integrate over thousands of
  hours (high biomass, time-averaged); active samplers (impinger/cyclone/filter)
  are short-term/low-biomass — record air volume sampled for the per-m³ denominator.
- **Optional BHI enrichment arm** (Plan 6 #4) to rescue biomass — run it ALONGSIDE
  the direct/environmental arm, never as a replacement (culture bias → Firmicutes).

### Minimum viable run (smallest design that still grades meaningfully)

Per batch: **1 extraction NTC + 1 field blank + internal spike in every sample +
1 Zymo positive control + duplicate the real samples.** That single batch NTC
enables Grade A (Tier 1); the spike gives absolute load; the field blank +
background separate air signal from contamination; the Zymo control proves the
workflow detects what it should.
