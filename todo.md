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

### 4-FollowUp — Validate the single-NTC NB dispersion prior

**What:** Quantify how the single-NTC negative-binomial background test's false-positive
rate varies with the chosen dispersion prior, and pick a defensible default.

**Why:** With one NTC, per-taxon dispersion is not estimable and the test leans on a shared
prior. That prior is the weakest link in the method's FPR claim — if FPR is highly sensitive
to it, the single-NTC (Tier 2) path's controlled-α guarantee is hollow.

**Depends on:** the multi-community truth set (the P1 data-acquisition task for the benchmark).
Cannot run until labeled data exists.

### T10 — Benchmark (CORE BUILT; full run blocked on data)

**Built + tested:** `pathogeniq/benchmark.py` (Kraken2 report parser,
kraken_to_grading_inputs adapter with a Wilson-CI stand-in for the bootstrap,
and PR scoring: average_precision / precision_at_recall) + `scripts/06_benchmark.py`
(kraken-raw vs kraken+grading on a Kraken2 report vs a truth-taxid set).
Smoke-tested on synthetic data.

**Metric finding (baked into the harness output):** grading's value is
false-positive removal, which moves *operating-point precision*; PR-AUC stays
flat when those FPs are already bottom-ranked by read count. Report both deltas.

**Blockers for the real run:**
- Kraken2 DB + run on CAMI/Zymo reads to produce `.kreport` files (CAMI reads
  are portal-gated at data.cami-challenge.org; Zymo reads are on disk).
- Truth taxid lists from the CAMI gold-standard profiles / known Zymo composition.
- taxid<->GCF map so the NB background filter (GCF-keyed) can join Kraken taxids
  (or rebuild the kitome background through Kraken2, taxid-keyed).
- A CZID-export parser for the third config.
- Held-out split once there are >=3-5 labeled communities (one Zymo can't anchor
  a PR curve — outside-voice c).

**Where to start:** sweep the dispersion prior over a plausible range, measure FPR at fixed α
on the held-out split, report sensitivity.
