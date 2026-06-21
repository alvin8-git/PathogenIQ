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

### 3D — False Positive Deduplication (Shigella/E.coli cross-mapping)

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

### 4-FollowUp — Validate the single-NTC NB dispersion prior

**What:** Quantify how the single-NTC negative-binomial background test's false-positive
rate varies with the chosen dispersion prior, and pick a defensible default.

**Why:** With one NTC, per-taxon dispersion is not estimable and the test leans on a shared
prior. That prior is the weakest link in the method's FPR claim — if FPR is highly sensitive
to it, the single-NTC (Tier 2) path's controlled-α guarantee is hollow.

**Depends on:** the multi-community truth set (the P1 data-acquisition task for the benchmark).
Cannot run until labeled data exists.

**Where to start:** sweep the dispersion prior over a plausible range, measure FPR at fixed α
on the held-out split, report sensitivity.
