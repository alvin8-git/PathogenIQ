# Changelog

All notable changes to PathogenIQ are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Planned
- Nextflow orchestration (local Docker, SLURM, AWS Batch)
- Batch-matched NTC (Tier 1) prospective evaluation for Grade-A claims
- Calibrate provisional cutoffs on labeled air data (breadth ratio, open-world completeness bands, air read floors)
- Viral-contig pathogenicity triage (ARG-carrying phage)
- Web dashboard for clinical users; FHIR-compatible report export

---

## [0.4.2] — 2026-07-02 — MAG Toolchain, Novelty Calibration & Report Rendering

### Added
- **MAG toolchain installer** `scripts/15_setup_mag_env.sh` — provisions metabat2 + CheckM (isolated `mag-bin` conda env) and GTDB-Tk (`gtdbtk` env) with thin `PATH` wrappers into the core env, same isolation discipline as the viral tools. Fetches the CheckM DB (~1.4 GB) and GTDB-Tk **r232** DB (~94 GB) onto a data volume (not the root partition), with loud FATAL verification on every stage.
- **Air assembly runner** `scripts/16_run_air_assemble.sh` — runs the `--assemble` MAG arm on the 5 PRJNA1228129 aircraft-filter samples (fail-fast toolchain check, sequential so GTDB-Tk's pplacer has the box's RAM).
- **MAG table + absolute-copies column in the PDF/HTML reports** (`pdf_report.py`, `html_report.py`) — the assembly arm's MAGs (bin, GTDB lineage, completeness, contamination, contigs, size, open-world grade) now render alongside the findings; an "Abs. copies" column appears only when a spike-in anchored quantification (so unspiked runs are not cluttered). Wired through `cli.py`; unit-tested.

### Validated
- **R1 novelty threshold** calibrated on real air data (`docs/novelty-threshold-validation-2026-07-01.md`). Kept `_DEFAULT_FLAG_THRESHOLD = 0.5`: it sits just above the air-NTC ceiling (0.472) so it never false-fires on kitome, but filters (≤0.478) and NTCs overlap, so the bare fraction (and its rank-resolution refinement) cannot separate assemble-worthy air from kitome. The gate is therefore **advisory** (–`--assemble` runs regardless); reliable air novelty is a post-assembly, per-MAG judgement.
- **MAG toolchain installs clean** — metabat2/CheckM/GTDB-Tk wrappers resolve on the core `PATH`; GTDB-Tk r232 (94 GB) + CheckM (1.4 GB) DBs verified on disk.

### Changed
- **AMR + VFDB now screen assembled contigs, not raw reads** — ~100× fewer sequences to align and clean full-length gene calls (measured: abricate-on-reads was ~16 min/DB on a 6.5 M-read air sample). The same contig set is shared with the viral arm, so assembly is paid for once.
- **EM bootstrap CIs parallelized across processes** — the resampling loop now fans out over cores.

---

## [0.4.1] — 2026-06-24 — AMR/VFDB Organism Attribution

### Fixed
- **Contig-based AMR/VFDB hits now attribute to findings.** With contig screening the ABRicate `SEQUENCE` is the assembly contig id (`k141_42`), which carries no organism name, so the old substring match always returned `unknown` and every finding showed empty `amr_genes`/`virulence_factors`. `map_contigs_to_organisms()` (`amr.py`) now aligns the contigs back to the screened genomes (minimap2) and attributes each hit to the organism(s) its contig maps to.
- **Novelty `flag_threshold` surfaced** in the report (`novelty.py`, `report.py`) — `flagged: false` is now interpretable against the bar it was tested against (the flag logic was already correct: 0.42 unclassified < 0.50).

### Added
- **Co-attribution across near-identical siblings** — a co-assembly collapses *E. coli* / *Shigella* into one shared contig set, so a contig is attributed to **every** genome within `_COMAP_MARGIN` (0.70) of its best aligned-bp, not just the single best. A shared core gene co-attributes to each sibling; an organism-specific region stays single. `AMRHit`/`VirulenceHit` gain `organism_matches` (best-first); the report groups by it.
- Tests: contig→organism mapping, co-mapping of near-ties vs specific-stays-single, the no-map "unknown" regression, and non-blocking fallback (`test_amr.py`).

### Validated
- **Air demo (3 Enterobacteriaceae):** *Shigella sonnei* 1→34 AMR / 0→16 VFDB, *Shigella flexneri* 0→30 / 0→12 — both 100% shared with *E. coli*'s core set (correct); *Pseudomonas* keeps its distinct `MexA`/`ParR` (0 shared across the family boundary); *E. coli*-specific `marA`/`emrE`/`kdpE`/`ugd` did not leak to *Shigella*.

---

## [0.4.0] — 2026-06-23 — Air Surveillance & Open-World Detection

This release extends PathogenIQ from clinical fluids to **air / bioaerosol surveillance** and adds an **open-world** arm that detects pathogens not in the reference database — including viruses and novel organisms. Design and validation: `docs/air-open-world-detection-2026-06-23.md`, `docs/air-pathogen-wedge-design-2026-06-22.md`, `docs/air-concordance-validation-2026-06-22.md`.

### Added
- **AIR specimen type** (`config.py`) — air/bioaerosol read floors + an air kitome contaminant prior (*Cutibacterium*, *Sphingomonas*, *Methylobacterium*, *Ralstonia*…), deliberately **not** genus-broad on *Pseudomonas*/*Acinetobacter* so real pathogens aren't demoted
- **Flag-not-subtract for dual-use pathogens** (`background.py`) — a contaminated NTC flags (not erases) *E. coli*, *Pseudomonas*, *Klebsiella*, *Shigella*, *Salmonella*, *S. aureus*, *Acinetobacter*, *Enterobacter*; pure contaminants still subtracted
- **Breadth-of-coverage gate wired into grading** — per-organism breadth computed from the alignment PAF (`align.py` → `coverage.py`); ratio < 0.25 → Grade X (clumped artifact). Surfaced as `breadth_ratio` in the report
- **PhiX removal** (`host_remove.py`) — minimap2 vs the bundled Illumina PhiX spike-in (NC_001422), non-blocking
- **VFDB virulence screen** + **spike-in absolute quantification** (`amr.py`, `quantify.py`) — copies / copies-per-volume anchored on a known spike
- **Open-world arms**:
  - R1 **novelty trigger** (`novelty.py`, `--novelty`) — broad-Kraken2 unclassified "dark-matter" fraction
  - R2 **viral arm** (`viral.py`, `--viral`) — geNomad (ID + ICTV taxonomy) → CheckV (completeness)
  - R3 **assembly/MAG recovery** (`assembly.py`, `--assemble`) — MEGAHIT → MetaBAT2 → CheckM → GTDB-Tk
  - R4 **pathogenicity triage** (`pathogenicity.py`) — VFDB/CARD markers + GTDB phylo-proximity → PATHOGEN_CANDIDATE / PATHOGEN_ADJACENT / ENVIRONMENTAL
  - R5 **open-world grading** (`report.py::grade_open_world`) — A/B/C/X for MAGs/viral from completeness + contamination + marker signal, **capped at B** (no same-run NTC)
- **CLI flags** `--novelty`, `--viral`, `--assemble`, `--spike-taxon/--spike-copies/--sample-volume`; JSON `novelty`/`viral`/`mags`/`pathogenicity` blocks; the open-world arms run even when the targeted screen finds zero hits
- **Tooling isolation** — `scripts/13_setup_viral_env.sh` installs geNomad/CheckV (and ABRicate) in their own conda envs with `PATH` symlinks/wrappers, so the core NumPy-2 stack is never perturbed
- **Validation harness & data** — `scripts/12_viral_insilico_spikein.py` (in-silico viral spike, offline `--selfcheck`); air datasets wired into `scripts/04` (Jeilu et al. 2025, PRJNA1228129)
- Unit tests for every new module (`test_novelty.py`, `test_viral.py`, `test_pathogenicity.py`, plus open-world grade + breadth + dual-use cases in `test_report.py`)

### Validated
- **Viral arm: 100% recall (3/3)** — in-silico T4 / Lambda / SARS-CoV-2 recovered with correct ICTV lineage from a real aircraft-filter background
- **Air concordance** — Zymo spike recall 6/10 under the air NTC; the 4 misses were **background over-subtraction**, the empirical origin of flag-not-subtract
- **Breadth gate** — real Zymo organisms score 0.79–0.96 (pass); the 0.25 cutoff catches clumped artifacts without touching real low-abundance hits

### Fixed
- `run_megahit` silently produced nothing when `assembly/` did not pre-exist (megahit's `os.mkdir` is single-level) — had made `--assemble` **and** `--viral` no-op on every run; caught by the viral harness
- geNomad/CheckV install had downgraded the core env's NumPy and broken SciPy — resolved by env isolation (above); pipeline env restored

---

## [0.3.0] — 2026-06-22 — NTC-Aware Grading & Multi-Community Benchmark

The headline of this release is **Plan 4**: a statistically principled, provenance-aware grading layer that subtracts reagent/kitome background and is validated against a strong classifier across nine held-out communities.

### Added
- **NTC background subtraction** (`pathogeniq/background.py`)
  - `BackgroundModel` (per-taxon RPM rates, control count, dispersion, tier) built from one or more no-template controls
  - Negative-binomial upper-tail test (`background_pvalue` / `is_background`) — depth-invariant via RPM normalisation, with a pseudocount floor (0.5 RPM) and a min-support floor (≥2 reads) so thin blanks don't over-suppress
  - `build_background`, `load_background_table`, `write_background_table`, `load_default_background` — pooled Tier-2 table format with a `# tier=N` header
- **Tiered grading** — Grade A requires a batch-matched NTC (Tier 1); pooled/foreign (Tier 2) and no-NTC (Tier 3) cap at Grade B
- **Pure `grade(GradingInput)` function** (`pathogeniq/report.py`) — read count, abundance, CI width, contaminant flag, specimen, tier, cross-map flag → grade; degenerate stats (NaN/inf/negative) → X. Auditable and classifier-agnostic.
- **Cross-mapping dedup** (`pathogeniq/crossmap.py`) — flags the minor member of a known >95%-identical group (E. coli / Shigella spp.) when outnumbered ≥10× by a relative, demoting it to Grade X (zero recall cost)
- **`taxon_id` threading** — GCF/taxid join key flows through `sketch.py` → `align.py` → `report.py` so the background and grading layers share one identifier space; `build_entries()` consolidates construction, dedup, contaminant flagging, and background filtering
- **Breadth-of-coverage** (`pathogeniq/coverage.py`) — Lander-Waterman expected breadth + depth-normalised breadth ratio (catches clumped/PCR artifacts)
- **HTML clinical report** (`pathogeniq/html_report.py`) rendering the same `ReportEntry` rows
- **Benchmark toolkit** (`pathogeniq/benchmark.py`) — Kraken2 report parser, `kraken_to_grading_inputs` adapter (Wilson-CI stand-in), CAMI `.profile` + `reads_mapping.tsv` truth loaders (per-sample / per-rank), and PR scoring (`average_precision`, `precision_at_recall`, `precision_recall`)
- **Validation & data scripts** — `04` ENA/NTC + blank downloader (title filter, run limit), `05` kitome-control selection, `06` single-community benchmark, `07` custom Kraken2 DB builder, `08` multi-community held-out PR-AUC (per-rank), `09` reads_mapping→truth rollup, `10` dispersion-prior validation, `11` blank-pool background builder
- **CLI flags** `--ntc`, `--background`, `--no-background`; shipped Tier-2 `pathogeniq/data/background_default.tsv`
- **`environment.yml`** conda environment (adds kraken2, bracken for benchmarking)
- Unit tests for every new module (`tests/test_background.py`, `test_benchmark.py`, `test_crossmap`-equivalents, `test_html_report.py`, `test_download_validation.py`, `test_select_kitome.py`, plus `--selfcheck` harnesses in scripts 09/10)

### Validated
- **Nine held-out communities, two ranks.** A read-count floor calibrated once on three Zymo D6300 runs lifts operating-point precision ~45× at preserved recall, generalising to CAMI II species communities (mouse-gut/marine/strain, mean P 0.011→0.494) and HMP genus body sites (skin/airway/oral, mean P 0.008→0.352). Details: `docs/benchmark-results-2026-06-21.md`.
- **Dispersion-prior study.** Leave-one-out FPR on kitome blanks is 30–65× α at *every* dispersion — coverage, not the prior, is the limit; the empirical justification for capping Tier-2 at Grade B. Details: `docs/dispersion-validation-2026-06-22.md`.

### Changed
- `ReportEntry`: added `taxon_id`, `tier`, `crossmap_of` fields; `grade` / `invalid_stats` now delegate to the pure `grade()` function via `as_input()`
- `AlignmentResult` gained a required `taxon_ids` field; `SketchHit` gained `taxon_id`
- `write_report()` now renders pre-built entries from `build_entries()`; JSON/TSV include `taxon_id`, `tier`, `contaminant_risk`, `invalid_stats`, `crossmap_of`
- Tier-1 clinical DB expanded to ~110 pathogen genomes; `CLAUDE.md` architecture + context-mode routing refreshed

### Fixed
- E. coli/Shigella reference cross-mapping inflating phantom polymicrobial calls (now deduped to Grade X)
- CAMI gold-standard profiles pooled all samples; `parse_cami_profile` now selects one sample's `@SampleID` block
- `load_truth` auto-detects CAMI `.profile` vs plain taxid lists (benchmark TP=0 bug)
- Kitome background provenance corrected: Salter Enterobacteriaceae are genuine kitome (stable across spike thresholds), not spike cross-mapping

---

## [0.2.0] — 2026-04-23 — Clinical Interpretation Engine

### Added
- **Specimen-aware contaminant suppression** (`pathogeniq/contaminants.py`)
  - Per-specimen contaminant priors for blood, CSF, BAL, and tissue
  - Cutibacterium acnes, Staphylococcus epidermidis, Bacillus spp., and others flagged automatically in blood and CSF
  - Oral flora (Streptococcus salivarius, Prevotella melaninogenica, etc.) flagged in BAL
  - Grade downgrade: flagged contaminants cannot reach Grade A; Grade B entries demote to Grade C
  - `contaminant_risk` boolean field in JSON and TSV report output
- **Antimicrobial resistance (AMR) gene detection** (`pathogeniq/amr.py`)
  - Wrapper around ABRicate + CARD database for AMR gene screening
  - Detects resistance genes (e.g., mecA, blaCTX-M-15, vanA) with identity/coverage thresholds (default: 90% / 80%)
  - Links AMR hits to detected organisms via genus/species substring matching
  - Graceful fallback: returns empty list if ABRicate is not installed
  - `AMRHit` dataclass with gene, drug class, identity, coverage, organism match, and database fields
- **PDF clinical report** (`pathogeniq/pdf_report.py`)
  - Single-page formatted report using ReportLab Platypus
  - Header with sample name, specimen type, read type, and timestamp
  - Findings table with colour-coded grades (A=green, B=orange, C=grey, X=red)
  - AMR gene summary table (when hits are present)
  - Grade legend and research-use disclaimer footer
  - Output path: `report/pathogeniq_report.pdf`
- **`--no-pdf` CLI flag** to skip PDF generation
- **`--amr-db` CLI flag** to select ABRicate database (default: `card`; supports `resfinder`, `vfdb`, etc.)
- **`amr_db` field** added to `PipelineConfig`
- Unit tests for all new modules:
  - `tests/test_contaminants.py` — 7 tests covering flagging logic, grade downgrade, and registry completeness
  - `tests/test_amr.py` — 6 tests covering TSV parsing, organism matching, unknown organism handling, and mock ABRicate runs
  - `tests/test_pdf_report.py` — 4 tests covering file creation, AMR table inclusion, output path correctness, and empty entries
- `reportlab>=4.0` dependency in `pyproject.toml`

### Changed
- `write_report()` signature extended: accepts `amr_hits` and `entries_override` parameters
- `ReportEntry` dataclass: added `contaminant_risk: bool = False` field
- `ReportEntry.grade` property: contaminant flag now blocks Grade A and demotes Grade B to C
- CLI pipeline (`pathogeniq/cli.py`): contaminant flagging, AMR screening, and PDF generation now wired into the `run` command
- `README.md`: pipeline diagram now shows AMR and PDF stages

### Fixed
- **`UnicodeDecodeError` in targeted alignment** (`pathogeniq/align.py`): minimap2 stderr contains ANSI color codes with invalid UTF-8 bytes. Removed `capture_output=True, text=True` from `subprocess.run()` call since minimap2 writes output to a PAF file via `-o`; stderr no longer needs decoding.
- **Defensive encoding in AMR module** (`pathogeniq/amr.py`): changed `text=True` to `encoding="utf-8", errors="replace"` in the abricate `subprocess.run()` call to prevent similar decode failures.
- `reportlab` import error when running tests outside the `pathogeniq` conda environment (resolved by installing pytest and dependencies in the correct environment)
- Circular import between `report.py` and `contaminants.py` resolved with `TYPE_CHECKING` guard

---

## [0.1.0] — 2026-04-22 — Core Pipeline

### Added
- **Project scaffold**: Python package structure, `pyproject.toml`, test harness, and CI/CD workflow
- **PipelineConfig dataclass** (`pathogeniq/config.py`)
  - `ReadType` enum: `SHORT` (Illumina) / `LONG` (Nanopore)
  - `SpecimenType` enum: `BLOOD`, `CSF`, `BAL`, `TISSUE`
  - Configurable parameters: threads, sketch threshold, bootstrap iterations, host reference path, tier-1 DB path
- **QC module** (`pathogeniq/qc.py`)
  - Short reads: fastp wrapper for quality trimming, adapter removal, and length filtering
  - Long reads: Chopper wrapper for Nanopore read quality filtering
  - Returns `QCMetrics` dataclass with total reads, passing reads, and adapter contamination rate
- **Host removal module** (`pathogeniq/host_remove.py`)
  - Short reads: BWA-MEM2 alignment against GRCh38 + decoy sequences
  - Long reads: minimap2 alignment against GRCh38
  - samtools-based unmapped read extraction → microbial-only FASTQ
  - Returns `HostRemovalMetrics` with microbial fraction and host reads removed
- **Sketch screening module** (`pathogeniq/sketch.py`)
  - sourmash MinHash sketching of non-human reads
  - Sequence Bloom Tree (SBT) containment search against 20,000+ reference genomes
  - Threshold filtering (default: 1% containment) → 50–200 candidate organisms
  - Typical runtime: < 60 seconds for full database
  - Returns `SketchHit` list with organism name, accession, and containment fraction
- **Targeted alignment module** (`pathogeniq/align.py`)
  - Per-candidate minimap2 alignment of non-human reads against shortlisted organism genomes
  - PAF (Pairwise mApping Format) parsing into binary read–organism assignment matrix
  - Returns `AlignmentResult` with alignment matrix and organism/accession lists
- **EM abundance estimation module** (`pathogeniq/em.py`)
  - Expectation-Maximization algorithm for resolving multi-mapping reads across organisms
  - Bootstrap resampling (n=100) for 95% Bayesian confidence intervals
  - Returns `EMResult` with abundance vector (summing to 1.0), read count, convergence flag, and iteration count
- **Clinical report module** (`pathogeniq/report.py`)
  - Evidence grading system: Grade A (high confidence), B (moderate), C (low / uncertain), X (insufficient reads)
  - Specimen-aware minimum read thresholds: blood=3, CSF=2, BAL=10, tissue=10
  - Grade A requires: sufficient reads + 95% CI width ≤ 0.15 (15 percentage points)
  - JSON and TSV output with abundance %, CI bounds, read counts, and grades
- **Click CLI** (`pathogeniq/cli.py`)
  - `pathogeniq run` command with `--input`, `--output`, `--db`, `--host-ref`, `--specimen`, `--read-type`, `--threads`, `--sketch-threshold`, `--n-bootstrap` flags
  - Progress messages for each pipeline stage
- **Database download scripts**
  - `scripts/01_download_zymo_db.sh` — ZymoBIOMICS D6300 10-genome validation database (~5 min)
  - `scripts/02_download_tier1_db.py` — Tier-1 clinical pathogen database (~66 genomes; WHO Priority Pathogens, CDC Select Agents; ~30–90 min)
  - `scripts/03_download_human_ref.sh` — GRCh38 + decoy + BWA-MEM2 index (~90 min, ~64 GB RAM)
- **ZymoBIOMICS integration validation** (`tests/integration/test_zymo_validation.py`)
  - End-to-end pipeline test against D6300 standard (8 bacteria at 12% + 2 yeasts at 2%)
  - Checks: all 10 organisms detected, abundance within 20% relative error, no Grade A false positives
  - Skips gracefully when `HUMAN_REF` and `TIER1_DB` environment variables are unset
- Unit tests for all pipeline stages (`tests/test_*.py`) — mock all external tools, run without databases

### Changed
- N/A (initial release)

### Fixed
- Removed invalid `--quiet` flag from sourmash sketch command (not supported in v4.9.4)
- Fixed `ncbi-genome-download` short flags to use correct syntax (`-A`, `-F`, `-o`, `-p`, `-r`)
- Fallback to GenBank format for eukaryotic genomes; tolerance for corrupt gzip files in download script
- Fixed PATH issues with `xargs` by using full sourmash path and simple for-loop

---

## Release Tag Conventions

| Tag | Meaning |
|-----|---------|
| `feat:` | New feature or functionality |
| `fix:` | Bug fix |
| `test:` | Test addition or improvement |
| `docs:` | Documentation update |
| `refactor:` | Code refactoring without behaviour change |
| `chore:` | Build, dependency, or tooling change |

---

## Version Compatibility

| Version | Python | Key Dependencies |
|---------|--------|------------------|
| 0.4.x | ≥3.11 | core: numpy ≥1.26, scipy ≥1.12, pysam ≥0.22, pandas ≥2.2, click ≥8.1, reportlab ≥4.0; external (core env): kraken2, minimap2, megahit; external (isolated `numpy<2` envs, `PATH`-wrapped): geNomad, CheckV, ABRicate, metabat2, CheckM, GTDB-Tk |
| 0.3.0 | ≥3.11 | numpy ≥1.26, scipy ≥1.12, pysam ≥0.22, pandas ≥2.2, click ≥8.1, reportlab ≥4.0; external: kraken2, bracken (benchmarking only) |
| 0.2.0 | ≥3.11 | numpy ≥1.26, scipy ≥1.12, pysam ≥0.22, pandas ≥2.2, click ≥8.1, reportlab ≥4.0 |
| 0.1.0 | ≥3.11 | numpy ≥1.26, scipy ≥1.12, pysam ≥0.22, pandas ≥2.2, click ≥8.1 |
