# Changelog

All notable changes to PathogenIQ are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Planned
- Nextflow orchestration (local Docker, SLURM, AWS Batch)
- Validation framework with ROC curves and benchmark suite
- Web dashboard for clinical users
- FHIR-compatible report export

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
| 0.2.0 | ≥3.11 | numpy ≥1.26, scipy ≥1.12, pysam ≥0.22, pandas ≥2.2, click ≥8.1, reportlab ≥4.0 |
| 0.1.0 | ≥3.11 | numpy ≥1.26, scipy ≥1.12, pysam ≥0.22, pandas ≥2.2, click ≥8.1 |
