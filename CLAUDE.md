# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Commands

```bash
# Install (editable, with dev deps)
pip install -e ".[dev]"

# Unit tests (fast, no external tools needed)
pytest tests/ -v

# Single test
pytest tests/test_report.py::test_grade_a -v

# Integration tests (requires databases + external tools)
pytest tests/integration/ -v -m integration

# Lint
ruff check pathogeniq/ tests/

# Type check
mypy pathogeniq/
```

## Architecture

PathogenIQ is a 5-stage clinical metagenomics pipeline. Each stage is a separate module that wraps an external bioinformatics tool and returns a typed dataclass.

**Pipeline flow** (`cli.py` orchestrates):
```
QC (fastp/Chopper) → Host removal (BWA-MEM2/minimap2) → Sketch screen (sourmash)
  → Targeted alignment (minimap2) → EM abundance → Report (JSON + TSV)
```

**Key data contracts between stages:**
- `PipelineConfig` (`config.py`) — single config dataclass passed to every stage; holds paths, `ReadType` (short/long), `SpecimenType` (blood/csf/bal/tissue), and tuning params
- `EMResult` (`em.py`) — `abundances` array (sums to 1) + `n_reads`; produced by `em_abundance()` and consumed by `write_report()`
- `AlignResult` (`align.py`) — `alignment_matrix` (reads × organisms binary array) + `organism_names`; bridges alignment → EM
- `ReportEntry` (`report.py`) — one row per organism with lazy `grade` property; contaminant flag set externally by `flag_contaminants()`

**Evidence grading** (`report.py` + `contaminants.py`):
- Grade A/B/C/X is a property on `ReportEntry`, not a stored field
- Grading uses specimen-specific minimum read thresholds (`_MIN_READS`) + CI width ≤ 0.15 for Grade A
- `contaminants.py` holds `CONTAMINANT_PRIORS` — a per-`SpecimenType` list of common contaminant species; `flag_contaminants()` sets `contaminant_risk=True` via case-insensitive substring match, which demotes Grade A → B/C

**Read-type branching:**
- Short reads: fastp for QC, BWA-MEM2 for host removal
- Long reads: Chopper for QC, minimap2 for host removal
- Both paths produce identical intermediate files; downstream stages are read-type agnostic

**External tool dependencies:** fastp, Chopper, BWA-MEM2, minimap2, samtools, sourmash — must be on `PATH` (install via conda bioconda channel)

## Test structure

- `tests/test_*.py` — unit tests; mock all subprocess/tool calls; run without any databases or external tools
- `tests/integration/test_zymo_validation.py` — end-to-end tests against ZymoBIOMICS D6300 standard; marked `@pytest.mark.integration`; require `HUMAN_REF` and `TIER1_DB` env vars pointing to built databases

## Database setup

Three databases must be built before running the pipeline or integration tests:

```bash
bash scripts/01_download_zymo_db.sh 16        # ZymoBIOMICS validation DB (~5 min)
python scripts/02_download_tier1_db.py --threads 16  # 66-genome clinical DB (~30-90 min)
bash scripts/03_download_human_ref.sh 16       # GRCh38 + index (~90 min, needs ~64 GB RAM)
```

Output lands in `databases/`. See `scripts/README.md` for accession tables.
