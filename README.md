# PathogenIQ

[![CI](https://github.com/alvin8-git/PathogenIQ/actions/workflows/ci.yml/badge.svg)](https://github.com/alvin8-git/PathogenIQ/actions/workflows/ci.yml)
[![Python 3.11+](https://img.shields.io/badge/python-3.11+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Code style: ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://github.com/astral-sh/ruff)
[![Platforms](https://img.shields.io/badge/platform-linux%20%7C%20macOS-lightgrey)](https://github.com/alvin8-git/PathogenIQ)

A lightweight, clinical-grade metagenomic pathogen detection pipeline for Illumina and Nanopore sequencing data. PathogenIQ uses a **sketch-first** architecture to screen 20,000+ genomes in under 60 seconds before targeted alignment and EM-based abundance estimation — delivering clinical reports in under one hour.

---

## Key Features

- **Sketch-first screening** — sourmash MinHash pre-filters to 50–200 candidates; no wasted alignment cycles
- **EM abundance estimation** — Expectation-Maximization resolves multi-mapping reads with bootstrap Bayesian confidence intervals
- **Evidence grading (A/B/C/X)** — specimen-aware thresholds suppress common contaminants automatically
- **Dual read-type support** — Illumina short reads (BWA-MEM2) and Nanopore long reads (minimap2)
- **All specimen types** — blood, CSF, BAL, tissue, with specimen-specific interpretation priors
- **Tiered pathogen database** — WHO Priority Pathogens, CDC Select Agents, common clinical viruses, fungi, parasites
- **Novel pathogen flagging** — hits outside the curated database are surfaced for review
- **AMR gene detection** — ABRicate+CARD scans non-human reads for antimicrobial resistance genes (e.g., mecA, blaCTX-M-15)
- **PDF clinical report** — formatted single-page report with grade colour-coding, AMR summary, and research-use disclaimer

---

## How PathogenIQ Compares

| Feature | CZID | Kraken2 + Bracken | **PathogenIQ** |
|---|---|---|---|
| **Architecture** | Cloud-only (AWS) | Local only | Local + cloud burst |
| **Read types** | Short + Long | Short read optimised | Short + Long |
| **Typical turnaround** | 4–8 hours | < 30 min | **< 1 hour** |
| **Screening method** | BLAST (NT/NR) | k-mer exact match | **MinHash sketch → targeted align** |
| **Abundance estimation** | Read count / RPM | Bracken re-estimation | **EM with bootstrap CI** |
| **Confidence scoring** | None | None | **Bayesian CI + Grade A/B/C/X** |
| **Specimen-aware priors** | No | No | **Yes** |
| **Host removal** | BWA | Bowtie2 / KneadData | BWA-MEM2 / minimap2 |
| **False positive control** | Minimum read filter | None | **Graded evidence system** |
| **Novel pathogen alert** | No | No | **Yes** |
| **Clinical report** | Web UI only | None | **JSON + TSV + PDF** |
| **AMR detection** | No | None | **ABRicate + CARD** |
| **Deployment** | Cloud-only ($$) | Local | **Hybrid (free self-hosted)** |
| **ZymoBIOMICS validated** | Yes | Partial | **Yes (integration tests)** |

> **Bottom line:** CZID is excellent for research but slow and cloud-locked. Kraken2+Bracken is fast but lacks confidence quantification and clinical context. PathogenIQ delivers sub-hour turnaround with statistically principled abundance estimates and a specimen-aware reporting layer designed for clinical decision support.

---

## Pipeline Overview

```
FASTQ (SR or LR)
     │
     ▼
 ┌─────────┐    fastp (SR) / Chopper (LR)
 │  QC     │──→ quality trim, adapter removal, length filter
 └─────────┘
     │
     ▼
 ┌──────────────┐    BWA-MEM2 (SR) / minimap2 (LR)
 │ Host removal │──→ subtract GRCh38 + PhiX reads
 └──────────────┘
     │
     ▼
 ┌─────────────────┐    sourmash sketch + SBT search
 │ Sketch screening│──→ 20k+ genomes → 50–200 candidates  (< 60 s)
 └─────────────────┘
     │
     ▼
 ┌──────────────────┐    minimap2 per candidate
 │ Targeted alignment│──→ multi-mapping read matrix
 └──────────────────┘
     │
     ▼
 ┌──────────────────┐    EM + bootstrap (n=100)
 │ EM abundance     │──→ θ estimates + 95% CI per organism
 └──────────────────┘
     │
     ▼
 ┌─────────────────┐    ABRicate + CARD
 │ AMR screening   │──→ resistance genes per organism
 └─────────────────┘
     │
     ▼
 ┌─────────────────┐    Grade A/B/C/X per organism
 │ Clinical report │──→ JSON + TSV + PDF report
 └─────────────────┘
```

---

## Installation

```bash
# 1. Clone
git clone https://github.com/alvin8-git/PathogenIQ.git
cd PathogenIQ

# 2. Create environment
conda create -n pathogeniq python=3.11
conda activate pathogeniq

# 3. Install bioinformatics tools
conda install -c bioconda bwa-mem2 minimap2 samtools fastp sourmash ncbi-genome-download

# 4. Install PathogenIQ
pip install -e .
```

**Verify installation:**
```bash
pathogeniq --help
```

---

## Database Setup

Run all three scripts concurrently in separate terminals (~90 min total on a 1 Gbps connection):

```bash
# Terminal 1 — ZymoBIOMICS validation DB (~5 min)
bash scripts/01_download_zymo_db.sh 16

# Terminal 2 — Full Tier-1 pathogen DB (~30–90 min, ~66 genomes)
python scripts/02_download_tier1_db.py --threads 16

# Terminal 3 — Human reference GRCh38 + BWA-MEM2 index (~90 min, needs ~64 GB RAM)
bash scripts/03_download_human_ref.sh 16
```

See [`scripts/README.md`](scripts/README.md) for full details, accession tables, and alternative pre-built database sources.

---

## Usage

### Run the full pipeline

```bash
export HUMAN_REF=$(pwd)/databases/human_ref/GRCh38_decoy.fa
export TIER1_DB=$(pwd)/databases/zymo_tier1/zymo_tier1.sbt.zip

pathogeniq run \
  --input sample.fq.gz \
  --output results/sample \
  --db "${TIER1_DB}" \
  --host-ref "${HUMAN_REF}" \
  --specimen blood \
  --read-type short \
  --threads 16
```

### Output files

```
results/sample/
├── report/
│   ├── pathogeniq_report.json   # Structured clinical report
│   ├── pathogeniq_report.tsv    # Tabular summary
│   └── pathogeniq_report.pdf    # Formatted clinical PDF
├── qc/                          # fastp HTML + JSON
├── host_removed.fq.gz           # Microbial reads only
└── alignments/                  # Per-organism PAF files
```

### Report fields

| Field | Description |
|---|---|
| `organism` | Species name |
| `accession` | Reference genome accession |
| `abundance_pct` | EM-estimated relative abundance |
| `ci_lower` / `ci_upper` | Bootstrap 95% confidence interval |
| `mapped_reads` | Reads assigned by EM |
| `grade` | Evidence grade: A (high), B (moderate), C (low), X (insufficient) |
| `contaminant_risk` | Flagged if organism is a known specimen-specific contaminant |
| `amr_genes` | Detected antimicrobial resistance genes linked to this organism |
| `specimen` | Specimen type used for grading thresholds |

---

## Validation

PathogenIQ ships with integration tests against the ZymoBIOMICS D6300 community standard (8 bacteria at 12% + 2 yeasts at 2% gDNA abundance):

```bash
export HUMAN_REF=$(pwd)/databases/human_ref/GRCh38_decoy.fa
export TIER1_DB=$(pwd)/databases/zymo_tier1/zymo_tier1.sbt.zip

pytest tests/integration/ -v -m integration
```

**Validation checks:**
- All 10 community members detected
- Measured abundance within 20% relative error of expected
- No Grade A false positives outside the ZymoBIOMICS panel

---

## Development

```bash
# Install dev dependencies
pip install -e ".[dev]"

# Run unit tests
pytest tests/ -v

# Lint
ruff check pathogeniq/ tests/

# Type check
mypy pathogeniq/
```

### Project structure

```
pathogeniq/
├── config.py        # ReadType, SpecimenType enums; PipelineConfig dataclass
├── qc.py            # fastp / Chopper wrapper
├── host_remove.py   # BWA-MEM2 / minimap2 host subtraction
├── sketch.py        # sourmash sketch + SBT search
├── align.py         # Targeted alignment (minimap2)
├── em.py            # EM abundance estimation + bootstrap CI
├── report.py        # Evidence grading + JSON/TSV output
├── contaminants.py  # Specimen-aware contaminant prior registry
├── amr.py           # ABRicate + CARD resistance gene detection
├── pdf_report.py    # ReportLab clinical PDF generator
└── cli.py           # Click CLI entry point

scripts/
├── 01_download_zymo_db.sh       # ZymoBIOMICS 10-genome DB
├── 02_download_tier1_db.py      # Full ~66-genome clinical DB
├── 03_download_human_ref.sh     # GRCh38 + BWA-MEM2 / minimap2 index
└── README.md                    # Database setup guide

tests/
├── integration/
│   └── test_zymo_validation.py  # ZymoBIOMICS D6300 end-to-end tests
└── test_*.py                    # Unit tests (all stages)
```

---

## Roadmap

- [x] Plan 2: Clinical interpretation engine (AMR overlay via ABRicate+CARD, PDF report)
- [ ] Plan 3: Nextflow orchestration (local Docker, SLURM, AWS Batch)
- [ ] Plan 4: Validation framework (ROC curves, benchmark suite)
- [ ] Web dashboard for clinical users
- [ ] FHIR-compatible report export

---

## Citation

If you use PathogenIQ in your research, please cite:

```
PathogenIQ: A sketch-first clinical metagenomics pipeline with EM abundance estimation.
[Manuscript in preparation]
```

---

## License

MIT License — see [LICENSE](LICENSE) for details.
