# PathogenIQ

[![CI](https://github.com/alvin8-git/PathogenIQ/actions/workflows/ci.yml/badge.svg)](https://github.com/alvin8-git/PathogenIQ/actions/workflows/ci.yml)
[![Python 3.11+](https://img.shields.io/badge/python-3.11+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Code style: ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://github.com/astral-sh/ruff)

A local clinical-metagenomics pipeline for detecting pathogens directly from sequencing reads — Illumina or Nanopore. PathogenIQ pairs a fast **sketch-first** targeted screen with a **statistically honest interpretation layer**: EM abundance with bootstrap confidence intervals, no-template-control (NTC) background subtraction, breadth-of-coverage gating, and a specimen-aware A/B/C/X evidence grade. On top of the targeted core it adds an **open-world** arm for **air / bioaerosol surveillance** that catches pathogens not in the reference database — including viruses and novel organisms.

The design goal is not "name every organism" — that is commodity. It is to turn a pile of reads into a **short, interpretable, controlled-confidence pathogen list** a clinician can act on, and to be explicit about what the evidence does and does not support.

---

## What it does

**Targeted core (always on)**

- **Sketch-first screening** — sourmash MinHash/SBT shortlists candidate genomes before any heavy alignment
- **EM abundance + bootstrap CI** — Expectation-Maximization resolves multi-mapping reads; 95% CIs quantify uncertainty
- **NTC background subtraction** — negative-binomial upper-tail test against a reagent/kitome background, the load-bearing control for low-biomass samples ([`background.py`](pathogeniq/background.py))
- **Flag-not-subtract for dual-use pathogens** — a contaminated control can never *erase* a real treatable pathogen (E. coli, Pseudomonas, Klebsiella, Salmonella, S. aureus…); it is flagged, not removed
- **Breadth-of-coverage gate** — reads clumped at one locus (PCR-dup / low-complexity / cross-map artifact) are graded out using the Lander–Waterman expected breadth ([`coverage.py`](pathogeniq/coverage.py))
- **Evidence grading (A/B/C/X)** — specimen-aware read floors, CI width, NTC tier cap (Grade A requires a same-run control), contaminant demotion, cross-mapping dedup
- **AMR + virulence** — ABRicate vs CARD (resistance) and VFDB (virulence factors)
- **Spike-in absolute quantification** — anchor reads to a known spike to report copies / copies-per-volume ([`quantify.py`](pathogeniq/quantify.py))
- **Specimen types** — blood, CSF, BAL, tissue, and **air** (bioaerosol), each with its own contaminant priors and read floors

**Open-world arm (opt-in, for air surveillance and novel pathogens)**

- **Novelty trigger** (`--novelty`) — classify reads against a broad Kraken2 DB; the unclassified "dark-matter" fraction flags content with no reference
- **Viral arm** (`--viral`) — assemble → geNomad (viral ID + ICTV taxonomy) → CheckV (completeness). Air pathogens are viral-heavy and the targeted/MAG arms are blind to them
- **Assembly / MAG recovery** (`--assemble`) — MEGAHIT → MetaBAT2 → CheckM → GTDB-Tk recovers reference-free genomes
- **Pathogenicity triage** — separates a novel *pathogen* from a novel *environmental* microbe via VFDB/CARD markers on the contigs + phylo-proximity to known pathogens
- **Open-world grading** — MAGs and viral contigs get the same A/B/C/X grade (from completeness + contamination + marker signal), capped at B (no same-run NTC for an assembled genome)

**Reports** — JSON + TSV always; PDF + HTML clinical reports.

---

## Pipeline

```
FASTQ (short or long)
   │
   ▼  QC ........................ fastp (short) / Chopper (long)
   ▼  Host removal .............. BWA-MEM2 (short) / minimap2 (long) vs GRCh38
   ▼  PhiX removal .............. minimap2 vs the Illumina PhiX spike-in (NC_001422)
   ▼  Novelty screen [--novelty]  Kraken2 broad DB → unclassified fraction
   ▼  Sketch screen ............ sourmash MinHash + SBT → candidate shortlist
   ▼  Targeted alignment ....... minimap2 per candidate → read×organism matrix
   │                             + breadth-of-coverage per organism
   ▼  EM abundance ............. EM + bootstrap (n=100) → θ + 95% CI
   ▼  AMR + virulence .......... ABRicate vs CARD + VFDB
   ▼  NTC background ........... NB upper-tail test; flag-not-subtract dual-use taxa
   ▼  Spike quant [--spike-*] .. absolute copies / copies-per-volume
   ▼  Grading ................. A/B/C/X: read floor · CI · NTC tier cap
   │                            · breadth gate · cross-map · contaminant prior
   │
   ├─ Open-world arm [--assemble] : MEGAHIT → MetaBAT2 → CheckM → GTDB-Tk
   │                                → pathogenicity triage → open-world grade
   ├─ Open-world arm [--viral]    : MEGAHIT → geNomad → CheckV → open-world grade
   │
   ▼  Report .................. JSON + TSV + PDF + HTML
```

Every external-tool stage is **non-blocking**: a missing tool or database degrades gracefully (the stage is skipped) rather than failing the run.

---

## Installation

```bash
git clone https://github.com/alvin8-git/PathogenIQ.git
cd PathogenIQ

# Core environment (pins the always-on tools + Python deps)
conda env create -f environment.yml
conda activate pathogeniq
pip install -e ".[dev]"

pathogeniq --help
```

The **open-world arm tools** (geNomad, CheckV, ABRicate) require `numpy<2` and conflict with the core env's NumPy-2 stack, so they live in **isolated conda envs** and are symlinked/wrapped onto `PATH`. One script sets them up:

```bash
bash scripts/13_setup_viral_env.sh     # geNomad + CheckV env + DBs (~3 GB), CLIs linked in
```

(ABRicate is set up the same way — its own env + a `PATH` wrapper — because it is a Perl/BLAST tool. See [`environment.yml`](environment.yml) for the exact commands.)

---

## Database setup

```bash
bash   scripts/01_download_zymo_db.sh 16          # ZymoBIOMICS validation DB (~5 min)
python scripts/02_download_tier1_db.py --threads 16  # Tier-1 clinical DB (~110 pathogen genomes)
bash   scripts/03_download_human_ref.sh 16          # GRCh38 + index (~90 min, ~64 GB RAM)
```

Optional databases for the open-world arm: a broad **Kraken2** DB (novelty), **geNomad** + **CheckV** DBs (viral, fetched by `scripts/13`). The heavy **GTDB-Tk** DB (~110 GB, for MAG taxonomy) is deliberately **not** part of the core install — it is a triggered Tier-2 extra. See [`scripts/README.md`](scripts/README.md).

---

## Usage

### Clinical sample (blood)

```bash
pathogeniq run \
  --input sample.fq.gz --output results/sample \
  --db databases/tier1/tier1_pathogens.sbt.zip \
  --host-ref /path/to/GRCh38/hg38.fa \
  --specimen blood --read-type short --threads 16
```

With a batch-matched NTC (the only way to reach Grade A):

```bash
pathogeniq run ... --ntc ntc_blank.fq.gz
```

### Air / bioaerosol surveillance (full open-world)

```bash
pathogeniq run \
  --input air_filter.fq.gz --output results/air \
  --db databases/tier1/tier1_pathogens.sbt.zip \
  --host-ref /path/to/GRCh38/hg38.fa \
  --specimen air --read-type short \
  --background air_ntc.tsv \
  --novelty --viral --threads 16
```

### Key options

| Option | Effect |
|---|---|
| `--specimen blood\|csf\|bal\|tissue\|air` | Selects read floors + contaminant priors |
| `--read-type short\|long` | fastp+BWA-MEM2 vs Chopper+minimap2 |
| `--ntc FASTQ` | Batch-matched NTC → Tier-1 background (enables Grade A) |
| `--background TSV` | Pooled/precomputed background → Tier-2 (capped at B) |
| `--no-background` | Disable background correction (Tier-3, uncorrected) |
| `--novelty` | Open-world Kraken2 unclassified-fraction trigger |
| `--viral` | Viral arm (geNomad + CheckV) |
| `--assemble` | Assembly/MAG arm + pathogenicity triage |
| `--spike-taxon / --spike-copies / --sample-volume` | Absolute quantification |
| `--no-pdf` | Skip the PDF (JSON/TSV/HTML still written) |

### Report

`results/<sample>/report/` holds `pathogeniq_report.{json,tsv,pdf,html}`. The JSON has a `findings` array (targeted hits with grade, `breadth_ratio`, NTC tier, AMR/virulence, absolute copies) plus optional `novelty`, `viral`, `mags`, and `pathogenicity` blocks — every hit, targeted or open-world, carries an A/B/C/X grade.

---

## Validation

PathogenIQ's design decisions are driven by real validation runs, not intuition. The full record is in [`Documentation.md`](Documentation.md) §6 and the `docs/` notes.

| Validation | Data | Result | Influenced |
|---|---|---|---|
| Mock standard | ZymoBIOMICS D6300 (8 bacteria + 2 yeast) | all members detected, abundance within tolerance | core grading thresholds |
| Held-out grading | CAMI II (mouse-gut/marine/strain) + HMP (skin/airway/oral) | precision **~45×** over raw Kraken2 at preserved recall, 9 unseen communities | grading wedge generalises |
| Dispersion prior | spike-free kitome blanks, leave-one-out | FPR 30–65× α at **every** dispersion → coverage, not the prior, is the limit | Tier-2 → Grade B cap |
| Air concordance | Jeilu et al. 2025 aircraft filters (PRJNA1228129) | Zymo spike recall 6/10; air NTC over-subtracted real E. coli/Pseudomonas | **flag-not-subtract** for dual-use pathogens |
| Viral arm | in-silico spike (T4, Lambda, SARS-CoV-2 → real air background) | **100% recall (3/3)**, correct ICTV lineage | confirmed the viral arm end-to-end; caught a silent `run_megahit` bug |

Run the mock-standard integration tests:

```bash
HUMAN_REF=... TIER1_DB=... pytest tests/integration/ -v -m integration
```

---

## Development

```bash
pip install -e ".[dev]"
pytest tests/ -v            # ~180 unit tests, no external tools needed (subprocess mocked)
ruff check pathogeniq/ tests/
mypy pathogeniq/
```

### Module map

| Module | Role |
|---|---|
| `config.py` | `ReadType` / `SpecimenType` enums; `PipelineConfig` |
| `qc.py`, `host_remove.py` | QC (fastp/Chopper) + host & PhiX removal |
| `sketch.py`, `align.py`, `em.py` | sketch screen → targeted alignment (+coverage) → EM |
| `coverage.py` | breadth-of-coverage (Lander–Waterman) |
| `background.py` | NB NTC background + dual-use flag-not-subtract |
| `contaminants.py`, `crossmap.py` | specimen contaminant priors; cross-mapping dedup |
| `report.py` | grading (targeted + open-world), JSON/TSV |
| `amr.py`, `quantify.py` | AMR/virulence; spike-in absolute quantification |
| `novelty.py` | open-world novelty trigger (Kraken2) |
| `assembly.py` | MEGAHIT→MetaBAT2→CheckM→GTDB-Tk MAG arm |
| `viral.py` | geNomad + CheckV viral arm |
| `pathogenicity.py` | pathogen-vs-environmental triage of MAGs |
| `pdf_report.py`, `html_report.py` | clinical report renderers |
| `benchmark.py` | Kraken2 grading adapter for the held-out benchmark |
| `cli.py` | Click entry point orchestrating all stages |

---

## Status & roadmap

Implemented: targeted core, NTC tiered grading, cross-mapping dedup, breadth gate, AIR specimen + air kitome priors, AMR + VFDB virulence, spike-in absolute quantification, and the full open-world arc (novelty → viral → MAG → pathogenicity triage → open-world grading). The viral arm is validated end-to-end (100% in-silico recall); the bacterial-MAG path is code-complete and unit-tested but needs the GTDB-Tk DB installed to run on real data.

Next: calibrate provisional cutoffs (breadth ratio, open-world completeness bands, air read floors) on labeled data; prospective batch-matched-NTC studies for clinical Grade-A; optional GTDB arm install.

---

## Citation & license

PathogenIQ — a sketch-first clinical-metagenomics pipeline with NTC-aware grading and an open-world air-surveillance arm. *Manuscript in preparation.* MIT License — see [LICENSE](LICENSE).

> **Research use only.** Not a validated in-vitro diagnostic. Findings must be corroborated with orthogonal testing and clinical context.
