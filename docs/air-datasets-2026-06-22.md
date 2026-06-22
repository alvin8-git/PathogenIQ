# Air / bioaerosol shotgun datasets (Plan 5 prep) — 2026-06-22

Candidate **ungated, downloadable** environmental-air shotgun metagenomics datasets
for adapting the pipeline to bioaerosol surveillance. All verified live against the
ENA portal API; air is even lower-biomass + more kitome-dominated than clinical
specimens, so datasets that ship **negative controls** and/or **spike-ins** are
prioritised. Download via `scripts/04_download_validation_data.py` (ENA filereport).

Filereport URL (swap `ACCESSION`):
`https://www.ebi.ac.uk/ena/portal/api/filereport?accession=ACCESSION&result=read_run&fields=run_accession,fastq_ftp,read_count,library_strategy,library_selection,sample_title&format=tsv`

## Tier 1 — true shotgun air (recommended)

| # | Accession | Study / context | Runs | Platform | Size | Controls / spike-ins |
|---|-----------|-----------------|------|----------|------|----------------------|
| 1 | **PRJNA1228129** | Aircraft filters + face masks (cabin/transit air) | 55 WGS Illumina | ~5.6M reads/run | ~20 GB | **Both:** 6 neg controls (filter `Control`), positive spike-ins (filter `Spike`) — ideal for kitome-subtraction + grading validation |
| 2 | **PRJNA561080** | MetaSUB subway & urban air (gCSD17 air arm) | 294 WGS Illumina | ~4.6M reads/run | ~109 GB | none surfaced (titles are library IDs) |
| 3 | **PRJEB59175** | Air microbiome via long-read bioaerosols | 285 WGS/RANDOM | PacBio Sequel IIe, ~136k reads/run | ~106 GB | **4 explicit blanks** ERR10900207–210 (`blank1`..`blank4`, 295–1200 reads). ⚠️ long-read branch |
| 4 | **PRJNA791118** | Lower-troposphere air, vertical stratification | 234 WGS Illumina | ~2.7M reads/run | ~94 GB | none surfaced |
| 5 | **PRJNA638794** | Air-sampling **method optimization** (low-biomass) | 98 WGS (+13 amplicon) | Illumina | ~48 GB | filter `library_strategy=WGS` |

**Best starting point: PRJNA1228129** — small, clean, and the only one with *both* negative controls and positive spike-ins, so it exercises background subtraction (controls) and detection sensitivity (spikes) in one dataset. Total ~42 GB (55 runs, 309M reads); the controls + spike-ins alone are only ~8 GB.

> **Reference:** Jeilu O, Sumner JT, Moghadam AA, Thompson KN, Huttenhower C, et al. (2025). *Metagenomic profiling of airborne microbial communities from aircraft filters and face masks.* **Microbiome.** doi:10.1186/s40168-025-02276-7 (PMID 41340070; medRxiv preprint doi:10.1101/2025.02.26.25322977). SRA study SRP566804.

## Supplementary verified Illumina shotgun air (taxid 655179)

PRJNA858396 "Global air" (71 runs, ~357 GB, deep ~40M/run) · PRJNA1045528 E. Mediterranean atmosphere (52) · PRJNA1049581 indoor day-night (130) · PRJNA486429 megacity PM (106, ~32M/run) · PRJNA643946 W. Siberia (61) · PRJEB67959 confined-environment exposome (52) · PRJNA436039 DVE time-series (783, shallow).

## Amplicon-only (NOT shotgun — lower priority)

PRJEB96751, PRJEB38899 (Tara aerobiome), PRJDB7249/PRJDB8920, PRJEB29770/PRJEB38776 (pig-farm), PRJNA1054501, PRJNA1079228, PRJNA1204957, PRJNA1145031, PRJNA1164146, PRJNA1424403.

## Excluded

- **PRJEB109029** (Aspergillus hospital bioaerosols) — cultured single-isolate WGS, not bulk metagenomes.
- **PRJNA732392** (MetaSUB 2016-17) — surface swabs, no air (air arm = PRJNA561080).
- CAMI II HMP "airways" — human body site, not environmental air.

## How this maps to the pipeline (Plan 5)

- **Controls → Tier-1 NTC path.** PRJNA1228129's negative controls are batch-matched air blanks — the right input for `build_background` to model the air kitome, and to test whether a real air NTC (Tier 1) lifts grades vs the foreign pooled prior.
- **Spike-ins → detection validation.** PRJNA1228129 spikes give known positives to score precision/recall, analogous to the Zymo/CAMI held-out benchmark but in the air domain.
- **Breadth-of-coverage gate** (`coverage.py`) is the air-specific addition: separate genome-wide real signal from clumped contamination artifacts at the very low depths typical of air.
