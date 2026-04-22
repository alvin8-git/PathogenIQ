# PathogenIQ — Design Specification

**Date:** 2026-04-22  
**Status:** Approved  
**Author:** Brainstorming session with Alvin  

---

## 1. Overview

PathogenIQ is a lightweight clinical metagenomics analysis platform that ingests NGS sequencing data (Illumina short reads and Oxford Nanopore long reads), performs QC, human host removal, and rapid pathogen identification against a tiered curated database, then generates a structured clinical report with evidence grading and antimicrobial resistance (AMR) profiling.

### Core Novelty

Three capabilities differentiate PathogenIQ from existing tools (CZID, Kraken2+Bracken):

1. **Speed** — sketch-first screening reduces alignment target from 20,000+ to 50–200 genomes, achieving <1 hour turnaround
2. **Confidence scoring** — Bayesian confidence intervals on every abundance estimate; evidence grades on every pathogen call
3. **Clinical intelligence** — specimen-aware thresholds, AMR genes linked to source organism, negative result interpretation, and human-readable evidence grading

---

## 2. Scope

### In Scope
- Short reads: Illumina paired-end 150–300 bp
- Long reads: Oxford Nanopore (Guppy/Dorado basecalled FASTQ)
- Specimen types: blood/plasma cfDNA, CSF, BAL/sputum, tissue biopsy
- Pathogen classes: bacteria, viruses, fungi, parasites
- AMR gene detection linked to source organism
- Novel/divergent pathogen flagging via de novo assembly
- Hybrid deployment: local server + cloud burst (AWS Batch / GCP)
- HPC support: SLURM via Nextflow

### Out of Scope (v1.0)
- 16S/ITS amplicon sequencing (whole-genome shotgun only)
- PacBio HiFi (added in v2.0)
- Adaptive Nanopore sampling (Nanopore streaming is read-only in v1.0)
- EHR integration

---

## 3. Pipeline Architecture

### High-Level Flow

```
FASTQ (Illumina SR / Nanopore LR)
         │
    ┌────▼────┐
    │  Stage 1 │  QC & Adapter Trimming
    │          │  fastp (SR) | Chopper + Porechop (LR)
    │          │  → quality report, microbial fraction metric
    └────┬────┘
         │
    ┌────▼────┐
    │  Stage 2 │  Host Removal
    │          │  BWA-MEM2 (SR) | minimap2 -x map-ont (LR)
    │          │  Reference: GRCh38 + decoys + lab contaminants
    │          │  → non-human reads only proceed
    └────┬────┘
         │
    ┌────▼────┐
    │  Stage 3 │  Sketch-First Screening  ← NOVEL CORE
    │          │  sourmash / Mash MinHash sketches (k=31, scaled=1000)
    │          │  against Tier-1 pathogen index (~2,000 genomes)
    │          │  → candidate shortlist in <60 seconds
    └────┬────┘
         │
    ┌────▼────┐
    │  Stage 4 │  Targeted Alignment + EM Abundance
    │          │  minimap2 (LR) | BWA-MEM2 (SR) → candidate genomes only
    │          │  EM algorithm resolves multi-mappers
    │          │  → abundance + 95% Bayesian CI per organism
    │          │  → de novo assembly trigger for unclassified reads
    └────┬────┘
         │
    ┌────▼────┐
    │  Stage 5 │  Clinical Interpretation Engine
    │          │  Specimen-aware thresholds, AMR overlay (CARD + ResFinder)
    │          │  Evidence grading, structured PDF + JSON report
    └──────────┘
```

### Read Type Branching

Stages 1–2 run identically for SR and LR. Stage 3 sketching is read-type agnostic. Stage 4 alignment branches (BWA-MEM2 vs minimap2 preset) then converges at the EM step. Stage 5 is unified.

---

## 4. Component Specifications

### Stage 1 — QC & Adapter Trimming

**Short reads:**
- Tool: fastp
- Parameters: Q20 quality filter, min length 50bp, poly-G tail trim, auto adapter detection, deduplication
- Output: filtered FASTQ, MultiQC-compatible JSON

**Long reads:**
- Tools: Chopper (quality/length filter, min Q8, min 200bp) + Porechop (adapter removal) + NanoStat (metrics)
- Output: filtered FASTQ, NanoStat report

**Key QC metrics reported:**
- `total_reads` — raw input
- `passing_reads` — post-QC
- `microbial_fraction` — % non-human (computed after Stage 2, reported with Stage 1 output)
- `effective_microbial_depth` — passing_reads × microbial_fraction

### Stage 2 — Host Removal

- **Short reads:** BWA-MEM2 against GRCh38 + alt contigs + decoy sequences + PhiX + common cloning vectors
- **Long reads:** minimap2 with `-x map-ont` preset, same reference
- Unmapped reads proceed; mapped reads discarded
- **Human homology masking:** pathogen database regions with >90% nucleotide identity to human are masked to prevent false negatives at host-pathogen sequence boundaries (e.g. *Neisseria* surface proteins)
- Output: non-human FASTQ + host alignment stats

### Stage 3 — Sketch-First Screening

- All Tier-1 pathogen genomes have pre-computed MinHash sketches (sourmash, k=31, scaled=1000)
- Incoming non-human reads are sketched in streaming batches (60s windows for Nanopore; full dataset for Illumina)
- Containment score computed per organism
- Organisms above threshold (default: containment ≥0.01, configurable) enter candidate shortlist
- Typical shortlist: 50–200 genomes from 2,000 Tier-1 genomes
- **Speed:** <60 seconds on 16-core machine for 10M reads
- **RAM:** ~8 GB for full Tier-1 sketch index (vs ~200 GB for Kraken2 full NT)

### Stage 4 — Targeted Alignment + EM Abundance

**Alignment:**
- minimap2 (LR) or BWA-MEM2 (SR) against candidate genomes only
- Multi-mapping reads retained with all alignment records

**EM Abundance Estimation:**
- Expectation-Maximization resolves multi-mapping reads (analogous to Salmon/kallisto for RNA-seq)
- Outputs per-organism: estimated read counts, genome-length-normalized abundance, 95% Bayesian CI, posterior probability of true presence

**Novel Pathogen Trigger:**
- Reads unmapped after Tier-1 candidate alignment → secondary pass against Tier-2 extended database (optional, +20–25 min)
- Reads remaining unmapped → de novo assembly (metaSPAdes for SR, Flye for LR)
- Assembled contigs → BLASTn vs NCBI NT (cloud burst)
- Contigs with <90% identity to any known sequence → Grade X flag: "Novel/Divergent Pathogen — Manual Review Required"

### Stage 5 — Clinical Interpretation Engine

**Specimen-aware detection thresholds:**

| Specimen Type | Min Reads | Min Confidence | Rationale |
|---|---|---|---|
| Blood / plasma cfDNA | 3 | 80% | Sterile site — any detection significant |
| CSF | 2 | 70% | Sterile site — very low biomass expected |
| BAL / sputum | Relative abundance | 60% | Normal flora present — relative deviation used |
| Tissue biopsy | 10 | 85% | Sterile site — contaminants common in processing |

**Evidence Grading:**

- **Grade A — Strong:** High read depth, narrow CI, Tier-1 organism, consistent with specimen type, AMR genes linkable
- **Grade B — Moderate:** Adequate depth, wider CI, or organism at specimen-type boundary
- **Grade C — Weak:** Low read depth, wide CI, or organism unexpected for specimen
- **Grade X — Unclassified:** Novel/divergent sequence, <90% identity to known sequence; manual review required

**AMR Overlay:**
- All non-human reads aligned against CARD + ResFinder via ABRicate (parallel to taxonomic pipeline)
- AMR genes linked to source organism via contig-of-origin when possible
- Output: per-organism resistance profile with gene name, mechanism, and drug class

---

## 5. Database Design

### Tier 1 — Clinical Core Panel
- ~2,000 genomes
- Sources: CZID pathogen list, WHO Priority Pathogens, CDC Select Agents, specimen-specific panels (CSF neurotropic viruses, blood gram-negatives/fungi, etc.)
- Versioned releases (SemVer), quarterly manual review
- Runs in <1 hour fast path

### Tier 2 — Extended Research Database
- ~50,000 genomes from NCBI RefSeq (bacterial complete + WGS assemblies, viral RefSeq, fungal RefSeq)
- Secondary pass on unclassified reads (+20–25 min)
- Optional — toggled by run mode (`--mode extended`)

### Tier 3 — Novel Pathogen (De Novo)
- De novo assembled contigs BLASTn'd against full NCBI NT
- Cloud burst only (AWS Batch)
- Async — report issued with "assembly in progress" flag, updated when complete

### Database Update Workflow

```
NCBI RefSeq nightly mirror (automated)
        │
   Validation pipeline:
   - CheckM quality filter (completeness >90%, contamination <5%)
   - Human contamination screen (remove sequences >90% human identity)
   - Taxonomy reconciliation (GTDB-Tk for bacteria/archaea)
        │
   Staging database → automated ZymoBIOMICS regression tests
        │
   Pass (10/10 organisms detected, abundance error <20%, FP rate <5%)?
   → Promote to Tier 1/2 with version bump + changelog
   Fail?
   → Hold for manual review, alert maintainers
```

Local deployments sync from signed, versioned release bundles — no live NCBI pulls in production. Eliminates false positives from contaminated new NCBI entries (a known Kraken2 failure mode).

---

## 6. Clinical Report Format

### Report Structure

```
PathogenIQ Report
Sample ID | Specimen type | Run date | DB version | Run time
─────────────────────────────────────────────────
SPECIMEN QUALITY
  Microbial fraction: X% (adequate / borderline / inadequate)
  Effective microbial depth: N reads
  Human depletion efficiency: X%
  QC status: PASS / FAIL (with reason)
─────────────────────────────────────────────────
TIER 1 — ACTIONABLE FINDINGS
  For each Grade A/B organism:
  ● Organism name          [Grade A/B]
    Abundance: X% of microbial reads
    95% CI: X–Y% | Read count: N
    AMR profile: [gene list + drug classes affected]
    Clinical note: [specimen-aware interpretive text]
─────────────────────────────────────────────────
TIER 2 — ORGANISMS OF UNCERTAIN SIGNIFICANCE
  For each Grade C organism:
  ● Organism name          [Grade C]
    Read count: N | 95% CI: X–Y%
    Note: [contaminant / colonizer / clinical correlation required]
─────────────────────────────────────────────────
TIER 3 — NOVEL / UNCLASSIFIED
  ● Contig description, length, best BLAST hit + % identity
    Manual review required
─────────────────────────────────────────────────
NEGATIVE RESULT INTERPRETATION
  Pathogens not detected: [bacterial / viral / fungal / parasitic]
  Sensitivity statement: [effective depth vs. specimen-type minimum]
  Pathogens specifically screened and not detected: [named list]
─────────────────────────────────────────────────
METHODS SUMMARY
  Pipeline version | DB version | Key parameters
```

### Output Formats
- **PDF** — formatted clinical report for medical record
- **JSON** — structured data for downstream systems / EHR integration
- **TSV** — full organism × abundance table for research use

---

## 7. Deployment Architecture

### Local Node (Minimum Specification)
- CPU: 32 cores (AMD EPYC or Intel Xeon)
- RAM: 64 GB
- Storage: 2 TB NVMe SSD
- GPU: Optional (NVIDIA A100/H100 for mm2-fast acceleration)
- OS: Ubuntu 22.04 LTS or RHEL 8+

### Containerization
Each stage ships as a signed Docker image (Singularity-compatible):

```
pathogeniq/qc:1.0           # fastp, Chopper, NanoStat, Porechop
pathogeniq/host-remove:1.0  # BWA-MEM2, minimap2, samtools
pathogeniq/sketch:1.0       # sourmash, Mash
pathogeniq/align-em:1.0     # minimap2, BWA-MEM2, custom EM, metaSPAdes, Flye
pathogeniq/interpret:1.0    # clinical engine, report generator
pathogeniq/db:2.1.3         # tiered pathogen database (versioned separately)
```

### Orchestration
- **Tool:** Nextflow
- Single `nextflow.config` switches between local Docker, SLURM (HPC), and AWS Batch (cloud)
- Cloud burst triggers: Tier 2/3 database pass, de novo assembly, >3 concurrent samples

### Data Sovereignty
- Raw FASTQ never transmitted to cloud
- Only non-human reads (post-Stage 2) transmitted for cloud burst jobs
- Encrypted in transit (TLS 1.3) and at rest (AES-256)

### Speed Budget

| Stage | Illumina SR (10M reads) | Nanopore LR (1M reads) |
|---|---|---|
| QC | ~2 min | ~3 min |
| Host removal | ~8 min | ~10 min |
| Sketch screening | ~1 min | ~1 min |
| Targeted alignment + EM | ~15 min | ~20 min |
| Clinical report | ~2 min | ~2 min |
| **Total (Tier 1)** | **~28 min** | **~36 min** |
| Tier 2 pass (optional) | +20–25 min | +20–25 min |
| De novo assembly (cloud, async) | +30–60 min | +30–60 min |

### Nanopore Streaming Mode (v1.0)
- Watch daemon monitors MinKNOW output directory
- 60-second batching of new reads through Stages 3–5
- Progressive confidence updates; live dashboard refresh
- Early-stop trigger: Grade A pathogen at >95% posterior probability AND effective depth exceeds specimen-type minimum

---

## 8. Validation Strategy

### ZymoBIOMICS Benchmark Dataset

Available at `/data/alvin/Metagenomics/ZymoStandardsFromMao/`:

| File | Description | Use |
|---|---|---|
| `V350234554_L01_UDB-32.Zymo_Std.fq.gz` | SR replicate 1 (1.15 GB) | Primary benchmark |
| `V350234554_L01_UDB-40.Zymo_Std-r2.fq.gz` | SR replicate 2 (987 MB) | Reproducibility |
| `V350234554_L01_UDB-48.Zymo_Std-r3.fq.gz` | SR replicate 3 (1.31 GB) | Reproducibility |
| `Zymobac_D4_3ng_1.fq.gz` | Low-input 3ng (818 MB) | Low-biomass sensitivity |
| `Zymobac_D4_6ng_2_1.fq.gz` | Higher-input 6ng (806 MB) | Input concentration effect |
| `ZymoM_10_1_1.fq.gz` | Mock community (779 MB) | Cross-validation |

**Ground truth:** 8 bacteria at 12% each (*P. aeruginosa*, *E. coli*, *S. enterica*, *L. fermentum*, *E. faecalis*, *S. aureus*, *L. monocytogenes*, *B. subtilis*) + 2 yeasts at 2% each (*S. cerevisiae*, *C. neoformans*).

### Validation Metrics

**Detection accuracy:**
- Sensitivity: organisms detected / 10 expected (target: 10/10)
- False positive rate: organisms detected not in ZymoBIOMICS panel (target: <5%)
- Species-level resolution: % of detections at species level (target: 100%)

**Abundance accuracy:**
- Pearson/Spearman correlation between measured and expected abundances
- Mean absolute error per organism (target: <20% relative error)

**Confidence interval calibration:**
- Coverage: do 95% CIs capture true abundance 95% of the time across replicates? (target: yes)
- Grade C call rate on low-input samples (3ng): expect higher uncertainty, lower grade

**Speed benchmarks:**
- Wall-clock time per stage on reference hardware
- Peak RAM and CPU utilization

### Comparative Benchmarks

Run identical FASTQ files through PathogenIQ, Kraken2+Bracken, and CZID. Compare:
- Detection accuracy
- Species-level resolution
- False positive rates
- Speed (wall-clock)
- Report completeness (CIs, contaminant flagging, negative interpretation)

### Automated Regression Suite

Every database version bump triggers:
1. ZymoBIOMICS benchmark (must pass: 10/10 detected, abundance error <20%, FP <5%)
2. Speed regression (must meet stage time budgets on reference hardware)
3. Report format validation (schema check on JSON output)

### Clinical Validation (Prospective)

Shadow-mode deployment alongside conventional culture in partner clinical laboratories. Compare PathogenIQ Grade A calls to final clinical microbiological diagnosis.

Targets:
- Sensitivity >90% for Grade A calls
- Specificity >95% for Grade A calls
- Time-to-result: median <45 minutes from FASTQ available

---

## 9. Blind Spots Addressed

| Gap in Current Tools | Root Cause | PathogenIQ Solution |
|---|---|---|
| No uncertainty quantification | Read counts only | Bayesian CI + posterior probability per organism |
| Contaminant blindness | No specimen model | Specimen-type contaminant priors and thresholds |
| AMR unlinked to organism | Separate analysis step | Contig-linked AMR overlay, per-organism resistance profile |
| Genus-level collapse | Kraken2 LCA algorithm | Targeted alignment resolves ambiguous reads with full genome context |
| Silent negatives | Tools report only positives | Explicit sensitivity statement + named pathogen negative screen |
| High RAM requirements | Full k-mer index in memory | Sketch index ~8 GB vs ~200 GB for Kraken2 full NT |
| Cloud-only data sovereignty | Architecture assumption | Patient data stays local; only non-human reads cloud-burst |
| Database contamination causes FPs | Live NCBI pulls | Signed versioned bundles, validated before release |
| No novel pathogen detection | Fixed database | De novo assembly + BLASTn flagging for divergent sequences |
| No clinical language | Bioinformatics output | Evidence grades + specimen-aware interpretive text |

---

## 10. Tool Summary

| Stage | Short Read Tool | Long Read Tool | Notes |
|---|---|---|---|
| QC | fastp | Chopper + Porechop + NanoStat | |
| Host removal | BWA-MEM2 | minimap2 -x map-ont | GRCh38 + decoys |
| Sketch screening | sourmash | sourmash | k=31, scaled=1000 |
| Alignment | BWA-MEM2 | minimap2 | Candidate genomes only |
| Abundance EM | Custom Python/Rust | Same | EM algorithm |
| De novo assembly | metaSPAdes | Flye | Cloud burst |
| AMR detection | ABRicate (CARD + ResFinder) | Same | Parallel to taxonomy |
| Orchestration | Nextflow | Same | SLURM + AWS Batch |
| Containerization | Docker / Singularity | Same | |

---

*Spec written from brainstorming session. Ready for implementation planning.*
