# PathogenIQ Documentation

## 1. Introduction

PathogenIQ is a clinical metagenomics pipeline designed to identify pathogens directly from patient samples using DNA sequencing. Unlike traditional culture-based diagnostics — which can take days and fail for fastidious or slow-growing organisms — PathogenIQ analyses the full genetic content of a clinical specimen and reports which microorganisms are present, how abundant they are, and what antimicrobial resistance genes they carry.

The name reflects two core design goals:

- **Pathogen**: clinical focus on detecting disease-causing organisms with sufficient confidence to inform treatment decisions.
- **IQ**: intelligent, statistically principled quantification — abundance estimates with confidence intervals, evidence grading, and specimen-aware interpretation rather than simple presence/absence calls.

### What is Metagenomics?

Metagenomics is the study of genetic material recovered directly from environmental or clinical samples without prior isolation or culture. In a clinical context, a single microbiome sample (blood, CSF, bronchoalveolar lavage, tissue biopsy) may contain DNA from hundreds of organisms: the human host, commensal bacteria, environmental contaminants, and — critically — the pathogen causing disease.

Metagenomic sequencing ("shotgun" sequencing) generates millions of short DNA fragments (reads) from all DNA in the sample. The computational challenge is to:

1. Remove host DNA (typically >90% of reads in blood, >99% in tissue)
2. Identify which microbial genomes the remaining reads originate from
3. Estimate the relative abundance of each organism
4. Distinguish true pathogens from commensals, contaminants, and sequencing noise
5. Detect antimicrobial resistance genes that may guide therapy

### Why PathogenIQ?

Existing tools fall into two categories, each with limitations:

- **Cloud platforms (e.g., CZ ID)** offer comprehensive analysis but require uploading patient data to third-party servers, have turnaround times of 4–8 hours, and lack statistical confidence measures.
- **Local k-mer tools (e.g., Kraken2 + Bracken)** are fast but provide no uncertainty quantification, no clinical context (a skin contaminant in blood is treated the same as a true pathogen), and no resistance gene detection.

PathogenIQ bridges this gap: it runs entirely locally, completes in under one hour, and produces clinically interpretable reports with graded evidence and specimen-aware contamination suppression.

---

## 2. The Clinical Problem

### Diagnostic Uncertainty in Sepsis and Meningitis

Bloodstream infections (sepsis) and central nervous system infections (meningitis, encephalitis) are medical emergencies where every hour of delay in appropriate antimicrobial therapy increases mortality. Yet traditional diagnostics face severe limitations:

| Challenge | Traditional Culture | Metagenomics |
|-----------|---------------------|--------------|
| Time to result | 24–72 hours | 4–24 hours |
| Sensitivity for fastidious organisms | Poor | High |
| Polymicrobial infections | Often missed | Detected |
| Antimicrobial resistance | Requires additional testing | Detectable concurrently |
| Prior antibiotic exposure | Culture frequently negative | DNA still detectable |

The core scientific problem is **discriminating signal from noise**. A blood sample may contain 10,000 microbial reads out of 10 million total reads. Of those 10,000, perhaps 8,000 come from skin contaminants introduced during venipuncture. Only 2,000 may represent the true bloodstream pathogen. Distinguishing these without prior knowledge is the central challenge clinical metagenomics must solve.

### The Contaminant Problem

Every specimen type has characteristic contaminants:

- **Blood**: skin flora (*Cutibacterium acnes*, *Staphylococcus epidermidis*) introduced during collection
- **CSF**: similar skin flora, plus *Corynebacterium* species from collection kit surfaces
- **BAL (bronchoalveolar lavage)**: oral flora (*Streptococcus salivarius*, *Prevotella* spp.) aspirated past the vocal cords during bronchoscopy
- **Tissue**: fewer systematic contaminants, but environmental organisms possible

A naive pipeline that reports every detected organism will generate false positives that mislead clinicians. PathogenIQ addresses this with a **specimen-aware contaminant prior registry** and a **graded evidence system** that prevents contaminants from receiving high-confidence grades.

### The Multi-Mapping Problem

When a microbial read aligns equally well to multiple closely related organisms (e.g., different *Staphylococcus* species, or multiple strains of *E. coli*), standard read-counting methods assign it arbitrarily to one. This leads to:

- Unstable abundance estimates (depends on which genome "won" the tie)
- Inflated confidence (the read was counted, but its origin is genuinely ambiguous)
- False polymicrobial calls (two related species both get "credit" for the same reads)

PathogenIQ addresses this with **Expectation-Maximization (EM) abundance estimation**, which probabilistically assigns multi-mapping reads based on the overall likelihood of each organism being present.

---

## 3. Scientific Architecture

PathogenIQ is built as a **sketch-first** pipeline: rather than aligning every read against every reference genome (prohibitively expensive), it first uses a fast probabilistic screening step to identify a shortlist of candidate organisms, then performs rigorous alignment and statistical abundance estimation only on those candidates.

### Stage-by-Stage Scientific Rationale

#### Stage 1: Quality Control (QC)

**Tools**: fastp (short reads) / Chopper (long reads)

**Why it matters**: Sequencing adapters, low-quality bases, and short fragments do not align reliably to reference genomes and add noise. Removing them improves alignment specificity and reduces false positives.

**Science**: fastp uses a sliding-window quality trimming algorithm that removes bases with Phred quality scores below a threshold, then filters reads shorter than a minimum length. Chopper applies similar logic for Nanopore reads, using mean read quality scores.

#### Stage 2: Host Removal

**Tools**: BWA-MEM2 (short reads) / minimap2 (long reads) + samtools

**Why it matters**: In clinical samples, >90% of reads typically originate from the human host. Aligning all reads against microbial databases would waste enormous compute and create alignment noise from low-complexity or repetitive human sequences.

**Science**: BWA-MEM2 is a Burrows-Wheeler transform (BWT) based aligner optimized for short reads. It indexes the reference genome into a suffix array that enables sub-linear time search. minimap2 uses a minimizer-based indexing scheme optimized for long, error-prone reads. Both identify reads that map confidently to the human reference (GRCh38 + decoy sequences) and exclude them from downstream analysis.

#### Stage 3: Sketch Screening

**Tool**: sourmash (MinHash + Sequence Bloom Tree)

**Why it matters**: A typical clinical database contains 20,000+ microbial genomes. Aligning all non-host reads against all of them would take hours. Sketch screening reduces this to ~50–200 candidates in under 60 seconds.

**Science**: MinHash is a locality-sensitive hashing technique originally developed for near-duplicate web page detection. For DNA, a MinHash sketch compresses a genome into a fixed-size set of the smallest k-mer hash values. Two genomes share similar MinHash values with probability proportional to their Jaccard similarity (the fraction of shared k-mers). sourmash extends this with a **Scaled MinHash** that preserves containment (the fraction of a small genome's k-mers found in a larger sample), which is more appropriate for metagenomics than symmetric similarity.

The **Sequence Bloom Tree (SBT)** is a hierarchical data structure where each node is a Bloom filter representing the union of k-mers in its child genomes. Searching proceeds top-down: if a query sketch's k-mers are not contained in a parent node's Bloom filter, all descendants are pruned. This enables sub-linear search across the entire database.

PathogenIQ uses a containment threshold (default: 1%) to shortlist candidate organisms. This threshold is calibrated to be permissive enough to catch low-abundance pathogens while filtering out organisms with only spurious k-mer matches.

#### Stage 4: Targeted Alignment

**Tool**: minimap2

**Why it matters**: Sketch screening tells us which organisms are *potentially* present, but MinHash containment is a coarse approximation. Targeted alignment maps every non-host read against the shortlisted genomes with base-level accuracy, producing the precise read–organism assignment matrix needed for abundance estimation.

**Science**: minimap2 is a long-read aligner that uses a chaining algorithm to find collinear sets of seed matches between reads and reference genomes. For short reads, it operates in sr (short-read) mode with tight parameters. The output is a PAF (Pairwise mApping Format) file containing alignment coordinates, mapping quality, and the number of matching bases.

PathogenIQ parses PAF records into a binary matrix where rows are reads and columns are organisms. A cell is 1 if the read maps to that organism with sufficient quality, 0 otherwise. This matrix captures multi-mapping: a single read may have 1s in multiple columns if it aligns to multiple closely related genomes.

#### Stage 5: EM Abundance Estimation

**Algorithm**: Expectation-Maximization with bootstrap confidence intervals

**Why it matters**: Simple read-counting (count reads per organism, divide by total) fails when reads multi-map. EM provides maximum-likelihood abundance estimates that account for multi-mapping ambiguity.

**Science**: The EM algorithm iterates between two steps:

1. **E-step (Expectation)**: Given current abundance estimates, compute the posterior probability that each read originated from each organism it maps to. For a read mapping to organisms A and B with abundances θ_A and θ_B, the probability it came from A is θ_A / (θ_A + θ_B).

2. **M-step (Maximization)**: Given the posterior probabilities, update abundance estimates as the fraction of total read-mass assigned to each organism.

This converges to maximum-likelihood estimates of the relative abundances θ, under the assumption that reads are drawn from a multinomial distribution over organisms. The estimates are guaranteed to sum to 1.0.

**Bootstrap Confidence Intervals**: PathogenIQ resamples reads with replacement 100 times, re-running EM on each resampled dataset. The 2.5th and 97.5th percentiles of the bootstrapped abundance estimates form a 95% confidence interval (CI). These CIs quantify uncertainty: a wide CI indicates low confidence in the abundance estimate (few reads, high multi-mapping), while a narrow CI indicates high confidence.

#### Stage 6: Evidence Grading

**Why it matters**: Clinicians need actionable information, not raw numbers. A report saying "Staphylococcus aureus: 0.8% abundance" is less useful than "Staphylococcus aureus: Grade A — high-confidence pathogen."

**Science**: PathogenIQ assigns grades based on three criteria:

| Grade | Criteria | Clinical Interpretation |
|-------|----------|------------------------|
| **A** | Sufficient reads for specimen type + narrow CI (≤15 pp) + not a known contaminant | High-confidence pathogen. Action likely warranted. |
| **B** | Sufficient reads + not a contaminant, but CI wider than 15 pp | Probable pathogen. Consider clinical context. |
| **C** | Sufficient reads, but known contaminant or uncertain origin | Possible pathogen / contaminant. Corroborate with culture or clinical picture. |
| **X** | Insufficient reads for specimen type | Insufficient evidence. Do not act on this finding alone. |

The specimen-dependent read thresholds reflect biological prior probability: blood and CSF normally have very few circulating organisms, so even 2–3 reads of a pathogen may be significant. BAL and tissue samples contain abundant commensal flora, so higher thresholds prevent false positives.

The 15 percentage-point CI width threshold for Grade A was chosen empirically: it corresponds roughly to a "well-estimated" abundance where the lower bound is within a factor of 2–3 of the point estimate, providing sufficient precision for clinical interpretation.

---

## 4. Clinical Interpretation Engine (Plan 2)

### Contaminant Suppression

The contaminant module (`pathogeniq/contaminants.py`) encodes domain knowledge about which organisms are commonly found as contaminants in which specimen types. This is not a heuristic — it is based on decades of clinical microbiology literature documenting the epidemiology of specimen contamination:

- *Cutibacterium acnes* and coagulase-negative staphylococci (*S. epidermidis*, *S. capitis*) are the most common blood culture contaminants, introduced during venipuncture.
- *Streptococcus salivarius* and related oral streptococci are frequently aspirated into BAL samples during bronchoscopy.
- *Corynebacterium striatum* is a skin commensal that contaminates both blood and CSF collections.

When an organism is flagged as a contaminant, its grade is capped: it cannot reach Grade A, and if it would otherwise be Grade B, it is downgraded to Grade C. This reflects the clinical principle that contaminant findings require corroboration (e.g., matching culture result, clinical syndrome) before action is taken.

Critically, true pathogens are **never** suppressed. *Staphylococcus aureus* is not in the contaminant registry for any specimen type because it is a well-documented bloodstream pathogen. The substring-matching design also allows for strain-level specificity: "Staphylococcus aureus subsp. aureus" will not match "Staphylococcus epidermidis."

### AMR Gene Detection

Antimicrobial resistance is a defining feature of modern clinical microbiology. Knowing that *Staphylococcus aureus* is present is useful; knowing that it carries **mecA** (conferring methicillin resistance) is critical for antibiotic selection.

PathogenIQ uses ABRicate, a command-line tool that scans sequences against curated resistance databases. The default database is **CARD** (Comprehensive Antibiotic Resistance Database), which catalogues experimentally validated resistance genes with associated drug classes and mechanisms.

The AMR module:
1. Converts non-human reads from FASTQ to FASTA format
2. Runs ABRicate with identity and coverage thresholds (default: 90% identity, 80% coverage)
3. Parses the TSV output into structured `AMRHit` records
4. Links each resistance gene to a detected organism by matching the sequence header against organism names

AMR detection is intentionally **non-blocking**: if ABRicate is not installed, the pipeline proceeds without resistance information. This ensures the pipeline remains usable in environments where resistance databases have not yet been set up.

### PDF Clinical Report

The PDF report (`pathogeniq/pdf_report.py`) translates the computational output into a format accessible to clinicians and laboratory directors. It includes:

- **Header**: sample identifier, specimen type, date/time
- **Findings table**: organism name, abundance, confidence interval, read count, grade (colour-coded), contaminant flag
- **AMR table**: resistance genes, drug classes, identity/coverage, organism match
- **Footer**: grade definitions and research-use disclaimer

The colour coding follows clinical conventions: green for high confidence, orange for moderate, grey for low, red for insufficient. The report is designed for single-page printing and integration into electronic health record systems.

---

## 5. Validation Philosophy

PathogenIQ is validated against the **ZymoBIOMICS D6300** microbial community standard, a commercially available mock community with known composition:

- 8 bacterial species at 12% gDNA abundance each
- 2 fungal species at 2% gDNA abundance each
- Total: 10 organisms at known relative abundances

This standard is the gold benchmark for metagenomics pipelines because:
1. The ground truth is known precisely
2. It includes both high-abundance (12%) and low-abundance (2%) organisms, testing sensitivity
3. It includes Gram-positive, Gram-negative, and fungal organisms, testing taxonomic breadth
4. It is commercially available and reproducible across labs

PathogenIQ's integration tests verify:
- **Completeness**: all 10 organisms are detected
- **Accuracy**: measured abundance within 20% relative error of expected
- **Specificity**: no Grade A false positives outside the known panel

These criteria are deliberately stringent: a pipeline that passes all three tests has demonstrated both sensitivity (catches real organisms) and specificity (does not invent organisms).

---

## 6. Future Directions

### Plan 3: Nextflow Orchestration

The current Python CLI is suitable for single-machine execution. For clinical deployment, PathogenIQ will be containerized with Docker and orchestrated via Nextflow, enabling:
- **Reproducibility**: identical execution across HPC clusters, cloud VMs, and local workstations
- **Scalability**: automatic parallelization of alignment and sketching steps
- **Portability**: run on SLURM, AWS Batch, or local Docker without code changes

### Plan 4: Validation Framework

Beyond ZymoBIOMICS, clinical validation requires:
- **Spike-in studies**: known pathogens added at varying concentrations to clinical matrices (blood, CSF)
- **Clinical sample benchmarking**: comparison against culture and multiplex PCR on retrospective patient samples
- **ROC analysis**: receiver operating characteristic curves quantifying sensitivity/specificity trade-offs at different read thresholds
- **Inter-laboratory reproducibility**: same samples sequenced and analysed at multiple sites

### Long-Term Vision

PathogenIQ aims to become a **reference implementation** for clinical metagenomics — not just a tool, but a demonstration of how principled statistical methods, specimen-aware clinical knowledge, and rigorous validation can combine to produce actionable diagnostic information from sequencing data.

The ultimate goal is integration into clinical microbiology laboratory workflows, with regulatory-grade validation (e.g., CE-IVD, FDA) and health information exchange compatibility (FHIR, HL7) enabling seamless communication with electronic health records and clinical decision support systems.

---

## 7. Glossary

| Term | Definition |
|------|------------|
| **Metagenomics** | Sequencing of DNA from mixed microbial communities without prior culturing |
| **Shotgun sequencing** | Random fragmentation and sequencing of all DNA in a sample |
| **Host removal** | Computational subtraction of human DNA reads from a clinical sample |
| **MinHash** | Probabilistic sketching technique for estimating Jaccard similarity between sets |
| **SBT (Sequence Bloom Tree)** | Hierarchical index for fast containment search of genomic sketches |
| **k-mer** | A substring of length k from a DNA sequence; fundamental unit of many bioinformatics algorithms |
| **PAF (Pairwise mApping Format)** | Text format for approximate mapping coordinates between reads and references |
| **EM (Expectation-Maximization)** | Iterative statistical algorithm for maximum-likelihood estimation with latent variables |
| **Bootstrap** | Resampling technique for estimating confidence intervals by repeated sampling with replacement |
| **CI (Confidence Interval)** | Range of values within which the true parameter is expected to lie with a given probability (95% for PathogenIQ) |
| **ABRicate** | Command-line tool for scanning sequences against curated antimicrobial resistance databases |
| **CARD** | Comprehensive Antibiotic Resistance Database; curated collection of resistance genes and associated phenotypes |
| **Commensal** | Microorganism that lives on or in a host without causing disease |
| **Contaminant** | Microorganism present due to sample collection or processing, not reflecting true infection |
| **BAL** | Bronchoalveolar lavage; fluid sample from the lungs |
| **CSF** | Cerebrospinal fluid; sample from the central nervous system |
| **GRCh38** | Reference genome assembly for Homo sapiens (Genome Reference Consortium Human build 38) |
| **ZymoBIOMICS** | Commercial mock microbial community standard with known composition |

---

## 8. References

1. **MinHash**: Broder, A. Z. (1997). On the resemblance and containment of documents. *Proceedings of Compression and Complexity of Sequences*.
2. **Scaled MinHash for genomics**: Brown, C. T., & Irber, L. (2016). sourmash: a library for MinHash sketching of DNA. *Journal of Open Source Software*.
3. **Sequence Bloom Tree**: Solomon, B., & Kingsford, C. (2016). Fast search of thousands of short-read sequencing experiments. *Nature Biotechnology*.
4. **BWA-MEM2**: Vasimuddin, M., et al. (2019). Efficient Architecture-Aware Acceleration of BWA-MEM for Multicore Systems. *IPDPS*.
5. **minimap2**: Li, H. (2018). Minimap2: pairwise alignment for nucleotide sequences. *Bioinformatics*.
6. **EM for metagenomics**: using the same probabilistic framework as Bracken (Lu, J., et al. 2017. *Nature Communications*) but applied to per-samplealignment rather than k-mer counts.
7. **CARD**: Alcock, B. P., et al. (2020). CARD 2020: antibiotic resistome surveillance with the comprehensive antibiotic resistance database. *Nucleic Acids Research*.
8. **Clinical metagenomics**: Wilson, M. R., et al. (2019). Actionable diagnosis of neuroleptospirosis by next-generation sequencing. *New England Journal of Medicine*.
9. **ZymoBIOMICS standard**: https://www.zymoresearch.com/products/zymobiomics-microbial-community-dna-standard

---

*Last updated: 2026-04-23*
