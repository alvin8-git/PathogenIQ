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

**Science**: Grading is a **pure function** of a single `GradingInput` record (read count, abundance, CI width, contaminant flag, specimen type, NTC tier, cross-mapping flag) — see `pathogeniq/report.py::grade()`. Keeping it pure and side-effect-free makes the decision auditable and unit-testable, and lets the *same* function grade output from any upstream classifier (it is reused to score Kraken2 output in the benchmark, §6).

| Grade | Criteria | Clinical Interpretation |
|-------|----------|------------------------|
| **A** | Sufficient reads for specimen type + narrow CI (≤15 pp) + not a contaminant + **batch-matched NTC available (Tier 1)** | High-confidence pathogen. Action likely warranted. |
| **B** | Sufficient reads + not a contaminant, but CI wider than 15 pp **or no batch-matched NTC (Tier 2/3)** | Probable pathogen. Consider clinical context. |
| **C** | Sufficient reads, but known contaminant or uncertain origin | Possible pathogen / contaminant. Corroborate with culture or clinical picture. |
| **X** | Insufficient reads for specimen type, **degenerate statistics** (NaN/inf/negative), **or flagged as a cross-mapping artifact / NTC background** | Insufficient evidence. Do not act on this finding alone. |

The specimen-dependent read thresholds reflect biological prior probability: blood and CSF normally have very few circulating organisms, so even 2–3 reads of a pathogen may be significant. BAL and tissue samples contain abundant commensal flora, so higher thresholds prevent false positives.

The 15 percentage-point CI width threshold for Grade A was chosen empirically: it corresponds roughly to a "well-estimated" abundance where the lower bound is within a factor of 2–3 of the point estimate, providing sufficient precision for clinical interpretation.

Two additional gates were added in Plan 4 and are the subject of §5:

- **NTC tier cap.** Grade A is *only* reachable when a batch-matched no-template control (NTC) accompanies the run (Tier 1). Without one, the best attainable grade is B — the laboratory cannot claim the highest confidence for a finding it has not controlled for reagent/kitome contamination. This cap is the central honesty mechanism of the grading layer.
- **Cross-mapping dedup.** When two reference genomes are >95% identical (e.g. *Escherichia coli* and the four *Shigella* species), reads from the dominant organism shadow-map onto its near-twins, manufacturing phantom polymicrobial calls. `pathogeniq/crossmap.py` detects when a known cross-mapping group member is outnumbered ≥10× by a relative and demotes the minor member to Grade X.

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

## 5. NTC-Aware Background Subtraction & Tiered Grading (Plan 4)

This is the core scientific contribution that distinguishes PathogenIQ from a fast classifier with a pretty report. It directly attacks the central problem of §2: **separating a real low-abundance pathogen from reagent and environmental contamination** in a low-biomass clinical specimen.

### 5.1 Motivation — the kitome problem

DNA extraction kits, enzymes, water, and lab surfaces are not sterile. They carry trace bacterial DNA — collectively the **"kitome"** and **"splashome"**. In a high-biomass sample this is invisible noise. But clinical metagenomics operates at the opposite extreme: a blood or CSF specimen may yield only a few hundred to a few thousand microbial reads total. At that depth, the kitome is not noise — it is a *competing signal*, often larger than the pathogen it obscures. The literature is unambiguous that contaminants scale inversely with input biomass (Salter et al. 2014), and that the same handful of genera — *Pseudomonas*, *Cutibacterium*, *Ralstonia*, *Bradyrhizobium*, Enterobacteriaceae — recur across studies and labs.

A pipeline that reports any organism above a fixed read threshold therefore generates a flood of false positives in exactly the samples where diagnostic stakes are highest. The clinically correct null hypothesis is not "this organism is absent" but **"this organism is indistinguishable from what my reagents would have produced anyway."** Testing that hypothesis requires a model of the background.

### 5.2 The no-template control (NTC)

The instrument for measuring reagent background is a **no-template control**: the full extraction-and-sequencing protocol run on a blank (water/buffer, no clinical material). Whatever DNA appears in an NTC came from the process, not the patient. PathogenIQ models the per-organism background rate from one or more NTCs and tests each sample organism against it (`pathogeniq/background.py`).

**Why rates, not counts.** A deep sample and a shallow NTC are not comparable by raw count. The background is normalised to **reads per million classified reads (RPM)**, making the test depth-invariant; depth re-enters only when a background *rate* is scaled back to an expected *count* at the sample's own sequencing depth.

**The statistical test.** For each organism, the expected background count at the sample's depth is `λ = rate_RPM × sample_depth / 1e6`. The observed sample count is compared against a **negative-binomial upper tail** with mean λ and a dispersion parameter `r` (overdispersion captures the fact that contamination is lumpier than a Poisson process). A low tail p-value means the observed count exceeds what background plausibly explains → real signal; a high p-value (≥ α, default 0.01) means the organism is indistinguishable from background and is zeroed out before grading. Two guards make this robust at low biomass:

- A **pseudocount** (0.5 RPM) gives every organism — even those absent from the NTC — a finite background floor, so a genuinely novel pathogen still has to clear a small bar rather than dividing by zero.
- A **minimum-support floor** (≥2 reads summed across controls) prevents a single spurious read in a thin blank from inflating to a huge RPM that would suppress real low-level pathogens.

### 5.3 The tier system — honesty about provenance

Not all background models are equal, and the grading layer is explicit about it:

| Tier | Background source | Max grade | Rationale |
|------|-------------------|-----------|-----------|
| **Tier 1** | Batch-matched NTC (run alongside this specimen, same reagent lot) | **A** | The only configuration that actually controls this run's contamination. |
| **Tier 2** | Pooled/foreign background (shipped default, built from public blanks) | **B** | Catches universal kitome genera, but a frozen prior cannot capture this batch's reagent lot. Capped. |
| **Tier 3** | No background available | **B** | Graded uncorrected; cannot claim contamination was controlled. |

This tiering is why a batch-matched NTC is *required* for any Grade A call (§3, Stage 6). It encodes a clinical-laboratory principle: you may not assert the highest confidence in a finding whose contamination you did not control.

### 5.4 What validation taught us about the limits

Two empirical studies (reproducible via `scripts/10_validate_dispersion.py` and `scripts/06`/`08`) shaped the final design and are documented honestly because they bound what the method can claim:

- **The dispersion prior is not the lever; coverage is.** A leave-one-out experiment across spike-free kitome blanks showed the single-NTC false-positive rate is 30–65× the nominal α at *every* dispersion value — it never approaches α. The cause is coverage, not the prior: ~10 of 18 kitome taxa appear in only one blank, so a contaminant absent from the (small) control pool hits the pseudocount floor and reads as "real" regardless of `r`. **This is the empirical justification for the Tier-2 → Grade B cap:** a single pooled background genuinely cannot deliver controlled-α detection, and no parameter tuning rescues it — only more/batch-matched NTCs do.
- **The kitome background must not be applied blindly.** On the Zymo benchmark, the pooled Salter background *dropped recall* by suppressing real Enterobacteriaceae (genuine Zymo members that are also common kitome). Consequently the background filter is not auto-applied to samples likely to contain Enterobacteriaceae; the cross-mapping dedup (§3) — which removes only phantom relatives of a dominant organism at zero recall cost — is the safe, always-on artifact filter.

### 5.5 Building the shipped Tier-2 background

The packaged `pathogeniq/data/background_default.tsv` is pooled from **genuine non-spiked negative-control shotgun datasets** (`scripts/11_pool_blank_background.py`), classified against the same Tier-1 reference DB so the per-organism rates are directly comparable. Because the dispersion study showed coverage is the binding constraint, the build deliberately pools blanks from **multiple independent labs and reagent types** to widen taxon coverage and dilute any single study's idiosyncrasy. The datasets are listed in §6.

---

## 6. Validation & Benchmarking

PathogenIQ's validation evolved from a single mock standard into a multi-community, cross-environment benchmark with a leakage-free held-out split. Three classes of data are used, for three different purposes.

### 6.1 Mock-standard validation — ZymoBIOMICS D6300

The **ZymoBIOMICS D6300** microbial community standard is a commercially available mock community with known composition:

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

These criteria are deliberately stringent: a pipeline that passes all three tests has demonstrated both sensitivity (catches real organisms) and specificity (does not invent organisms). Three sequencing runs of D6300 (UDB-32/-40/-48) are used together to test run-to-run reproducibility.

### 6.2 Cross-community truth — CAMI II simulated communities

A mock standard is a single community; it cannot show whether the grading thresholds *generalise*. For that PathogenIQ uses the **CAMI II** (Critical Assessment of Metagenome Interpretation) challenge datasets — CAMISIM-simulated short-read communities shipped with a **gold-standard taxonomic profile** (`.profile`), parsed per-sample by `pathogeniq/benchmark.py::parse_cami_profile`. Three structurally different environments are used as held-out communities:

| Community | Character | Truth (per-sample, species) |
|-----------|-----------|------------------------------|
| **Mouse gut** | Host-associated anaerobes | 64 species |
| **Marine** | Open-environment diversity | 259 species |
| **Strain-madness** | Few species, many near-identical strains | 20 species |

Strain-madness is deliberately the cross-mapping stress test (§3) — it is where read-count signal alone struggles most.

### 6.3 Real low-biomass communities — HMP body sites

To test on genuinely host-associated, low-biomass material (the regime closest to clinical specimens), PathogenIQ adds three **CAMI II Human Microbiome Project** body-site samples — **skin, airway, oral**. These ship per-read truth (`reads_mapping.tsv`) rather than a profile; truth is reconstructed by `scripts/09_reads_mapping_truth.py`. Because their 97%-identity OTU truth resolves only to **genus**, they are scored at genus rank as a *separate* panel (mixing ranks in one average would be apples-to-oranges).

### 6.4 NTC background datasets

The Tier-2 background (§5.5) is built from public **negative-control** shotgun runs, downloaded via `scripts/04_download_validation_data.py` (ENA filereport API) and classified by `scripts/11_pool_blank_background.py`:

| Dataset | ENA accession | Role |
|---------|---------------|------|
| Salter et al. 2014 reagent controls | ERP006808 | Initial background; **deprecated** as a default source — it is a spiked *S. bongori* dilution series from one lab, so it is idiosyncratic and Enterobacteriaceae-heavy. |
| HUNT One Health reagent/extraction blanks | PRJEB66439 | Non-spiked, true-`RANDOM` shotgun blanks; deep enough to anchor per-organism RPM. |
| Qiita 14332 nucleic-acid pipeline blanks | PRJEB56784 | Non-spiked EtOH/IPA/swab/buffer reagent blanks across reagent types — coverage *breadth*. |

The shift from Salter-only to a pool of HUNT + Qiita blanks across independent labs is a direct response to the §5.4 finding that **NTC coverage**, not the statistical prior, limits the background's quality.

### 6.5 Held-out benchmark results

The grading layer is scored against raw Kraken2 (a strong, widely-used classifier) on each labeled community. Crucially, the read-count operating point (floor) is **calibrated once on the three Zymo runs and then frozen**, so every CAMI/HMP community is genuinely out-of-sample (`scripts/08_heldout_pr_auc.py`). The headline metric is **operating-point precision** at preserved recall — grading's value is removing the low-support false-positive tail, which a ranking metric like PR-AUC is insensitive to when those false positives are already bottom-ranked.

| Held-out panel | Communities | Raw precision | Graded precision | Recall |
|----------------|-------------|---------------|------------------|--------|
| CAMI II (species rank) | mouse-gut, marine, strain | 0.011 | **0.494** (≈45×) | 0.684 |
| HMP (genus rank) | skin, airway, oral | 0.008 | **0.352** (≈44×) | 0.952 |

A single floor learned on a mock standard lifts precision ~45× across **nine held-out communities** spanning two taxonomic ranks and four environment types, none seen during calibration — strong evidence the grading wedge is community- and rank-agnostic rather than overfit to Zymo. Full per-community results: `docs/benchmark-results-2026-06-21.md`.

---

## 7. Future Directions

### Plan 3: Nextflow Orchestration

The current Python CLI is suitable for single-machine execution. For clinical deployment, PathogenIQ will be containerized with Docker and orchestrated via Nextflow, enabling:
- **Reproducibility**: identical execution across HPC clusters, cloud VMs, and local workstations
- **Scalability**: automatic parallelization of alignment and sketching steps
- **Portability**: run on SLURM, AWS Batch, or local Docker without code changes

### Plan 4: NTC-Aware Grading & Validation — largely delivered

The NTC-aware tiered grading, NB background subtraction, cross-mapping dedup, and the multi-community held-out benchmark (§5, §6) are implemented and tested. What remains to reach clinical-grade validation:

- **Batch-matched NTC studies**: the §5.4 result shows a single pooled background cannot deliver controlled-α detection — Tier 1 (a per-batch NTC) is the path to Grade A, and needs prospective evaluation on paired sample+NTC runs.
- **Spike-in studies**: known pathogens added at varying concentrations to clinical matrices (blood, CSF) to characterise the limit of detection per specimen type.
- **Clinical sample benchmarking**: comparison against culture and multiplex PCR on retrospective patient samples.
- **Inter-laboratory reproducibility**: same samples sequenced and analysed at multiple sites.

### Plan 5: Air / bioaerosol surveillance (exploratory)

The same NTC backbone, applied to environmental air samples (even lower biomass, even more kitome-dominated), with breadth-of-coverage gating to separate clumped contamination artifacts from genome-wide real signal. Parked for future exploration; see `todo.md`.

### Long-Term Vision

PathogenIQ aims to become a **reference implementation** for clinical metagenomics — not just a tool, but a demonstration of how principled statistical methods, specimen-aware clinical knowledge, and rigorous validation can combine to produce actionable diagnostic information from sequencing data.

The ultimate goal is integration into clinical microbiology laboratory workflows, with regulatory-grade validation (e.g., CE-IVD, FDA) and health information exchange compatibility (FHIR, HL7) enabling seamless communication with electronic health records and clinical decision support systems.

---

## 8. Glossary

| Term | Definition |
|------|------------|
| **Metagenomics** | Sequencing of DNA from mixed microbial communities without prior culturing |
| **NTC (No-Template Control)** | A blank (water/buffer, no clinical material) run through the full protocol to measure reagent/process background |
| **Kitome** | Background microbial DNA contributed by extraction kits, reagents, and lab surfaces; dominant noise in low-biomass samples |
| **RPM (Reads Per Million)** | Depth-normalised abundance unit (reads per million classified reads) making samples and NTCs comparable |
| **Negative binomial** | Overdispersed count distribution used to model background; tail p-value tests whether an observed count exceeds background |
| **Dispersion (r)** | Negative-binomial parameter controlling background variance; smaller r = more overdispersion |
| **Tier (1/2/3)** | Provenance level of the background model (batch-matched / pooled / none); caps the maximum attainable grade |
| **Cross-mapping** | Reads from one organism mapping onto a >95%-identical relative, manufacturing phantom detections (e.g. *E. coli* → *Shigella*) |
| **CAMI** | Critical Assessment of Metagenome Interpretation; simulated communities with gold-standard truth profiles |
| **OTU** | Operational Taxonomic Unit; a clustered group of sequences (97% identity here), often resolving only to genus |
| **PR-AUC / precision@recall** | Ranking and operating-point summaries of a precision–recall curve; the benchmark metrics |
| **Held-out evaluation** | Calibrating a threshold on one set of communities and reporting on disjoint, unseen communities |
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

## 9. References

1. **MinHash**: Broder, A. Z. (1997). On the resemblance and containment of documents. *Proceedings of Compression and Complexity of Sequences*.
2. **Scaled MinHash for genomics**: Brown, C. T., & Irber, L. (2016). sourmash: a library for MinHash sketching of DNA. *Journal of Open Source Software*.
3. **Sequence Bloom Tree**: Solomon, B., & Kingsford, C. (2016). Fast search of thousands of short-read sequencing experiments. *Nature Biotechnology*.
4. **BWA-MEM2**: Vasimuddin, M., et al. (2019). Efficient Architecture-Aware Acceleration of BWA-MEM for Multicore Systems. *IPDPS*.
5. **minimap2**: Li, H. (2018). Minimap2: pairwise alignment for nucleotide sequences. *Bioinformatics*.
6. **EM for metagenomics**: using the same probabilistic framework as Bracken (Lu, J., et al. 2017. *Nature Communications*) but applied to per-samplealignment rather than k-mer counts.
7. **CARD**: Alcock, B. P., et al. (2020). CARD 2020: antibiotic resistome surveillance with the comprehensive antibiotic resistance database. *Nucleic Acids Research*.
8. **Clinical metagenomics**: Wilson, M. R., et al. (2019). Actionable diagnosis of neuroleptospirosis by next-generation sequencing. *New England Journal of Medicine*.
9. **ZymoBIOMICS standard**: https://www.zymoresearch.com/products/zymobiomics-microbial-community-dna-standard
10. **Reagent contamination / kitome**: Salter, S. J., et al. (2014). Reagent and laboratory contamination can critically impact sequence-based microbiome analyses. *BMC Biology*.
11. **CAMI benchmark**: Meyer, F., et al. (2022). Critical Assessment of Metagenome Interpretation — the second round of challenges. *Nature Methods*.
12. **Kraken2 (benchmark baseline)**: Wood, D. E., Lu, J., & Langmead, B. (2019). Improved metagenomic analysis with Kraken 2. *Genome Biology*.
13. **NTC datasets**: HUNT One Health (ENA PRJEB66439); Qiita study 14332 (ENA PRJEB56784) — non-spiked negative-control shotgun runs pooled for the Tier-2 background.

---

*Last updated: 2026-06-22*

*Design changes reflected: Plan 4 NTC-aware tiered grading, negative-binomial background subtraction (`background.py`), cross-mapping dedup (`crossmap.py`), the Kraken2 grading adapter and multi-community held-out benchmark (`benchmark.py`, `scripts/06`/`08`/`09`/`10`/`11`), and the dispersion-prior validation. See also `docs/benchmark-results-2026-06-21.md` and `docs/dispersion-validation-2026-06-22.md`.*
