# PathogenIQ — Database Setup Scripts

Run scripts 01–03 **concurrently in separate terminals** to save time.
Total setup time: ~90 min on a 1 Gbps connection with 16-core machine.

---

## Quick Start (Concurrent Setup)

Open three terminals and run each block simultaneously:

**Terminal 1 — ZymoBIOMICS validation DB** (~5 min):
```bash
cd /data/alvin/Metagenomics
bash scripts/01_download_zymo_db.sh 16
```

**Terminal 2 — Full Tier-1 pathogen DB** (~30–90 min):
```bash
cd /data/alvin/Metagenomics
pip install ncbi-genome-download sourmash
python scripts/02_download_tier1_db.py --threads 16
```

**Terminal 3 — Human reference** (~90 min total: download + BWA-MEM2 index):
```bash
cd /data/alvin/Metagenomics
bash scripts/03_download_human_ref.sh 16
```

---

## Prerequisites

```bash
# Bioinformatics tools (choose one)
conda install -c bioconda bwa-mem2 minimap2 samtools fastp sourmash ncbi-genome-download

# Or pip (Python tools only)
pip install ncbi-genome-download sourmash
```

---

## What Each Script Does

### `01_download_zymo_db.sh` — ZymoBIOMICS Validation DB

Downloads exactly the 10 reference genomes in the ZymoBIOMICS D6300 standard:

| Organism | NCBI Accession | Expected Abundance |
|---|---|---|
| *Pseudomonas aeruginosa* PAO1 | GCF_000006765.1 | 12% |
| *Escherichia coli* K-12 MG1655 | GCF_000005845.2 | 12% |
| *Salmonella enterica* LT2 | GCF_000006945.2 | 12% |
| *Lactobacillus fermentum* ATCC14931 | GCF_000759415.1 | 12% |
| *Enterococcus faecalis* V583 | GCF_000007785.1 | 12% |
| *Staphylococcus aureus* NCTC8325 | GCF_000013425.1 | 12% |
| *Listeria monocytogenes* EGD-e | GCF_000196035.1 | 12% |
| *Bacillus subtilis* 168 | GCF_000009045.1 | 12% |
| *Saccharomyces cerevisiae* S288C | GCF_000146045.2 | 2% |
| *Cryptococcus neoformans* H99 | GCF_000149245.2 | 2% |

**Output:** `databases/zymo_tier1/zymo_tier1.sbt.zip`

### `02_download_tier1_db.py` — Full Tier-1 Pathogen DB

Downloads ~110 reference genomes covering:
- WHO Priority Pathogens (A/B/C tiers)
- CDC Select Agents
- Common clinical viruses (HSV, CMV, EBV, SARS-CoV-2, influenza)
- Pathogenic fungi (*Candida*, *Aspergillus*, *Cryptococcus*)
- Key parasites (*Plasmodium*, *Toxoplasma*, *Cryptosporidium*)

**ZymoBIOMICS-only mode** (for quick start before full DB is ready):
```bash
python scripts/02_download_tier1_db.py --zymo-only
```

**Output:** `databases/tier1/tier1_pathogens.sbt.zip`

### `03_download_human_ref.sh` — Human Reference for Host Removal

Downloads GRCh38 primary assembly + PhiX decoy, builds BWA-MEM2 and minimap2 indexes.

> **Note:** BWA-MEM2 indexing requires ~64 GB RAM and ~60 min.
> If your machine has <64 GB RAM, use `minimap2` for host removal instead
> (set `--read-type long` or modify `host_remove.py` to use minimap2 for SR).

**Output:** `databases/human_ref/GRCh38_decoy.fa` + index files

---

## After Download — Run Integration Tests

```bash
export HUMAN_REF=$(pwd)/databases/human_ref/GRCh38_decoy.fa
export TIER1_DB=$(pwd)/databases/zymo_tier1/zymo_tier1.sbt.zip  # start with zymo

pytest tests/integration/ -v -m integration
```

---

## Run the Full Pipeline on ZymoBIOMICS Data

```bash
export HUMAN_REF=$(pwd)/databases/human_ref/GRCh38_decoy.fa
export TIER1_DB=$(pwd)/databases/zymo_tier1/zymo_tier1.sbt.zip

pathogeniq run \
  --input ZymoStandardsFromMao/V350234554_L01_UDB-32.Zymo_Std.fq.gz \
  --output results/zymo_rep1 \
  --db "${TIER1_DB}" \
  --host-ref "${HUMAN_REF}" \
  --specimen blood \
  --read-type short \
  --threads 16
```

Report will be at `results/zymo_rep1/report/pathogeniq_report.json`.

---

## Validation, Benchmark & Open-World Scripts (04–16)

Scripts 01–03 build the databases the pipeline needs. The rest support validation,
benchmarking, the NTC background, and the open-world / air-surveillance arm.

| Script | Purpose |
|--------|---------|
| `04_download_validation_data.py` | Download labelled data from ENA (CAMI/HMP, NTC blanks, and the Jeilu et al. aircraft datasets — `blanks-hunt`, `blanks-qiita`, `air-aircraft-filter`, `air-aircraft-ntc`, `air-aircraft-spike`) with title filters + `--limit` |
| `05_select_kitome_controls.py` | Select kitome/negative controls for the background pool |
| `06_benchmark.py` | Single-community grading benchmark vs raw Kraken2 |
| `07_build_kraken_db.py` | Build a custom Kraken2 DB from the tier-1 genomes (benchmark baseline) |
| `08_heldout_pr_auc.py` | Multi-community held-out PR-AUC / precision@recall (floor frozen on Zymo) |
| `09_reads_mapping_truth.py` | Roll CAMI `reads_mapping.tsv` up into per-sample truth |
| `10_validate_dispersion.py` | Leave-one-out dispersion-prior validation on kitome blanks (`--selfcheck`) |
| `11_pool_blank_background.py` | Classify + pool negative-control blanks into a Tier-2 background table |
| `12_viral_insilico_spikein.py` | In-silico viral spike-in validation (wgsim → assembly → geNomad/CheckV → score); offline `--selfcheck` |
| `13_setup_viral_env.sh` | Set up the viral arm: isolated `genomad` env (geNomad + CheckV), fetch DBs, symlink CLIs into the pathogeniq env |
| `15_setup_mag_env.sh` | Set up the MAG arm: isolated `mag-bin` env (metabat2 + CheckM) and `gtdbtk` env, fetch the CheckM DB (~1.4 GB) + GTDB-Tk r232 DB (~94 GB) onto a data volume, wrap the CLIs onto the core `PATH`. Loud FATAL verification per stage |
| `16_run_air_assemble.sh` | Run the `--assemble` MAG arm on the 5 PRJNA1228129 aircraft-filter samples (fail-fast toolchain check, one sample at a time) |

(`build_background_default.py` is a helper that rebuilds the shipped Tier-2
`pathogeniq/data/background_default.tsv` from a pooled blank classification.)

**Open-world DBs (optional).** The novelty trigger needs a broad Kraken2 DB
(`databases/kraken2`, ~14 GB); the viral arm needs geNomad + CheckV DBs (`scripts/13`,
~3 GB). The MAG arm's **GTDB-Tk DB (~94 GB, release r232)** plus the CheckM DB
(~1.4 GB) are installed by `scripts/15` — deliberately not part of the core install,
they are a triggered Tier-2 extra. geNomad/CheckV/ABRicate and metabat2/CheckM/GTDB-Tk
are all installed in their own conda envs (they require `numpy<2`) and exposed on
`PATH`; never install them into the core `pathogeniq` env.

---

## Alternative: Pre-built sourmash Databases

If you want to skip genome downloads entirely, sourmash provides pre-built databases
covering all RefSeq genomes. These are larger but require no custom build step:

```bash
# GTDB representative genomes (bacteria/archaea) — ~15 GB
curl -L https://farm.cse.ucdavis.edu/~ctbrown/sourmash-db/gtdb-rs214/gtdb-rs214-reps.k31.sbt.zip \
  -o databases/gtdb_reps.k31.sbt.zip

# GenBank viral — ~2 GB
curl -L https://farm.cse.ucdavis.edu/~ctbrown/sourmash-db/genbank-2022.03/genbank-2022.03-viral.k31.sbt.zip \
  -o databases/genbank_viral.k31.sbt.zip
```

> These pre-built DBs are ideal for research. For clinical use, the custom Tier-1 DB
> (script 02) is preferable because it is validated and versioned.
