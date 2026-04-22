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

Downloads ~60+ reference genomes covering:
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
