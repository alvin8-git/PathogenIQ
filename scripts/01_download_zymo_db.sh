#!/usr/bin/env bash
# Download the 10 ZymoBIOMICS D6300 reference genomes from NCBI RefSeq
# and build a sourmash SBT for immediate pipeline validation.
#
# Runtime: ~5 min on a decent connection
# Output:  databases/zymo_tier1/zymo_tier1.sbt.zip
#
# Prerequisites:
#   pip install ncbi-genome-download sourmash
#   (or: conda install -c bioconda ncbi-genome-download sourmash)

set -euo pipefail

OUTDIR="databases/zymo_tier1"
GENOME_DIR="${OUTDIR}/genomes"
SIG_DIR="${OUTDIR}/sigs"
DB_OUT="${OUTDIR}/zymo_tier1.sbt.zip"
THREADS=${1:-8}

mkdir -p "${GENOME_DIR}" "${SIG_DIR}"

echo "=== Step 1: Download ZymoBIOMICS reference genomes in parallel ==="

# NCBI RefSeq accessions for ZymoBIOMICS D6300 community members
# These are the representative/reference strains closest to the ZymoBIOMICS material
ACCESSION_LIST="GCF_000006765.1,GCF_000005845.2,GCF_000006945.2,GCF_000759415.1,GCF_000007785.1,GCF_000013425.1,GCF_000196035.1,GCF_000009045.1,GCF_000146045.2,GCF_000149245.2"

# Organism legend:
#   GCF_000006765.1  Pseudomonas aeruginosa PAO1
#   GCF_000005845.2  Escherichia coli K-12 MG1655
#   GCF_000006945.2  Salmonella enterica LT2
#   GCF_000759415.1  Lactobacillus fermentum ATCC14931
#   GCF_000007785.1  Enterococcus faecalis V583
#   GCF_000013425.1  Staphylococcus aureus NCTC8325
#   GCF_000196035.1  Listeria monocytogenes EGD-e
#   GCF_000009045.1  Bacillus subtilis 168
#   GCF_000146045.2  Saccharomyces cerevisiae S288C
#   GCF_000149245.2  Cryptococcus neoformans H99

# Download all 10 genomes in parallel (ncbi-genome-download handles concurrency)
# Correct short flags: -A accessions, -F formats, -o output, -p parallel, -r retries
ncbi-genome-download \
  -A "${ACCESSION_LIST}" \
  -F fasta \
  -o "${GENOME_DIR}" \
  -p ${THREADS} \
  -r 3 \
  all

echo "=== Step 2: Decompress and flatten genome files ==="

find "${GENOME_DIR}" -name "*.fna.gz" | while read -r gz; do
  name=$(basename "$(dirname "${gz}")")
  gunzip -c "${gz}" > "${GENOME_DIR}/${name}.fna"
done

echo "=== Step 3: Sketch genomes (sourmash, k=31, scaled=1000) ==="

SOURMASH=$(command -v sourmash) || { echo "ERROR: sourmash not found. Run: pip install sourmash"; exit 1; }
echo "Using sourmash: ${SOURMASH}"

while IFS= read -r fna; do
  base=$(basename "${fna}" .fna)
  sig="${SIG_DIR}/${base}.sig"
  if [ ! -f "${sig}" ]; then
    echo "  Sketching ${base}..."
    "${SOURMASH}" sketch dna -p k=31,scaled=1000 "${fna}" -o "${sig}" --quiet
  else
    echo "  Skipping ${base} (already sketched)"
  fi
done < <(find "${GENOME_DIR}" -name "*.fna")

echo "=== Step 4: Build sourmash SBT index ==="

"${SOURMASH}" index \
  --ksize 31 \
  "${DB_OUT}" \
  "${SIG_DIR}"/*.sig

echo ""
echo "=== Done ==="
echo "ZymoBIOMICS validation DB: ${DB_OUT}"
echo ""
echo "Run integration tests with:"
echo "  TIER1_DB=$(pwd)/${DB_OUT} HUMAN_REF=/path/to/GRCh38.fa \\"
echo "    pytest tests/integration/ -v -m integration"
