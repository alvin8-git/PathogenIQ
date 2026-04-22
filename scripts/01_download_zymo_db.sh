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
declare -A ACCESSIONS=(
  ["Pseudomonas_aeruginosa_PAO1"]="GCF_000006765.1"
  ["Escherichia_coli_K12_MG1655"]="GCF_000005845.2"
  ["Salmonella_enterica_LT2"]="GCF_000006945.2"
  ["Lactobacillus_fermentum_ATCC14931"]="GCF_000759415.1"
  ["Enterococcus_faecalis_V583"]="GCF_000007785.1"
  ["Staphylococcus_aureus_NCTC8325"]="GCF_000013425.1"
  ["Listeria_monocytogenes_EGDe"]="GCF_000196035.1"
  ["Bacillus_subtilis_168"]="GCF_000009045.1"
  ["Saccharomyces_cerevisiae_S288C"]="GCF_000146045.2"
  ["Cryptococcus_neoformans_H99"]="GCF_000149245.2"
)

# Build comma-separated accession list for ncbi-genome-download
ACCESSION_LIST=$(IFS=,; echo "${ACCESSIONS[*]}")

# Download all 10 genomes in parallel (ncbi-genome-download handles concurrency)
ncbi-genome-download \
  --accessions "${ACCESSION_LIST}" \
  --formats fasta \
  --output-folder "${GENOME_DIR}" \
  --parallel ${THREADS} \
  --retries 3 \
  all

echo "=== Step 2: Decompress and flatten genome files ==="

find "${GENOME_DIR}" -name "*.fna.gz" | while read -r gz; do
  name=$(basename "$(dirname "${gz}")")
  gunzip -c "${gz}" > "${GENOME_DIR}/${name}.fna"
done

echo "=== Step 3: Sketch genomes in parallel (sourmash, k=31, scaled=1000) ==="

# Use GNU parallel if available, otherwise xargs -P
FASTA_LIST=$(find "${GENOME_DIR}" -name "*.fna" | tr '\n' ' ')

if command -v parallel &>/dev/null; then
  find "${GENOME_DIR}" -name "*.fna" | \
    parallel -j ${THREADS} \
    sourmash sketch dna -p k=31,scaled=1000 {} -o "${SIG_DIR}/{/.}.sig"
else
  find "${GENOME_DIR}" -name "*.fna" | \
    xargs -P ${THREADS} -I{} bash -c \
    'sourmash sketch dna -p k=31,scaled=1000 "$1" -o "'"${SIG_DIR}"'/$(basename "$1" .fna).sig"' _ {}
fi

echo "=== Step 4: Build sourmash SBT index ==="

sourmash index \
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
