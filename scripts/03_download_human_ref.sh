#!/usr/bin/env bash
# Download GRCh38 + decoy sequences and index for BWA-MEM2 and minimap2.
#
# Runtime: ~30-60 min (download) + ~60 min (BWA-MEM2 index)
# Output:  databases/human_ref/GRCh38_decoy.fa + index files
#
# Prerequisites:
#   conda install -c bioconda bwa-mem2 minimap2 samtools
#   (bwa-mem2 index requires ~64 GB RAM for GRCh38)

set -euo pipefail

OUTDIR="databases/human_ref"
THREADS=${1:-16}
mkdir -p "${OUTDIR}"

echo "=== Step 1: Download GRCh38 primary assembly + decoy sequences ==="

# Primary assembly (no alt contigs — cleaner for metagenomics)
if [ ! -f "${OUTDIR}/GRCh38_no_alt.fa.gz" ]; then
  curl -C - --retry 5 --retry-delay 10 -L \
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/GCA_000001405.15_GRCh38_assembly_structure/Primary_Assembly/assembled_chromosomes/FASTA/chroms/" \
    -o "${OUTDIR}/GRCh38_primary.html"

  # Download chromosomes 1-22, X, Y, M in parallel
  echo "Downloading chromosomes in parallel..."
  for CHR in {1..22} X Y M; do
    URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/GCA_000001405.15_GRCh38_assembly_structure/Primary_Assembly/assembled_chromosomes/FASTA/chr${CHR}.fna.gz"
    curl -C - --retry 5 -L -o "${OUTDIR}/chr${CHR}.fna.gz" "${URL}" &
  done
  wait
  echo "All chromosomes downloaded."
else
  echo "GRCh38 already downloaded, skipping."
fi

echo "=== Step 2: Download decoy sequences (viral, unmapped) ==="

# hs37d5 decoy — commonly used set of known non-human contaminants
DECOY_URL="https://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz"
if [ ! -f "${OUTDIR}/hs37d5_decoy.fa.gz" ]; then
  curl -C - --retry 5 -L -o "${OUTDIR}/hs37d5_decoy.fa.gz" "${DECOY_URL}" &
fi

# PhiX control genome
PHIX_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/819/615/GCF_000819615.1_ViralProj14015/GCF_000819615.1_ViralProj14015_genomic.fna.gz"
if [ ! -f "${OUTDIR}/phix.fna.gz" ]; then
  curl -C - --retry 5 -L -o "${OUTDIR}/phix.fna.gz" "${PHIX_URL}" &
fi

wait
echo "Decoy sequences downloaded."

echo "=== Step 3: Concatenate into single reference ==="

COMBINED="${OUTDIR}/GRCh38_decoy.fa"
if [ ! -f "${COMBINED}" ]; then
  gunzip -c "${OUTDIR}"/chr*.fna.gz > "${COMBINED}"
  gunzip -c "${OUTDIR}/phix.fna.gz" >> "${COMBINED}"
  echo "Combined reference created: ${COMBINED}"
else
  echo "Combined reference already exists."
fi

echo "=== Step 4: Build BWA-MEM2 index (requires ~64 GB RAM, ~60 min) ==="

if [ ! -f "${COMBINED}.bwt.2bit.64" ]; then
  bwa-mem2 index -p "${COMBINED}" "${COMBINED}"
  echo "BWA-MEM2 index built."
else
  echo "BWA-MEM2 index already exists."
fi

echo "=== Step 5: Build minimap2 index for Nanopore (optional, ~5 min) ==="

MMAP_IDX="${OUTDIR}/GRCh38_decoy.mmi"
if [ ! -f "${MMAP_IDX}" ]; then
  minimap2 -x map-ont -d "${MMAP_IDX}" "${COMBINED}"
  echo "minimap2 index built: ${MMAP_IDX}"
else
  echo "minimap2 index already exists."
fi

echo ""
echo "=== Done ==="
echo "Human reference: ${COMBINED}"
echo ""
echo "Set environment variable:"
echo "  export HUMAN_REF=$(pwd)/${COMBINED}"
