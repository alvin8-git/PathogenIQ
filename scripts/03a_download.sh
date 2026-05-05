#!/usr/bin/env bash
# Download GRCh38 + decoy sequences and index for BWA-MEM2 and minimap2.
# Optimized for your 4x RTX 4090 / 754GB RAM rig.

set -euo pipefail

OUTDIR="databases/human_ref"
THREADS=${1:-32} # Bumped default threads given your HPC specs
mkdir -p "${OUTDIR}"

echo "=== Step 1: Download GRCh38 Primary Assembly ==="

# Using the Analysis Set version of GRCh38 is generally cleaner for bioinformatics pipelines
# It handles the naming conventions (chr1 vs 1) consistently.
REF_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz"
COMBINED_GZ="${OUTDIR}/GRCh38_full.fna.gz"
COMBINED="${OUTDIR}/GRCh38_decoy.fa"

if [ ! -f "${COMBINED_GZ}" ]; then
  echo "Downloading GRCh38 analysis set..."
  curl -C - --retry 5 -L -o "${COMBINED_GZ}" "${REF_URL}"
fi

echo "=== Step 2: Download PhiX Control ==="
PHIX_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/819/615/GCF_000819615.1_ViralProj14015/GCF_000819615.1_ViralProj14015_genomic.fna.gz"
if [ ! -f "${OUTDIR}/phix.fna.gz" ]; then
  curl -C - --retry 5 -L -o "${OUTDIR}/phix.fna.gz" "${PHIX_URL}"
fi

echo "=== Step 3: Concatenate into single reference ==="
if [ ! -f "${COMBINED}" ]; then
  gunzip -c "${COMBINED_GZ}" > "${COMBINED}"
  gunzip -c "${OUTDIR}/phix.fna.gz" >> "${COMBINED}"
  echo "Combined reference created: ${COMBINED}"
else
  echo "Combined reference already exists."
fi

echo "=== Step 4: Build BWA-MEM2 index ==="
# Since you have 754GB RAM, BWA-MEM2 will run smoothly.
if [ ! -f "${COMBINED}.bwt.2bit.64" ]; then
  bwa-mem2 index -p "${COMBINED}" "${COMBINED}"
  echo "BWA-MEM2 index built."
fi

echo "=== Step 5: Build minimap2 index for Nanopore ==="
# Since you're targeting Oxford Nanopore FAS roles, this is critical.
MMAP_IDX="${OUTDIR}/GRCh38_decoy.mmi"
if [ ! -f "${MMAP_IDX}" ]; then
  minimap2 -t "${THREADS}" -x map-ont -d "${MMAP_IDX}" "${COMBINED}"
  echo "minimap2 index built: ${MMAP_IDX}"
fi

echo -e "\n=== Done ===\nHuman reference: ${COMBINED}"
