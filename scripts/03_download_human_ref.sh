#!/usr/bin/env bash
# Download GRCh38.p14 (single-file) + PhiX decoy and index for BWA-MEM2 / minimap2.
#
# Runtime: ~20-40 min download + ~60 min BWA-MEM2 index
# Output:  databases/human_ref/GRCh38_decoy.fa + index files
#
# RAM requirement: BWA-MEM2 index needs ~64 GB RAM.
#                  If RAM < 64 GB, only minimap2 index is built (usable for both SR and LR).
#
# Prerequisites (install before running):
#   conda install -c bioconda bwa-mem2 minimap2 samtools

set -euo pipefail

OUTDIR="databases/human_ref"
THREADS=${1:-16}
mkdir -p "${OUTDIR}"

# ── Check required tools ───────────────────────────────────────────────────
for tool in bwa-mem2 minimap2 samtools; do
  if ! command -v "${tool}" &>/dev/null; then
    echo "ERROR: ${tool} not found."
    echo "Install with: conda install -c bioconda bwa-mem2 minimap2 samtools"
    exit 1
  fi
done

COMBINED="${OUTDIR}/GRCh38_decoy.fa"

echo "=== Step 1: Download GRCh38.p14 primary assembly (single file, ~900 MB) ==="

GRCh38_GZ="${OUTDIR}/GRCh38.p14_genomic.fna.gz"
if [ ! -f "${GRCh38_GZ}" ]; then
  curl -C - --retry 5 --retry-delay 10 -L \
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz" \
    -o "${GRCh38_GZ}"
  echo "GRCh38.p14 downloaded."
else
  echo "GRCh38.p14 already downloaded, skipping."
fi

echo "=== Step 2: Download PhiX control genome (spike-in contaminant) ==="

PHIX_GZ="${OUTDIR}/phix.fna.gz"
if [ ! -f "${PHIX_GZ}" ]; then
  curl -C - --retry 5 -L \
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/819/615/GCF_000819615.1_ViralProj14015/GCF_000819615.1_ViralProj14015_genomic.fna.gz" \
    -o "${PHIX_GZ}"
  echo "PhiX downloaded."
else
  echo "PhiX already downloaded."
fi

echo "=== Step 3: Concatenate into single reference ==="

if [ ! -f "${COMBINED}" ]; then
  gunzip -c "${GRCh38_GZ}" > "${COMBINED}"
  gunzip -c "${PHIX_GZ}" >> "${COMBINED}"
  echo "Combined reference created: ${COMBINED} ($(du -sh "${COMBINED}" | cut -f1))"
else
  echo "Combined reference already exists: ${COMBINED}"
fi

echo "=== Step 4: Build BWA-MEM2 index (requires ~64 GB RAM, ~60 min) ==="

BWA_IDX="${COMBINED}.bwt.2bit.64"
if [ ! -f "${BWA_IDX}" ]; then
  AVAIL_GB=$(awk '/MemAvailable/ {printf "%.0f", $2/1024/1024}' /proc/meminfo)
  echo "Available RAM: ~${AVAIL_GB} GB"
  if [ "${AVAIL_GB}" -ge 60 ]; then
    bwa-mem2 index "${COMBINED}"
    echo "BWA-MEM2 index built."
  else
    echo "WARNING: Only ${AVAIL_GB} GB RAM available — skipping BWA-MEM2 index."
    echo "         Pipeline will fall back to minimap2 for short-read host removal."
  fi
else
  echo "BWA-MEM2 index already exists."
fi

echo "=== Step 5: Build minimap2 index (~5 min, works for both SR and LR) ==="

MMAP_IDX="${OUTDIR}/GRCh38_decoy.sr.mmi"
if [ ! -f "${MMAP_IDX}" ]; then
  minimap2 -x sr -t "${THREADS}" -d "${MMAP_IDX}" "${COMBINED}"
  echo "minimap2 SR index built: ${MMAP_IDX}"
else
  echo "minimap2 SR index already exists."
fi

MMAP_LR_IDX="${OUTDIR}/GRCh38_decoy.ont.mmi"
if [ ! -f "${MMAP_LR_IDX}" ]; then
  minimap2 -x map-ont -t "${THREADS}" -d "${MMAP_LR_IDX}" "${COMBINED}"
  echo "minimap2 ONT index built: ${MMAP_LR_IDX}"
else
  echo "minimap2 ONT index already exists."
fi

echo ""
echo "=== Done ==="
echo "Human reference : ${COMBINED}"
echo "BWA-MEM2 index  : ${BWA_IDX}"
echo "minimap2 SR idx : ${MMAP_IDX}"
echo "minimap2 ONT idx: ${MMAP_LR_IDX}"
echo ""
echo "Set environment variable:"
echo "  export HUMAN_REF=$(pwd)/${COMBINED}"
