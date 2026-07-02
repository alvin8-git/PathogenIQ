#!/usr/bin/env bash
# Re-run the open-world assembly arm (--assemble) on the 5 environmental aircraft
# filters, now that the MAG toolchain (metabat2 + CheckM + GTDB-Tk) is installed.
# Recovers MAGs (megahit -> metabat2 -> CheckM QC -> GTDB-Tk taxonomy) and runs R4
# pathogenicity triage on them. Goal: recover the environmental dominants
# (Sphingomonas / Methylobacterium) the read-based tier-1 screen can't see.
#
# Sequential (one sample at a time) so gtdbtk's pplacer has the box's RAM to itself.
# Non-blocking arms: a gtdbtk OOM/failure leaves MAGs recovered but unnamed.
set -uo pipefail

source /home/alvin/miniconda3/etc/profile.d/conda.sh
conda activate pathogeniq

REPO=/data/alvin/Metagenomics
DB=$REPO/databases/tier1/tier1_pathogens.sbt.zip
HOST=/data/alvin/ref/GRCh38/hg38.fa
READS=$REPO/databases/aircraft/air-aircraft-filter
OUT=/data/alvin/tmp/air_assemble
THREADS=12
FILTERS="SRR32514297 SRR32514308 SRR32514319 SRR32514330 SRR32514331"

# fail fast if the toolchain isn't actually wired in
for t in megahit metabat2 jgi_summarize_bam_contig_depths checkm gtdbtk; do
  command -v "$t" >/dev/null || { echo "FATAL: $t not on PATH — run scripts/15 first"; exit 1; }
done
echo "### toolchain OK: $(gtdbtk --version 2>&1 | head -1)"

mkdir -p "$OUT"
for acc in $FILTERS; do
  r1=$READS/${acc}_1.fastq.gz
  [ -f "$r1" ] || { echo "SKIP $acc: $r1 missing"; continue; }
  echo "### [$acc] assemble start $(date)"
  pathogeniq run \
    --input "$r1" --output "$OUT/$acc" \
    --db "$DB" --host-ref "$HOST" \
    --specimen air --read-type short \
    --assemble --threads "$THREADS" \
    > "$OUT/${acc}.log" 2>&1 \
    && echo "### [$acc] DONE $(date)" \
    || echo "### [$acc] FAILED (exit $?) — see $OUT/${acc}.log $(date)"
done
echo "### AIR_ASSEMBLE_COMPLETE $(date)"
