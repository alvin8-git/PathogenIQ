#!/usr/bin/env bash
# Install the MAG-recovery toolchain (metabat2 + CheckM + GTDB-Tk) in ISOLATED
# conda envs and expose them to the core `pathogeniq` env via thin PATH wrappers.
# Same discipline as scripts/13 (viral): never install these into the core env
# (CheckM/GTDB-Tk pin old numpy/pplacer and would break the core SciPy stack).
#
# Two isolated envs so their solves don't fight (gtdbtk vs checkm numpy pins):
#   mag-bin : metabat2 (+ jgi_summarize_bam_contig_depths) + checkm-genome (CheckM v1)
#   gtdbtk  : gtdbtk
#
# DBs go on /data (NOT the root/home partition): CheckM ~1.4 GB, GTDB-Tk ~110 GB.
# Idempotent: skips env creation / downloads that are already complete, so it can
# be re-run after a partial failure. Downloads use the pathogeniq env's aria2c by
# ABSOLUTE path (multi-connection, resumable) and every stage VERIFIES its output
# and exits loudly on failure — no silent no-ops.
set -uo pipefail

PQ_BIN=/home/alvin/miniconda3/envs/pathogeniq/bin
ARIA2="$PQ_BIN/aria2c"
MAG_ENV=mag-bin
GTDB_ENV=gtdbtk
CHECKM_DB=/data/alvin/ref/checkm
GTDB_DB=/data/alvin/ref/gtdbtk
CHECKM_URL=https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz
GTDB_URL=https://data.gtdb.aau.ecogenomic.org/releases/release232/232.0/auxillary_files/gtdbtk_package/full_package/gtdbtk_r232_data.tar.gz

source /home/alvin/miniconda3/etc/profile.d/conda.sh
mkdir -p "$CHECKM_DB" "$GTDB_DB"
[ -x "$ARIA2" ] || { echo "FATAL: aria2c not at $ARIA2"; exit 1; }
DL() { "$ARIA2" -x8 -s8 -j1 --file-allocation=none --auto-file-renaming=false -c -d "$1" -o "$2" "$3"; }

echo "### [1/5] $MAG_ENV (metabat2 + checkm-genome) $(date)"
if [ ! -x /home/alvin/miniconda3/envs/$MAG_ENV/bin/metabat2 ]; then
  mamba create -y -n "$MAG_ENV" -c conda-forge -c bioconda metabat2 checkm-genome \
    || { echo "FATAL: $MAG_ENV solve failed"; exit 1; }
else echo "  exists, skip"; fi

echo "### [2/5] $GTDB_ENV (gtdbtk) $(date)"
if [ ! -x /home/alvin/miniconda3/envs/$GTDB_ENV/bin/gtdbtk ]; then
  mamba create -y -n "$GTDB_ENV" -c conda-forge -c bioconda gtdbtk \
    || { echo "FATAL: $GTDB_ENV solve failed"; exit 1; }
else echo "  exists, skip"; fi

MAG_BIN=/home/alvin/miniconda3/envs/$MAG_ENV/bin
GTDB_BIN=/home/alvin/miniconda3/envs/$GTDB_ENV/bin

echo "### [3/5] CheckM DB -> $CHECKM_DB $(date)"
if [ ! -f "$CHECKM_DB/taxon_marker_sets.tsv" ]; then
  DL "$CHECKM_DB" checkm_db.tar.gz "$CHECKM_URL" || { echo "FATAL: CheckM DB download failed"; exit 1; }
  tar xzf "$CHECKM_DB/checkm_db.tar.gz" -C "$CHECKM_DB" || { echo "FATAL: CheckM DB extract failed"; exit 1; }
  rm -f "$CHECKM_DB/checkm_db.tar.gz"
fi
[ -f "$CHECKM_DB/taxon_marker_sets.tsv" ] || { echo "FATAL: CheckM DB missing taxon_marker_sets.tsv after extract"; exit 1; }
conda activate "$MAG_ENV"; export CHECKM_DATA_PATH="$CHECKM_DB"
checkm data setRoot "$CHECKM_DB" || echo "WARN: setRoot non-zero (wrapper exports CHECKM_DATA_PATH anyway)"
conda deactivate
echo "  CheckM DB OK"

echo "### [4/5] GTDB-Tk DB -> $GTDB_DB (~66 GB download, ~110 GB extracted) $(date)"
# Download the tarball directly to /data (download-db.sh ignores GTDBTK_DATA_PATH
# and writes to the root partition single-threaded). Extract flat so $GTDB_DB
# directly holds taxonomy/ markers/ ... (what GTDBTK_DATA_PATH must point at).
if [ ! -d "$GTDB_DB/taxonomy" ]; then
  DL "$GTDB_DB" gtdbtk_data.tar.gz "$GTDB_URL" || { echo "FATAL: GTDB DB download failed"; exit 1; }
  tar xzf "$GTDB_DB/gtdbtk_data.tar.gz" -C "$GTDB_DB" --strip-components=1 \
    || { echo "FATAL: GTDB DB extract failed"; exit 1; }
  rm -f "$GTDB_DB/gtdbtk_data.tar.gz"
fi
[ -d "$GTDB_DB/taxonomy" ] && [ -d "$GTDB_DB/markers" ] \
  || { echo "FATAL: GTDB DB missing taxonomy/ or markers/ after extract"; exit 1; }
echo "  GTDB DB OK"

echo "### [5/5] wrappers -> $PQ_BIN $(date)"
cat > "$PQ_BIN/metabat2" <<EOF
#!/usr/bin/env bash
exec "$MAG_BIN/metabat2" "\$@"
EOF
cat > "$PQ_BIN/jgi_summarize_bam_contig_depths" <<EOF
#!/usr/bin/env bash
exec "$MAG_BIN/jgi_summarize_bam_contig_depths" "\$@"
EOF
cat > "$PQ_BIN/checkm" <<EOF
#!/usr/bin/env bash
export PATH="$MAG_BIN:\$PATH"
export CHECKM_DATA_PATH="$CHECKM_DB"
exec "$MAG_BIN/checkm" "\$@"
EOF
cat > "$PQ_BIN/gtdbtk" <<EOF
#!/usr/bin/env bash
export PATH="$GTDB_BIN:\$PATH"
export GTDBTK_DATA_PATH="$GTDB_DB"
exec "$GTDB_BIN/gtdbtk" "\$@"
EOF
chmod +x "$PQ_BIN/metabat2" "$PQ_BIN/jgi_summarize_bam_contig_depths" "$PQ_BIN/checkm" "$PQ_BIN/gtdbtk"

echo "### verify $(date)"
source /home/alvin/miniconda3/etc/profile.d/conda.sh; conda activate pathogeniq
FAIL=0
# Presence on PATH is the real check (the pipeline gates on shutil.which). Do NOT
# use --version: metabat2/jgi/checkm don't support it and exit non-zero even when
# perfectly installed.
for t in metabat2 jgi_summarize_bam_contig_depths checkm gtdbtk; do
  if command -v "$t" >/dev/null 2>&1; then echo "  $t -> OK ($(command -v "$t"))"; else echo "  $t -> MISSING"; FAIL=1; fi
done
[ "$FAIL" = 0 ] || { echo "FATAL: a wrapper is not on PATH"; exit 1; }
echo "### MAG_SETUP_COMPLETE $(date)"
