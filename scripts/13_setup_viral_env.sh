#!/usr/bin/env bash
# Set up the viral arm (--viral) tools. geNomad + CheckV require numpy<2 and
# conflict with the pathogeniq env's numpy-2.x stack (scipy/pandas), so they live
# in their OWN conda env and their CLIs are symlinked into the pathogeniq env —
# the shebangs point at the genomad env's python, so they run self-consistently
# regardless of which env is active (they are external binaries, like abricate).
#
# Idempotent. Run headless/logout-surviving for the ~1.4 GB geNomad DB:
#   setsid nohup bash scripts/13_setup_viral_env.sh > setup.log 2>&1 < /dev/null &
set -uo pipefail
source "$(conda info --base)/etc/profile.d/conda.sh"
BASE=$(conda info --base)
PENV="$BASE/envs/pathogeniq/bin"
GENV="$BASE/envs/genomad/bin"
DB=databases     # run from the repo root

echo "VIRAL_ENV_START $(date)"

echo "## [1/3] isolated 'genomad' env (genomad + checkv, self-consistent numpy<2) ##"
if [ -x "$GENV/genomad" ] && [ -x "$GENV/checkv" ]; then
  echo "  already present, skipping create"
else
  mamba create -y -n genomad -c conda-forge -c bioconda genomad checkv \
    || { echo "ENV CREATE FAILED"; exit 1; }
fi

echo "## [2/3] geNomad DB -> $DB/genomad_db ##"
if [ -d "$DB/genomad_db" ] && [ -n "$(ls -A "$DB/genomad_db" 2>/dev/null)" ]; then
  echo "  already present, skipping"
else
  conda run -n genomad genomad download-database "$DB" || echo "GENOMAD DB DOWNLOAD FAILED"
fi

echo "## [2b] CheckV DB -> $DB/checkv_db (versioned dir symlinked) ##"
if [ -e "$DB/checkv_db" ]; then
  echo "  already present, skipping"
else
  conda run -n genomad checkv download_database "$DB" || echo "CHECKV DB DOWNLOAD FAILED"
  ver=$(ls -d "$DB"/checkv-db-* 2>/dev/null | sort | tail -1)
  [ -n "$ver" ] && ln -sfn "$(basename "$ver")" "$DB/checkv_db"
fi

echo "## [3/3] symlink viral CLIs into the pathogeniq env ##"
for t in genomad checkv; do
  [ -x "$GENV/$t" ] && ln -sfn "$GENV/$t" "$PENV/$t" && echo "  linked $t" || echo "  MISSING $GENV/$t"
done

echo "## verify (from pathogeniq env) ##"
conda activate pathogeniq
genomad --version >/dev/null 2>&1 && echo "  OK genomad" || echo "  genomad broken"
checkv -h          >/dev/null 2>&1 && echo "  OK checkv"  || echo "  checkv broken"
[ -d "$DB/genomad_db" ] && echo "  OK genomad_db" || echo "  MISSING genomad_db"
[ -e "$DB/checkv_db" ]  && echo "  OK checkv_db"  || echo "  MISSING checkv_db"
echo "### VIRAL_ENV_DONE $(date) ###"
