"""Open-world novelty trigger for airborne-pathogen surveillance.

The targeted arm (sketch → align → EM) is closed-world: it can only see the ~110
genomes in the tier-1 DB and is blind to anything else — including novel or simply
uncatalogued pathogens, which is exactly what air surveillance must not miss. This
module quantifies the "dark matter": it classifies the non-host / non-PhiX reads
against a BROAD database (Kraken2 Standard, NOT the tier-1 DB — a narrow DB would
mark everything off-target as unclassified and falsely inflate novelty) and reports
the unclassified fraction. A high unclassified fraction is the cheap gate that says
"something here has no reference" and should trigger the expensive discovery arms
(assembly/MAG, viral).

    non-host reads ──► kraken2 (Standard DB) ──► report ──► unclassified fraction
                                                              │
                                              flagged if >= threshold ──► run discovery arms

Non-blocking like the AMR/assembly arms: a missing kraken2 binary or DB degrades to
``None`` rather than failing the run.
"""
from __future__ import annotations

import os
import shutil
import subprocess
from dataclasses import dataclass, field
from pathlib import Path

from .config import PipelineConfig

# Fraction of reads with NO classification (against a broad DB) above which the
# sample is flagged as harbouring novel/uncatalogued content worth assembling.
# ponytail: crude single global threshold — refine with rank-resolution (reads
# stuck above species level) if the bare unclassified fraction proves too blunt.
_DEFAULT_FLAG_THRESHOLD = 0.5


@dataclass
class NoveltyResult:
    total_reads: int
    classified_reads: int
    unclassified_reads: int
    unclassified_fraction: float
    n_species: int                              # distinct species-rank taxa seen
    top_taxa: list[tuple[str, int]] = field(default_factory=list)  # (species, reads)
    flagged: bool = False                       # unclassified_fraction >= threshold


def parse_kraken_report(text: str, *, flag_threshold: float = _DEFAULT_FLAG_THRESHOLD) -> NoveltyResult:
    """Parse a Kraken2 ``--report`` table into a NoveltyResult.

    Columns (tab-separated): pct, clade_reads, taxon_reads, rank_code, taxid,
    name. The ``U`` row (taxid 0) is unclassified; the ``root`` row (``R``) holds
    all classified reads; ``S`` rows are species.
    """
    unclassified = classified = 0
    classified_from_sum = 0
    species: list[tuple[str, int]] = []
    for line in text.splitlines():
        cols = line.split("\t")
        if len(cols) < 6:
            continue
        try:
            clade_reads = int(cols[1])
            taxon_reads = int(cols[2])
        except ValueError:
            continue
        rank = cols[3].strip()
        name = cols[5].strip()
        if rank == "U":
            unclassified = clade_reads
        else:
            classified_from_sum += taxon_reads
            if name == "root":
                classified = clade_reads
        if rank == "S" and taxon_reads > 0:
            species.append((name, taxon_reads))
    if classified == 0:            # no explicit root row -> fall back to the sum
        classified = classified_from_sum
    total = classified + unclassified
    fraction = unclassified / total if total else 0.0
    species.sort(key=lambda x: x[1], reverse=True)
    return NoveltyResult(
        total_reads=total,
        classified_reads=classified,
        unclassified_reads=unclassified,
        unclassified_fraction=fraction,
        n_species=len(species),
        top_taxa=species[:5],
        flagged=fraction >= flag_threshold,
    )


def _resolve_kraken_db(path: Path) -> Path | None:
    """A valid Kraken2 DB is a directory containing ``taxo.k2d`` (+ hash/opts).
    Kraken DBs ship as tarballs that extract into a named subdir, so if ``path``
    itself isn't the DB, descend one level to find the subdir that is."""
    if not path.exists():
        return None
    if (path / "taxo.k2d").exists():
        return path
    for sub in sorted(p for p in path.iterdir() if p.is_dir()):
        if (sub / "taxo.k2d").exists():
            return sub
    return None


def kraken2_db_path() -> Path | None:
    """Resolve the broad Kraken2 DB: ``$KRAKEN2_DB`` if set, else the conventional
    ``databases/kraken2`` (the Standard build, NOT ``kraken2_tier1``). Descends into
    a single extracted subdir if needed. None if no valid DB (with ``taxo.k2d``) is
    found."""
    env = os.environ.get("KRAKEN2_DB")
    return _resolve_kraken_db(Path(env) if env else Path("databases/kraken2"))


def run_kraken2(reads: Path, db: Path, threads: int, out_dir: Path) -> Path | None:
    """Classify ``reads`` against the broad DB; return the report path. None if
    kraken2 or the DB is missing, or the run fails. Per-read output is discarded
    (only the aggregate report is needed)."""
    if not shutil.which("kraken2") or not Path(db).exists():
        return None
    out_dir.mkdir(parents=True, exist_ok=True)
    report = out_dir / "kraken2.report"
    cmd = ["kraken2", "--db", str(db), "--threads", str(threads),
           "--report", str(report), "--output", os.devnull]
    if str(reads).endswith(".gz"):
        cmd.append("--gzip-compressed")
    cmd.append(str(reads))
    try:
        subprocess.run(cmd, capture_output=True, check=True)
    except subprocess.CalledProcessError:
        return None
    return report if report.exists() else None


def assess_novelty(
    cfg: PipelineConfig,
    reads: Path,
    *,
    db: Path | None = None,
    flag_threshold: float = _DEFAULT_FLAG_THRESHOLD,
) -> NoveltyResult | None:
    """Open-world novelty trigger: classify non-host reads against the broad DB
    and quantify the unclassified fraction. None (non-blocking) when kraken2 or
    the DB is unavailable."""
    db = db or kraken2_db_path()
    if db is None:
        return None
    report = run_kraken2(reads, db, cfg.threads, cfg.output_dir / "novelty")
    if report is None:
        return None
    return parse_kraken_report(report.read_text(), flag_threshold=flag_threshold)
