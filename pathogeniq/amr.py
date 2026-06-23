from __future__ import annotations

import csv
import gzip
import io
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path

from .config import PipelineConfig


@dataclass
class AMRHit:
    gene: str
    drug_class: str
    identity_pct: float
    coverage_pct: float
    organism_match: str  # matched organism name or "unknown"
    database: str


@dataclass
class VirulenceHit:
    gene: str
    factor: str  # VFDB PRODUCT — the virulence factor description
    identity_pct: float
    coverage_pct: float
    organism_match: str  # matched organism name or "unknown"
    database: str


def _fastq_to_fasta(fastq_path: Path, fasta_path: Path) -> None:
    """Convert FASTQ to FASTA for ABRicate input.

    Handles gzipped input (the pipeline hands AMR the gzipped non-host reads) and
    decodes leniently (``errors="replace"``) so a stray non-ASCII byte degrades a
    single base rather than aborting the whole run."""
    opener = gzip.open if str(fastq_path).endswith(".gz") else open
    with opener(fastq_path, "rt", encoding="ascii", errors="replace") as fq, \
            open(fasta_path, "w") as fa:
        while True:
            header = fq.readline()
            if not header:
                break
            seq = fq.readline().strip()
            fq.readline()  # +
            fq.readline()  # quality
            fa.write(f">{header.lstrip('@').strip()}\n{seq}\n")


def _match_organism(sequence: str, organism_names: list[str]) -> str:
    """Match an ABRicate sequence header to a known organism by substring
    (underscore-normalised); 'unknown' if none match."""
    seq_norm = sequence.replace("_", " ").lower()
    for org in organism_names:
        if org.lower() in seq_norm or org.replace(" ", "_").lower() in seq_norm:
            return org
    return "unknown"


def _parse_abricate_tsv(tsv_text: str, organism_names: list[str]) -> list[AMRHit]:
    """Parse ABRicate TSV output into AMRHit objects."""
    hits: list[AMRHit] = []
    reader = csv.DictReader(io.StringIO(tsv_text), delimiter="\t")
    for row in reader:
        if row.get("#FILE", "").startswith("#"):
            continue
        hits.append(
            AMRHit(
                gene=row.get("GENE", ""),
                drug_class=row.get("RESISTANCE", ""),
                identity_pct=float(row.get("%IDENTITY", 0)),
                coverage_pct=float(row.get("%COVERAGE", 0)),
                organism_match=_match_organism(row.get("SEQUENCE", ""), organism_names),
                database=row.get("DATABASE", ""),
            )
        )
    return hits


def _run_abricate(
    cfg: PipelineConfig, reads_path: Path, db: str, min_identity: float, min_coverage: float,
) -> str | None:
    """Convert reads to FASTA and run ABRicate against ``db``. Returns the TSV
    stdout, or None if ABRicate is not on PATH or the run failed (non-blocking)."""
    if not shutil.which("abricate"):
        return None
    # Non-blocking: any failure (bad input, IO, abricate crash) skips the screen
    # rather than aborting the pipeline — AMR/virulence is an overlay, not core.
    try:
        fasta_path = cfg.output_dir / "amr_reads.fa"
        _fastq_to_fasta(reads_path, fasta_path)
        result = subprocess.run(
            ["abricate", "--db", db, "--minid", str(min_identity),
             "--mincov", str(min_coverage), str(fasta_path)],
            capture_output=True, encoding="utf-8", errors="replace",
        )
    except (OSError, ValueError, subprocess.SubprocessError):
        return None
    return result.stdout if result.returncode == 0 else None


def run_amr_screen(
    cfg: PipelineConfig,
    reads_path: Path,
    organism_names: list[str],
    db: str = "card",
    min_identity: float = 90.0,
    min_coverage: float = 80.0,
) -> list[AMRHit]:
    """Run ABRicate (default CARD) for resistance genes. Empty list if ABRicate
    not on PATH (non-blocking)."""
    tsv = _run_abricate(cfg, reads_path, db, min_identity, min_coverage)
    return _parse_abricate_tsv(tsv, organism_names) if tsv is not None else []


def run_virulence_screen(
    cfg: PipelineConfig,
    reads_path: Path,
    organism_names: list[str],
    db: str = "vfdb",
    min_identity: float = 90.0,
    min_coverage: float = 80.0,
) -> list[VirulenceHit]:
    """Run ABRicate against VFDB (virulence factor database) alongside the AMR
    screen. Same machinery as run_amr_screen, but keeps the PRODUCT column (the
    virulence factor description) rather than RESISTANCE. Non-blocking."""
    tsv = _run_abricate(cfg, reads_path, db, min_identity, min_coverage)
    if tsv is None:
        return []
    hits: list[VirulenceHit] = []
    reader = csv.DictReader(io.StringIO(tsv), delimiter="\t")
    for row in reader:
        if row.get("#FILE", "").startswith("#"):
            continue
        hits.append(
            VirulenceHit(
                gene=row.get("GENE", ""),
                factor=row.get("PRODUCT", "") or row.get("RESISTANCE", ""),
                identity_pct=float(row.get("%IDENTITY", 0)),
                coverage_pct=float(row.get("%COVERAGE", 0)),
                organism_match=_match_organism(row.get("SEQUENCE", ""), organism_names),
                database=row.get("DATABASE", ""),
            )
        )
    return hits
