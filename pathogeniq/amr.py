from __future__ import annotations

import csv
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
    cfg: PipelineConfig, contigs: Path | None, db: str, min_identity: float, min_coverage: float,
) -> str | None:
    """Run ABRicate against ``db`` on an assembled-contig FASTA. Returns the TSV
    stdout, or None if there are no contigs, abricate is absent, or the run fails.

    Contigs (not raw reads) are the right input: ABRicate is a BLAST-over-assembly
    tool. Screening ~tens of thousands of contigs instead of millions of reads is
    ~100x less work and yields clean full-length gene hits rather than fragmented
    per-read partials. Non-blocking — any failure skips the overlay."""
    if contigs is None or not shutil.which("abricate"):
        return None
    try:
        result = subprocess.run(
            ["abricate", "--db", db, "--minid", str(min_identity),
             "--mincov", str(min_coverage), str(contigs)],
            capture_output=True, encoding="utf-8", errors="replace",
        )
    except (OSError, subprocess.SubprocessError):
        return None
    return result.stdout if result.returncode == 0 else None


def run_amr_screen(
    cfg: PipelineConfig,
    contigs: Path | None,
    organism_names: list[str],
    db: str = "card",
    min_identity: float = 90.0,
    min_coverage: float = 80.0,
) -> list[AMRHit]:
    """Run ABRicate (default CARD) for resistance genes on assembled ``contigs``.
    Empty list if there are no contigs or ABRicate is absent (non-blocking)."""
    tsv = _run_abricate(cfg, contigs, db, min_identity, min_coverage)
    return _parse_abricate_tsv(tsv, organism_names) if tsv is not None else []


def run_virulence_screen(
    cfg: PipelineConfig,
    contigs: Path | None,
    organism_names: list[str],
    db: str = "vfdb",
    min_identity: float = 90.0,
    min_coverage: float = 80.0,
) -> list[VirulenceHit]:
    """Run ABRicate against VFDB (virulence factor database) on assembled ``contigs``,
    alongside the AMR screen. Same machinery as run_amr_screen, but keeps the PRODUCT
    column (the virulence factor description) rather than RESISTANCE. Non-blocking."""
    tsv = _run_abricate(cfg, contigs, db, min_identity, min_coverage)
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
