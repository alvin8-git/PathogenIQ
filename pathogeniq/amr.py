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


def _fastq_to_fasta(fastq_path: Path, fasta_path: Path) -> None:
    """Convert FASTQ to FASTA for ABRicate input."""
    with open(fastq_path, "rb") as fq, open(fasta_path, "w") as fa:
        while True:
            header = fq.readline()
            if not header:
                break
            seq = fq.readline().decode("ascii").strip()
            fq.readline()  # +
            fq.readline()  # quality
            fa.write(f">{header.decode('ascii').lstrip('@').strip()}\n{seq}\n")


def _parse_abricate_tsv(tsv_text: str, organism_names: list[str]) -> list[AMRHit]:
    """Parse ABRicate TSV output into AMRHit objects."""
    hits: list[AMRHit] = []
    reader = csv.DictReader(io.StringIO(tsv_text), delimiter="\t")
    for row in reader:
        if row.get("#FILE", "").startswith("#"):
            continue
        sequence = row.get("SEQUENCE", "")
        # Match sequence header to a known organism by substring (underscore-normalised)
        seq_norm = sequence.replace("_", " ").lower()
        matched = "unknown"
        for org in organism_names:
            if org.lower() in seq_norm or org.replace(" ", "_").lower() in seq_norm:
                matched = org
                break
        hits.append(
            AMRHit(
                gene=row.get("GENE", ""),
                drug_class=row.get("RESISTANCE", ""),
                identity_pct=float(row.get("%IDENTITY", 0)),
                coverage_pct=float(row.get("%COVERAGE", 0)),
                organism_match=matched,
                database=row.get("DATABASE", ""),
            )
        )
    return hits


def run_amr_screen(
    cfg: PipelineConfig,
    reads_path: Path,
    organism_names: list[str],
    db: str = "card",
    min_identity: float = 90.0,
    min_coverage: float = 80.0,
) -> list[AMRHit]:
    """Run ABRicate on non-human reads. Returns empty list if ABRicate not on PATH."""
    if not shutil.which("abricate"):
        return []

    # Convert FASTQ → FASTA in a temp file for abricate
    fasta_path = cfg.output_dir / "amr_reads.fa"
    _fastq_to_fasta(reads_path, fasta_path)

    result = subprocess.run(
        [
            "abricate",
            "--db", db,
            "--minid", str(min_identity),
            "--mincov", str(min_coverage),
            str(fasta_path),
        ],
        capture_output=True,
        encoding="utf-8",
        errors="replace",
    )
    if result.returncode != 0:
        return []

    return _parse_abricate_tsv(result.stdout, organism_names)
