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


def _match_organism(
    sequence: str,
    organism_names: list[str],
    contig_to_org: dict[str, str] | None = None,
) -> str:
    """Resolve an ABRicate SEQUENCE header to a known organism.

    With contig-based input the header is an assembly contig id (e.g. ``k141_42``)
    that carries no organism name, so substring matching always fails. When a
    ``contig_to_org`` map (built by aligning contigs to the targeted genomes) is
    supplied, look the contig up there first. Fall back to underscore-normalised
    substring matching (covers the legacy read-based path and tests). 'unknown'
    if nothing matches."""
    if contig_to_org:
        org = contig_to_org.get(sequence)
        if org:
            return org
    seq_norm = sequence.replace("_", " ").lower()
    for org in organism_names:
        if org.lower() in seq_norm or org.replace(" ", "_").lower() in seq_norm:
            return org
    return "unknown"


def map_contigs_to_organisms(cfg: PipelineConfig, contigs: Path | None, hits) -> dict[str, str]:
    """Map each assembled contig to the targeted organism whose reference genome
    it best aligns to (most aligned bases via minimap2). This is what lets a
    contig-based AMR/VFDB hit be attributed to a finding — abricate only knows the
    contig id, not which organism it came from.

    ``hits`` are the sketch hits (``.name`` + ``.genome_path``), so this aligns the
    contigs against just the handful of screened genomes, not the whole DB. Contigs
    with no alignment are simply absent (-> 'unknown'). Non-blocking: returns ``{}``
    if minimap2 or the contigs are missing.

    ponytail: best-single-organism by aligned bp — a chimeric contig is assigned
    whole to its dominant source. Per-region binning only if mixed contigs matter."""
    if contigs is None or not shutil.which("minimap2"):
        return {}
    out = cfg.output_dir / "amr"
    out.mkdir(parents=True, exist_ok=True)
    best: dict[str, tuple[int, str]] = {}   # contig -> (aligned_bp, organism)
    for idx, hit in enumerate(hits):
        genome = getattr(hit, "genome_path", None)
        if genome is None or not Path(genome).exists():
            continue
        paf = out / f"contigmap_{idx}.paf"
        try:
            subprocess.run(
                ["minimap2", "-x", "asm5", "-t", str(cfg.threads),
                 "-o", str(paf), str(genome), str(contigs)],
                capture_output=True, check=True,
            )
        except (OSError, subprocess.SubprocessError):
            continue
        if not paf.exists():
            continue
        for line in paf.read_text().splitlines():
            cols = line.split("\t")
            if len(cols) < 10:
                continue
            contig = cols[0]
            try:
                matches = int(cols[9])   # PAF col 10: residue matches
            except ValueError:
                continue
            if matches > best.get(contig, (0, ""))[0]:
                best[contig] = (matches, hit.name)
    return {contig: org for contig, (_, org) in best.items()}


def _parse_abricate_tsv(
    tsv_text: str,
    organism_names: list[str],
    contig_to_org: dict[str, str] | None = None,
) -> list[AMRHit]:
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
                organism_match=_match_organism(
                    row.get("SEQUENCE", ""), organism_names, contig_to_org),
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
    contig_to_org: dict[str, str] | None = None,
) -> list[AMRHit]:
    """Run ABRicate (default CARD) for resistance genes on assembled ``contigs``.
    Empty list if there are no contigs or ABRicate is absent (non-blocking)."""
    tsv = _run_abricate(cfg, contigs, db, min_identity, min_coverage)
    return _parse_abricate_tsv(tsv, organism_names, contig_to_org) if tsv is not None else []


def run_virulence_screen(
    cfg: PipelineConfig,
    contigs: Path | None,
    organism_names: list[str],
    db: str = "vfdb",
    min_identity: float = 90.0,
    min_coverage: float = 80.0,
    contig_to_org: dict[str, str] | None = None,
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
                organism_match=_match_organism(
                    row.get("SEQUENCE", ""), organism_names, contig_to_org),
                database=row.get("DATABASE", ""),
            )
        )
    return hits
