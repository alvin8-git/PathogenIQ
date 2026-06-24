from __future__ import annotations

import csv
import io
import shutil
import subprocess
from dataclasses import dataclass, field
from pathlib import Path

from .config import PipelineConfig


# A contig is attributed to every reference whose aligned-bp is within this
# fraction of its best match. Sibling species that share a gene's region (E. coli
# / Shigella) co-attribute it; a contig carrying organism-specific sequence aligns
# more fully to one genome and stays single. ponytail: 0.95 margin on aligned bp —
# loosen if real shared genes are being dropped from a sibling.
_COMAP_MARGIN = 0.95


@dataclass
class AMRHit:
    gene: str
    drug_class: str
    identity_pct: float
    coverage_pct: float
    organism_match: str  # single best-matching organism or "unknown"
    database: str
    organism_matches: list[str] = field(default_factory=list)  # all co-mapped orgs


@dataclass
class VirulenceHit:
    gene: str
    factor: str  # VFDB PRODUCT — the virulence factor description
    identity_pct: float
    coverage_pct: float
    organism_match: str  # single best-matching organism or "unknown"
    database: str
    organism_matches: list[str] = field(default_factory=list)  # all co-mapped orgs


def _attribute(
    sequence: str,
    organism_names: list[str],
    contig_to_org: dict[str, list[str]] | None = None,
) -> list[str]:
    """Resolve an ABRicate SEQUENCE header to the organism(s) it belongs to,
    best-first.

    With contig-based input the header is an assembly contig id (e.g. ``k141_42``)
    that carries no organism name, so substring matching fails. When a
    ``contig_to_org`` map (contig -> co-mapped organisms, built by aligning contigs
    to the targeted genomes) is supplied, use it. Fall back to underscore-normalised
    substring matching (legacy read-based path and tests). ``["unknown"]`` if
    nothing matches."""
    if contig_to_org:
        orgs = contig_to_org.get(sequence)
        if orgs:
            return orgs
    seq_norm = sequence.replace("_", " ").lower()
    for org in organism_names:
        if org.lower() in seq_norm or org.replace(" ", "_").lower() in seq_norm:
            return [org]
    return ["unknown"]


def map_contigs_to_organisms(
    cfg: PipelineConfig, contigs: Path | None, hits,
) -> dict[str, list[str]]:
    """Map each assembled contig to the targeted organism(s) it aligns to, so a
    contig-based AMR/VFDB hit can be attributed to a finding (abricate only knows
    the contig id, not the organism).

    For near-identical siblings (E. coli / Shigella) the co-assembly collapses
    into one shared contig set, so forcing a single winner buried every Shigella
    gene under E. coli. Instead each contig is attributed to **all** references
    whose total aligned bases are within ``_COMAP_MARGIN`` of its best match:
    a gene's region present in several siblings co-attributes to each, while an
    organism-specific region (more aligned bases to one genome) stays single.
    Returned lists are best-first. Non-blocking: ``{}`` if minimap2 or contigs are
    missing.

    ``hits`` are the sketch hits (``.name`` + ``.genome_path``), so this aligns the
    contigs against just the handful of screened genomes, not the whole DB."""
    if contigs is None or not shutil.which("minimap2"):
        return {}
    out = cfg.output_dir / "amr"
    out.mkdir(parents=True, exist_ok=True)
    # contig -> {organism: total aligned bp}
    aligned: dict[str, dict[str, int]] = {}
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
            per_org = aligned.setdefault(contig, {})
            per_org[hit.name] = per_org.get(hit.name, 0) + matches
    mapping: dict[str, list[str]] = {}
    for contig, per_org in aligned.items():
        best = max(per_org.values())
        if best <= 0:
            continue
        cutoff = best * _COMAP_MARGIN
        orgs = sorted((o for o, bp in per_org.items() if bp >= cutoff),
                      key=lambda o: per_org[o], reverse=True)
        mapping[contig] = orgs
    return mapping


def _parse_abricate_tsv(
    tsv_text: str,
    organism_names: list[str],
    contig_to_org: dict[str, list[str]] | None = None,
) -> list[AMRHit]:
    """Parse ABRicate TSV output into AMRHit objects."""
    hits: list[AMRHit] = []
    reader = csv.DictReader(io.StringIO(tsv_text), delimiter="\t")
    for row in reader:
        if row.get("#FILE", "").startswith("#"):
            continue
        orgs = _attribute(row.get("SEQUENCE", ""), organism_names, contig_to_org)
        hits.append(
            AMRHit(
                gene=row.get("GENE", ""),
                drug_class=row.get("RESISTANCE", ""),
                identity_pct=float(row.get("%IDENTITY", 0)),
                coverage_pct=float(row.get("%COVERAGE", 0)),
                organism_match=orgs[0],
                database=row.get("DATABASE", ""),
                organism_matches=orgs,
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
    contig_to_org: dict[str, list[str]] | None = None,
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
    contig_to_org: dict[str, list[str]] | None = None,
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
        orgs = _attribute(row.get("SEQUENCE", ""), organism_names, contig_to_org)
        hits.append(
            VirulenceHit(
                gene=row.get("GENE", ""),
                factor=row.get("PRODUCT", "") or row.get("RESISTANCE", ""),
                identity_pct=float(row.get("%IDENTITY", 0)),
                coverage_pct=float(row.get("%COVERAGE", 0)),
                organism_match=orgs[0],
                database=row.get("DATABASE", ""),
                organism_matches=orgs,
            )
        )
    return hits
