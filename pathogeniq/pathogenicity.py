"""R4: pathogenicity triage for reference-free hits (assembled MAGs).

A novel assembled organism only matters as a PATHOGEN if it (a) carries
pathogenicity markers — virulence factors (VFDB) or resistance genes (CARD) on its
own contigs — or (b) sits phylogenetically next to a known pathogen (its GTDB
genus/species matches one in the tier-1 clinical DB). Air is full of novel
ENVIRONMENTAL microbes (Sphingomonas, Methylobacterium, ...); this is the
discriminator that separates a novel pathogen from a novel soil/dust bug, so the
open-world arms don't bury the clinician in benign novelty.

    MAG FASTA ──► abricate vfdb + card  (marker genes on the genome)
             └──► GTDB lineage ──► genus/species match vs known pathogens
                                        │
                          assess ──► verdict: PATHOGEN_CANDIDATE (markers) /
                                     PATHOGEN_ADJACENT (lineage only) / ENVIRONMENTAL

Markers are the stronger signal: a novel genome carrying virulence/AMR genes is
concerning regardless of taxonomy. Lineage proximity is the fallback (DBs miss
genes). Non-blocking: abricate missing -> no markers, verdict falls back to phylo.
"""
from __future__ import annotations

import csv
import io
import json
import shutil
import subprocess
from dataclasses import dataclass, field
from pathlib import Path

from .config import PipelineConfig

PATHOGEN_CANDIDATE = "PATHOGEN_CANDIDATE"
PATHOGEN_ADJACENT = "PATHOGEN_ADJACENT"
ENVIRONMENTAL = "ENVIRONMENTAL"


@dataclass
class PathogenicityAssessment:
    name: str                                   # MAG bin_id
    taxonomy: str | None
    n_virulence: int
    n_amr: int
    phylo_match: str | None                     # matched known-pathogen taxon, or None
    verdict: str
    score: int
    virulence_genes: list[str] = field(default_factory=list)
    amr_genes: list[str] = field(default_factory=list)


def pathogen_taxa_from_name_map(name_map_path: Path) -> set[str]:
    """Build the known-pathogen reference set (lowercased genus + 'genus species'
    tokens) from the tier-1 DB's name_map.json. Empty set if the file is missing."""
    p = Path(name_map_path)
    if not p.exists():
        return set()
    taxa: set[str] = set()
    for name in json.loads(p.read_text()).values():
        words = str(name).split()
        if not words:
            continue
        taxa.add(words[0].lower())                       # genus
        if len(words) >= 2:
            taxa.add(f"{words[0]} {words[1]}".lower())   # genus species
    return taxa


def _lineage_taxa(taxonomy: str | None) -> tuple[str | None, str | None]:
    """Extract (genus, species) from a GTDB lineage string, lowercased and prefix-
    stripped. e.g. 'd__Bacteria;...;g__Bacillus;s__Bacillus cereus' -> (bacillus,
    bacillus cereus). Returns (None, None) if absent."""
    if not taxonomy:
        return None, None
    genus = species = None
    for tok in taxonomy.split(";"):
        tok = tok.strip()
        if tok.startswith("g__") and tok[3:]:
            genus = tok[3:].strip().lower()
        elif tok.startswith("s__") and tok[3:]:
            species = tok[3:].strip().lower()
    return genus, species


def phylo_match(taxonomy: str | None, pathogen_taxa: set[str]) -> str | None:
    """Return the matched known-pathogen taxon (species preferred over genus), or
    None. Species match is stronger evidence than a shared genus."""
    genus, species = _lineage_taxa(taxonomy)
    if species and species in pathogen_taxa:
        return species
    if genus and genus in pathogen_taxa:
        return genus
    return None


def assess_pathogenicity(
    taxonomy: str | None, n_virulence: int, n_amr: int, matched: str | None,
) -> tuple[str, int]:
    """Pure verdict + score. Markers dominate (virulence/AMR genes = concerning
    regardless of taxonomy); lineage proximity is the fallback signal."""
    markers = n_virulence + n_amr
    score = markers * 2 + (3 if matched else 0)
    if markers >= 1:
        return PATHOGEN_CANDIDATE, score
    if matched:
        return PATHOGEN_ADJACENT, score
    return ENVIRONMENTAL, score


def _abricate_genes(fasta: Path, db: str, min_identity: float, min_coverage: float) -> list[str]:
    """Run abricate on a contig FASTA (genome, not reads) -> list of GENE names.
    Empty if abricate is absent or the run fails (non-blocking)."""
    if not shutil.which("abricate"):
        return []
    result = subprocess.run(
        ["abricate", "--db", db, "--minid", str(min_identity),
         "--mincov", str(min_coverage), str(fasta)],
        capture_output=True, encoding="utf-8", errors="replace",
    )
    if result.returncode != 0:
        return []
    genes = []
    for row in csv.DictReader(io.StringIO(result.stdout), delimiter="\t"):
        gene = (row.get("GENE") or "").strip()
        if gene and not gene.startswith("#"):
            genes.append(gene)
    return genes


def triage_mags(
    cfg: PipelineConfig,
    mags: list,
    pathogen_taxa: set[str] | None = None,
    *,
    min_identity: float = 90.0,
    min_coverage: float = 80.0,
) -> list[PathogenicityAssessment]:
    """Triage each MAG: VFDB + CARD markers on its contigs + GTDB phylo-proximity to
    known pathogens -> a pathogen/environmental verdict. ``pathogen_taxa`` defaults
    to the set derived from the tier-1 DB's name_map.json."""
    if pathogen_taxa is None:
        pathogen_taxa = pathogen_taxa_from_name_map(cfg.db_tier1.parent / "name_map.json")
    out: list[PathogenicityAssessment] = []
    for m in mags:
        vir = _abricate_genes(m.fasta_path, "vfdb", min_identity, min_coverage)
        amr = _abricate_genes(m.fasta_path, "card", min_identity, min_coverage)
        matched = phylo_match(m.taxonomy, pathogen_taxa)
        verdict, score = assess_pathogenicity(m.taxonomy, len(vir), len(amr), matched)
        out.append(PathogenicityAssessment(
            name=m.bin_id, taxonomy=m.taxonomy, n_virulence=len(vir), n_amr=len(amr),
            phylo_match=matched, verdict=verdict, score=score,
            virulence_genes=vir, amr_genes=amr,
        ))
    out.sort(key=lambda a: a.score, reverse=True)
    return out
