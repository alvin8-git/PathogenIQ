import csv
import json
import math
from dataclasses import dataclass
from enum import Enum
from pathlib import Path

import numpy as np

from .background import BackgroundModel, is_background, is_dual_use
from .config import PipelineConfig, SpecimenType
from .contaminants import flag_contaminants
from .crossmap import deduplicate_closely_related
from .em import EMResult


class EvidenceGrade(str, Enum):
    A = "A"
    B = "B"
    C = "C"
    X = "X"


_MIN_READS: dict[SpecimenType, int] = {
    SpecimenType.BLOOD: 3,
    SpecimenType.CSF: 2,
    SpecimenType.BAL: 10,
    SpecimenType.TISSUE: 10,
    # ponytail: provisional air floor — low for sensitivity (air is ultra-low
    # biomass), since NTC subtraction + breadth gate are the real FP controls.
    # Calibrate on labeled air data (cf. the held-out floor work) before clinical use.
    SpecimenType.AIR: 5,
}

_MAX_CI_WIDTH_A = 0.15
# Breadth-of-coverage gate: reads covering far less of the genome than their depth
# predicts (Lander-Waterman) are clumped at a few loci -> PCR-dup / low-complexity /
# cross-mapping artifact, not a real organism. Below this observed/expected ratio the
# finding is non-gradeable (Grade X). ponytail: provisional cutoff — calibrate on
# labeled air data (cf. the read-floor work) before clinical use.
_MIN_BREADTH_RATIO = 0.25


@dataclass(frozen=True)
class GradingInput:
    """Normalized inputs the grading rule reads. Built from a ReportEntry today;
    the Kraken2 benchmark adapter (T10) will build the same struct from a
    classifier table. Keeping grading a pure function of this struct gives one
    grade path shared by the pipeline and the benchmark (DRY)."""
    read_count: int
    abundance: float
    ci_width: float
    contaminant_risk: bool
    specimen_type: SpecimenType
    tier: int = 3       # NTC background tier: 1 batch-matched, 2 pooled, 3 none
    crossmap: bool = False   # likely a cross-mapping artifact of a dominant relative
    breadth_ratio: float | None = None   # observed/expected coverage breadth; None = not computed


def _has_invalid_stats(g: GradingInput) -> bool:
    """Degenerate, non-gradeable inputs: NaN/inf abundance or CI width, or a
    negative read count (INT64_MIN underflow from astype(int) on a NaN
    abundance). A non-finite ci_width also catches a NaN/inf in either CI bound."""
    return (
        g.read_count < 0
        or not math.isfinite(g.abundance)
        or not math.isfinite(g.ci_width)
    )


def grade(g: GradingInput) -> EvidenceGrade:
    """Pure grading rule — the single source of truth for A/B/C/X, shared by
    ReportEntry.grade and the benchmark adapter.

    Grade A requires a batch-matched NTC (tier 1): without one, the highest a
    finding can reach is B, even with ample reads and a tight CI. This is the
    tier cap — you cannot claim top confidence without a same-run control.
    """
    if _has_invalid_stats(g):
        return EvidenceGrade.X
    if g.crossmap:
        return EvidenceGrade.X   # cross-mapping artifact of a dominant relative
    min_reads = _MIN_READS.get(g.specimen_type, 5)
    if g.read_count < min_reads:
        return EvidenceGrade.X
    if g.breadth_ratio is not None and g.breadth_ratio < _MIN_BREADTH_RATIO:
        return EvidenceGrade.X   # reads clumped far below expected breadth -> artifact
    if g.contaminant_risk:
        return EvidenceGrade.C
    if g.ci_width <= _MAX_CI_WIDTH_A and g.tier == 1:
        return EvidenceGrade.A
    return EvidenceGrade.B


# --- Open-world grading (R5): one A/B/C/X scale for reference-free hits ---------
# Reference-free hits (assembled MAGs, viral contigs) have no read-count / NTC /
# CI, so they cannot use grade(). They are graded instead on genome-quality
# evidence: completeness (CheckM bacterial / CheckV viral), contamination, and a
# supporting signal (viral hallmark genes / pathogenicity markers). They are
# CAPPED AT B — Grade A means NTC-controlled detection (a same-run control), which
# an assembled genome by definition lacks. Thresholds mirror the paper's MAG
# retention (>=50% completeness, <=10% contamination) and high-quality (>=90%/<=5%).
_OW_COMPLETE_MIN = 50.0   # below this = fragmentary, non-gradeable
_OW_COMPLETE_B = 70.0     # at/above this = a solid genome -> B (the open-world ceiling)
_OW_MAX_CONTAM = 10.0     # above this = likely chimeric bin -> X


@dataclass(frozen=True)
class OpenWorldGradingInput:
    completeness: float | None   # CheckM (bacterial) or CheckV (viral) % completeness
    contamination: float | None  # CheckM contamination % (None for viral / not assessed)
    supporting_signal: int = 0   # viral hallmark genes / pathogenicity markers


def grade_open_world(g: OpenWorldGradingInput) -> EvidenceGrade:
    """Grade a reference-free hit on the shared A/B/C/X scale. Capped at B: Grade A
    is reserved for NTC-controlled targeted detection. Completeness is the spine; a
    supporting signal (hallmarks/markers) rescues a hit with no completeness QC."""
    comp = g.completeness
    if comp is None:
        # no CheckM/CheckV estimate (DB absent / not assessed): a recovered genome
        # with a supporting signal is real-but-unverified (C); nothing else gradeable.
        return EvidenceGrade.C if g.supporting_signal >= 1 else EvidenceGrade.X
    if comp < _OW_COMPLETE_MIN:
        return EvidenceGrade.X
    if g.contamination is not None and g.contamination > _OW_MAX_CONTAM:
        return EvidenceGrade.X
    if comp >= _OW_COMPLETE_B:
        return EvidenceGrade.B
    return EvidenceGrade.C


def grade_mag(mag, *, n_markers: int = 0) -> EvidenceGrade:
    """Open-world grade for an assembly MAG (CheckM completeness/contamination +
    pathogenicity marker count as supporting signal)."""
    return grade_open_world(OpenWorldGradingInput(
        completeness=mag.completeness, contamination=mag.contamination,
        supporting_signal=n_markers,
    ))


def grade_viral(vc) -> EvidenceGrade:
    """Open-world grade for a viral contig (CheckV completeness + geNomad hallmark
    gene count as supporting signal)."""
    return grade_open_world(OpenWorldGradingInput(
        completeness=vc.completeness, contamination=None,
        supporting_signal=(vc.n_hallmarks or 0),
    ))


@dataclass
class ReportEntry:
    organism: str
    abundance: float
    ci_lower: float
    ci_upper: float
    read_count: int
    specimen_type: SpecimenType
    contaminant_risk: bool = False
    taxon_id: str = ""    # stable GCF/GCA accession; join key for NTC background
    tier: int = 3         # NTC background tier applied to this finding (1/2/3)
    crossmap_of: str = ""  # dominant relative this is likely cross-mapping from
    absolute_copies: float | None = None  # spike-in absolute load (None = no spike)
    breadth_ratio: float | None = None  # observed/expected coverage breadth (None = not computed)

    def as_input(self) -> GradingInput:
        """Normalize this entry into the struct the grading rule reads."""
        return GradingInput(
            read_count=self.read_count,
            abundance=self.abundance,
            ci_width=self.ci_upper - self.ci_lower,
            contaminant_risk=self.contaminant_risk,
            specimen_type=self.specimen_type,
            tier=self.tier,
            crossmap=bool(self.crossmap_of),
            breadth_ratio=self.breadth_ratio,
        )

    @property
    def invalid_stats(self) -> bool:
        """Surfaced in the report so a degenerate finding is visibly flagged
        rather than silently scored Grade X."""
        return _has_invalid_stats(self.as_input())

    @property
    def grade(self) -> EvidenceGrade:
        return grade(self.as_input())


def build_entries(
    cfg: PipelineConfig,
    organism_names: list[str],
    em_result: EMResult,
    ci_lower: np.ndarray,
    ci_upper: np.ndarray,
    taxon_ids: list[str],
    *,
    background: BackgroundModel | None = None,
    coverage: list | None = None,   # CoverageStats parallel to organism_names
) -> list[ReportEntry]:
    """Build the report entries once, for every renderer (JSON/TSV/PDF/HTML).

    Applies the NTC tier, the contaminant blocklist, and the background filter,
    then sorts by abundance:

        construct ──► set tier ──► (tier != 1) static blocklist   [CQ3]
                  ──► (background) drop taxa where p >= alpha
                  ──► sort by abundance desc
    """
    tier = background.tier if background is not None else 3
    read_counts = (em_result.abundances * em_result.n_reads).astype(int)
    entries = [
        ReportEntry(
            organism=name,
            abundance=float(em_result.abundances[i]),
            ci_lower=float(ci_lower[i]),
            ci_upper=float(ci_upper[i]),
            read_count=int(read_counts[i]),
            specimen_type=cfg.specimen_type,
            taxon_id=taxon_ids[i],
            tier=tier,
            breadth_ratio=(coverage[i].breadth_ratio if coverage else None),
        )
        for i, name in enumerate(organism_names)
    ]
    # Plan-3D: flag cross-mapping artifacts (e.g. Shigella from a dominant E. coli)
    # so grade() drops them as X.
    entries = deduplicate_closely_related(entries)
    # CQ3: a batch-matched NTC (Tier 1) supersedes the static contaminant
    # blocklist; the blocklist still applies at Tier 2/3.
    if tier != 1:
        entries = flag_contaminants(entries)
    if background is not None:
        kept = []
        for e in entries:
            if is_background(e.taxon_id, e.read_count, em_result.n_reads, background):
                # Dual-use pathogen + background overlap: flag, don't erase. A
                # contaminated NTC must not be allowed to subtract away a real
                # treatable pathogen (validated on the air NTC, 2026-06-23).
                if is_dual_use(e.organism):
                    e.contaminant_risk = True
                    kept.append(e)
                # else: pure contaminant -> subtract (drop) as before
            else:
                kept.append(e)
        entries = kept
    entries.sort(key=lambda e: e.abundance, reverse=True)
    return entries


def write_report(
    cfg: PipelineConfig,
    entries: list[ReportEntry],
    em_result: EMResult,
    amr_hits: list | None = None,
    virulence_hits: list | None = None,
    spike_info=None,
    mags: list | None = None,
    novelty=None,
    viral: list | None = None,
    pathogenicity: list | None = None,
    timings: dict | None = None,
) -> Path:
    out = cfg.output_dir / "report"
    out.mkdir(parents=True, exist_ok=True)

    entries = sorted(entries, key=lambda e: e.abundance, reverse=True)

    # A contig shared by near-identical siblings (E. coli / Shigella) attributes its
    # hit to each co-mapped organism, so the gene shows under every finding whose
    # genome carries that region — not just the single best.
    amr_by_org: dict[str, list[dict]] = {}
    for hit in (amr_hits or []):
        for org in (hit.organism_matches or [hit.organism_match]):
            amr_by_org.setdefault(org, []).append({
                "gene": hit.gene,
                "drug_class": hit.drug_class,
                "identity_pct": hit.identity_pct,
                "coverage_pct": hit.coverage_pct,
            })

    vir_by_org: dict[str, list[dict]] = {}
    for hit in (virulence_hits or []):
        for org in (hit.organism_matches or [hit.organism_match]):
            vir_by_org.setdefault(org, []).append({
                "gene": hit.gene,
                "factor": hit.factor,
                "identity_pct": hit.identity_pct,
                "coverage_pct": hit.coverage_pct,
            })

    tsv_path = out / "pathogeniq_report.tsv"
    with open(tsv_path, "w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=["organism", "taxon_id", "abundance_pct", "ci_lower_pct", "ci_upper_pct",
                        "read_count", "grade", "tier", "contaminant_risk", "invalid_stats",
                        "crossmap_of", "absolute_copies"],
            delimiter="\t",
        )
        writer.writeheader()
        for e in entries:
            writer.writerow({
                "organism": e.organism,
                "taxon_id": e.taxon_id,
                "abundance_pct": f"{e.abundance * 100:.2f}",
                "ci_lower_pct": f"{e.ci_lower * 100:.2f}",
                "ci_upper_pct": f"{e.ci_upper * 100:.2f}",
                "read_count": e.read_count,
                "grade": e.grade.value,
                "tier": e.tier,
                "contaminant_risk": e.contaminant_risk,
                "invalid_stats": e.invalid_stats,
                "crossmap_of": e.crossmap_of,
                "absolute_copies": "" if e.absolute_copies is None else f"{e.absolute_copies:.6g}",
            })

    vol = getattr(spike_info, "sample_volume", None) if spike_info is not None else None

    def _per_volume(copies: float | None) -> float | None:
        return round(copies / vol, 4) if (copies is not None and vol) else None

    json_path = out / "pathogeniq_report.json"
    payload = {
        "sample": str(cfg.input_fastq.name),
        "specimen_type": cfg.specimen_type.value,
        "read_type": cfg.read_type.value,
        "total_classified_reads": em_result.n_reads,
        "findings": [
            {
                "organism": e.organism,
                "taxon_id": e.taxon_id,
                "abundance_pct": round(e.abundance * 100, 2),
                "ci_lower_pct": round(e.ci_lower * 100, 2),
                "ci_upper_pct": round(e.ci_upper * 100, 2),
                "read_count": e.read_count,
                "grade": e.grade.value,
                "tier": e.tier,
                "contaminant_risk": e.contaminant_risk,
                "invalid_stats": e.invalid_stats,
                "crossmap_of": e.crossmap_of,
                "breadth_ratio": (None if e.breadth_ratio is None
                                  else round(e.breadth_ratio, 3)),
                "absolute_copies": (None if e.absolute_copies is None
                                    else round(e.absolute_copies, 4)),
                "copies_per_volume": _per_volume(e.absolute_copies),
                "amr_genes": amr_by_org.get(e.organism, []),
                "virulence_factors": vir_by_org.get(e.organism, []),
            }
            for e in entries
        ],
    }
    if spike_info is not None:
        payload["spike_in"] = {
            "taxon_id": spike_info.taxon_id,
            "copies_added": spike_info.copies_added,
            "spike_reads": spike_info.spike_reads,
            "sample_volume": spike_info.sample_volume,
            "found": spike_info.found,
        }
    # R5: pathogenicity marker count per MAG feeds its open-world grade.
    markers_by_mag = {a.name: a.n_virulence + a.n_amr for a in (pathogenicity or [])}
    if mags is not None:
        payload["mags"] = [
            {
                "bin_id": m.bin_id,
                "taxonomy": m.taxonomy,
                "completeness": m.completeness,
                "contamination": m.contamination,
                "n_contigs": m.n_contigs,
                "total_bp": m.total_bp,
                "grade": grade_mag(m, n_markers=markers_by_mag.get(m.bin_id, 0)).value,
                "fasta_path": str(m.fasta_path),
            }
            for m in mags
        ]
    if viral is not None:
        payload["viral"] = [
            {
                "contig_id": v.contig_id,
                "length": v.length,
                "taxonomy": v.taxonomy,
                "topology": v.topology,
                "virus_score": v.virus_score,
                "n_hallmarks": v.n_hallmarks,
                "completeness": v.completeness,
                "checkv_quality": v.checkv_quality,
                "grade": grade_viral(v).value,
            }
            for v in viral
        ]
    if pathogenicity is not None:
        payload["pathogenicity"] = [
            {
                "mag": a.name,
                "taxonomy": a.taxonomy,
                "verdict": a.verdict,
                "score": a.score,
                "n_virulence": a.n_virulence,
                "n_amr": a.n_amr,
                "phylo_match": a.phylo_match,
                "virulence_genes": a.virulence_genes,
                "amr_genes": a.amr_genes,
            }
            for a in pathogenicity
        ]
    if timings:
        payload["stage_seconds"] = timings
    if novelty is not None:
        payload["novelty"] = {
            "total_reads": novelty.total_reads,
            "classified_reads": novelty.classified_reads,
            "unclassified_reads": novelty.unclassified_reads,
            "unclassified_fraction": round(novelty.unclassified_fraction, 4),
            "n_species": novelty.n_species,
            "top_taxa": [{"species": s, "reads": n} for s, n in novelty.top_taxa],
            "flagged": novelty.flagged,
            "flag_threshold": novelty.flag_threshold,
        }
    with open(json_path, "w") as f:
        json.dump(payload, f, indent=2)

    return out
