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
    if g.contaminant_risk:
        return EvidenceGrade.C
    if g.ci_width <= _MAX_CI_WIDTH_A and g.tier == 1:
        return EvidenceGrade.A
    return EvidenceGrade.B


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
) -> Path:
    out = cfg.output_dir / "report"
    out.mkdir(parents=True, exist_ok=True)

    entries = sorted(entries, key=lambda e: e.abundance, reverse=True)

    amr_by_org: dict[str, list[dict]] = {}
    for hit in (amr_hits or []):
        amr_by_org.setdefault(hit.organism_match, []).append({
            "gene": hit.gene,
            "drug_class": hit.drug_class,
            "identity_pct": hit.identity_pct,
            "coverage_pct": hit.coverage_pct,
        })

    vir_by_org: dict[str, list[dict]] = {}
    for hit in (virulence_hits or []):
        vir_by_org.setdefault(hit.organism_match, []).append({
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
    if mags is not None:
        payload["mags"] = [
            {
                "bin_id": m.bin_id,
                "taxonomy": m.taxonomy,
                "completeness": m.completeness,
                "contamination": m.contamination,
                "n_contigs": m.n_contigs,
                "total_bp": m.total_bp,
                "fasta_path": str(m.fasta_path),
            }
            for m in mags
        ]
    with open(json_path, "w") as f:
        json.dump(payload, f, indent=2)

    return out
