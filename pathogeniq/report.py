import csv
import json
import math
from dataclasses import dataclass
from enum import Enum
from pathlib import Path

import numpy as np

from .config import PipelineConfig, SpecimenType
from .contaminants import flag_contaminants
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
    ReportEntry.grade and the benchmark adapter."""
    if _has_invalid_stats(g):
        return EvidenceGrade.X
    min_reads = _MIN_READS.get(g.specimen_type, 5)
    if g.read_count >= min_reads and g.ci_width <= _MAX_CI_WIDTH_A and not g.contaminant_risk:
        return EvidenceGrade.A
    if g.read_count >= min_reads and not g.contaminant_risk:
        return EvidenceGrade.B
    if g.read_count >= min_reads:
        return EvidenceGrade.C
    return EvidenceGrade.X


@dataclass
class ReportEntry:
    organism: str
    abundance: float
    ci_lower: float
    ci_upper: float
    read_count: int
    specimen_type: SpecimenType
    contaminant_risk: bool = False
    taxon_id: str = ""   # stable GCF/GCA accession; join key for NTC background

    def as_input(self) -> GradingInput:
        """Normalize this entry into the struct the grading rule reads."""
        return GradingInput(
            read_count=self.read_count,
            abundance=self.abundance,
            ci_width=self.ci_upper - self.ci_lower,
            contaminant_risk=self.contaminant_risk,
            specimen_type=self.specimen_type,
        )

    @property
    def invalid_stats(self) -> bool:
        """Surfaced in the report so a degenerate finding is visibly flagged
        rather than silently scored Grade X."""
        return _has_invalid_stats(self.as_input())

    @property
    def grade(self) -> EvidenceGrade:
        return grade(self.as_input())


def write_report(
    cfg: PipelineConfig,
    organism_names: list[str],
    em_result: EMResult,
    ci_lower: np.ndarray,
    ci_upper: np.ndarray,
    amr_hits: list | None = None,
    entries_override: list[ReportEntry] | None = None,
) -> Path:
    out = cfg.output_dir / "report"
    out.mkdir(parents=True, exist_ok=True)

    if entries_override is not None:
        entries = entries_override
    else:
        read_counts = (em_result.abundances * em_result.n_reads).astype(int)
        entries = [
            ReportEntry(
                organism=name,
                abundance=float(em_result.abundances[i]),
                ci_lower=float(ci_lower[i]),
                ci_upper=float(ci_upper[i]),
                read_count=int(read_counts[i]),
                specimen_type=cfg.specimen_type,
            )
            for i, name in enumerate(organism_names)
        ]
        entries = flag_contaminants(entries)
    entries.sort(key=lambda e: e.abundance, reverse=True)

    amr_by_org: dict[str, list[dict]] = {}
    for hit in (amr_hits or []):
        amr_by_org.setdefault(hit.organism_match, []).append({
            "gene": hit.gene,
            "drug_class": hit.drug_class,
            "identity_pct": hit.identity_pct,
            "coverage_pct": hit.coverage_pct,
        })

    tsv_path = out / "pathogeniq_report.tsv"
    with open(tsv_path, "w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=["organism", "taxon_id", "abundance_pct", "ci_lower_pct", "ci_upper_pct",
                        "read_count", "grade", "contaminant_risk", "invalid_stats"],
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
                "contaminant_risk": e.contaminant_risk,
                "invalid_stats": e.invalid_stats,
            })

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
                "contaminant_risk": e.contaminant_risk,
                "invalid_stats": e.invalid_stats,
                "amr_genes": amr_by_org.get(e.organism, []),
            }
            for e in entries
        ],
    }
    with open(json_path, "w") as f:
        json.dump(payload, f, indent=2)

    return out
