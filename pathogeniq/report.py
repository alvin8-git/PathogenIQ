import csv
import json
from dataclasses import dataclass
from enum import Enum
from pathlib import Path

import numpy as np

from .config import PipelineConfig, SpecimenType
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


@dataclass
class ReportEntry:
    organism: str
    abundance: float
    ci_lower: float
    ci_upper: float
    read_count: int
    specimen_type: SpecimenType
    contaminant_risk: bool = False

    @property
    def grade(self) -> EvidenceGrade:
        min_reads = _MIN_READS.get(self.specimen_type, 5)
        ci_width = self.ci_upper - self.ci_lower
        if self.read_count >= min_reads and ci_width <= _MAX_CI_WIDTH_A and not self.contaminant_risk:
            return EvidenceGrade.A
        if self.read_count >= min_reads and not self.contaminant_risk:
            return EvidenceGrade.B
        if self.read_count >= min_reads:
            return EvidenceGrade.C
        return EvidenceGrade.X


def write_report(
    cfg: PipelineConfig,
    organism_names: list[str],
    em_result: EMResult,
    ci_lower: np.ndarray,
    ci_upper: np.ndarray,
) -> Path:
    out = cfg.output_dir / "report"
    out.mkdir(parents=True, exist_ok=True)

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
    entries.sort(key=lambda e: e.abundance, reverse=True)

    tsv_path = out / "pathogeniq_report.tsv"
    with open(tsv_path, "w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=["organism", "abundance_pct", "ci_lower_pct", "ci_upper_pct", "read_count", "grade"],
            delimiter="\t",
        )
        writer.writeheader()
        for e in entries:
            writer.writerow({
                "organism": e.organism,
                "abundance_pct": f"{e.abundance * 100:.2f}",
                "ci_lower_pct": f"{e.ci_lower * 100:.2f}",
                "ci_upper_pct": f"{e.ci_upper * 100:.2f}",
                "read_count": e.read_count,
                "grade": e.grade.value,
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
                "abundance_pct": round(e.abundance * 100, 2),
                "ci_lower_pct": round(e.ci_lower * 100, 2),
                "ci_upper_pct": round(e.ci_upper * 100, 2),
                "read_count": e.read_count,
                "grade": e.grade.value,
            }
            for e in entries
        ],
    }
    with open(json_path, "w") as f:
        json.dump(payload, f, indent=2)

    return out
