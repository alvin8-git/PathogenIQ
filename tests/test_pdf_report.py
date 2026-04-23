from pathlib import Path

import pytest

from pathogeniq.amr import AMRHit
from pathogeniq.config import PipelineConfig, ReadType, SpecimenType
from pathogeniq.pdf_report import write_pdf_report
from pathogeniq.report import ReportEntry


def _cfg(tmp_path: Path) -> PipelineConfig:
    return PipelineConfig(
        input_fastq=tmp_path / "sample.fq.gz",
        read_type=ReadType.SHORT,
        specimen_type=SpecimenType.BLOOD,
        output_dir=tmp_path,
        db_tier1=tmp_path / "db.sbt.zip",
        host_reference=tmp_path / "ref.fa",
    )


def _entries() -> list[ReportEntry]:
    return [
        ReportEntry("Staphylococcus aureus", 0.80, 0.72, 0.88, 80, SpecimenType.BLOOD),
        ReportEntry("Cutibacterium acnes", 0.15, 0.08, 0.22, 15, SpecimenType.BLOOD, contaminant_risk=True),
        ReportEntry("Escherichia coli", 0.05, 0.01, 0.09, 5, SpecimenType.BLOOD),
    ]


def test_pdf_report_creates_file(tmp_path):
    cfg = _cfg(tmp_path)
    pdf_path = write_pdf_report(cfg, _entries(), amr_hits=[])
    assert pdf_path.exists()
    assert pdf_path.suffix == ".pdf"
    assert pdf_path.stat().st_size > 1000  # non-trivial PDF


def test_pdf_report_with_amr_hits(tmp_path):
    cfg = _cfg(tmp_path)
    hits = [
        AMRHit("mecA", "beta-lactam", 99.5, 100.0, "Staphylococcus aureus", "card"),
        AMRHit("vanA", "glycopeptide", 98.1, 95.0, "unknown", "card"),
    ]
    pdf_path = write_pdf_report(cfg, _entries(), amr_hits=hits)
    assert pdf_path.exists()
    assert pdf_path.stat().st_size > 1000


def test_pdf_report_output_path(tmp_path):
    cfg = _cfg(tmp_path)
    pdf_path = write_pdf_report(cfg, _entries(), amr_hits=[])
    assert pdf_path == tmp_path / "report" / "pathogeniq_report.pdf"


def test_pdf_report_empty_entries(tmp_path):
    cfg = _cfg(tmp_path)
    pdf_path = write_pdf_report(cfg, [], amr_hits=[])
    assert pdf_path.exists()
    assert pdf_path.stat().st_size > 1000
