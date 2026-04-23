import json
import csv
import numpy as np
from pathogeniq.config import PipelineConfig, ReadType, SpecimenType
from pathogeniq.em import EMResult
from pathogeniq.report import write_report, ReportEntry, EvidenceGrade


def _cfg(tmp_path):
    return PipelineConfig(
        input_fastq=tmp_path / "in.fastq.gz",
        read_type=ReadType.SHORT,
        specimen_type=SpecimenType.BLOOD,
        output_dir=tmp_path,
        db_tier1=tmp_path / "db",
        host_reference=tmp_path / "h.fa",
    )


def _em_result():
    return EMResult(
        abundances=np.array([0.70, 0.28, 0.02]),
        n_reads=1000,
        n_organisms=3,
        iterations=45,
    )


def test_write_report_creates_tsv(tmp_path):
    cfg = _cfg(tmp_path)
    em = _em_result()
    names = ["Escherichia coli", "Klebsiella pneumoniae", "Cutibacterium acnes"]
    lower = np.array([0.65, 0.23, 0.005])
    upper = np.array([0.75, 0.33, 0.04])
    write_report(cfg, names, em, lower, upper)
    tsv = tmp_path / "report" / "pathogeniq_report.tsv"
    assert tsv.exists()


def test_write_report_creates_json(tmp_path):
    cfg = _cfg(tmp_path)
    em = _em_result()
    names = ["Escherichia coli", "Klebsiella pneumoniae", "Cutibacterium acnes"]
    lower = np.array([0.65, 0.23, 0.005])
    upper = np.array([0.75, 0.33, 0.04])
    write_report(cfg, names, em, lower, upper)
    jsn = tmp_path / "report" / "pathogeniq_report.json"
    assert jsn.exists()
    with open(jsn) as f:
        data = json.load(f)
    assert "findings" in data
    assert len(data["findings"]) == 3


def test_evidence_grade_blood_high_confidence():
    entry = ReportEntry(
        organism="Escherichia coli",
        abundance=0.70,
        ci_lower=0.65,
        ci_upper=0.75,
        read_count=700,
        specimen_type=SpecimenType.BLOOD,
    )
    assert entry.grade == EvidenceGrade.A


def test_evidence_grade_blood_low_reads():
    entry = ReportEntry(
        organism="Cutibacterium acnes",
        abundance=0.001,
        ci_lower=0.0,
        ci_upper=0.01,
        read_count=1,
        specimen_type=SpecimenType.BLOOD,
    )
    assert entry.grade == EvidenceGrade.X


def test_tsv_columns(tmp_path):
    cfg = _cfg(tmp_path)
    em = _em_result()
    names = ["Organism A", "Organism B", "Organism C"]
    lower = np.array([0.60, 0.20, 0.00])
    upper = np.array([0.80, 0.36, 0.04])
    write_report(cfg, names, em, lower, upper)
    tsv = tmp_path / "report" / "pathogeniq_report.tsv"
    with open(tsv) as f:
        reader = csv.DictReader(f, delimiter="\t")
        rows = list(reader)
    assert "organism" in rows[0]
    assert "abundance_pct" in rows[0]
    assert "ci_lower_pct" in rows[0]
    assert "ci_upper_pct" in rows[0]
    assert "grade" in rows[0]
    assert "contaminant_risk" in rows[0]


def test_report_json_includes_contaminant_risk(tmp_path):
    cfg = _cfg(tmp_path)
    em = _em_result()
    names = ["Staphylococcus aureus", "Cutibacterium acnes", "Escherichia coli"]
    lower = np.array([0.60, 0.01, 0.20])
    upper = np.array([0.80, 0.10, 0.36])
    write_report(cfg, names, em, lower, upper)
    json_path = tmp_path / "report" / "pathogeniq_report.json"
    with open(json_path) as f:
        data = json.load(f)

    findings = {f["organism"]: f for f in data["findings"]}
    assert findings["Staphylococcus aureus"]["contaminant_risk"] is False
    assert findings["Cutibacterium acnes"]["contaminant_risk"] is True
    assert findings["Escherichia coli"]["contaminant_risk"] is False


def test_report_tsv_includes_contaminant_risk(tmp_path):
    cfg = _cfg(tmp_path)
    em = _em_result()
    names = ["Staphylococcus aureus", "Cutibacterium acnes"]
    lower = np.array([0.60, 0.01])
    upper = np.array([0.80, 0.10])
    write_report(cfg, names, em, lower, upper)
    tsv = tmp_path / "report" / "pathogeniq_report.tsv"
    with open(tsv) as f:
        reader = csv.DictReader(f, delimiter="\t")
        rows = {r["organism"]: r for r in reader}
    assert rows["Staphylococcus aureus"]["contaminant_risk"] == "False"
    assert rows["Cutibacterium acnes"]["contaminant_risk"] == "True"
