import json
import csv
from dataclasses import replace
import numpy as np
from pathogeniq.config import PipelineConfig, ReadType, SpecimenType
from pathogeniq.em import EMResult
from pathogeniq.background import build_background
from pathogeniq.report import (
    write_report,
    build_entries,
    ReportEntry,
    EvidenceGrade,
    GradingInput,
    grade,
)


def _write(cfg, names, em, lower, upper):
    """Build entries (Tier 3, no background) and write the report — mirrors the
    pipeline path now that write_report renders pre-built entries."""
    entries = build_entries(cfg, names, em, lower, upper, [""] * len(names))
    return write_report(cfg, entries, em)


def _gi(**kw):
    base = dict(
        read_count=700,
        abundance=0.70,
        ci_width=0.10,
        contaminant_risk=False,
        specimen_type=SpecimenType.BLOOD,
    )
    base.update(kw)
    return GradingInput(**base)


def test_grade_fn_grade_a():
    # Grade A requires a batch-matched NTC (tier 1)
    assert grade(_gi(tier=1)) == EvidenceGrade.A


def test_grade_fn_tier_caps_a_to_b():
    # without a batch-matched NTC, a would-be-A finding caps at B
    assert grade(_gi(tier=2)) == EvidenceGrade.B
    assert grade(_gi(tier=3)) == EvidenceGrade.B


def test_grade_fn_b_when_ci_wide():
    # wide CI, not a contaminant -> drops from A to B
    assert grade(_gi(ci_width=0.50)) == EvidenceGrade.B


def test_grade_fn_c_when_contaminant():
    # contaminant flag fails both A and B (which require not-contaminant) -> C
    assert grade(_gi(contaminant_risk=True)) == EvidenceGrade.C


def test_grade_fn_x_low_reads():
    assert grade(_gi(read_count=1)) == EvidenceGrade.X  # blood min_reads=3


def test_grade_fn_x_on_invalid():
    assert grade(_gi(abundance=float("nan"))) == EvidenceGrade.X
    assert grade(_gi(ci_width=float("nan"))) == EvidenceGrade.X
    assert grade(_gi(read_count=-5)) == EvidenceGrade.X


def test_report_entry_grade_delegates_to_fn():
    entry = ReportEntry(
        organism="Escherichia coli",
        abundance=0.70,
        ci_lower=0.65,
        ci_upper=0.75,
        read_count=700,
        specimen_type=SpecimenType.BLOOD,
    )
    assert entry.grade == grade(entry.as_input())


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
    _write(cfg, names, em, lower, upper)
    tsv = tmp_path / "report" / "pathogeniq_report.tsv"
    assert tsv.exists()


def test_write_report_creates_json(tmp_path):
    cfg = _cfg(tmp_path)
    em = _em_result()
    names = ["Escherichia coli", "Klebsiella pneumoniae", "Cutibacterium acnes"]
    lower = np.array([0.65, 0.23, 0.005])
    upper = np.array([0.75, 0.33, 0.04])
    _write(cfg, names, em, lower, upper)
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
        tier=1,   # batch-matched NTC required for Grade A
    )
    assert entry.grade == EvidenceGrade.A
    # same finding without a batch-matched NTC caps at B
    assert replace(entry, tier=3).grade == EvidenceGrade.B


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
    _write(cfg, names, em, lower, upper)
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
    _write(cfg, names, em, lower, upper)
    json_path = tmp_path / "report" / "pathogeniq_report.json"
    with open(json_path) as f:
        data = json.load(f)

    findings = {f["organism"]: f for f in data["findings"]}
    assert findings["Staphylococcus aureus"]["contaminant_risk"] is False
    assert findings["Cutibacterium acnes"]["contaminant_risk"] is True
    assert findings["Escherichia coli"]["contaminant_risk"] is False


def test_grade_x_on_nan_ci():
    # CRITICAL regression: NaN CI bounds (degenerate EM output) must not grade A/B/C
    entry = ReportEntry(
        organism="Escherichia coli",
        abundance=0.70,
        ci_lower=float("nan"),
        ci_upper=float("nan"),
        read_count=700,
        specimen_type=SpecimenType.BLOOD,
    )
    assert entry.invalid_stats is True
    assert entry.grade == EvidenceGrade.X


def test_grade_x_on_int64min_read_count():
    # CRITICAL regression: NaN abundance -> astype(int) underflows to INT64_MIN.
    # A negative read count is impossible and must be flagged, not graded.
    with np.errstate(invalid="ignore"):   # the cast IS the bug we're reproducing
        int64min = int(np.array([np.nan]).astype(int)[0])
    assert int64min < 0
    entry = ReportEntry(
        organism="Klebsiella pneumoniae",
        abundance=float("nan"),
        ci_lower=0.0,
        ci_upper=0.0,
        read_count=int64min,
        specimen_type=SpecimenType.BLOOD,
    )
    assert entry.invalid_stats is True
    assert entry.grade == EvidenceGrade.X


def test_valid_entry_not_flagged():
    # Regression guard: the NaN check must not trip on healthy values
    entry = ReportEntry(
        organism="Escherichia coli",
        abundance=0.70,
        ci_lower=0.65,
        ci_upper=0.75,
        read_count=700,
        specimen_type=SpecimenType.BLOOD,
        tier=1,
    )
    assert entry.invalid_stats is False
    assert entry.grade == EvidenceGrade.A


def test_report_flags_invalid_stats(tmp_path):
    cfg = _cfg(tmp_path)
    em = _em_result()
    with np.errstate(invalid="ignore"):   # the cast IS the bug we're reproducing
        bad_count = int(np.array([np.nan]).astype(int)[0])
    entries = [
        ReportEntry(
            organism="Degenerate org",
            abundance=float("nan"),
            ci_lower=float("nan"),
            ci_upper=float("nan"),
            read_count=bad_count,
            specimen_type=SpecimenType.BLOOD,
        ),
    ]
    write_report(cfg, entries, em)
    with open(tmp_path / "report" / "pathogeniq_report.json") as f:
        finding = json.load(f)["findings"][0]
    assert finding["invalid_stats"] is True
    assert finding["grade"] == "X"


def test_build_entries_sets_tier_and_filters_background(tmp_path):
    # the consolidated builder is the single path all renderers use: it must set
    # the tier AND drop background taxa (previously HTML skipped both)
    cfg = _cfg(tmp_path)
    em = EMResult(abundances=np.array([0.9, 0.1]), n_reads=1_000_000, n_organisms=2, iterations=5)
    names = ["Real pathogen", "Reagent bug"]
    lower = np.array([0.88, 0.05])
    upper = np.array([0.92, 0.15])
    taxon_ids = ["GCF_real", "GCF_bg"]
    # GCF_bg sits at the sample's level in the NTC (filtered); GCF_real is absent (kept)
    bg = build_background([({"GCF_bg": 100_000}, 1_000_000)], tier=2)
    entries = build_entries(cfg, names, em, lower, upper, taxon_ids, background=bg)
    orgs = [e.organism for e in entries]
    assert "Real pathogen" in orgs
    assert "Reagent bug" not in orgs
    assert all(e.tier == 2 for e in entries)


def test_report_includes_taxon_id(tmp_path):
    cfg = _cfg(tmp_path)
    em = _em_result()
    entries = [
        ReportEntry(
            organism="Escherichia coli",
            abundance=0.70,
            ci_lower=0.65,
            ci_upper=0.75,
            read_count=700,
            specimen_type=SpecimenType.BLOOD,
            taxon_id="GCF_000005845.2",
        ),
    ]
    write_report(cfg, entries, em)
    with open(tmp_path / "report" / "pathogeniq_report.json") as f:
        data = json.load(f)
    assert data["findings"][0]["taxon_id"] == "GCF_000005845.2"
    with open(tmp_path / "report" / "pathogeniq_report.tsv") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    assert rows[0]["taxon_id"] == "GCF_000005845.2"


def test_report_tsv_includes_contaminant_risk(tmp_path):
    cfg = _cfg(tmp_path)
    em = _em_result()
    names = ["Staphylococcus aureus", "Cutibacterium acnes"]
    lower = np.array([0.60, 0.01])
    upper = np.array([0.80, 0.10])
    _write(cfg, names, em, lower, upper)
    tsv = tmp_path / "report" / "pathogeniq_report.tsv"
    with open(tsv) as f:
        reader = csv.DictReader(f, delimiter="\t")
        rows = {r["organism"]: r for r in reader}
    assert rows["Staphylococcus aureus"]["contaminant_risk"] == "False"
    assert rows["Cutibacterium acnes"]["contaminant_risk"] == "True"
