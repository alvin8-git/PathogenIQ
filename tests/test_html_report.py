from pathlib import Path
import numpy as np
from pathogeniq.config import PipelineConfig, ReadType, SpecimenType
from pathogeniq.qc import QCMetrics
from pathogeniq.host_remove import HostRemovalMetrics
from pathogeniq.sketch import SketchHit
from pathogeniq.em import EMResult
from pathogeniq.amr import AMRHit
from pathogeniq.report import build_entries
from pathogeniq.html_report import write_html_report


def _cfg(tmp_path: Path) -> PipelineConfig:
    return PipelineConfig(
        input_fastq=tmp_path / "sample.fq.gz",
        read_type=ReadType.SHORT,
        specimen_type=SpecimenType.BLOOD,
        output_dir=tmp_path,
        db_tier1=tmp_path / "db.sbt.zip",
        host_reference=tmp_path / "human.fa",
    )


def _inputs(tmp_path):
    qc = QCMetrics(total_reads=1_000_000, passing_reads=950_000)
    hr = HostRemovalMetrics(total_reads=950_000, human_reads=5_000, nonhuman_reads=945_000)
    hits = [
        SketchHit("Staphylococcus aureus", 0.85, tmp_path / "staph.fna"),
        SketchHit("Escherichia coli", 0.42, tmp_path / "ecoli.fna"),
    ]
    organisms = ["Staphylococcus aureus", "Escherichia coli"]
    em = EMResult(abundances=np.array([0.80, 0.20]), n_reads=50000, n_organisms=2, iterations=12)
    ci_lower = np.array([0.75, 0.15])
    ci_upper = np.array([0.85, 0.25])
    entries = build_entries(_cfg(tmp_path), organisms, em, ci_lower, ci_upper, ["", ""])
    return qc, hr, hits, entries, em


def test_html_report_creates_file(tmp_path):
    cfg = _cfg(tmp_path)
    qc, hr, hits, entries, em = _inputs(tmp_path)
    path = write_html_report(cfg, qc, hr, hits, entries, em)
    assert path.exists()
    assert path.suffix == ".html"
    assert path.stat().st_size > 2000


def test_html_report_output_path(tmp_path):
    cfg = _cfg(tmp_path)
    qc, hr, hits, entries, em = _inputs(tmp_path)
    path = write_html_report(cfg, qc, hr, hits, entries, em)
    assert path == tmp_path / "report" / "pathogeniq_report.html"


def test_html_report_contains_key_sections(tmp_path):
    cfg = _cfg(tmp_path)
    qc, hr, hits, entries, em = _inputs(tmp_path)
    path = write_html_report(cfg, qc, hr, hits, entries, em)
    content = path.read_text()
    assert "Quality Control" in content
    assert "Host Removal" in content
    assert "Sketch Screening" in content
    assert "Microbial Findings" in content
    assert "Antimicrobial Resistance" in content


def test_html_report_organisms_present(tmp_path):
    cfg = _cfg(tmp_path)
    qc, hr, hits, entries, em = _inputs(tmp_path)
    path = write_html_report(cfg, qc, hr, hits, entries, em)
    content = path.read_text()
    assert "Staphylococcus aureus" in content
    assert "Escherichia coli" in content


def test_html_report_with_amr_hits(tmp_path):
    cfg = _cfg(tmp_path)
    qc, hr, hits, entries, em = _inputs(tmp_path)
    amr = [AMRHit("mecA", "beta-lactam", 99.5, 100.0, "Staphylococcus aureus", "card")]
    path = write_html_report(cfg, qc, hr, hits, entries, em, amr_hits=amr)
    content = path.read_text()
    assert "mecA" in content
    assert "beta-lactam" in content


def test_html_report_no_amr_message(tmp_path):
    cfg = _cfg(tmp_path)
    qc, hr, hits, entries, em = _inputs(tmp_path)
    path = write_html_report(cfg, qc, hr, hits, entries, em)
    content = path.read_text()
    assert "No AMR genes detected" in content


def test_html_report_amr_collapsible_sortable(tmp_path):
    cfg = _cfg(tmp_path)
    qc, hr, hits, entries, em = _inputs(tmp_path)
    amr = [AMRHit(f"gene{i}", "beta-lactam", 99.0, 95.0, "Escherichia coli", "card")
           for i in range(8)]
    content = write_html_report(cfg, qc, hr, hits, entries, em, amr_hits=amr).read_text()
    assert 'class="table-scroll"' in content   # wide table stays in the page
    assert 'class="sortable"' in content        # click-to-sort headers
    assert 'class="extra"' in content           # rows past the first 5 collapsed
    assert "Show 3 more" in content             # 8 hits -> 3 behind the toggle
    assert "<script>" in content                # sort/toggle JS present


def test_html_report_with_mags(tmp_path):
    from pathogeniq.assembly import MAG
    cfg = _cfg(tmp_path)
    qc, hr, hits, entries, em = _inputs(tmp_path)
    mags = [MAG("bin.1", tmp_path / "bin.1.fa", 92.0, 2.0,
                "d__Bacteria;g__Sphingomonas", 45, 3_800_000)]
    path = write_html_report(cfg, qc, hr, hits, entries, em, mags=mags)
    content = path.read_text()
    assert "Metagenome-Assembled Genomes" in content
    assert "Sphingomonas" in content


def test_html_report_absolute_copies_column(tmp_path):
    cfg = _cfg(tmp_path)
    qc, hr, hits, entries, em = _inputs(tmp_path)
    entries[0].absolute_copies = 1234.5
    path = write_html_report(cfg, qc, hr, hits, entries, em)
    content = path.read_text()
    assert "Abs. copies" in content
    # No column when no entry carries a spike-derived copy number.
    _, _, _, plain, _ = _inputs(tmp_path)
    plain_content = write_html_report(cfg, qc, hr, hits, plain, em).read_text()
    assert "Abs. copies" not in plain_content


def test_html_report_contaminant_flag(tmp_path):
    cfg = _cfg(tmp_path)
    qc = QCMetrics(total_reads=100_000, passing_reads=95_000)
    hr = HostRemovalMetrics(total_reads=95_000, human_reads=1_000, nonhuman_reads=94_000)
    hits = [SketchHit("Cutibacterium acnes", 0.30, tmp_path / "cutibac.fna")]
    em = EMResult(abundances=np.array([1.0]), n_reads=500, n_organisms=1, iterations=5)
    cil = np.array([0.90])
    ciu = np.array([1.00])
    entries = build_entries(cfg, ["Cutibacterium acnes"], em, cil, ciu, [""])
    path = write_html_report(cfg, qc, hr, hits, entries, em)
    content = path.read_text()
    assert "contaminant" in content
