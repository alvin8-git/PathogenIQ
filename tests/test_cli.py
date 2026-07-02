import numpy as np
from click.testing import CliRunner
from unittest.mock import patch
from pathogeniq.cli import cli, _resolve_background
from pathogeniq.em import EMResult
from pathogeniq.qc import QCMetrics
from pathogeniq.host_remove import HostRemovalMetrics
from pathogeniq.sketch import SketchHit
from pathogeniq.align import AlignmentResult


def test_cli_run_requires_input():
    runner = CliRunner()
    result = runner.invoke(cli, ["run"])
    assert result.exit_code != 0
    assert "Missing" in result.output or "Error" in result.output


def test_cli_run_help():
    runner = CliRunner()
    result = runner.invoke(cli, ["run", "--help"])
    assert result.exit_code == 0
    assert "--input" in result.output
    assert "--specimen" in result.output
    assert "--read-type" in result.output
    assert "--ntc" in result.output
    assert "--background" in result.output
    assert "--no-background" in result.output


def test_resolve_background_no_background_returns_none():
    assert _resolve_background(None, None, None, True) is None


def test_resolve_background_none_when_no_flags(monkeypatch):
    # no flags + no curated default -> Tier 3 (None). Monkeypatched so the test
    # does not depend on whether the shipped default table is populated.
    monkeypatch.setattr("pathogeniq.cli.load_default_background", lambda: None)
    assert _resolve_background(None, None, None, False) is None


def test_resolve_background_uses_populated_default(tmp_path, monkeypatch):
    from pathogeniq import background as bg
    out = tmp_path / "bg.tsv"
    bg.write_background_table(out, {"GCF_a": 5.0}, tier=2)
    monkeypatch.setattr("pathogeniq.cli.load_default_background",
                        lambda: bg.load_background_table(out))
    model = _resolve_background(None, None, None, False)
    assert model is not None and model.tier == 2


def test_no_background_overrides_populated_default(monkeypatch):
    # --no-background must force Tier 3 even when a default table exists
    monkeypatch.setattr("pathogeniq.cli.load_default_background",
                        lambda: (_ for _ in ()).throw(AssertionError("default must not load")))
    assert _resolve_background(None, None, None, True) is None


def test_resolve_background_table_path(tmp_path):
    table = tmp_path / "bg.tsv"
    table.write_text("# tier=2\ntaxon_id\trpm\nGCF_x\t1.0\n")
    model = _resolve_background(None, None, table, False)
    assert model is not None and model.tier == 2


def test_resolve_background_table_beats_ntc(tmp_path):
    # --background wins over --ntc: the NTC must NOT be classified
    table = tmp_path / "bg.tsv"
    table.write_text("# tier=2\ntaxon_id\trpm\nGCF_x\t1.0\n")
    ntc = tmp_path / "ntc.fastq.gz"
    ntc.touch()
    with patch("pathogeniq.cli._classify_taxon_counts",
               side_effect=AssertionError("should not classify when --background given")):
        model = _resolve_background(None, ntc, table, False)
    assert model.tier == 2


def test_resolve_background_ntc_builds_tier1(tmp_path):
    ntc = tmp_path / "ntc.fastq.gz"
    ntc.touch()
    with patch("pathogeniq.cli._classify_taxon_counts",
               return_value=({"GCF_bg": 100}, 1_000_000)):
        model = _resolve_background(None, ntc, None, False)
    assert model is not None and model.tier == 1
    assert "GCF_bg" in model.rates


def test_cli_run_invokes_pipeline(tmp_path):
    runner = CliRunner()
    fake_fastq = tmp_path / "sample.fastq.gz"
    fake_fastq.touch()
    fake_db = tmp_path / "db.sbt.zip"
    fake_db.touch()
    fake_ref = tmp_path / "human.fa"
    fake_ref.touch()

    em_result = EMResult(
        abundances=np.array([1.0]),
        n_reads=100,
        n_organisms=1,
        iterations=10,
    )

    with patch("pathogeniq.cli.run_qc") as mock_qc, \
         patch("pathogeniq.cli.run_host_removal") as mock_hr, \
         patch("pathogeniq.cli.run_phix_removal", side_effect=lambda cfg, fq: (fq, 0)), \
         patch("pathogeniq.cli.run_sketch_screen") as mock_sketch, \
         patch("pathogeniq.cli.run_targeted_alignment") as mock_align, \
         patch("pathogeniq.cli.em_abundance") as mock_em, \
         patch("pathogeniq.cli.bootstrap_ci") as mock_ci, \
         patch("pathogeniq.cli.run_megahit", return_value=None), \
         patch("pathogeniq.cli.run_amr_screen", return_value=[]), \
         patch("pathogeniq.cli.run_virulence_screen", return_value=[]), \
         patch("pathogeniq.cli.write_report") as mock_report, \
         patch("pathogeniq.cli.write_html_report") as mock_html:

        mock_qc.return_value = (tmp_path / "filtered.fastq.gz", QCMetrics(100, 95))
        mock_hr.return_value = (tmp_path / "nonhuman.fastq.gz", HostRemovalMetrics(95, 80, 15))
        mock_sketch.return_value = [SketchHit("org1", 0.5, tmp_path / "org1.fa")]
        mock_align.return_value = AlignmentResult(
            alignment_matrix=np.array([[1.0]]),
            organism_names=["org1"],
            taxon_ids=["GCF_000005845.2"],
            read_ids=["r1"],
        )
        mock_em.return_value = em_result
        mock_ci.return_value = (np.array([0.9]), np.array([1.0]))
        mock_report.return_value = tmp_path / "report"
        mock_html.return_value = tmp_path / "report" / "pathogeniq_report.html"

        result = runner.invoke(cli, [
            "run",
            "--input", str(fake_fastq),
            "--output", str(tmp_path / "out"),
            "--db", str(fake_db),
            "--host-ref", str(fake_ref),
            "--specimen", "blood",
            "--read-type", "short",
        ])

    assert result.exit_code == 0, result.output


def test_cli_run_skip_host_removal(tmp_path):
    """--skip-host-removal (F4): host removal is bypassed for host-free specimens;
    run_host_removal must NOT be called and the pipeline still completes."""
    runner = CliRunner()
    for name in ("s.fq.gz", "db.sbt.zip", "ref.fa"):
        (tmp_path / name).touch()
    fake_em = EMResult(abundances=np.array([1.0]), n_reads=10, n_organisms=1, iterations=3)
    fake_align = AlignmentResult(alignment_matrix=np.ones((10, 1)), organism_names=["org1"],
                                 read_ids=["r1"], taxon_ids=["GCF_000005845.2"])
    with patch("pathogeniq.cli.run_qc", return_value=(tmp_path / "filt.fq.gz", QCMetrics(100, 90))), \
         patch("pathogeniq.cli.run_host_removal",
               side_effect=AssertionError("host removal must be skipped")), \
         patch("pathogeniq.cli.run_phix_removal", side_effect=lambda cfg, fq: (fq, 0)), \
         patch("pathogeniq.cli.run_sketch_screen", return_value=[SketchHit("org1", 0.05, tmp_path / "g.fa")]), \
         patch("pathogeniq.cli.run_targeted_alignment", return_value=fake_align), \
         patch("pathogeniq.cli.em_abundance", return_value=fake_em), \
         patch("pathogeniq.cli.bootstrap_ci", return_value=(np.array([0.9]), np.array([1.0]))), \
         patch("pathogeniq.cli.run_megahit", return_value=None), \
         patch("pathogeniq.cli.run_amr_screen", return_value=[]), \
         patch("pathogeniq.cli.run_virulence_screen", return_value=[]), \
         patch("pathogeniq.cli.write_report", return_value=tmp_path / "report"), \
         patch("pathogeniq.cli.write_pdf_report", return_value=tmp_path / "r.pdf"), \
         patch("pathogeniq.cli.write_html_report", return_value=tmp_path / "r.html"):
        result = runner.invoke(cli, [
            "run", "--input", str(tmp_path / "s.fq.gz"), "--output", str(tmp_path / "out"),
            "--db", str(tmp_path / "db.sbt.zip"), "--host-ref", str(tmp_path / "ref.fa"),
            "--specimen", "air", "--read-type", "short", "--skip-host-removal",
        ])
    assert result.exit_code == 0, result.output
    assert "Skipped (--skip-host-removal)" in result.output


def test_cli_run_produces_pdf(tmp_path):
    """End-to-end CLI smoke test: PDF is produced when pipeline succeeds."""
    from unittest.mock import patch, MagicMock

    fq = tmp_path / "s.fq.gz"
    fq.touch()
    db = tmp_path / "db.sbt.zip"
    db.touch()
    ref = tmp_path / "ref.fa"
    ref.touch()

    fake_em = EMResult(abundances=np.array([1.0]), n_reads=10, n_organisms=1, iterations=3)
    fake_align = AlignmentResult(
        alignment_matrix=np.ones((10, 1)),
        organism_names=["Staphylococcus aureus"],
        read_ids=["r1"],
        taxon_ids=["GCF_000005845.2"],
    )

    runner = CliRunner()
    with patch("pathogeniq.cli.run_qc", return_value=(tmp_path / "filt.fq.gz", QCMetrics(100, 90))), \
         patch("pathogeniq.cli.run_host_removal", return_value=(tmp_path / "nonhuman.fq.gz", HostRemovalMetrics(90, 80, 10))), \
         patch("pathogeniq.cli.run_phix_removal", side_effect=lambda cfg, fq: (fq, 0)), \
         patch("pathogeniq.cli.run_sketch_screen", return_value=[SketchHit("Staphylococcus aureus", 0.05, tmp_path / "genome.fa")]), \
         patch("pathogeniq.cli.run_targeted_alignment", return_value=fake_align), \
         patch("pathogeniq.cli.em_abundance", return_value=fake_em), \
         patch("pathogeniq.cli.bootstrap_ci", return_value=(np.array([0.9]), np.array([1.0]))), \
         patch("pathogeniq.cli.run_megahit", return_value=None), \
         patch("pathogeniq.cli.run_amr_screen", return_value=[]), \
         patch("pathogeniq.cli.run_virulence_screen", return_value=[]), \
         patch("pathogeniq.cli.write_report", return_value=tmp_path / "report"), \
         patch("pathogeniq.cli.write_pdf_report") as mock_pdf, \
         patch("pathogeniq.cli.write_html_report") as mock_html:
        mock_pdf.return_value = tmp_path / "report" / "pathogeniq_report.pdf"
        mock_html.return_value = tmp_path / "report" / "pathogeniq_report.html"
        result = runner.invoke(cli, [
            "run",
            "--input", str(fq),
            "--output", str(tmp_path / "out"),
            "--db", str(db),
            "--host-ref", str(ref),
            "--specimen", "blood",
            "--read-type", "short",
        ])

    assert result.exit_code == 0, result.output
    mock_pdf.assert_called_once()
    mock_html.assert_called_once()


def test_cli_run_no_pdf_flag(tmp_path):
    """--no-pdf suppresses PDF generation."""
    from unittest.mock import patch

    fq = tmp_path / "s.fq.gz"
    fq.touch()
    db = tmp_path / "db.sbt.zip"
    db.touch()
    ref = tmp_path / "ref.fa"
    ref.touch()

    fake_em = EMResult(abundances=np.array([1.0]), n_reads=10, n_organisms=1, iterations=3)
    fake_align = AlignmentResult(
        alignment_matrix=np.ones((10, 1)),
        organism_names=["Staphylococcus aureus"],
        read_ids=["r1"],
        taxon_ids=["GCF_000005845.2"],
    )

    runner = CliRunner()
    with patch("pathogeniq.cli.run_qc", return_value=(tmp_path / "filt.fq.gz", QCMetrics(100, 90))), \
         patch("pathogeniq.cli.run_host_removal", return_value=(tmp_path / "nonhuman.fq.gz", HostRemovalMetrics(90, 80, 10))), \
         patch("pathogeniq.cli.run_phix_removal", side_effect=lambda cfg, fq: (fq, 0)), \
         patch("pathogeniq.cli.run_sketch_screen", return_value=[SketchHit("Staphylococcus aureus", 0.05, tmp_path / "genome.fa")]), \
         patch("pathogeniq.cli.run_targeted_alignment", return_value=fake_align), \
         patch("pathogeniq.cli.em_abundance", return_value=fake_em), \
         patch("pathogeniq.cli.bootstrap_ci", return_value=(np.array([0.9]), np.array([1.0]))), \
         patch("pathogeniq.cli.run_megahit", return_value=None), \
         patch("pathogeniq.cli.run_amr_screen", return_value=[]), \
         patch("pathogeniq.cli.run_virulence_screen", return_value=[]), \
         patch("pathogeniq.cli.write_report", return_value=tmp_path / "report"), \
         patch("pathogeniq.cli.write_pdf_report") as mock_pdf, \
         patch("pathogeniq.cli.write_html_report") as mock_html:
        mock_html.return_value = tmp_path / "report" / "pathogeniq_report.html"
        result = runner.invoke(cli, [
            "run",
            "--input", str(fq),
            "--output", str(tmp_path / "out"),
            "--db", str(db),
            "--host-ref", str(ref),
            "--specimen", "blood",
            "--no-pdf",
        ])

    assert result.exit_code == 0, result.output
    mock_pdf.assert_not_called()
    mock_html.assert_called_once()


def test_cli_run_no_candidates_exits_cleanly(tmp_path):
    runner = CliRunner()
    fake_fastq = tmp_path / "sample.fastq.gz"
    fake_fastq.touch()
    fake_db = tmp_path / "db.sbt.zip"
    fake_db.touch()
    fake_ref = tmp_path / "human.fa"
    fake_ref.touch()

    with patch("pathogeniq.cli.run_qc") as mock_qc, \
         patch("pathogeniq.cli.run_host_removal") as mock_hr, \
         patch("pathogeniq.cli.run_phix_removal", side_effect=lambda cfg, fq: (fq, 0)), \
         patch("pathogeniq.cli.run_sketch_screen") as mock_sketch:

        mock_qc.return_value = (tmp_path / "filtered.fastq.gz", QCMetrics(100, 95))
        mock_hr.return_value = (tmp_path / "nonhuman.fastq.gz", HostRemovalMetrics(95, 80, 15))
        mock_sketch.return_value = []  # no candidates

        result = runner.invoke(cli, [
            "run",
            "--input", str(fake_fastq),
            "--output", str(tmp_path / "out"),
            "--db", str(fake_db),
            "--host-ref", str(fake_ref),
            "--specimen", "csf",
            "--read-type", "short",
        ])

    assert result.exit_code == 0
    assert "No targeted pathogens detected" in result.output
