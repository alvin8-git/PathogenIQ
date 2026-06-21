import numpy as np
from click.testing import CliRunner
from unittest.mock import patch
from pathogeniq.cli import cli
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
         patch("pathogeniq.cli.run_sketch_screen") as mock_sketch, \
         patch("pathogeniq.cli.run_targeted_alignment") as mock_align, \
         patch("pathogeniq.cli.em_abundance") as mock_em, \
         patch("pathogeniq.cli.bootstrap_ci") as mock_ci, \
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
         patch("pathogeniq.cli.run_sketch_screen", return_value=[SketchHit("Staphylococcus aureus", 0.05, tmp_path / "genome.fa")]), \
         patch("pathogeniq.cli.run_targeted_alignment", return_value=fake_align), \
         patch("pathogeniq.cli.em_abundance", return_value=fake_em), \
         patch("pathogeniq.cli.bootstrap_ci", return_value=(np.array([0.9]), np.array([1.0]))), \
         patch("pathogeniq.cli.run_amr_screen", return_value=[]), \
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
         patch("pathogeniq.cli.run_sketch_screen", return_value=[SketchHit("Staphylococcus aureus", 0.05, tmp_path / "genome.fa")]), \
         patch("pathogeniq.cli.run_targeted_alignment", return_value=fake_align), \
         patch("pathogeniq.cli.em_abundance", return_value=fake_em), \
         patch("pathogeniq.cli.bootstrap_ci", return_value=(np.array([0.9]), np.array([1.0]))), \
         patch("pathogeniq.cli.run_amr_screen", return_value=[]), \
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
    assert "No pathogens detected" in result.output
