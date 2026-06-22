from unittest.mock import MagicMock, patch
import pytest
from pathogeniq.config import PipelineConfig, ReadType, SpecimenType
from pathogeniq.host_remove import (
    run_host_removal,
    run_phix_removal,
    phix_reference_path,
    HostRemovalMetrics,
)


def _cfg(tmp_path, read_type=ReadType.SHORT):
    return PipelineConfig(
        input_fastq=tmp_path / "filtered.fastq.gz",
        read_type=read_type,
        specimen_type=SpecimenType.BLOOD,
        output_dir=tmp_path,
        db_tier1=tmp_path / "db",
        host_reference=tmp_path / "human.fa",
        threads=4,
    )


def test_host_removal_sr_uses_bwa(tmp_path):
    cfg = _cfg(tmp_path, ReadType.SHORT)
    filtered = tmp_path / "filtered.fastq.gz"
    with patch("subprocess.run") as mock_run:
        mock_run.return_value = MagicMock(returncode=0)
        nonhuman, metrics = run_host_removal(cfg, filtered)
    calls = [" ".join(c[0][0]) for c in mock_run.call_args_list]
    assert any("bwa-mem2" in c or "minimap2" in c for c in calls)


def test_host_removal_lr_uses_minimap2(tmp_path):
    cfg = _cfg(tmp_path, ReadType.LONG)
    filtered = tmp_path / "filtered.fastq.gz"
    with patch("subprocess.run") as mock_run:
        mock_run.return_value = MagicMock(returncode=0)
        nonhuman, metrics = run_host_removal(cfg, filtered)
    calls = [" ".join(c[0][0]) for c in mock_run.call_args_list]
    assert any("minimap2" in c for c in calls)


def test_host_removal_creates_output_dir(tmp_path):
    cfg = _cfg(tmp_path)
    filtered = tmp_path / "filtered.fastq.gz"
    with patch("subprocess.run"):
        nonhuman, _ = run_host_removal(cfg, filtered)
    assert (tmp_path / "host_removal").is_dir()


def test_host_removal_metrics_fields():
    m = HostRemovalMetrics(total_reads=1000, human_reads=800, nonhuman_reads=200)
    assert m.microbial_fraction == pytest.approx(0.2)
    assert m.human_fraction == pytest.approx(0.8)


def test_phix_reference_is_bundled():
    # the PhiX-174 reference (NC_001422.1) must ship with the package
    assert phix_reference_path().exists()


def test_phix_removal_short_uses_minimap2_sr(tmp_path):
    cfg = _cfg(tmp_path, ReadType.SHORT)
    with patch("subprocess.run") as mock_run:
        mock_run.return_value = MagicMock(returncode=0)
        out, n_phix = run_phix_removal(cfg, tmp_path / "nonhuman.fastq.gz")
    calls = [" ".join(c[0][0]) for c in mock_run.call_args_list]
    assert any("minimap2" in c and "-ax sr" in c for c in calls)   # short-read preset
    assert out.name == "nophix.fastq.gz"
    assert n_phix == 0   # mocked flagstat -> nothing parsed


def test_phix_removal_long_uses_map_ont(tmp_path):
    cfg = _cfg(tmp_path, ReadType.LONG)
    with patch("subprocess.run") as mock_run:
        mock_run.return_value = MagicMock(returncode=0)
        run_phix_removal(cfg, tmp_path / "nonhuman.fastq.gz")
    calls = [" ".join(c[0][0]) for c in mock_run.call_args_list]
    assert any("map-ont" in c for c in calls)


def test_phix_removal_graceful_when_reference_missing(tmp_path):
    # non-blocking: a missing reference returns the input unchanged, runs nothing
    cfg = _cfg(tmp_path)
    inp = tmp_path / "nonhuman.fastq.gz"
    with patch("pathogeniq.host_remove.phix_reference_path", return_value=tmp_path / "nope.fasta"), \
         patch("subprocess.run") as mock_run:
        out, n_phix = run_phix_removal(cfg, inp)
    assert out == inp and n_phix == 0
    mock_run.assert_not_called()
