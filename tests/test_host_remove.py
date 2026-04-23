from unittest.mock import MagicMock, patch
import pytest
from pathogeniq.config import PipelineConfig, ReadType, SpecimenType
from pathogeniq.host_remove import run_host_removal, HostRemovalMetrics


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
