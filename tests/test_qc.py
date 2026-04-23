import json
from unittest.mock import MagicMock, patch, mock_open
from pathogeniq.config import PipelineConfig, ReadType, SpecimenType
from pathogeniq.qc import run_qc


def _cfg(tmp_path, read_type=ReadType.SHORT):
    return PipelineConfig(
        input_fastq=tmp_path / "in.fastq.gz",
        read_type=read_type,
        specimen_type=SpecimenType.BLOOD,
        output_dir=tmp_path,
        db_tier1=tmp_path / "db",
        host_reference=tmp_path / "h.fa",
        threads=4,
    )


FASTP_JSON = {
    "summary": {
        "before_filtering": {"total_reads": 1000},
        "after_filtering": {"total_reads": 950},
    }
}


def test_run_qc_sr_calls_fastp(tmp_path):
    cfg = _cfg(tmp_path, ReadType.SHORT)
    with patch("subprocess.run") as mock_run, \
         patch("builtins.open", mock_open(read_data=json.dumps(FASTP_JSON))):
        mock_run.return_value = MagicMock(returncode=0)
        filtered, metrics = run_qc(cfg)
    cmd = mock_run.call_args[0][0]
    assert cmd[0] == "fastp"
    assert "--qualified_quality_phred" in cmd
    assert "--low_complexity_filter" in cmd


def test_run_qc_sr_returns_metrics(tmp_path):
    cfg = _cfg(tmp_path, ReadType.SHORT)
    with patch("subprocess.run"), \
         patch("builtins.open", mock_open(read_data=json.dumps(FASTP_JSON))):
        filtered, metrics = run_qc(cfg)
    assert metrics.total_reads == 1000
    assert metrics.passing_reads == 950


def test_run_qc_lr_calls_chopper(tmp_path):
    cfg = _cfg(tmp_path, ReadType.LONG)
    with patch("subprocess.run") as mock_run:
        mock_run.return_value = MagicMock(returncode=0)
        filtered, metrics = run_qc(cfg)
    cmd = mock_run.call_args[0][0]
    assert "chopper" in " ".join(cmd).lower()


def test_run_qc_creates_output_dir(tmp_path):
    cfg = _cfg(tmp_path, ReadType.SHORT)
    with patch("subprocess.run"), \
         patch("builtins.open", mock_open(read_data=json.dumps(FASTP_JSON))):
        filtered, _ = run_qc(cfg)
    assert (tmp_path / "qc").is_dir()
