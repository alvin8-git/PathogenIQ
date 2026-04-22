import numpy as np
from pathlib import Path
from unittest.mock import MagicMock, patch
import pytest
from pathogeniq.config import PipelineConfig, ReadType, SpecimenType
from pathogeniq.sketch import SketchHit
from pathogeniq.align import run_targeted_alignment, AlignmentResult


def _cfg(tmp_path):
    return PipelineConfig(
        input_fastq=tmp_path / "nonhuman.fastq.gz",
        read_type=ReadType.SHORT,
        specimen_type=SpecimenType.BLOOD,
        output_dir=tmp_path,
        db_tier1=tmp_path / "db",
        host_reference=tmp_path / "h.fa",
        threads=4,
    )


def _hits(tmp_path):
    return [
        SketchHit("ecoli", 0.35, tmp_path / "ecoli.fa"),
        SketchHit("staph", 0.12, tmp_path / "staph.fa"),
    ]


def test_targeted_alignment_calls_minimap2(tmp_path):
    cfg = _cfg(tmp_path)
    nonhuman = tmp_path / "nonhuman.fastq.gz"
    hits = _hits(tmp_path)
    with patch("subprocess.run") as mock_run, \
         patch("pathogeniq.align._parse_paf", return_value=set()):
        mock_run.return_value = MagicMock(returncode=0)
        run_targeted_alignment(cfg, nonhuman, hits)
    calls = [c[0][0] for c in mock_run.call_args_list]
    assert any(c[0] == "minimap2" for c in calls)


def test_targeted_alignment_result_shape(tmp_path):
    cfg = _cfg(tmp_path)
    nonhuman = tmp_path / "nonhuman.fastq.gz"
    hits = _hits(tmp_path)
    # org 0 gets read_a, org 1 gets read_b, org 0 also gets read_c
    def fake_parse(paf_file):
        if "org_0" in str(paf_file):
            return {"read_a", "read_c"}
        return {"read_b"}

    with patch("subprocess.run"), \
         patch("pathogeniq.align._parse_paf", side_effect=fake_parse):
        result = run_targeted_alignment(cfg, nonhuman, hits)
    assert result.alignment_matrix.shape[1] == len(hits)
    assert result.organism_names == ["ecoli", "staph"]


def test_targeted_alignment_result_fields(tmp_path):
    cfg = _cfg(tmp_path)
    nonhuman = tmp_path / "nonhuman.fastq.gz"
    hits = _hits(tmp_path)
    with patch("subprocess.run"), \
         patch("pathogeniq.align._parse_paf", return_value=set()):
        result = run_targeted_alignment(cfg, nonhuman, hits)
    assert hasattr(result, "alignment_matrix")
    assert hasattr(result, "organism_names")
    assert hasattr(result, "read_ids")


def test_targeted_alignment_matrix_values(tmp_path):
    cfg = _cfg(tmp_path)
    nonhuman = tmp_path / "nonhuman.fastq.gz"
    hits = _hits(tmp_path)

    def fake_parse(paf_file):
        if "org_0" in str(paf_file):
            return {"read_a"}
        return {"read_b"}

    with patch("subprocess.run"), \
         patch("pathogeniq.align._parse_paf", side_effect=fake_parse):
        result = run_targeted_alignment(cfg, nonhuman, hits)

    assert result.alignment_matrix.shape == (2, 2)
    # read_a maps to org 0 only
    read_a_idx = result.read_ids.index("read_a")
    assert result.alignment_matrix[read_a_idx, 0] == 1.0
    assert result.alignment_matrix[read_a_idx, 1] == 0.0
