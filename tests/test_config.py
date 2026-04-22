from pathlib import Path
import pytest
from pathogeniq.config import PipelineConfig, ReadType, SpecimenType


def test_pipeline_config_defaults(tmp_path):
    cfg = PipelineConfig(
        input_fastq=tmp_path / "sample.fastq.gz",
        read_type=ReadType.SHORT,
        specimen_type=SpecimenType.BLOOD,
        output_dir=tmp_path / "out",
        db_tier1=tmp_path / "db",
        host_reference=tmp_path / "human.fa",
    )
    assert cfg.threads == 8
    assert cfg.sketch_threshold == 0.01


def test_read_type_enum_values():
    assert ReadType.SHORT == "short"
    assert ReadType.LONG == "long"


def test_specimen_type_enum_values():
    assert SpecimenType.BLOOD == "blood"
    assert SpecimenType.CSF == "csf"
    assert SpecimenType.BAL == "bal"
    assert SpecimenType.TISSUE == "tissue"


def test_pipeline_config_custom_threads(tmp_path):
    cfg = PipelineConfig(
        input_fastq=tmp_path / "s.fastq.gz",
        read_type=ReadType.LONG,
        specimen_type=SpecimenType.CSF,
        output_dir=tmp_path,
        db_tier1=tmp_path,
        host_reference=tmp_path / "h.fa",
        threads=32,
    )
    assert cfg.threads == 32
