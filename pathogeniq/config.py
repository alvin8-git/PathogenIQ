from dataclasses import dataclass
from enum import Enum
from pathlib import Path


class ReadType(str, Enum):
    SHORT = "short"
    LONG = "long"


class SpecimenType(str, Enum):
    BLOOD = "blood"
    CSF = "csf"
    BAL = "bal"
    TISSUE = "tissue"
    AIR = "air"      # environmental air / bioaerosol (filters, impingers, masks)


@dataclass
class PipelineConfig:
    input_fastq: Path
    read_type: ReadType
    specimen_type: SpecimenType
    output_dir: Path
    db_tier1: Path
    host_reference: Path
    threads: int = 8
    sketch_threshold: float = 0.003
    n_bootstrap: int = 100
    amr_db: str = "card"
