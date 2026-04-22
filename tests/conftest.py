import gzip
import pytest
from pathlib import Path


@pytest.fixture
def synthetic_sr_fastq(tmp_path) -> Path:
    """10 synthetic 150bp Illumina reads, all high quality."""
    reads = []
    for i in range(10):
        reads.append(f"@read_{i}")
        reads.append("ACGTACGTACGTACGT" * 9 + "ACGTACGTAC")  # 150 bp
        reads.append("+")
        reads.append("I" * 150)  # Q=40 throughout
    fastq = tmp_path / "synthetic_sr.fastq.gz"
    with gzip.open(fastq, "wt") as f:
        f.write("\n".join(reads) + "\n")
    return fastq


@pytest.fixture
def synthetic_lr_fastq(tmp_path) -> Path:
    """5 synthetic 1000bp Nanopore reads."""
    reads = []
    for i in range(5):
        reads.append(f"@read_lr_{i}")
        reads.append("ACGTACGTACGT" * 83 + "ACGT")  # 1000 bp
        reads.append("+")
        reads.append("I" * 1000)
    fastq = tmp_path / "synthetic_lr.fastq.gz"
    with gzip.open(fastq, "wt") as f:
        f.write("\n".join(reads) + "\n")
    return fastq
