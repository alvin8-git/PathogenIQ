import subprocess
from dataclasses import dataclass
from pathlib import Path

from .config import PipelineConfig, ReadType


@dataclass
class HostRemovalMetrics:
    total_reads: int
    human_reads: int
    nonhuman_reads: int

    @property
    def microbial_fraction(self) -> float:
        return self.nonhuman_reads / self.total_reads if self.total_reads else 0.0

    @property
    def human_fraction(self) -> float:
        return self.human_reads / self.total_reads if self.total_reads else 0.0


def run_host_removal(cfg: PipelineConfig, filtered_fastq: Path) -> tuple[Path, HostRemovalMetrics]:
    out = cfg.output_dir / "host_removal"
    out.mkdir(parents=True, exist_ok=True)
    nonhuman_fastq = out / "nonhuman.fastq.gz"
    sam_file = out / "host_aligned.sam"

    if cfg.read_type == ReadType.SHORT:
        align_cmd = [
            "bwa-mem2", "mem",
            "-t", str(cfg.threads),
            str(cfg.host_reference),
            str(filtered_fastq),
            "-o", str(sam_file),
        ]
    else:
        align_cmd = [
            "minimap2",
            "-ax", "map-ont",
            "-t", str(cfg.threads),
            str(cfg.host_reference),
            str(filtered_fastq),
            "-o", str(sam_file),
        ]

    subprocess.run(align_cmd, capture_output=True, text=True, check=True)

    extract_cmd = [
        "bash", "-c",
        f"samtools view -f 4 -b {sam_file} "
        f"| samtools fastq - "
        f"| gzip > {nonhuman_fastq}",
    ]
    subprocess.run(extract_cmd, capture_output=True, text=True, check=True)

    flagstat = subprocess.run(
        ["samtools", "flagstat", str(sam_file)],
        capture_output=True, text=True, check=True,
    ).stdout
    total_reads = 0
    mapped_reads = 0
    for line in flagstat.splitlines():
        if "in total" in line:
            total_reads = int(line.split()[0])
        elif "mapped" in line and "primary mapped" not in line and "%" in line:
            mapped_reads = int(line.split()[0])
    nonhuman_reads = total_reads - mapped_reads
    metrics = HostRemovalMetrics(
        total_reads=total_reads,
        human_reads=mapped_reads,
        nonhuman_reads=nonhuman_reads,
    )
    return nonhuman_fastq, metrics
