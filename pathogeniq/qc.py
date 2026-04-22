import json
import subprocess
from dataclasses import dataclass
from pathlib import Path

from .config import PipelineConfig, ReadType


@dataclass
class QCMetrics:
    total_reads: int
    passing_reads: int

    @property
    def pass_rate(self) -> float:
        return self.passing_reads / self.total_reads if self.total_reads else 0.0


def run_qc(cfg: PipelineConfig) -> tuple[Path, QCMetrics]:
    out = cfg.output_dir / "qc"
    out.mkdir(parents=True, exist_ok=True)
    filtered = out / "filtered.fastq.gz"
    json_report = out / "fastp.json"

    if cfg.read_type == ReadType.SHORT:
        cmd = [
            "fastp",
            "--in1", str(cfg.input_fastq),
            "--out1", str(filtered),
            "--json", str(json_report),
            "--qualified_quality_phred", "20",
            "--length_required", "50",
            "--low_complexity_filter",
            "--thread", str(cfg.threads),
            "--disable_adapter_trimming",
        ]
        subprocess.run(cmd, capture_output=True, text=True, check=True)
        with open(json_report) as f:
            data = json.load(f)
        metrics = QCMetrics(
            total_reads=data["summary"]["before_filtering"]["total_reads"],
            passing_reads=data["summary"]["after_filtering"]["total_reads"],
        )
    else:
        cmd = [
            "bash", "-c",
            f"gunzip -c {cfg.input_fastq} "
            f"| chopper -q 8 -l 200 --threads {cfg.threads} "
            f"| gzip > {filtered}",
        ]
        subprocess.run(cmd, capture_output=True, text=True, check=True)
        metrics = QCMetrics(total_reads=0, passing_reads=0)

    return filtered, metrics
