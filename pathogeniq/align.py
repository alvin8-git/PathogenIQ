import subprocess
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from .config import PipelineConfig, ReadType
from .sketch import SketchHit


@dataclass
class AlignmentResult:
    alignment_matrix: np.ndarray   # shape (n_reads, n_organisms)
    organism_names: list[str]
    read_ids: list[str]


def run_targeted_alignment(
    cfg: PipelineConfig,
    nonhuman_fastq: Path,
    hits: list[SketchHit],
) -> AlignmentResult:
    out = cfg.output_dir / "align"
    out.mkdir(parents=True, exist_ok=True)

    all_read_ids: list[str] = []
    org_reads: dict[int, set[str]] = {i: set() for i in range(len(hits))}

    preset = "sr" if cfg.read_type == ReadType.SHORT else "map-ont"

    for org_idx, hit in enumerate(hits):
        paf_file = out / f"org_{org_idx}.paf"
        cmd = [
            "minimap2",
            "-c",
            "-x", preset,
            "-t", str(cfg.threads),
            str(hit.genome_path),
            str(nonhuman_fastq),
            "-o", str(paf_file),
        ]
        subprocess.run(cmd, capture_output=True, text=True, check=True)
        read_set = _parse_paf(paf_file)
        org_reads[org_idx] = read_set
        for r in read_set:
            if r not in all_read_ids:
                all_read_ids.append(r)

    n_reads = len(all_read_ids)
    n_orgs = len(hits)
    matrix = np.zeros((n_reads, n_orgs), dtype=float)
    read_idx = {rid: i for i, rid in enumerate(all_read_ids)}

    for org_idx, read_set in org_reads.items():
        for rid in read_set:
            if rid in read_idx:
                matrix[read_idx[rid], org_idx] = 1.0

    return AlignmentResult(
        alignment_matrix=matrix,
        organism_names=[h.name for h in hits],
        read_ids=all_read_ids,
    )


def _parse_paf(paf_file: Path) -> set[str]:
    """Return set of read IDs that aligned (mapq > 0)."""
    read_ids: set[str] = set()
    if not paf_file.exists():
        return read_ids
    with open(paf_file) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 12:
                continue
            mapq = int(parts[11])
            if mapq > 0:
                read_ids.add(parts[0])
    return read_ids
