import csv
import subprocess
from dataclasses import dataclass
from pathlib import Path

from .config import PipelineConfig


@dataclass
class SketchHit:
    name: str
    containment: float
    genome_path: Path


def run_sketch_screen(cfg: PipelineConfig, nonhuman_fastq: Path) -> list[SketchHit]:
    out = cfg.output_dir / "sketch"
    out.mkdir(parents=True, exist_ok=True)
    sample_sig = out / "sample.sig"
    results_csv = out / "results.csv"

    sketch_cmd = [
        "sourmash", "sketch", "dna",
        "-p", "k=31,scaled=1000",
        str(nonhuman_fastq),
        "-o", str(sample_sig),
    ]
    subprocess.run(sketch_cmd, capture_output=True, text=True, check=True)

    search_cmd = [
        "sourmash", "search",
        str(sample_sig),
        str(cfg.db_tier1),
        "--containment",
        "--threshold", str(cfg.sketch_threshold),
        "-o", str(results_csv),
    ]
    subprocess.run(search_cmd, capture_output=True, text=True)

    hits: list[SketchHit] = []
    try:
        with open(results_csv) as f:
            reader = csv.DictReader(f)
            for row in reader:
                # sourmash may output 'containment' or 'similarity' column depending on version
                if "containment" in row:
                    containment = float(row["containment"])
                elif "similarity" in row:
                    containment = float(row["similarity"])
                else:
                    continue
                if containment >= cfg.sketch_threshold:
                    hits.append(SketchHit(
                        name=row["name"],
                        containment=containment,
                        genome_path=Path(row["filename"]),
                    ))
    except FileNotFoundError:
        pass
    return hits
