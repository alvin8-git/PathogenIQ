import csv
import json
import re
import subprocess
from dataclasses import dataclass
from pathlib import Path

from .config import PipelineConfig


@dataclass
class SketchHit:
    name: str
    containment: float
    genome_path: Path
    taxon_id: str = ""   # stable GCF/GCA accession; join key for NTC background


_ACCESSION_RE = re.compile(r"(GC[AF]_\d+\.\d+)")


def _build_md5_map(db_path: Path) -> dict[str, tuple[str, str, Path]]:
    """Return md5 → (species_name, taxon_id, genome_path) from .sig files beside the SBT zip.

    taxon_id is the GCF/GCA accession parsed from the sig filename stem; it is the
    stable join key for NTC background subtraction (display names are not stable).
    """
    sig_dir = db_path.parent / "sigs"
    genome_dir = db_path.parent / "genomes"

    name_map: dict[str, str] = {}
    name_map_path = db_path.parent / "name_map.json"
    if name_map_path.exists():
        name_map = json.loads(name_map_path.read_text())

    md5_map: dict[str, tuple[str, Path]] = {}
    if not sig_dir.exists():
        return md5_map
    for sig_file in sig_dir.glob("*.sig"):
        try:
            data = json.loads(sig_file.read_text())
            for entry in data:
                for sig in entry.get("signatures", []):
                    md5 = sig.get("md5sum", "")
                    if not md5:
                        continue
                    stem = sig_file.stem
                    genome: Path = genome_dir / f"{stem}.fna"
                    if not genome.exists():
                        genome = genome_dir / f"{stem}.fna.gz"
                    if not genome.exists():
                        m_acc = _ACCESSION_RE.search(stem)
                        if m_acc:
                            acc = m_acc.group(1)
                            for pat in (f"{acc}*.fna", f"{acc}*.fna.gz"):
                                found = sorted(genome_dir.glob(pat))
                                if found:
                                    genome = found[0]
                                    break
                    m = _ACCESSION_RE.search(stem)
                    taxon_id = m.group(1) if m else ""
                    name = name_map.get(taxon_id, stem) if m else stem
                    md5_map[md5] = (name, taxon_id, genome)
        except Exception:
            continue
    return md5_map


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
    subprocess.run(sketch_cmd, capture_output=True, check=True)

    search_cmd = [
        "sourmash", "search",
        str(sample_sig),
        str(cfg.db_tier1),
        "--containment",
        "--threshold", str(cfg.sketch_threshold),
        "-o", str(results_csv),
    ]
    subprocess.run(search_cmd, capture_output=True)

    md5_map = _build_md5_map(cfg.db_tier1)

    hits: list[SketchHit] = []
    try:
        with open(results_csv) as f:
            reader = csv.DictReader(f)
            for row in reader:
                if "containment" in row:
                    containment = float(row["containment"])
                elif "similarity" in row:
                    containment = float(row["similarity"])
                else:
                    continue
                if containment < cfg.sketch_threshold:
                    continue
                md5 = row.get("md5", "")
                if md5 and md5 in md5_map:
                    name, taxon_id, genome_path = md5_map[md5]
                else:
                    # fallback: use CSV values as-is (e.g. pre-named signatures)
                    name = row.get("name", "")
                    genome_path = Path(row.get("filename", ""))
                    m_acc = _ACCESSION_RE.search(name)
                    taxon_id = m_acc.group(1) if m_acc else ""
                hits.append(SketchHit(
                    name=name,
                    containment=containment,
                    genome_path=genome_path,
                    taxon_id=taxon_id,
                ))
    except FileNotFoundError:
        pass
    return hits
