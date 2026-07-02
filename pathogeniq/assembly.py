"""Optional assembly + MAG recovery arm (Plan 6 #3).

The targeted read-based pipeline (sketch → align → EM) can only find organisms
already in the reference DB. Air and other environmental samples are full of
novel/uncultured taxa that have no reference. This arm assembles the non-host
reads de novo and bins them into metagenome-assembled genomes (MAGs), so genuinely
unknown organisms can be recovered and classified by genome rather than by read
match:

    MEGAHIT (assemble) → minimap2+jgi depth → MetaBAT2 (bin)
        → CheckM (completeness/contamination QC) → GTDB-Tk (genome taxonomy)

Every step is **non-blocking**: a missing tool (or its reference database — CheckM
and GTDB-Tk need multi-GB DBs) degrades gracefully rather than failing the run.
CheckM/GTDB-Tk are optional refinements — bins are still returned without them.
This arm is expensive and off by default (CLI ``--assemble``). Prokka gene
annotation and MAG-level ARG screening are deferred (see todo Plan 6).
"""
from __future__ import annotations

import csv
import io
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

from .config import PipelineConfig, ReadType


@dataclass
class MAG:
    bin_id: str
    fasta_path: Path
    completeness: float | None   # None if CheckM unavailable
    contamination: float | None
    taxonomy: str | None         # GTDB lineage string, None if GTDB-Tk unavailable
    n_contigs: int
    total_bp: int


def _fasta_stats(path: Path) -> tuple[int, int]:
    """(#contigs, total bp) of a bin FASTA."""
    n_contigs = total_bp = 0
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                n_contigs += 1
            else:
                total_bp += len(line.strip())
    return n_contigs, total_bp


def _parse_checkm_table(text: str) -> dict[str, tuple[float, float]]:
    """Parse CheckM ``--tab_table`` output -> {bin_id: (completeness, contamination)}."""
    out: dict[str, tuple[float, float]] = {}
    for row in csv.DictReader(io.StringIO(text), delimiter="\t"):
        bid = (row.get("Bin Id") or row.get("Bin ID") or "").strip()
        if not bid:
            continue
        try:
            out[bid] = (float(row.get("Completeness", 0)), float(row.get("Contamination", 0)))
        except ValueError:
            continue
    return out


def _parse_gtdbtk_summary(text: str) -> dict[str, str]:
    """Parse a GTDB-Tk summary.tsv -> {bin_id: classification lineage}."""
    out: dict[str, str] = {}
    for row in csv.DictReader(io.StringIO(text), delimiter="\t"):
        genome = (row.get("user_genome") or "").strip()
        if genome:
            out[genome] = (row.get("classification") or "").strip()
    return out


def run_megahit(cfg: PipelineConfig, reads: Path) -> Path | None:
    """Assemble non-host reads into contigs. None if MEGAHIT is absent or fails."""
    if not shutil.which("megahit"):
        return None
    out = cfg.output_dir / "assembly" / "megahit"
    # megahit's own mkdir is single-level (os.mkdir) and dies if the parent is
    # missing; it also wants to create `out` itself under --force, so make the
    # PARENT but not `out`.
    out.parent.mkdir(parents=True, exist_ok=True)
    try:
        subprocess.run(
            ["megahit", "-r", str(reads), "-o", str(out), "-t", str(cfg.threads), "--force"],
            capture_output=True, check=True,
        )
    except subprocess.CalledProcessError:
        return None
    contigs = out / "final.contigs.fa"
    return contigs if contigs.exists() else None


def _contig_depths(cfg: PipelineConfig, contigs: Path, reads: Path) -> Path | None:
    """Map reads back to contigs and summarise per-contig depth for MetaBAT2."""
    if not all(shutil.which(t) for t in ("minimap2", "samtools", "jgi_summarize_bam_contig_depths")):
        return None
    out = cfg.output_dir / "assembly"
    out.mkdir(parents=True, exist_ok=True)
    bam = out / "aln.bam"
    depth = out / "depth.txt"
    preset = "sr" if cfg.read_type == ReadType.SHORT else "map-ont"
    try:
        subprocess.run(
            ["bash", "-c",
             f"minimap2 -ax {preset} -t {cfg.threads} {contigs} {reads} "
             f"| samtools sort -@ {cfg.threads} -o {bam}"],
            capture_output=True, check=True,
        )
        subprocess.run(["samtools", "index", str(bam)], capture_output=True, check=True)
        subprocess.run(["jgi_summarize_bam_contig_depths", "--outputDepth", str(depth), str(bam)],
                       capture_output=True, check=True)
    except subprocess.CalledProcessError:
        return None
    return depth if depth.exists() else None


def run_metabat2(cfg: PipelineConfig, contigs: Path, depth: Path) -> list[Path]:
    """Bin contigs into MAGs. Empty list if MetaBAT2 absent or no bins formed."""
    if not shutil.which("metabat2"):
        return []
    bindir = cfg.output_dir / "assembly" / "bins"
    bindir.mkdir(parents=True, exist_ok=True)
    try:
        subprocess.run(
            ["metabat2", "-i", str(contigs), "-a", str(depth),
             "-o", str(bindir / "bin"), "-t", str(cfg.threads)],
            capture_output=True, check=True,
        )
    except subprocess.CalledProcessError:
        return []
    return sorted(bindir.glob("bin.*.fa"))


def run_checkm(cfg: PipelineConfig, bins: list[Path]) -> dict[str, tuple[float, float]]:
    """CheckM completeness/contamination per bin. Empty if CheckM/DB unavailable."""
    if not bins or not shutil.which("checkm"):
        return {}
    bindir = bins[0].parent
    out = cfg.output_dir / "assembly" / "checkm"
    tsv = out / "checkm.tsv"
    try:
        subprocess.run(
            ["checkm", "lineage_wf", "-x", "fa", "--tab_table", "-f", str(tsv),
             "-t", str(cfg.threads), str(bindir), str(out)],
            capture_output=True, check=True,
        )
    except subprocess.CalledProcessError:
        return {}
    return _parse_checkm_table(tsv.read_text()) if tsv.exists() else {}


def run_gtdbtk(cfg: PipelineConfig, bins: list[Path]) -> dict[str, str]:
    """GTDB-Tk genome taxonomy per bin. Empty if GTDB-Tk/DB unavailable."""
    if not bins or not shutil.which("gtdbtk"):
        return {}
    bindir = bins[0].parent
    out = cfg.output_dir / "assembly" / "gtdbtk"
    try:
        # ponytail: no --skip_ani_screen — gtdbtk 2.7.2 removed it (older versions
        # required it); 2.7.2 runs the skani ANI screen by default (skani DB present).
        subprocess.run(
            ["gtdbtk", "classify_wf", "--genome_dir", str(bindir), "--out_dir", str(out),
             "-x", "fa", "--cpus", str(cfg.threads)],
            capture_output=True, check=True,
        )
    except subprocess.CalledProcessError as e:
        # Non-blocking, but do NOT swallow silently: a bad flag / OOM here would
        # otherwise leave every MAG unnamed with no trace (this exact silent swallow
        # hid the removed-flag bug for a whole air run).
        err = (e.stderr or b"").decode("utf-8", "replace").strip().splitlines()
        print(f"WARN: gtdbtk classify_wf failed, MAGs left unnamed: "
              f"{err[-1] if err else f'exit {e.returncode}'}", file=sys.stderr)
        return {}
    taxa: dict[str, str] = {}
    for name in ("gtdbtk.bac120.summary.tsv", "gtdbtk.ar53.summary.tsv"):
        p = out / name
        if p.exists():
            taxa.update(_parse_gtdbtk_summary(p.read_text()))
    return taxa


def run_assembly_stage(
    cfg: PipelineConfig,
    reads: Path,
    *,
    min_completeness: float = 50.0,
    max_contamination: float = 10.0,
) -> list[MAG]:
    """Full assembly arm: assemble → bin → QC → classify. Returns high-quality
    MAGs (CheckM completeness ≥ ``min_completeness`` and contamination ≤
    ``max_contamination`` — the paper's thresholds). When CheckM is unavailable
    the quality filter is skipped (every bin is kept, completeness=None) rather
    than dropping all bins. Returns [] if assembly/binning produced nothing."""
    contigs = run_megahit(cfg, reads)
    if contigs is None:
        return []
    depth = _contig_depths(cfg, contigs, reads)
    bins = run_metabat2(cfg, contigs, depth) if depth is not None else []
    if not bins:
        return []
    quality = run_checkm(cfg, bins)
    taxa = run_gtdbtk(cfg, bins)
    mags: list[MAG] = []
    for b in bins:
        bid = b.stem
        completeness, contamination = quality.get(bid, (None, None))
        if completeness is not None and (
            completeness < min_completeness or contamination > max_contamination
        ):
            continue   # CheckM ran and this bin failed QC
        n_contigs, total_bp = _fasta_stats(b)
        mags.append(MAG(
            bin_id=bid, fasta_path=b, completeness=completeness, contamination=contamination,
            taxonomy=taxa.get(bid), n_contigs=n_contigs, total_bp=total_bp,
        ))
    return mags
