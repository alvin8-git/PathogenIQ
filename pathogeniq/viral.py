"""Viral identification arm for airborne-pathogen surveillance.

Airborne transmission is dominated by viruses, yet the read-based bacterial arm
(sketch → align → EM) and the GTDB-Tk MAG arm are both blind to them. This arm
identifies viral sequences in the de novo assembly and classifies them:

    contigs ──► geNomad (viral identification + ICTV taxonomy) ──► viral contigs
                                                                       │
                                          CheckV (completeness / quality QC) ◄┘

geNomad scores each contig for viral origin (hallmark genes, marker enrichment)
and assigns ICTV lineage; CheckV estimates genome completeness so a near-complete
viral genome can be told apart from a short fragment. Both steps are non-blocking
(missing tool or its DB → degrades gracefully), like the AMR/assembly arms.

Source caveat: DNA metagenomics captures DNA viruses and integrated proviruses,
NOT RNA viruses (influenza, SARS-CoV-2, RSV) — those need an RNA-seq library. This
arm is source-agnostic: it classifies whatever contigs it is given, so RNA-virus
contigs are identified too IF the wet-lab produced a metatranscriptome.
"""
from __future__ import annotations

import csv
import io
import os
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path

from .config import PipelineConfig


@dataclass
class ViralContig:
    contig_id: str
    length: int
    taxonomy: str | None          # ICTV lineage from geNomad
    topology: str | None          # e.g. "DTR", "Provirus", "No terminal repeats"
    virus_score: float | None     # geNomad viral score (0-1)
    n_hallmarks: int | None       # count of viral hallmark genes
    completeness: float | None    # CheckV % completeness, None if CheckV absent
    checkv_quality: str | None    # CheckV quality tier, None if absent


def _to_float(v: str | None) -> float | None:
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def _to_int(v: str | None) -> int | None:
    try:
        return int(float(v))
    except (TypeError, ValueError):
        return None


def _parse_genomad_summary(text: str) -> dict[str, dict]:
    """Parse a geNomad ``*_virus_summary.tsv`` -> {seq_name: row fields}."""
    out: dict[str, dict] = {}
    for row in csv.DictReader(io.StringIO(text), delimiter="\t"):
        seq = (row.get("seq_name") or "").strip()
        if seq:
            out[seq] = row
    return out


def _parse_checkv_summary(text: str) -> dict[str, dict]:
    """Parse a CheckV ``quality_summary.tsv`` -> {contig_id: row fields}."""
    out: dict[str, dict] = {}
    for row in csv.DictReader(io.StringIO(text), delimiter="\t"):
        cid = (row.get("contig_id") or "").strip()
        if cid:
            out[cid] = row
    return out


def genomad_db_path() -> Path | None:
    """``$GENOMAD_DB`` if set, else ``databases/genomad_db``; None if missing."""
    env = os.environ.get("GENOMAD_DB")
    p = Path(env) if env else Path("databases/genomad_db")
    return p if p.exists() else None


def checkv_db_path() -> Path | None:
    """``$CHECKV_DB`` if set, else ``databases/checkv_db``; None if missing."""
    env = os.environ.get("CHECKV_DB")
    p = Path(env) if env else Path("databases/checkv_db")
    return p if p.exists() else None


def run_genomad(cfg: PipelineConfig, contigs: Path) -> Path | None:
    """Identify + classify viral contigs. Returns the virus_summary.tsv path, or
    None if geNomad/its DB is missing or the run fails."""
    db = genomad_db_path()
    if not shutil.which("genomad") or db is None:
        return None
    out = cfg.output_dir / "viral" / "genomad"
    out.mkdir(parents=True, exist_ok=True)
    try:
        subprocess.run(
            ["genomad", "end-to-end", "--cleanup", "-t", str(cfg.threads),
             str(contigs), str(out), str(db)],
            capture_output=True, check=True,
        )
    except subprocess.CalledProcessError:
        return None
    summaries = sorted(out.glob("*_summary/*_virus_summary.tsv"))
    return summaries[0] if summaries else None


def run_checkv(cfg: PipelineConfig, viral_fasta: Path) -> dict[str, dict]:
    """CheckV completeness/quality per viral contig. Empty if CheckV/DB absent."""
    db = checkv_db_path()
    if not shutil.which("checkv") or db is None:
        return {}
    out = cfg.output_dir / "viral" / "checkv"
    out.mkdir(parents=True, exist_ok=True)
    try:
        subprocess.run(
            ["checkv", "end_to_end", str(viral_fasta), str(out),
             "-d", str(db), "-t", str(cfg.threads)],
            capture_output=True, check=True,
        )
    except subprocess.CalledProcessError:
        return {}
    qs = out / "quality_summary.tsv"
    return _parse_checkv_summary(qs.read_text()) if qs.exists() else {}


def run_viral_stage(cfg: PipelineConfig, contigs: Path) -> list[ViralContig]:
    """Full viral arm: geNomad identify+classify → CheckV QC. Returns the viral
    contigs found (empty if geNomad/DB unavailable or none found). CheckV
    completeness is attached when available, else left None."""
    summary = run_genomad(cfg, contigs)
    if summary is None:
        return []
    records = _parse_genomad_summary(summary.read_text())
    if not records:
        return []
    viral_fasta = next(iter(summary.parent.glob("*_virus.fna")), None)
    quality = run_checkv(cfg, viral_fasta) if viral_fasta is not None else {}
    out: list[ViralContig] = []
    for seq, r in records.items():
        q = quality.get(seq, {})
        out.append(ViralContig(
            contig_id=seq,
            length=_to_int(r.get("length")) or 0,
            taxonomy=(r.get("taxonomy") or "").strip() or None,
            topology=(r.get("topology") or "").strip() or None,
            virus_score=_to_float(r.get("virus_score")),
            n_hallmarks=_to_int(r.get("n_hallmarks")),
            completeness=_to_float(q.get("completeness")),
            checkv_quality=(q.get("checkv_quality") or "").strip() or None,
        ))
    out.sort(key=lambda v: (v.completeness or -1, v.virus_score or -1), reverse=True)
    return out
