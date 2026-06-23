#!/usr/bin/env python
"""In-silico viral spike-in validation for the viral arm (geNomad + CheckV).

Our real air data cannot validate the viral arm: the aircraft filters are DNA-seq
(no RNA respiratory viruses) with no labelled viral truth, and the Zymo D6300 spike
has no viral members. So we generate gold truth synthetically: simulate reads from
known viral genomes, optionally mix them into a real aircraft-filter background, run
the viral arm, and score whether geNomad recovers each genome with the right ICTV
lineage (and CheckV's completeness vs the known genome).

    viral genomes (RefSeq) ──► wgsim reads ──(+ optional filter background)──►
        megahit ──► run_viral_stage (geNomad + CheckV) ──► score vs truth

This validates the *bioinformatic* arm (what we built). End-to-end air-virus
capability still needs real RNA-seq air virome — a wet-lab decision. See
docs/air-open-world-detection-2026-06-23.md.

Tools: wgsim (reads), megahit + genomad + checkv (+ DBs) for the arm. Non-blocking:
missing tools/DBs are reported and the run stops cleanly. The scoring logic runs
offline via ``--selfcheck`` (no external tools needed).

Usage:
    python scripts/12_viral_insilico_spikein.py --selfcheck
    python scripts/12_viral_insilico_spikein.py --out /data/alvin/tmp/viral_val \\
        [--background databases/aircraft/air-aircraft-filter/SRR32514319_1.fastq.gz] \\
        [--depth 30] [--threads 16]
"""
from __future__ import annotations

import argparse
import gzip
import shutil
import subprocess
import sys
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from pathogeniq.config import PipelineConfig, ReadType, SpecimenType  # noqa: E402
from pathogeniq.assembly import run_megahit  # noqa: E402
from pathogeniq.viral import run_viral_stage, ViralContig  # noqa: E402

# Truth set: (label, RefSeq accession, expected ICTV lineage token geNomad should
# emit). Tokens are at the most stable level available; note Lambda and T4 both sit
# in class Caudoviricetes, so the class-level token can alias between tailed phages
# (a per-genome contig tag is not recoverable after co-assembly).
# ponytail: 3 genomes spanning a plausible-in-air dsDNA phage, a second tailed
# phage, and a respiratory RNA virus (DNA reads from its genome — proves the arm
# would classify it IF present in contigs). Extend the list to widen coverage.
TRUTH = [
    ("Escherichia phage T4",      "NC_000866.4", "Straboviridae"),   # dsDNA tailed phage
    ("Escherichia phage Lambda",  "NC_001416.1", "Caudoviricetes"),  # dsDNA tailed phage (class)
    ("SARS-CoV-2",                "NC_045512.2", "Coronaviridae"),   # respiratory RNA virus
]

_EFETCH = ("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
           "?db=nuccore&id={acc}&rettype=fasta&retmode=text")


def fetch_genome(acc: str, dest: Path) -> Path:
    """Download a RefSeq genome FASTA by accession (cached)."""
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    with urllib.request.urlopen(_EFETCH.format(acc=acc)) as r:  # noqa: S310
        dest.write_bytes(r.read())
    return dest


def _genome_len(fasta: Path) -> int:
    return sum(len(line.strip()) for line in fasta.read_text().splitlines()
               if not line.startswith(">"))


def simulate_reads(genome: Path, depth: int, out_r1: Path, read_len: int = 150) -> int:
    """wgsim single-end-style: keep R1 only. Returns #reads written."""
    n_reads = max(1, _genome_len(genome) * depth // read_len)
    r1 = out_r1.with_suffix(".r1.fq")
    r2 = out_r1.with_suffix(".r2.fq")
    subprocess.run(
        ["wgsim", "-N", str(n_reads), "-1", str(read_len), "-2", str(read_len),
         "-S", "11", str(genome), str(r1), str(r2)],
        capture_output=True, check=True,
    )
    shutil.move(str(r1), str(out_r1))
    r2.unlink(missing_ok=True)
    return n_reads


def build_spiked_fastq(out_dir: Path, depth: int, background: Path | None) -> tuple[Path, int]:
    """Fetch genomes, simulate reads, concatenate (+ optional background) into one
    FASTQ. Returns (fastq path, total spiked reads)."""
    genomes_dir = out_dir / "genomes"
    reads_dir = out_dir / "reads"
    reads_dir.mkdir(parents=True, exist_ok=True)
    combined = out_dir / "spiked.fastq"
    total = 0
    with open(combined, "wb") as out:
        if background is not None:
            opener = gzip.open if str(background).endswith(".gz") else open
            with opener(background, "rb") as bg:
                shutil.copyfileobj(bg, out)
        for label, acc, _ in TRUTH:
            g = fetch_genome(acc, genomes_dir / f"{acc}.fasta")
            r1 = reads_dir / f"{acc}.fq"
            n = simulate_reads(g, depth, r1)
            total += n
            with open(r1, "rb") as fr:
                shutil.copyfileobj(fr, out)
            print(f"  simulated {n:>8,} reads from {label} ({acc})")
    return combined, total


def score_recovery(
    viral_contigs: list[ViralContig],
    truth: list[tuple[str, str, str]] = TRUTH,
) -> tuple[list[dict], float]:
    """For each truth genome, did geNomad recover a viral contig carrying the
    expected lineage token? Returns (per-genome results, recall)."""
    results = []
    for label, acc, token in truth:
        match = next(
            (v for v in viral_contigs
             if v.taxonomy and token.lower() in v.taxonomy.lower()),
            None,
        )
        results.append({
            "label": label,
            "accession": acc,
            "expected": token,
            "recovered": match is not None,
            "taxonomy": match.taxonomy if match else None,
            "completeness": match.completeness if match else None,
        })
    recall = sum(r["recovered"] for r in results) / len(results) if results else 0.0
    return results, recall


def _selfcheck() -> None:
    """Offline check of the scoring logic (no external tools)."""
    def vc(tax):
        return ViralContig(contig_id="c", length=1000, taxonomy=tax, topology=None,
                           virus_score=0.9, n_hallmarks=3, completeness=95.0,
                           checkv_quality="High-quality")
    contigs = [vc("Viruses;...;Straboviridae"), vc("Viruses;...;Coronaviridae")]
    results, recall = score_recovery(contigs)
    assert results[0]["recovered"] is True            # T4 -> Straboviridae
    assert results[2]["recovered"] is True            # SARS-CoV-2 -> Coronaviridae
    # Lambda token Caudoviricetes not present in either taxonomy string -> missed
    assert results[1]["recovered"] is False
    assert abs(recall - 2 / 3) < 1e-9
    # empty input -> zero recall, no crash
    _, r0 = score_recovery([])
    assert r0 == 0.0
    print("selfcheck OK: scoring logic verified offline")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=Path, help="Working/output dir")
    ap.add_argument("--background", type=Path, default=None,
                    help="Real aircraft-filter FASTQ to mix the spike into (optional)")
    ap.add_argument("--depth", type=int, default=30, help="Per-genome simulated coverage")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--selfcheck", action="store_true",
                    help="Verify scoring logic offline and exit")
    args = ap.parse_args()

    if args.selfcheck:
        _selfcheck()
        return
    if not args.out:
        ap.error("--out is required (or use --selfcheck)")

    missing = [t for t in ("wgsim", "megahit", "genomad", "checkv") if not shutil.which(t)]
    if missing:
        print(f"Missing tools: {', '.join(missing)}. Install them (+ geNomad/CheckV DBs) "
              f"to run the spike-in. Scoring logic is verifiable now with --selfcheck.")
        return

    args.out.mkdir(parents=True, exist_ok=True)
    print(f"[1/3] Building spiked FASTQ (depth={args.depth}x)"
          + (f", background={args.background.name}" if args.background else ", no background"))
    fastq, n_spike = build_spiked_fastq(args.out, args.depth, args.background)
    print(f"      {n_spike:,} spiked viral reads -> {fastq}")

    print("[2/3] Assembly + viral arm (geNomad + CheckV)...")
    cfg = PipelineConfig(
        input_fastq=fastq, read_type=ReadType.SHORT, specimen_type=SpecimenType.AIR,
        output_dir=args.out, db_tier1=Path("databases/tier1/tier1_pathogens.sbt.zip"),
        host_reference=Path("/data/alvin/ref/GRCh38/hg38.fa"), threads=args.threads,
    )
    contigs = run_megahit(cfg, fastq)
    viral = run_viral_stage(cfg, contigs) if contigs is not None else []
    print(f"      {len(viral)} viral contig(s) identified")

    print("[3/3] Scoring vs truth...")
    results, recall = score_recovery(viral)
    print(f"\n  {'genome':28s} {'expected':16s} {'recovered':10s} {'completeness':12s} taxonomy")
    for r in results:
        comp = "-" if r["completeness"] is None else f"{r['completeness']:.1f}"
        print(f"  {r['label']:28s} {r['expected']:16s} "
              f"{'YES' if r['recovered'] else 'no':10s} {comp:12s} {r['taxonomy'] or ''}")
    print(f"\n  RECALL: {recall:.0%} ({sum(r['recovered'] for r in results)}/{len(results)})")


if __name__ == "__main__":
    main()
