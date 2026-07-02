#!/usr/bin/env python
"""Pool genuine non-spiked NTC/blank datasets into the Tier-2 background.

The shipped background_default.tsv was built from one spiked dilution series
(Salter 2014), which makes it idiosyncratic and Enterobacteriaceae-heavy, and
the dispersion validation (docs/dispersion-validation-2026-06-22.md) showed the
real limit is NTC *coverage* — run-unique contaminants hitting the pseudocount
floor. This pools several genuinely non-spiked blank studies (more labs/kits ->
broader coverage, no spike confound) into a cleaner per-taxon background.

Each blank FASTQ is classified against the SAME Tier-1 DB as the existing table
(GCF-keyed), so rates pool correctly. Per-run counts are cached as counts.json
so re-pooling with different thresholds is cheap. No spike filtering — these are
real blanks; every classified read is background.

    python scripts/11_pool_blank_background.py \\
        --study hunt:databases/validation/blanks-hunt \\
        --study qiita:databases/validation/blanks-qiita \\
        --db databases/tier1/tier1_pathogens.sbt.zip \\
        --host-ref /data/alvin/ref/GRCh38/hg38.fa \\
        --out pathogeniq/data/background_default.tsv
"""
import argparse
import json
from pathlib import Path

from pathogeniq.background import build_background, write_background_table
from pathogeniq.cli import _classify_taxon_counts
from pathogeniq.config import PipelineConfig, ReadType, SpecimenType


def classify_run(fq: Path, *, db: Path, host_ref: Path, workdir: Path,
                 read_type: str, threads: int) -> tuple[dict[str, int], int]:
    """Classify one blank (R1) to GCF-keyed counts, caching counts.json."""
    out = workdir / fq.stem
    out.mkdir(parents=True, exist_ok=True)
    cache = out / "counts.json"
    if cache.exists():
        d = json.loads(cache.read_text())
        return {k: int(v) for k, v in d["counts"].items()}, int(d["total"])
    cfg = PipelineConfig(
        input_fastq=fq, read_type=ReadType(read_type), specimen_type=SpecimenType.BLOOD,
        output_dir=out, db_tier1=db, host_reference=host_ref, threads=threads,
    )
    counts, total = _classify_taxon_counts(cfg, fq)
    cache.write_text(json.dumps({"counts": counts, "total": total}))
    return counts, total


def study_runs(spec: str) -> tuple[str, list[Path]]:
    """'label:dir' -> (label, sorted R1 fastqs in dir)."""
    label, _, d = spec.partition(":")
    fqs = sorted(Path(d).glob("*_1.fastq.gz")) or sorted(Path(d).glob("*_1.fastq"))
    return label, fqs


def main() -> None:
    ap = argparse.ArgumentParser(description="Pool non-spiked blanks into the Tier-2 background.")
    ap.add_argument("--one", type=Path, default=None,
                    help="classify a single R1 fastq into its counts.json cache and exit "
                         "(for parallel pre-population via xargs -P, then pool with --study)")
    ap.add_argument("--label", default="misc", help="study label for --one's cache subdir")
    ap.add_argument("--study", action="append", default=[],
                    help="label:dir (R1 *_1.fastq[.gz] in dir are the blank runs); repeatable")
    ap.add_argument("--include-counts", action="append", default=[],
                    help="extra dir of cached <run>/counts.json to pool in (e.g. Salter spike-free); repeatable")
    ap.add_argument("--db", required=True, type=Path)
    ap.add_argument("--host-ref", required=True, type=Path)
    ap.add_argument("--workdir", default=Path("/data/alvin/tmp/blank_pool"), type=Path)
    ap.add_argument("--out", default=Path("pathogeniq/data/background_default.tsv"), type=Path)
    ap.add_argument("--read-type", choices=[r.value for r in ReadType], default="short")
    ap.add_argument("--min-reads", type=int, default=2)
    ap.add_argument("--threads", type=int, default=16)
    ap.add_argument("--dry-run", action="store_true", help="list runs per study, classify nothing")
    args = ap.parse_args()

    if args.one is not None:
        counts, total = classify_run(args.one, db=args.db, host_ref=args.host_ref,
                                     workdir=args.workdir / args.label, read_type=args.read_type,
                                     threads=args.threads)
        print(f"{args.one.name}: total={total}, nonzero taxa={sum(1 for v in counts.values() if v > 0)}")
        return

    pooled: list[tuple[dict[str, int], int]] = []
    prov: list[str] = []
    for spec in args.study:
        label, fqs = study_runs(spec)
        print(f"[{label}] {len(fqs)} blank run(s)")
        if args.dry_run:
            continue
        usable = 0
        for fq in fqs:
            counts, total = classify_run(fq, db=args.db, host_ref=args.host_ref,
                                         workdir=args.workdir / label, read_type=args.read_type,
                                         threads=args.threads)
            if total > 0:
                pooled.append((counts, total)); usable += 1
        prov.append(f"{label}: {usable} blank runs")
    for cdir in args.include_counts:
        n = 0
        for cache in sorted(Path(cdir).glob("*/counts.json")):
            d = json.loads(cache.read_text())
            if int(d["total"]) > 0:
                pooled.append(({k: int(v) for k, v in d["counts"].items()}, int(d["total"]))); n += 1
        prov.append(f"{Path(cdir).name}: {n} cached runs")
    if args.dry_run:
        return

    model = build_background(pooled, tier=2, min_reads=args.min_reads)
    if model is None:
        raise SystemExit("no usable blanks — nothing pooled")
    notes = [
        "PROVENANCE: pooled genuine NON-SPIKED negative-control shotgun datasets",
        "  (cleaner/broader than the prior Salter-only spiked-series build):",
        *[f"    {p}" for p in prov],
        f"  {model.n_controls} usable controls, {len(model.rates)} taxa, min_reads={args.min_reads}.",
        "  Classified against the Tier-1 clinical DB (GCF-keyed). Tier 2 = pooled,",
        "  foreign prior -> grades capped at B. A per-batch NTC (Tier 1) is stronger;",
        "  see docs/dispersion-validation-2026-06-22.md (coverage, not dispersion, is the limit).",
    ]
    write_background_table(args.out, model.rates, tier=2, notes=notes)
    print(f"\nwrote {len(model.rates)} taxa from {model.n_controls} controls -> {args.out}")


if __name__ == "__main__":
    main()
