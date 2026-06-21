#!/usr/bin/env python
"""Build pathogeniq/data/background_default.tsv from public NTC FASTQs.

This is the curation recipe for the Tier-2 pooled background. It classifies each
no-template-control FASTQ through the same QC->host->sketch->align stages the
pipeline uses, pools the per-taxon RPM rates, and writes the default table so
NB background subtraction runs by default for users without a batch-matched NTC.

Suggested source data (download first; not bundled):
  - Salter et al. 2014, "Reagent and laboratory contamination can critically
    impact sequence-based microbiome analyses" — metagenomic reads at
    ENA ERP006808. The ultrapure-water negative controls are the kitome signal.
  - Additional NTC / reagent-blank samples from clinical-mNGS BioProjects.

Usage:
  python scripts/build_background_default.py \
      --ntc ntc1.fq.gz ntc2.fq.gz \
      --db databases/tier1/tier1_pathogens.sbt.zip \
      --host-ref databases/human_ref/GRCh38.fa \
      [--out pathogeniq/data/background_default.tsv] [--threads 16]

Requires the bioinformatics tools on PATH (fastp/Chopper, BWA-MEM2/minimap2,
samtools, sourmash) and the built Tier-1 DB + host reference.
"""
import argparse
from pathlib import Path

from pathogeniq.background import (
    build_background,
    default_background_path,
    write_background_table,
)
from pathogeniq.cli import _classify_taxon_counts
from pathogeniq.config import PipelineConfig, ReadType, SpecimenType


def main() -> None:
    ap = argparse.ArgumentParser(description="Build the Tier-2 pooled background table.")
    ap.add_argument("--ntc", nargs="+", required=True, help="NTC FASTQ file(s)")
    ap.add_argument("--db", required=True, help="Tier-1 sourmash SBT database")
    ap.add_argument("--host-ref", required=True, help="host reference (indexed)")
    ap.add_argument("--read-type", choices=[r.value for r in ReadType], default="short")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--out", default=None, help="output table (default: packaged path)")
    ap.add_argument("--workdir", default="/data/alvin/tmp/bg_build",
                    help="scratch dir for intermediate classification outputs")
    args = ap.parse_args()

    ntc_counts: list[tuple[dict[str, int], int]] = []
    for fq in args.ntc:
        cfg = PipelineConfig(
            input_fastq=Path(fq),
            read_type=ReadType(args.read_type),
            specimen_type=SpecimenType.BLOOD,   # specimen is irrelevant for NTC classification
            output_dir=Path(args.workdir) / Path(fq).stem,
            db_tier1=Path(args.db),
            host_reference=Path(args.host_ref),
            threads=args.threads,
        )
        counts, total = _classify_taxon_counts(cfg, Path(fq))
        print(f"{fq}: {len(counts)} taxa, {total} classified reads")
        ntc_counts.append((counts, total))

    model = build_background(ntc_counts, tier=2)
    if model is None:
        raise SystemExit("No NTC produced classified reads — nothing to pool.")

    out = Path(args.out) if args.out else default_background_path()
    write_background_table(out, model.rates, tier=2)
    print(f"Wrote {len(model.rates)} taxa to {out} (pooled from {model.n_controls} control(s))")


if __name__ == "__main__":
    main()
