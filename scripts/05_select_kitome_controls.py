#!/usr/bin/env python
"""Select reagent-dominated (kitome) controls from a spiked dilution series.

The Salter 2014 set (ENA ERP006808) is a *Salmonella bongori* serial dilution,
not a set of blanks: high-biomass runs are ~98% spike, but as biomass drops the
reagent/kitome contamination takes over. The low-spike runs are therefore real
reagent-background controls. This turns one water blank into ~10-20 usable
controls from data already on disk.

For each input FASTQ this script classifies it through the same QC->host->
sketch->align stages, computes the spike fraction (reads assigned to the spike
organism / total classified), keeps runs at or below --max-spike-frac, and pools
the survivors into a Tier-2 background table.

Usage:
  # report only (no build) — see the spike fraction per run
  python scripts/05_select_kitome_controls.py --ntc databases/validation/salter-ntc/*_1.fastq.gz \\
      --db databases/tier1/tier1_pathogens.sbt.zip --host-ref databases/human_ref/GRCh38.fa --dry-run

  # select low-spike runs and write the pooled background
  python scripts/05_select_kitome_controls.py --ntc databases/validation/salter-ntc/*_1.fastq.gz \\
      --db databases/tier1/tier1_pathogens.sbt.zip --host-ref databases/human_ref/GRCh38.fa \\
      --max-spike-frac 0.10 --min-reads 2 --threads 16

Note: this assumes a SINGLE known spike organism (the Salter design). For a
dataset with multiple spikes, pass each as a repeated --spike-taxon prefix.
"""
import argparse
from pathlib import Path

from pathogeniq.background import build_background, default_background_path, write_background_table
from pathogeniq.cli import _classify_taxon_counts
from pathogeniq.config import PipelineConfig, ReadType, SpecimenType

# Salmonella enterica reference (the S. bongori spike maps here). Prefix match is
# version-agnostic (GCF_000006945.2 etc.).
_DEFAULT_SPIKE = "GCF_000006945"


def spike_fraction(counts: dict[str, int], total: int, spike_prefixes: list[str]) -> float:
    """Fraction of classified reads assigned to any spike organism."""
    if total <= 0:
        return 0.0
    spike = sum(c for tx, c in counts.items()
                if any(tx.startswith(p) for p in spike_prefixes))
    return spike / total


def main() -> None:
    ap = argparse.ArgumentParser(description="Select kitome controls from a spiked dilution series.")
    ap.add_argument("--ntc", nargs="+", required=True, help="dilution-series FASTQ file(s)")
    ap.add_argument("--db", required=True)
    ap.add_argument("--host-ref", required=True)
    ap.add_argument("--spike-taxon", action="append", default=None,
                    help=f"spike taxon_id prefix (default {_DEFAULT_SPIKE}); repeatable")
    ap.add_argument("--max-spike-frac", type=float, default=0.10,
                    help="keep runs with spike fraction <= this (reagent-dominated)")
    ap.add_argument("--min-reads", type=int, default=2)
    ap.add_argument("--read-type", choices=[r.value for r in ReadType], default="short")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--out", default=None)
    ap.add_argument("--workdir", default="/data/alvin/tmp/kitome_select")
    ap.add_argument("--dry-run", action="store_true", help="report spike fractions, do not build")
    args = ap.parse_args()

    spikes = args.spike_taxon or [_DEFAULT_SPIKE]

    classified: list[tuple[str, dict[str, int], int, float]] = []
    for fq in args.ntc:
        cfg = PipelineConfig(
            input_fastq=Path(fq),
            read_type=ReadType(args.read_type),
            specimen_type=SpecimenType.BLOOD,
            output_dir=Path(args.workdir) / Path(fq).stem,
            db_tier1=Path(args.db),
            host_reference=Path(args.host_ref),
            threads=args.threads,
        )
        counts, total = _classify_taxon_counts(cfg, Path(fq))
        frac = spike_fraction(counts, total, spikes)
        classified.append((Path(fq).name, counts, total, frac))

    print(f"\n{'run':40s} {'reads':>8s} {'spike%':>7s}  keep")
    selected: list[tuple[dict[str, int], int]] = []
    for name, counts, total, frac in sorted(classified, key=lambda x: x[3]):
        keep = total > 0 and frac <= args.max_spike_frac
        print(f"{name:40s} {total:8d} {frac*100:6.1f}%  {'YES' if keep else 'no'}")
        if keep:
            selected.append((counts, total))

    print(f"\n{len(selected)}/{len(classified)} runs selected "
          f"(spike fraction <= {args.max_spike_frac:.0%})")
    if args.dry_run:
        return
    if not selected:
        raise SystemExit("No runs below the spike-fraction threshold — relax --max-spike-frac.")

    model = build_background(selected, tier=2, min_reads=args.min_reads)
    if model is None:
        raise SystemExit("Selected runs produced no taxa above the support floor.")
    out = Path(args.out) if args.out else default_background_path()
    write_background_table(out, model.rates, tier=2)
    print(f"Wrote {len(model.rates)} taxa to {out} (pooled from {model.n_controls} kitome control(s))")


if __name__ == "__main__":
    main()
