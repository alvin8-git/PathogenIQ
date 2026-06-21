#!/usr/bin/env python
"""Run the grading benchmark on a Kraken2 report against a truth set (T10).

Compares two configs on one labeled community:
  - kraken-raw : every Kraken2 species call, ranked by read count (no filtering)
  - kraken+grading : the same, but taxa graded X (below the specimen read floor /
                     degenerate stats) are dropped, and the tier cap applies

If grading removes low-support false positives, its precision-at-fixed-sensitivity
and PR-AUC (average precision) beat raw. That is the benchmark claim.

This script CONSUMES a Kraken2 report — it does not run Kraken2 (that needs a
Kraken2 DB + the reads, your step). Produce the report with:
    kraken2 --db <DB> --report sample.kreport --output - sample.fastq

Truth file: one NCBI taxid per line = the species actually present (from the CAMI
gold-standard profile, or the known ZymoBIOMICS composition).

Usage:
    python scripts/06_benchmark.py --kraken-report sample.kreport \\
        --truth truth_taxids.txt --specimen blood --recall 0.95

NOT YET WIRED (documented blockers):
  - CZID baseline: add a --czid-report config once a CZID export parser exists.
  - NB background filter: the shipped background is keyed by GCF accession but
    Kraken2 reports use NCBI taxids; the background join needs a taxid-keyed
    background (rebuild the kitome controls through Kraken2) or a GCF<->taxid map.
    Until then the graded config applies the read-floor / CI / tier rules only.
"""
import argparse
from pathlib import Path

from pathogeniq.background import is_background, load_background_table
from pathogeniq.benchmark import (
    average_precision,
    kraken_to_grading_inputs,
    parse_kraken2_report,
    precision_at_recall,
    precision_recall,
)
from pathogeniq.config import SpecimenType
from pathogeniq.report import EvidenceGrade, grade


def _score(label, scored, predicted, truth, recall):
    ap = average_precision(scored, truth)
    par = precision_at_recall(scored, truth, recall)
    p, r = precision_recall(predicted, truth)
    tp = len(predicted & truth)
    print(f"{label:18s} {ap:7.3f} {par:11.3f} {p:7.3f} {r:7.3f} "
          f"{tp:4d} {len(predicted - truth):4d} {len(truth - predicted):4d}")
    return ap, p


def main() -> None:
    ap = argparse.ArgumentParser(description="Grading benchmark on a Kraken2 report (T10).")
    ap.add_argument("--kraken-report", required=True, type=Path)
    ap.add_argument("--truth", required=True, type=Path, help="one NCBI taxid per line")
    ap.add_argument("--specimen", required=True, choices=[s.value for s in SpecimenType])
    ap.add_argument("--tier", type=int, default=3, help="NTC tier for the graded config (1/2/3)")
    ap.add_argument("--background", type=Path, default=None,
                    help="taxid-keyed background table; applies the NB filter to the graded config")
    ap.add_argument("--recall", type=float, default=0.95, help="fixed sensitivity for precision@recall")
    args = ap.parse_args()

    taxa = parse_kraken2_report(args.kraken_report.read_text())
    truth = {ln.strip() for ln in args.truth.read_text().splitlines() if ln.strip()}
    specimen = SpecimenType(args.specimen)
    print(f"{len(taxa)} species in report, {len(truth)} in truth set\n")
    print(f"{'config':18s} {'PR-AUC':>7s} {'P@%.0f%%R' % (args.recall * 100):>11s} "
          f"{'prec':>7s} {'rec':>7s} {'TP':>4s} {'FP':>4s} {'FN':>4s}")

    # config 1: raw Kraken2 — every species, ranked by reads, no filtering
    raw_scored = [(t.taxid, float(t.reads)) for t in taxa]
    raw_pred = {t.taxid for t in taxa}
    raw_ap, raw_p = _score("kraken-raw", raw_scored, raw_pred, truth, args.recall)

    # config 2: Kraken2 + grading — NB background filter (if given) + drop Grade-X,
    # tier cap applied
    bg = load_background_table(args.background) if args.background else None
    tier = bg.tier if bg is not None else args.tier
    total = sum(t.reads for t in taxa)
    gi = kraken_to_grading_inputs(taxa, specimen=specimen, tier=tier)
    kept = set()
    for t in taxa:
        if bg is not None and is_background(t.taxid, t.reads, total, bg):
            continue   # indistinguishable from NTC background -> drop
        if grade(gi[t.taxid]) != EvidenceGrade.X:
            kept.add(t.taxid)
    label = "kraken+grading+bg" if bg is not None else "kraken+grading"
    graded_scored = [(t.taxid, float(t.reads)) for t in taxa if t.taxid in kept]
    graded_ap, graded_p = _score(label, graded_scored, kept, truth, args.recall)

    # Grading's value is removing low-support false positives, which moves
    # operating-point precision; PR-AUC can be flat when those FPs are already
    # bottom-ranked by read count. Report both.
    print(f"\nprecision delta (grading - raw): {graded_p - raw_p:+.3f}  "
          f"<- grading's main effect (false-positive removal)")
    print(f"PR-AUC delta    (grading - raw): {graded_ap - raw_ap:+.3f}  "
          f"(insensitive to already-bottom-ranked FPs)")


if __name__ == "__main__":
    main()
