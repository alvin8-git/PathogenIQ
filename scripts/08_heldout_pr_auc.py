#!/usr/bin/env python
"""Multi-community held-out PR-AUC for the grading benchmark (T10).

Reads a manifest of communities (tab-separated: name, kreport, truth, specimen),
scores kraken-raw vs kraken+grading per community, and reports per-community plus
the mean PR-AUC and operating-point precision. The grading thresholds are fixed a
priori (not tuned on these data), so every community is an out-of-sample
evaluation and the mean is a held-out aggregate.

Optional calibrated split: --train N uses the first N communities to pick the
read-count floor that maximizes mean graded precision, then reports the held-out
(remaining) communities with that floor — a leakage-free train/test split.

Truth files: one taxid per line, OR a CAMI gold-standard .profile (auto-detected,
so CAMI communities drop in as soon as their reads are classified to a .kreport).

Manifest columns: name, kreport, truth, specimen, [rank]. The optional 5th
column is the Kraken2 rank code to score at (default S=species); use G for
communities whose truth is only genus-resolved (e.g. CAMI 97%-OTU HMP body
sites). Keep different ranks in separate manifests — the mean averages
PR/precision and mixing ranks is apples-to-oranges.

Manifest example:
    zymo_r1<TAB>r1.kreport<TAB>zymo_truth.txt<TAB>blood
    cami_mousegut<TAB>mg.kreport<TAB>cami2_mouse_gut_gs.profile<TAB>tissue
    hmp_skin<TAB>skin.kreport<TAB>skin_truth.txt<TAB>blood<TAB>G
"""
import argparse
import statistics
from pathlib import Path

from pathogeniq.benchmark import (
    average_precision,
    kraken_to_grading_inputs,
    load_truth,
    parse_kraken2_report,
    precision_recall,
)
from pathogeniq.config import SpecimenType
from pathogeniq.crossmap import find_crossmappers
from pathogeniq.report import EvidenceGrade, grade


def graded_kept(taxa, specimen, read_floor) -> set[str]:
    gi = kraken_to_grading_inputs(taxa, specimen=specimen, tier=3)
    cross = find_crossmappers([(t.taxid, t.name, t.reads) for t in taxa])
    kept = set()
    for t in taxa:
        if t.taxid in cross or t.reads < read_floor:
            continue
        if grade(gi[t.taxid]) != EvidenceGrade.X:
            kept.add(t.taxid)
    return kept


def _f1(p: float, r: float) -> float:
    return 2 * p * r / (p + r) if (p + r) else 0.0


def metrics(taxa, truth, specimen, read_floor):
    raw_scored = [(t.taxid, float(t.reads)) for t in taxa]
    kept = graded_kept(taxa, specimen, read_floor)
    g_scored = [(t.taxid, float(t.reads)) for t in taxa if t.taxid in kept]
    raw_p, _ = precision_recall({t.taxid for t in taxa}, truth)
    g_p, g_r = precision_recall(kept, truth)
    return (average_precision(raw_scored, truth), average_precision(g_scored, truth),
            raw_p, g_p, g_r)


def main() -> None:
    ap = argparse.ArgumentParser(description="Held-out PR-AUC across communities (T10).")
    ap.add_argument("--manifest", required=True, type=Path)
    ap.add_argument("--read-floor", type=int, default=0, help="extra global read floor on the graded config")
    ap.add_argument("--train", type=int, default=0, help="calibrate the floor on the first N communities, report the rest")
    args = ap.parse_args()

    comms = []
    for ln in args.manifest.read_text().splitlines():
        if not ln.strip() or ln.startswith("#"):
            continue
        cols = ln.split("\t")
        name, kr, tr, sp = cols[:4]
        rank = cols[4].strip() if len(cols) > 4 and cols[4].strip() else "S"  # kraken rank code; G for genus-resolved truth (HMP OTUs)
        taxa = parse_kraken2_report(Path(kr).read_text(), rank=rank)
        comms.append((name, taxa, load_truth(Path(tr)), SpecimenType(sp)))

    read_floor = args.read_floor
    test = comms
    if args.train > 0:
        train, test = comms[:args.train], comms[args.train:]
        best = (-1.0, read_floor)
        for f in (0, 2, 5, 10, 20, 50, 100, 200, 500):
            # calibrate on F1 (balanced) — precision alone picks a degenerate floor
            mean_f1 = statistics.mean(_f1(*metrics(t, tr, sp, f)[3:5]) for (_n, t, tr, sp) in train)
            if mean_f1 > best[0]:
                best = (mean_f1, f)
        read_floor = best[1]
        print(f"calibrated read floor = {read_floor} (train mean F1 {best[0]:.3f})\n")

    print(f"{'community':18s} {'raw PR':>7s} {'grd PR':>7s} {'raw P':>7s} {'grd P':>7s} {'grd R':>7s}")
    rows = []
    for name, taxa, truth, sp in test:
        m = metrics(taxa, truth, sp, read_floor)
        rows.append(m)
        print(f"{name:18s} {m[0]:7.3f} {m[1]:7.3f} {m[2]:7.3f} {m[3]:7.3f} {m[4]:7.3f}")
    if rows:
        mean = [statistics.mean(r[i] for r in rows) for i in range(5)]
        print(f"\n{'MEAN (held-out)':18s} {mean[0]:7.3f} {mean[1]:7.3f} {mean[2]:7.3f} {mean[3]:7.3f} {mean[4]:7.3f}")
        print(f"grading precision lift (mean): {mean[3] - mean[2]:+.3f}")


if __name__ == "__main__":
    main()
