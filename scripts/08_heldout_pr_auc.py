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

Manifest example:
    zymo_r1<TAB>r1.kreport<TAB>zymo_truth.txt<TAB>blood
    cami_mousegut<TAB>mg.kreport<TAB>cami2_mouse_gut_gs.profile<TAB>tissue
"""
import argparse
import statistics
from pathlib import Path

from pathogeniq.benchmark import (
    average_precision,
    kraken_to_grading_inputs,
    parse_cami_profile,
    parse_kraken2_report,
    precision_recall,
)
from pathogeniq.config import SpecimenType
from pathogeniq.crossmap import find_crossmappers
from pathogeniq.report import EvidenceGrade, grade


def load_truth(path: Path) -> set[str]:
    txt = Path(path).read_text()
    if txt.lstrip().startswith("@") or "@@TAXID" in txt:   # CAMI profile
        return parse_cami_profile(txt)
    return {ln.strip() for ln in txt.splitlines() if ln.strip()}


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


def metrics(taxa, truth, specimen, read_floor):
    raw_scored = [(t.taxid, float(t.reads)) for t in taxa]
    kept = graded_kept(taxa, specimen, read_floor)
    g_scored = [(t.taxid, float(t.reads)) for t in taxa if t.taxid in kept]
    raw_p, _ = precision_recall({t.taxid for t in taxa}, truth)
    g_p, _ = precision_recall(kept, truth)
    return average_precision(raw_scored, truth), average_precision(g_scored, truth), raw_p, g_p


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
        name, kr, tr, sp = ln.split("\t")[:4]
        taxa = parse_kraken2_report(Path(kr).read_text())
        comms.append((name, taxa, load_truth(Path(tr)), SpecimenType(sp)))

    read_floor = args.read_floor
    test = comms
    if args.train > 0:
        train, test = comms[:args.train], comms[args.train:]
        best = (-1.0, read_floor)
        for f in (0, 2, 5, 10, 20, 50, 100, 200, 500):
            mean_p = statistics.mean(metrics(t, tr, sp, f)[3] for (_n, t, tr, sp) in train)
            if mean_p > best[0]:
                best = (mean_p, f)
        read_floor = best[1]
        print(f"calibrated read floor = {read_floor} (train mean graded precision {best[0]:.3f})\n")

    print(f"{'community':18s} {'raw PR':>7s} {'grd PR':>7s} {'raw P':>7s} {'grd P':>7s}")
    rows = []
    for name, taxa, truth, sp in test:
        r_ap, g_ap, r_p, g_p = metrics(taxa, truth, sp, read_floor)
        rows.append((r_ap, g_ap, r_p, g_p))
        print(f"{name:18s} {r_ap:7.3f} {g_ap:7.3f} {r_p:7.3f} {g_p:7.3f}")
    if rows:
        mean = [statistics.mean(r[i] for r in rows) for i in range(4)]
        print(f"\n{'MEAN (held-out)':18s} {mean[0]:7.3f} {mean[1]:7.3f} {mean[2]:7.3f} {mean[3]:7.3f}")
        print(f"grading precision lift (mean): {mean[3] - mean[2]:+.3f}")


if __name__ == "__main__":
    main()
