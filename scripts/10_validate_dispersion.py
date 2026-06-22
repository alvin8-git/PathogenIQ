#!/usr/bin/env python
"""Validate the single-NTC NB dispersion prior (Plan-4 follow-up).

With one NTC the per-taxon NB dispersion `r` is not estimable, so background.py
uses a shared prior (_DEFAULT_DISPERSION). This script quantifies how the
background test's false-positive rate depends on `r` and picks a defensible
default, using the cached per-run kitome counts from script 05.

Two curves, same data, over a sweep of `r`:

1. FALSE-POSITIVE RATE (the calibration target). Leave-one-out across the
   spike-free control runs: build the background from the others, then test
   every taxon in the held-out blank. The held-out run is itself a kitome blank,
   so any taxon called "above background" (p < alpha) is a FALSE POSITIVE — there
   is no real signal in a blank. The spike taxon is excluded (it IS real). A
   well-chosen `r` keeps empirical FPR <= alpha; too-large `r` (too-narrow null)
   fails to absorb genuine run-to-run kitome variation and inflates FPR.

2. LIMIT OF DETECTION (the sensitivity cost). Minimum sample read count for a
   taxon ABSENT from the NTC (background = pseudocount floor) to clear the test
   at a reference depth. Smaller `r` (wider null) raises the LoD — the price of
   FPR control. Reported so the `r` choice is an explicit FPR/sensitivity trade.

Recommendation = the largest `r` whose LOO FPR <= alpha (lowest LoD subject to
FPR control).

    python scripts/10_validate_dispersion.py --workdir /data/alvin/tmp/kitome_select
"""
import argparse
import importlib.util
import json
from pathlib import Path

from pathogeniq.background import build_background, is_background

# reuse the spike selection from script 05 (module name starts with a digit)
_spec = importlib.util.spec_from_file_location(
    "kitome05", Path(__file__).resolve().parent / "05_select_kitome_controls.py")
_k05 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_k05)
spike_fraction, _DEFAULT_SPIKE = _k05.spike_fraction, _k05._DEFAULT_SPIKE

_DISPERSIONS = [0.5, 1.0, 2.0, 5.0, 10.0, 50.0, 1000.0]   # 1000 ~ Poisson limit


def load_runs(workdir: Path) -> list[tuple[str, dict[str, int], int]]:
    runs = []
    for cache in sorted(workdir.glob("*/counts.json")):
        d = json.loads(cache.read_text())
        runs.append((cache.parent.name, {k: int(v) for k, v in d["counts"].items()}, int(d["total"])))
    return runs


def loo_fpr(spike_free, dispersion, alpha, spike_prefixes):
    """Leave-one-out false-positive rate over the blank runs at one dispersion."""
    fp = tested = 0
    for i, (_name, counts, total) in enumerate(spike_free):
        others = [(c, t) for j, (_n, c, t) in enumerate(spike_free) if j != i]
        model = build_background(others, tier=2, dispersion=dispersion)
        if model is None:
            continue
        for tx, c in counts.items():
            if c <= 0 or any(tx.startswith(p) for p in spike_prefixes):
                continue   # spike taxon is real signal, not a null -> exclude
            tested += 1
            if not is_background(tx, c, total, model, alpha=alpha):
                fp += 1   # called "above background" in a blank => false positive
    return fp / tested if tested else 0.0, tested


def limit_of_detection(spike_free, dispersion, alpha, depth=1_000_000):
    """Min read count for an NTC-absent taxon (mu = pseudocount) to clear the test."""
    model = build_background([(c, t) for _n, c, t in spike_free], tier=2, dispersion=dispersion)
    absent = "__not_in_any_ntc__"
    assert model.rates.get(absent) is None
    for k in range(1, 100_000):
        if not is_background(absent, k, depth, model, alpha=alpha):
            return k
    return -1


def main() -> None:
    ap = argparse.ArgumentParser(description="Validate the NB dispersion prior (Plan-4).")
    ap.add_argument("--workdir", type=Path, default=Path("/data/alvin/tmp/kitome_select"))
    ap.add_argument("--max-spike-frac", type=float, default=0.10)
    ap.add_argument("--alpha", type=float, default=0.01)
    ap.add_argument("--depth", type=int, default=1_000_000, help="reference sample depth for LoD")
    args = ap.parse_args()

    spikes = [_DEFAULT_SPIKE]
    runs = load_runs(args.workdir)
    spike_free = [(n, c, t) for n, c, t in runs
                  if t > 0 and spike_fraction(c, t, spikes) <= args.max_spike_frac]
    print(f"{len(runs)} runs, {len(spike_free)} spike-free blanks "
          f"(spike <= {args.max_spike_frac:.0%}); alpha={args.alpha}, LoD depth={args.depth:,}\n")

    print(f"{'dispersion r':>12s} {'LOO FPR':>9s} {'FPR<=a':>7s} {'LoD reads':>10s}")
    ok_rs = []
    for r in _DISPERSIONS:
        fpr, n = loo_fpr(spike_free, r, args.alpha, spikes)
        lod = limit_of_detection(spike_free, r, args.alpha, args.depth)
        ctrl = fpr <= args.alpha
        ok_rs.append((r, fpr, lod, ctrl))
        print(f"{r:12.1f} {fpr:9.4f} {'yes' if ctrl else 'NO':>7s} {lod:10d}")

    controlled = [(r, lod) for r, _f, lod, c in ok_rs if c]
    print()
    if controlled:
        best_r, best_lod = max(controlled, key=lambda x: x[0])   # largest r that controls FPR
        print(f"RECOMMENDED r = {best_r} (largest dispersion with FPR <= {args.alpha}; "
              f"LoD {best_lod} reads at depth {args.depth:,})")
        print(f"current _DEFAULT_DISPERSION = 2.0 -> "
              f"{'controls FPR' if any(r == 2.0 and c for r, _f, _l, c in ok_rs) else 'does NOT control FPR'}")
    else:
        # diagnose WHY: a taxon absent from the training pool defaults to the
        # pseudocount floor and is called real at any r — so run-unique taxa set
        # an irreducible FPR floor that dispersion cannot touch.
        from collections import Counter
        seen = Counter(tx for _n, c, _t in spike_free for tx, k in c.items() if k > 0)
        uniq = sum(1 for v in seen.values() if v == 1)
        print(f"no swept dispersion controls FPR at alpha={args.alpha}. The limit is NTC "
              f"COVERAGE, not the prior: {uniq}/{len(seen)} kitome taxa occur in only one "
              f"blank, so they hit the pseudocount floor in leave-one-out and read as real "
              f"at every r. Fix = more/broader NTCs or a batch-matched NTC (Tier 1); the "
              f"shared dispersion prior is not the lever. Keep Tier-2 capped at Grade B.")


def _selfcheck() -> None:
    # two toy blanks; a taxon spuriously huge in one fold should register as FP at high r
    blanks = [("a", {"t1": 5, "t2": 1}, 1000), ("b", {"t1": 4, "t2": 800}, 1000)]
    fpr_wide, _ = loo_fpr(blanks, 0.5, 0.01, [_DEFAULT_SPIKE])
    fpr_narrow, _ = loo_fpr(blanks, 1000.0, 0.01, [_DEFAULT_SPIKE])
    assert 0.0 <= fpr_wide <= fpr_narrow <= 1.0, (fpr_wide, fpr_narrow)  # narrower null >= FPR
    assert limit_of_detection(blanks, 0.5, 0.01) >= limit_of_detection(blanks, 1000.0, 0.01)
    print("selfcheck OK")


if __name__ == "__main__":
    import sys
    if "--selfcheck" in sys.argv:
        _selfcheck()
    else:
        main()
