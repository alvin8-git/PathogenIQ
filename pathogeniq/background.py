"""Negative-binomial NTC background subtraction for clinical mNGS grading.

Reagent/kitome contamination dominates low-biomass mNGS false positives. This
module models the per-taxon background rate from no-template control(s) (NTCs)
and tests each sample taxon against it, so a taxon indistinguishable from
background can be zeroed out before grading.

Pipeline:

    NTC count(s) ──► build_background ──► BackgroundModel (per-taxon RPM rates)
                                               │
    sample count, depth ───────────────────────┤
                                               ▼
                                       background_pvalue  (NB upper tail)
                                               │
                                   is_background(p >= alpha) ──► zero out

Why RPM, not raw counts: a deep sample and a shallow NTC are not comparable by
count. Rates (reads per million classified reads) make the test depth-invariant;
depth re-enters only when the rate is scaled back to an expected count at the
sample's own depth.

Single-NTC limitation: with one control, per-taxon dispersion is not estimable,
so a shared weakly-informative dispersion prior is used. Quantifying the
sensitivity of the false-positive rate to that prior is tracked as the Plan-4
follow-up in todo.md.
"""
from __future__ import annotations

import csv
import re
from dataclasses import dataclass
from pathlib import Path

from scipy.stats import nbinom

# Pseudocount (RPM) added to every background rate so taxa absent from the NTC
# still have a finite background floor (no div-by-zero, no infinite confidence).
_PSEUDOCOUNT_RPM = 0.5
# Shared NB dispersion (the scipy `n`/size parameter r). Smaller r = more
# overdispersion. Weakly-informative default for the single-NTC case where the
# per-taxon dispersion cannot be estimated.
_DEFAULT_DISPERSION = 2.0
# Controlled false-positive rate: a taxon with tail p >= ALPHA is treated as
# indistinguishable from background.
_DEFAULT_ALPHA = 0.01


@dataclass(frozen=True)
class BackgroundModel:
    rates: dict[str, float]   # taxon_id -> background rate in RPM (incl. pseudocount)
    n_controls: int           # number of usable (non-empty) NTC samples
    dispersion: float         # shared NB dispersion (r); per-taxon estimation is future work
    tier: int                 # 1 = batch-matched NTC, 2 = pooled, 3 = none


def _rpm(count: int, total: int) -> float:
    return (count / total) * 1e6 if total > 0 else 0.0


def build_background(
    ntc_counts: list[tuple[dict[str, int], int]],
    *,
    tier: int,
    dispersion: float = _DEFAULT_DISPERSION,
    pseudocount_rpm: float = _PSEUDOCOUNT_RPM,
) -> BackgroundModel | None:
    """Build a background model from one or more NTC samples.

    ``ntc_counts`` is a list of ``(per-taxon counts keyed by taxon_id, total
    classified reads)``. The per-taxon background rate is the mean RPM across
    usable controls plus the pseudocount.

    Returns ``None`` when no control has any classified reads (the mandatory
    zero-NTC-reads edge case) — the caller must then grade uncorrected and flag
    the run as Tier 3.
    """
    usable = [(counts, total) for counts, total in ntc_counts if total > 0]
    if not usable:
        return None
    taxa: set[str] = set().union(*(set(counts) for counts, _ in usable))
    rates: dict[str, float] = {}
    for tx in taxa:
        per_ntc = [_rpm(counts.get(tx, 0), total) for counts, total in usable]
        rates[tx] = sum(per_ntc) / len(per_ntc) + pseudocount_rpm
    return BackgroundModel(
        rates=rates,
        n_controls=len(usable),
        dispersion=dispersion,
        tier=tier,
    )


def load_background_table(
    path: Path,
    *,
    dispersion: float = _DEFAULT_DISPERSION,
) -> BackgroundModel:
    """Load a precomputed pooled background table (the Tier-2 ``--background`` path).

    Format: an optional ``# tier=N`` comment line, then a ``taxon_id<TAB>rpm``
    header and one row per taxon. Values are final per-taxon background rates in
    RPM (any pseudocount is already baked in by whoever built the table). Tier
    defaults to 2 (pooled) when the comment is absent.
    """
    tier = 2
    data_lines: list[str] = []
    for line in Path(path).read_text().splitlines():
        stripped = line.strip()
        if not stripped:
            continue
        if stripped.startswith("#"):
            m = re.search(r"tier\s*=\s*(\d+)", stripped)
            if m:
                tier = int(m.group(1))
            continue
        data_lines.append(line)
    rates: dict[str, float] = {}
    for row in csv.DictReader(data_lines, delimiter="\t"):
        taxon_id = (row.get("taxon_id") or "").strip()
        if taxon_id:
            rates[taxon_id] = float(row["rpm"])
    return BackgroundModel(rates=rates, n_controls=0, dispersion=dispersion, tier=tier)


def default_background_path() -> Path:
    """Path to the packaged pooled background table (Tier-2 default)."""
    return Path(__file__).resolve().parent / "data" / "background_default.tsv"


def load_default_background() -> BackgroundModel | None:
    """Load the packaged pooled background (Tier-2 default path), or None when it
    is missing or still an empty placeholder — in which case the caller falls
    back to Tier 3 rather than capping grades with no background benefit."""
    path = default_background_path()
    if not path.exists():
        return None
    model = load_background_table(path)
    return model if model.rates else None


def write_background_table(path: Path, rates: dict[str, float], *, tier: int = 2) -> None:
    """Write per-taxon background rates as a ``# tier=N`` table that
    load_background_table can read back. Used by scripts/build_background_default.py."""
    lines = [f"# tier={tier}", "taxon_id\trpm"]
    for taxon_id, rpm in sorted(rates.items()):
        lines.append(f"{taxon_id}\t{rpm:.6g}")
    Path(path).write_text("\n".join(lines) + "\n")


def background_pvalue(
    taxon_id: str,
    sample_count: int,
    sample_total: int,
    model: BackgroundModel,
) -> float:
    """Upper-tail p-value ``P(count >= sample_count | NB background at sample depth)``.

    Low p  -> the observed count exceeds what background explains -> real signal.
    High p -> indistinguishable from background.

    A taxon absent from the NTC uses the pseudocount floor as its background rate,
    so genuine signal still passes.
    """
    if sample_count <= 0:
        return 1.0
    mu_rpm = model.rates.get(taxon_id, _PSEUDOCOUNT_RPM)
    lam = mu_rpm * sample_total / 1e6   # expected background count at the sample's depth
    if lam <= 0:
        return 0.0   # no background expected but reads observed -> definitively real
    r = model.dispersion
    p = r / (r + lam)
    return float(nbinom.sf(sample_count - 1, r, p))   # sf(k-1) = P(X >= k)


def is_background(
    taxon_id: str,
    sample_count: int,
    sample_total: int,
    model: BackgroundModel,
    alpha: float = _DEFAULT_ALPHA,
) -> bool:
    """True if the taxon is indistinguishable from background (tail p >= alpha)
    and should be zeroed out before grading."""
    return background_pvalue(taxon_id, sample_count, sample_total, model) >= alpha
