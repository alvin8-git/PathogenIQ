"""Benchmark support for the NTC-aware grading method (T10).

Two reusable pieces, both classifier-agnostic and unit-testable without running
any external tool:

1. A Kraken2 adapter — parse a Kraken2 report and turn each species into the
   same GradingInput the pipeline grades, so the grading layer can score Kraken2
   output (the whole point of the wedge: grading is independent of the upstream
   classifier).
2. Scoring — average precision (PR-AUC) and precision-at-fixed-recall of a set of
   predicted-present taxa against a truth set, to compare Kraken2-raw vs
   Kraken2+grading vs an external baseline on a labeled community.

The actual benchmark RUN (Kraken2 on CAMI/Zymo reads, CZID, the held-out split)
is orchestrated by scripts/06_benchmark.py and needs the data + tools; this
module is the pure core it builds on.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path

from .config import SpecimenType
from .report import GradingInput


@dataclass(frozen=True)
class KrakenTaxon:
    taxid: str
    name: str
    reads: int   # clade-level reads (Kraken2 report column 2)
    rank: str


def parse_kraken2_report(text: str, *, rank: str = "S") -> list[KrakenTaxon]:
    """Parse a Kraken2 report (TSV) and return taxa at the given rank.

    Kraken2 report columns: percent, clade_reads, taxon_reads, rank_code,
    ncbi_taxid, name (name is depth-indented). Default rank "S" = species.
    """
    taxa: list[KrakenTaxon] = []
    for line in text.splitlines():
        cols = line.split("\t")
        if len(cols) < 6:
            continue
        clade_reads, rank_code, taxid, name = cols[1].strip(), cols[3].strip(), cols[4].strip(), cols[5].strip()
        if rank_code != rank:
            continue
        try:
            reads = int(clade_reads)
        except ValueError:
            continue
        taxa.append(KrakenTaxon(taxid=taxid, name=name, reads=reads, rank=rank_code))
    return taxa


def _wilson_halfwidth(k: int, n: int, z: float = 1.96) -> float:
    """Half-width of the Wilson score interval for the proportion k/n. Kraken2
    gives no bootstrap CI, so this stands in: it narrows with more reads, like
    the pipeline's bootstrap CI, so grade()'s ci_width threshold behaves sensibly."""
    if n <= 0:
        return 1.0
    p = k / n
    return (z * math.sqrt((p * (1 - p) + z * z / (4 * n)) / n)) / (1 + z * z / n)


def kraken_to_grading_inputs(
    taxa: list[KrakenTaxon],
    *,
    specimen: SpecimenType,
    tier: int = 3,
) -> dict[str, GradingInput]:
    """Map Kraken2 taxa to GradingInput records keyed by NCBI taxid. ci_width is a
    Wilson interval on the read proportion (Kraken2 has no bootstrap CI).

    Note: the NB background join keys on taxon_id; for the benchmark the background
    must be keyed by the same id space (taxid) as the Kraken output — a taxid<->GCF
    reconciliation handled by the benchmark script, not here.
    """
    total = sum(t.reads for t in taxa)
    out: dict[str, GradingInput] = {}
    for t in taxa:
        abundance = t.reads / total if total else 0.0
        ci_width = 2 * _wilson_halfwidth(t.reads, total)
        out[t.taxid] = GradingInput(
            read_count=t.reads,
            abundance=abundance,
            ci_width=ci_width,
            contaminant_risk=False,
            specimen_type=specimen,
            tier=tier,
        )
    return out


def parse_cami_profile(
    text: str,
    *,
    rank: str = "species",
    min_pct: float = 0.0,
    sample_id: str | None = None,
) -> set[str]:
    """Parse a CAMI gold-standard profile into a truth set of taxids at ``rank``.

    CAMI profile rows (after the ``@@TAXID...`` header) are tab-separated:
    TAXID, RANK, TAXPATH, TAXPATHSN, PERCENTAGE. A file may hold many samples,
    each a ``@SampleID:...`` block; pass ``sample_id`` (exact match) to score one
    sample, or None to pool all (the union). Keeps taxa at ``rank`` above
    ``min_pct``.
    """
    truth: set[str] = set()
    in_scope = sample_id is None
    for line in text.splitlines():
        s = line.strip()
        if not s:
            continue
        if s.startswith("@SampleID:"):
            in_scope = sample_id is None or s.split(":", 1)[1].strip() == sample_id
            continue
        if s.startswith("@") or s.startswith("#") or not in_scope:
            continue
        cols = s.split("\t")
        if len(cols) < 5:
            continue
        taxid, row_rank, pct = cols[0].strip(), cols[1].strip(), cols[4].strip()
        if row_rank != rank:
            continue
        try:
            if float(pct) > min_pct:
                truth.add(taxid)
        except ValueError:
            continue
    return truth


def load_truth(path: Path) -> set[str]:
    """Load a truth taxid set: a CAMI gold-standard .profile (auto-detected by its
    @ header) or a plain one-taxid-per-line file."""
    txt = Path(path).read_text()
    if txt.lstrip().startswith("@") or "@@TAXID" in txt:
        return parse_cami_profile(txt)
    return {ln.strip() for ln in txt.splitlines() if ln.strip()}


def precision_recall(predicted: set[str], truth: set[str]) -> tuple[float, float]:
    """Precision and recall of a predicted-present set against the truth set."""
    tp = len(predicted & truth)
    fp = len(predicted - truth)
    fn = len(truth - predicted)
    precision = tp / (tp + fp) if (tp + fp) else 1.0
    recall = tp / (tp + fn) if (tp + fn) else 1.0
    return precision, recall


def average_precision(scored: list[tuple[str, float]], truth: set[str]) -> float:
    """Average precision (the standard PR-AUC summary) of a ranked list.

    ``scored`` is (taxon_id, score) with higher score = more confident present.
    Perfect ranking (all truth above all non-truth) -> 1.0; a false positive
    ranked above a true positive lowers it.
    """
    total_pos = len(truth)
    if total_pos == 0:
        return 1.0
    ranked = sorted(scored, key=lambda x: x[1], reverse=True)
    tp = fp = 0
    prev_recall = 0.0
    ap = 0.0
    for taxon_id, _score in ranked:
        if taxon_id in truth:
            tp += 1
        else:
            fp += 1
        recall = tp / total_pos
        precision = tp / (tp + fp)
        ap += (recall - prev_recall) * precision
        prev_recall = recall
    return ap


def precision_at_recall(scored: list[tuple[str, float]], truth: set[str], target_recall: float) -> float:
    """Precision at the first threshold where recall reaches ``target_recall``
    (the design's 'precision at fixed sensitivity' metric)."""
    total_pos = len(truth)
    if total_pos == 0:
        return 1.0
    ranked = sorted(scored, key=lambda x: x[1], reverse=True)
    tp = fp = 0
    for taxon_id, _score in ranked:
        if taxon_id in truth:
            tp += 1
        else:
            fp += 1
        if tp / total_pos >= target_recall:
            return tp / (tp + fp)
    return 0.0   # target recall never reached
