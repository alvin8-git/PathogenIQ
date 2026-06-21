"""Breadth-of-coverage validation for flagged taxa.

Read count alone can't separate a real organism from an artifact: 100 reads piled
at one locus (PCR duplicate / low-complexity / cross-mapping) look identical, by
count, to 100 reads spread across the genome (real signal). Breadth of coverage
distinguishes them, and the Lander-Waterman expected breadth makes it
depth-normalized so it works at any sequencing depth:

    observed_breadth = covered_bases / genome_length
    depth            = aligned_bases / genome_length
    expected_breadth = 1 - exp(-depth)          # uniform reads (Lander-Waterman)
    breadth_ratio    = observed_breadth / expected_breadth

    breadth_ratio ~ 1   -> reads spread as expected for the depth  -> real
    breadth_ratio << 1  -> reads clumped far below expectation     -> artifact

This is the alignment-path complement to the k-mer read count: it needs read
*positions*, which the targeted-alignment stage already produces in its minimap2
PAF (and which align.py currently discards). It is the natural home for the air-
surveillance design's "coverage breadth over read depth" rule.
"""
from __future__ import annotations

import math
from dataclasses import dataclass


@dataclass(frozen=True)
class CoverageStats:
    genome_length: int
    n_alignments: int
    aligned_bases: int
    covered_bases: int
    breadth: float          # covered / genome_length
    depth: float            # aligned_bases / genome_length
    expected_breadth: float  # 1 - exp(-depth), uniform-read expectation
    breadth_ratio: float    # observed / expected; ~1 real, <<1 clumped artifact


def _merged_covered(intervals: list[tuple[int, int]]) -> int:
    """Total bases covered by a set of [start, end) intervals (union)."""
    if not intervals:
        return 0
    covered = 0
    cur_s, cur_e = sorted(intervals)[0]
    for s, e in sorted(intervals)[1:]:
        if s <= cur_e:
            cur_e = max(cur_e, e)
        else:
            covered += cur_e - cur_s
            cur_s, cur_e = s, e
    return covered + (cur_e - cur_s)


def _stats(genome_length: int, n: int, aligned: int, covered: int) -> CoverageStats:
    breadth = covered / genome_length if genome_length else 0.0
    depth = aligned / genome_length if genome_length else 0.0
    expected = 1 - math.exp(-depth) if depth > 0 else 0.0
    ratio = breadth / expected if expected > 0 else 0.0
    return CoverageStats(genome_length, n, aligned, covered, breadth, depth, expected, ratio)


def coverage_from_paf(paf_text: str) -> CoverageStats:
    """Aggregate breadth-of-coverage across all reference contigs in one PAF
    (one organism, possibly multi-contig). Only mapped (mapq > 0) records count."""
    per_target: dict[str, tuple[list[tuple[int, int]], int]] = {}
    for line in paf_text.splitlines():
        cols = line.split("\t")
        if len(cols) < 12:
            continue
        try:
            tname, tlen = cols[5], int(cols[6])
            tstart, tend, mapq = int(cols[7]), int(cols[8]), int(cols[11])
        except ValueError:
            continue
        if mapq <= 0:
            continue
        ivs, _ = per_target.setdefault(tname, ([], tlen))
        ivs.append((tstart, tend))
    genome_length = sum(tlen for _ivs, tlen in per_target.values())
    covered = sum(_merged_covered(ivs) for ivs, _tlen in per_target.values())
    aligned = sum(e - s for ivs, _tlen in per_target.values() for s, e in ivs)
    n = sum(len(ivs) for ivs, _tlen in per_target.values())
    return _stats(genome_length, n, aligned, covered)
