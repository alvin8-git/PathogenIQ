"""Spike-in absolute quantification (Plan 6 #2).

Relative abundance (RPM / %) loses the absolute scale: "organism X is 2% of the
reads" says nothing about how many copies of X are actually present — two samples
with the same percentage can differ 1000-fold in true load. A known internal
spike-in (a fixed, accurately-quantified amount of an organism absent from the
sample, added before extraction) restores the scale. With ``spike_reads`` reads
from a spike of known input ``spike_copies``, every target's absolute load is

    copies(target) = reads(target) / reads(spike) * copies(spike)

and, given the volume the sample represents, ``copies(target) / sample_volume``
gives the clinically meaningful per-mL / per-m^3 load.

Requirements: the spike organism must be in the reference DB (so it earns a read
count) and is removed from the findings — it is a process control, not a
detection. See docs / todo wet-lab section for dose + protocol (add before lysis,
verify input by ddPCR, target ~0.5-5% of reads).
"""
from __future__ import annotations

from dataclasses import dataclass

from .report import ReportEntry


@dataclass
class SpikeInfo:
    taxon_id: str
    copies_added: float
    spike_reads: int
    sample_volume: float | None = None   # units the caller chose (mL, m^3, ...)
    found: bool = False                   # spike detected with >0 reads


def absolute_copies(target_reads: int, spike_reads: int, spike_copies: float) -> float | None:
    """copies(target) = reads(target)/reads(spike) * copies(spike). None when the
    spike has no reads — without a denominator the absolute scale is undefined."""
    if spike_reads <= 0:
        return None
    return target_reads / spike_reads * spike_copies


def quantify_entries(
    entries: list[ReportEntry],
    *,
    spike_taxon_id: str,
    spike_copies: float,
    sample_volume: float | None = None,
) -> tuple[list[ReportEntry], SpikeInfo]:
    """Anchor absolute load to the spike entry, drop the spike from the findings
    (it is a control), and set ``absolute_copies`` on each remaining entry.

    If the spike is absent or has zero reads, the spike row is still removed but
    ``absolute_copies`` stays None (``SpikeInfo.found`` is False) — the run should
    be flagged, since a missing spike means quantification failed, not that the
    sample was empty."""
    spike = next((e for e in entries if e.taxon_id == spike_taxon_id), None)
    spike_reads = spike.read_count if spike is not None else 0
    info = SpikeInfo(
        taxon_id=spike_taxon_id,
        copies_added=spike_copies,
        spike_reads=spike_reads,
        sample_volume=sample_volume,
        found=spike is not None and spike_reads > 0,
    )
    kept = [e for e in entries if e.taxon_id != spike_taxon_id]
    if info.found:
        for e in kept:
            e.absolute_copies = absolute_copies(e.read_count, spike_reads, spike_copies)
    return kept, info
