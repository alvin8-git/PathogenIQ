"""Collapse cross-mapping false positives from closely related references.

Some references are >95% genome-identical — Shigella *is* Escherichia coli
taxonomically — so reads from the dominant organism cross-map and produce phantom
calls for its relatives (the Zymo/benchmark E. coli -> Shigella artifact). When
one member of a known cross-mapping group vastly outnumbers another in the same
sample, the minor call is almost certainly an artifact: it's flagged with the
dominant organism's name, which grades it X (removed), not just demoted.

Scope is a hardcoded registry of known problematic groups, per the design — not
runtime ANI. Matching is case-insensitive substring against the organism name,
like the contaminant blocklist.
"""
from __future__ import annotations

_CROSSMAP_GROUPS: list[list[str]] = [
    # Shigella spp. are genomically E. coli; reads cross-map across all of these.
    ["Escherichia coli", "Shigella flexneri", "Shigella sonnei",
     "Shigella dysenteriae", "Shigella boydii"],
]

_DEFAULT_RATIO = 10.0


def find_crossmappers(
    items: list[tuple[str, str, int]],
    *,
    ratio: float = _DEFAULT_RATIO,
) -> dict[str, str]:
    """Identify cross-mapping artifacts.

    ``items`` is ``(id, organism_name, read_count)``. For each known cross-mapping
    group present in the sample, the member with the most reads is the likely true
    source; any other member it outnumbers by >= ``ratio`` is flagged. Returns
    ``{minor_id: major_organism_name}``. Pure — usable by both the report builder
    and the Kraken2 benchmark.
    """
    flagged: dict[str, str] = {}
    for group in _CROSSMAP_GROUPS:
        members = [
            (i, name, reads) for (i, name, reads) in items
            if any(g.lower() in name.lower() for g in group)
        ]
        if len(members) < 2:
            continue
        major = max(members, key=lambda m: m[2])
        for i, _name, reads in members:
            if i == major[0]:
                continue
            if reads > 0 and major[2] >= ratio * reads:
                flagged[i] = major[1]
    return flagged


def deduplicate_closely_related(entries, *, ratio: float = _DEFAULT_RATIO):
    """Set ``crossmap_of`` on report entries that are cross-mapping artifacts of a
    dominant relative in the same sample. Mutates and returns the entries."""
    items = [(str(i), e.organism, e.read_count) for i, e in enumerate(entries)]
    flagged = find_crossmappers(items, ratio=ratio)
    for i, e in enumerate(entries):
        major = flagged.get(str(i))
        if major:
            e.crossmap_of = major
    return entries
