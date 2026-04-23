from __future__ import annotations

from typing import TYPE_CHECKING

from pathogeniq.config import SpecimenType

if TYPE_CHECKING:
    from pathogeniq.report import ReportEntry


# Organisms that are routine contaminants in specific specimen types.
# Matching is substring (case-insensitive) against ReportEntry.organism.
CONTAMINANT_PRIORS: dict[SpecimenType, list[str]] = {
    SpecimenType.BLOOD: [
        "Cutibacterium acnes",
        "Staphylococcus epidermidis",
        "Staphylococcus capitis",
        "Staphylococcus hominis",
        "Bacillus cereus",
        "Bacillus subtilis",
        "Micrococcus luteus",
        "Corynebacterium striatum",
        "Corynebacterium jeikeium",
    ],
    SpecimenType.CSF: [
        "Cutibacterium acnes",
        "Staphylococcus epidermidis",
        "Staphylococcus capitis",
        "Corynebacterium striatum",
    ],
    SpecimenType.BAL: [
        "Streptococcus salivarius",
        "Streptococcus mitis",
        "Streptococcus oralis",
        "Prevotella melaninogenica",
        "Veillonella parvula",
        "Rothia mucilaginosa",
        "Neisseria sicca",
        "Fusobacterium nucleatum",
    ],
    SpecimenType.TISSUE: [],
}


def flag_contaminants(entries: list[ReportEntry]) -> list[ReportEntry]:
    """Set contaminant_risk=True on entries matching known contaminants for their specimen type."""
    for entry in entries:
        known = CONTAMINANT_PRIORS.get(entry.specimen_type, [])
        entry.contaminant_risk = any(
            c.lower() in entry.organism.lower() for c in known
        )
    return entries
