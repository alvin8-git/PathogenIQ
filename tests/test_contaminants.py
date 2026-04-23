from pathogeniq.contaminants import flag_contaminants, CONTAMINANT_PRIORS
from pathogeniq.config import SpecimenType
from pathogeniq.report import ReportEntry, EvidenceGrade


def _entry(organism: str, specimen: SpecimenType, read_count: int = 5,
           abundance: float = 0.05, ci_lower: float = 0.01, ci_upper: float = 0.10) -> ReportEntry:
    return ReportEntry(
        organism=organism,
        abundance=abundance,
        ci_lower=ci_lower,
        ci_upper=ci_upper,
        read_count=read_count,
        specimen_type=specimen,
        contaminant_risk=False,
    )


def test_cutibacterium_flagged_in_blood():
    entries = [_entry("Cutibacterium acnes", SpecimenType.BLOOD)]
    result = flag_contaminants(entries)
    assert result[0].contaminant_risk is True


def test_cutibacterium_not_flagged_in_tissue():
    entries = [_entry("Cutibacterium acnes", SpecimenType.TISSUE)]
    result = flag_contaminants(entries)
    assert result[0].contaminant_risk is False


def test_streptococcus_salivarius_flagged_in_bal():
    entries = [_entry("Streptococcus salivarius", SpecimenType.BAL)]
    result = flag_contaminants(entries)
    assert result[0].contaminant_risk is True


def test_staphylococcus_aureus_not_flagged_blood():
    # S. aureus is a true pathogen even in blood — must NOT be suppressed
    entries = [_entry("Staphylococcus aureus", SpecimenType.BLOOD)]
    result = flag_contaminants(entries)
    assert result[0].contaminant_risk is False


def test_staph_epidermidis_flagged_in_csf():
    entries = [_entry("Staphylococcus epidermidis", SpecimenType.CSF)]
    result = flag_contaminants(entries)
    assert result[0].contaminant_risk is True


def test_contaminant_registry_has_all_specimen_types():
    for stype in SpecimenType:
        assert stype in CONTAMINANT_PRIORS


def test_grade_b_contaminant_downgrades_to_c():
    # Grade B entry (enough reads, wide CI) that is a contaminant → reported as C
    entry = _entry("Cutibacterium acnes", SpecimenType.BLOOD,
                   read_count=10, ci_lower=0.00, ci_upper=0.20)
    flagged = flag_contaminants([entry])
    # Grade should fall through to C because contaminant_risk=True makes grade() return C
    assert flagged[0].grade == EvidenceGrade.C
