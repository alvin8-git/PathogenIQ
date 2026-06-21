from pathogeniq.crossmap import find_crossmappers, deduplicate_closely_related
from pathogeniq.config import SpecimenType
from pathogeniq.report import ReportEntry, EvidenceGrade


def _entry(name, reads, taxid=""):
    return ReportEntry(
        organism=name, abundance=0.1, ci_lower=0.05, ci_upper=0.15,
        read_count=reads, specimen_type=SpecimenType.BLOOD, taxon_id=taxid, tier=1,
    )


def test_find_crossmappers_flags_dominated_shigella():
    # E. coli dominates Shigella >10x -> Shigella flagged as cross-mapping
    items = [("562", "Escherichia coli", 100000),
             ("623", "Shigella flexneri", 800),
             ("624", "Shigella sonnei", 700)]
    flagged = find_crossmappers(items)
    assert flagged == {"623": "Escherichia coli", "624": "Escherichia coli"}


def test_find_crossmappers_no_flag_when_ratio_not_met():
    # a genuine co-detection (within 10x) is not flagged
    items = [("562", "Escherichia coli", 1000), ("623", "Shigella flexneri", 500)]
    assert find_crossmappers(items) == {}


def test_find_crossmappers_single_member_untouched():
    items = [("562", "Escherichia coli", 5000)]
    assert find_crossmappers(items) == {}


def test_deduplicate_sets_crossmap_of_and_grades_x():
    entries = [_entry("Escherichia coli", 100000, "GCF_ec"),
               _entry("Shigella flexneri", 500, "GCF_sf")]
    deduplicate_closely_related(entries)
    ec, sf = entries
    assert ec.crossmap_of == ""            # the dominant organism is untouched
    assert sf.crossmap_of == "Escherichia coli"
    assert ec.grade != EvidenceGrade.X     # E. coli still graded
    assert sf.grade == EvidenceGrade.X     # Shigella dropped as an artifact
