from pathogeniq.config import SpecimenType
from pathogeniq.report import grade, EvidenceGrade
from pathogeniq.benchmark import (
    parse_kraken2_report,
    kraken_to_grading_inputs,
    precision_recall,
    average_precision,
    precision_at_recall,
)


# Kraken2 report: percent, clade_reads, taxon_reads, rank, taxid, name (indented)
KRAKEN = (
    "50.00\t5000\t5000\tU\t0\tunclassified\n"
    "30.00\t3000\t10\tD\t2\t  Bacteria\n"
    "20.00\t2000\t2000\tS\t562\t    Escherichia coli\n"
    "8.00\t800\t800\tS\t1280\t    Staphylococcus aureus\n"
    "0.50\t5\t5\tS\t1773\t    Mycobacterium tuberculosis\n"
)


def test_parse_kraken2_report_species_only():
    taxa = parse_kraken2_report(KRAKEN)
    ids = {t.taxid for t in taxa}
    assert ids == {"562", "1280", "1773"}          # only rank S
    ecoli = next(t for t in taxa if t.taxid == "562")
    assert ecoli.reads == 2000
    assert ecoli.name == "Escherichia coli"


def test_kraken_to_grading_inputs_abundance_and_ci():
    taxa = parse_kraken2_report(KRAKEN)
    gi = kraken_to_grading_inputs(taxa, specimen=SpecimenType.BLOOD, tier=1)
    total = 2000 + 800 + 5
    assert abs(gi["562"].abundance - 2000 / total) < 1e-9
    # Wilson CI on the proportion approximates the bootstrap (resampling all
    # reads), so the DOMINANT taxon has the wider absolute CI (variance peaks at
    # p=0.5), not the rare one. Document that property.
    assert 0.0 < gi["562"].ci_width < 1.0
    assert gi["562"].ci_width > gi["1773"].ci_width


def test_grading_scores_kraken_output_and_tier_cap():
    # the wedge: grade() consumes Kraken-derived GradingInput unchanged, and the
    # tier cap still applies — Grade A only with a batch-matched NTC (tier 1)
    taxa = parse_kraken2_report(KRAKEN)
    gi1 = kraken_to_grading_inputs(taxa, specimen=SpecimenType.BLOOD, tier=1)
    gi3 = kraken_to_grading_inputs(taxa, specimen=SpecimenType.BLOOD, tier=3)
    assert grade(gi1["562"]) == EvidenceGrade.A   # 2000 reads, tight CI, tier 1
    assert grade(gi3["562"]) == EvidenceGrade.B   # same taxon, no NTC -> capped at B


def test_precision_recall():
    p, r = precision_recall({"a", "b", "c"}, {"a", "b", "d"})
    assert abs(p - 2 / 3) < 1e-9   # a,b correct; c wrong
    assert abs(r - 2 / 3) < 1e-9   # a,b found; d missed


def test_average_precision_perfect_ranking():
    truth = {"A", "B"}
    scored = [("A", 10.0), ("B", 8.0), ("C", 5.0), ("D", 2.0)]
    assert abs(average_precision(scored, truth) - 1.0) < 1e-9


def test_average_precision_false_positive_hurts():
    truth = {"A", "B"}
    # a false positive (C) ranked between the two true positives
    scored = [("A", 10.0), ("C", 9.0), ("B", 8.0)]
    ap = average_precision(scored, truth)
    assert ap < 1.0
    assert abs(ap - (0.5 * 1.0 + 0.5 * (2 / 3))) < 1e-9   # 0.5 + 0.333 = 0.833


def test_precision_at_recall():
    truth = {"A", "B"}
    scored = [("A", 10.0), ("C", 9.0), ("B", 8.0)]
    # reach 100% recall only after including C (the FP) -> precision 2/3
    assert abs(precision_at_recall(scored, truth, 1.0) - 2 / 3) < 1e-9
    # 50% recall reached at A alone -> precision 1.0
    assert precision_at_recall(scored, truth, 0.5) == 1.0
