import json

from pathogeniq.pathogenicity import (
    pathogen_taxa_from_name_map, _lineage_taxa, phylo_match, assess_pathogenicity,
    PATHOGEN_CANDIDATE, PATHOGEN_ADJACENT, ENVIRONMENTAL,
)


def test_pathogen_taxa_from_name_map(tmp_path):
    p = tmp_path / "name_map.json"
    p.write_text(json.dumps({"GCF_1": "Bacillus anthracis", "GCF_2": "Escherichia coli"}))
    taxa = pathogen_taxa_from_name_map(p)
    assert "bacillus" in taxa and "bacillus anthracis" in taxa
    assert "escherichia" in taxa and "escherichia coli" in taxa


def test_pathogen_taxa_missing_file(tmp_path):
    assert pathogen_taxa_from_name_map(tmp_path / "nope.json") == set()


def test_lineage_taxa_parses_gtdb():
    g, s = _lineage_taxa("d__Bacteria;p__Bacillota;g__Bacillus;s__Bacillus cereus")
    assert g == "bacillus" and s == "bacillus cereus"
    assert _lineage_taxa(None) == (None, None)
    assert _lineage_taxa("d__Bacteria;g__;s__") == (None, None)   # empty prefixes


def test_phylo_match_species_then_genus():
    taxa = {"bacillus", "bacillus anthracis"}
    # species-level match wins
    assert phylo_match("g__Bacillus;s__Bacillus anthracis", taxa) == "bacillus anthracis"
    # genus-only match when species not a known pathogen
    assert phylo_match("g__Bacillus;s__Bacillus subtilis", taxa) == "bacillus"
    # environmental genus -> no match
    assert phylo_match("g__Sphingomonas;s__Sphingomonas hankookensis", taxa) is None


def test_assess_markers_dominate():
    # any marker -> PATHOGEN_CANDIDATE regardless of taxonomy
    v, s = assess_pathogenicity("g__Sphingomonas", n_virulence=1, n_amr=0, matched=None)
    assert v == PATHOGEN_CANDIDATE and s == 2


def test_assess_phylo_only_is_adjacent():
    v, s = assess_pathogenicity("g__Bacillus", n_virulence=0, n_amr=0, matched="bacillus")
    assert v == PATHOGEN_ADJACENT and s == 3


def test_assess_environmental():
    v, s = assess_pathogenicity("g__Sphingomonas", 0, 0, None)
    assert v == ENVIRONMENTAL and s == 0


def test_assess_score_combines():
    v, s = assess_pathogenicity("g__Bacillus", n_virulence=2, n_amr=1, matched="bacillus")
    assert v == PATHOGEN_CANDIDATE and s == 3 * 2 + 3   # markers*2 + phylo
