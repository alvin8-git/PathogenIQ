from pathogeniq.novelty import parse_kraken_report

# Synthetic Kraken2 --report: 600 unclassified, 400 classified (300 E. coli, 100 S. aureus)
_REPORT = "\n".join([
    "\t".join(["60.00", "600", "600", "U", "0", "unclassified"]),
    "\t".join(["40.00", "400", "0", "R", "1", "root"]),
    "\t".join(["40.00", "400", "0", "D", "2", "  Bacteria"]),
    "\t".join(["30.00", "300", "300", "S", "562", "    Escherichia coli"]),
    "\t".join(["10.00", "100", "100", "S", "1280", "    Staphylococcus aureus"]),
])


def test_parse_counts_and_fraction():
    r = parse_kraken_report(_REPORT)
    assert r.total_reads == 1000
    assert r.classified_reads == 400
    assert r.unclassified_reads == 600
    assert abs(r.unclassified_fraction - 0.6) < 1e-9
    assert r.n_species == 2
    assert r.top_taxa[0] == ("Escherichia coli", 300)


def test_flag_threshold():
    # 60% unclassified: flagged at default 0.5, not at 0.7
    assert parse_kraken_report(_REPORT).flagged is True
    assert parse_kraken_report(_REPORT, flag_threshold=0.7).flagged is False


def test_no_root_row_falls_back_to_sum():
    # report without an explicit root row still totals classified from species
    rep = "\n".join([
        "\t".join(["50.00", "100", "100", "U", "0", "unclassified"]),
        "\t".join(["50.00", "100", "100", "S", "562", "Escherichia coli"]),
    ])
    r = parse_kraken_report(rep)
    assert r.classified_reads == 100
    assert r.total_reads == 200
    assert abs(r.unclassified_fraction - 0.5) < 1e-9


def test_resolve_kraken_db_descends_into_subdir(tmp_path):
    from pathogeniq.novelty import _resolve_kraken_db
    # DB directly in the dir
    (tmp_path / "taxo.k2d").write_bytes(b"x")
    assert _resolve_kraken_db(tmp_path) == tmp_path
    # DB one level down (the tarball-extracted layout that silently skipped novelty)
    nested = tmp_path / "nested"
    sub = nested / "k2_standard_08gb"
    sub.mkdir(parents=True)
    (sub / "taxo.k2d").write_bytes(b"x")
    assert _resolve_kraken_db(nested) == sub
    # no DB anywhere
    assert _resolve_kraken_db(tmp_path / "missing") is None
    empty = tmp_path / "empty"
    empty.mkdir()
    assert _resolve_kraken_db(empty) is None


def test_empty_report_no_crash():
    r = parse_kraken_report("")
    assert r.total_reads == 0
    assert r.unclassified_fraction == 0.0
    assert r.flagged is False
