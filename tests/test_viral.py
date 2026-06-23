from pathogeniq.viral import (
    _parse_genomad_summary, _parse_checkv_summary, _to_float, _to_int, ViralContig,
)

_GENOMAD = "\n".join([
    "\t".join(["seq_name", "length", "topology", "n_genes", "virus_score",
               "fdr", "n_hallmarks", "taxonomy"]),
    "\t".join(["k141_5", "42310", "DTR", "60", "0.99", "0.0", "8",
               "Viruses;Duplodnaviria;Caudoviricetes"]),
    "\t".join(["k141_9", "5120", "No terminal repeats", "7", "0.81", "0.01", "1",
               "Viruses;Monodnaviria"]),
])

_CHECKV = "\n".join([
    "\t".join(["contig_id", "contig_length", "checkv_quality", "completeness", "contamination"]),
    "\t".join(["k141_5", "42310", "High-quality", "98.5", "0.0"]),
    "\t".join(["k141_9", "5120", "Low-quality", "NA", "0.0"]),
])


def test_parse_genomad():
    r = _parse_genomad_summary(_GENOMAD)
    assert set(r) == {"k141_5", "k141_9"}
    assert r["k141_5"]["taxonomy"].endswith("Caudoviricetes")
    assert r["k141_5"]["virus_score"] == "0.99"


def test_parse_checkv():
    q = _parse_checkv_summary(_CHECKV)
    assert q["k141_5"]["completeness"] == "98.5"
    assert q["k141_9"]["checkv_quality"] == "Low-quality"


def test_numeric_coercion_handles_NA():
    assert _to_float("98.5") == 98.5
    assert _to_float("NA") is None
    assert _to_float(None) is None
    assert _to_int("60") == 60
    assert _to_int("NA") is None


def test_viral_contig_join_shape():
    # mimic run_viral_stage's join: genomad record + checkv quality -> ViralContig
    g = _parse_genomad_summary(_GENOMAD)["k141_5"]
    c = _parse_checkv_summary(_CHECKV)["k141_5"]
    v = ViralContig(
        contig_id="k141_5", length=_to_int(g["length"]) or 0,
        taxonomy=g["taxonomy"], topology=g["topology"],
        virus_score=_to_float(g["virus_score"]), n_hallmarks=_to_int(g["n_hallmarks"]),
        completeness=_to_float(c["completeness"]), checkv_quality=c["checkv_quality"],
    )
    assert v.completeness == 98.5 and v.n_hallmarks == 8 and v.virus_score == 0.99
