from pathlib import Path
from unittest.mock import patch


from types import SimpleNamespace

from pathogeniq.amr import (
    _parse_abricate_tsv,
    _run_abricate,
    map_contigs_to_organisms,
    run_amr_screen,
    run_virulence_screen,
)
from pathogeniq.config import PipelineConfig, ReadType, SpecimenType


def _cfg(tmp_path: Path) -> PipelineConfig:
    return PipelineConfig(
        input_fastq=tmp_path / "sample.fq.gz",
        read_type=ReadType.SHORT,
        specimen_type=SpecimenType.BLOOD,
        output_dir=tmp_path,
        db_tier1=tmp_path / "db.sbt.zip",
        host_reference=tmp_path / "ref.fa",
    )


def test_run_abricate_skips_without_contigs(tmp_path):
    # no contigs (assembly produced nothing) -> screen skips, even if abricate exists
    with patch("pathogeniq.amr.shutil.which", return_value="abricate"):
        assert _run_abricate(_cfg(tmp_path), None, "card", 90.0, 80.0) is None


def test_run_abricate_is_non_blocking_on_error(tmp_path):
    # abricate present, contigs present, but the subprocess fails -> None (run survives)
    contigs = tmp_path / "contigs.fa"
    contigs.write_text(">k1\nACGT\n")
    with patch("pathogeniq.amr.shutil.which", return_value="abricate"), \
         patch("pathogeniq.amr.subprocess.run", side_effect=OSError("boom")):
        assert _run_abricate(_cfg(tmp_path), contigs, "card", 90.0, 80.0) is None


def test_parse_abricate_tsv_basic():
    tsv = (
        "#FILE\tSEQUENCE\tSTART\tEND\tGENE\tCOVERAGE\tGAPS\t%COVERAGE\t%IDENTITY\t"
        "DATABASE\tACCESSION\tPRODUCT\tRESISTANCE\n"
        "amr_reads.fa\tStaphylococcus_aureus_read1\t1\t100\tmecA\t1-100/100\t0\t100\t99.5\t"
        "card\tNC_000001\tmethicillin resistance protein\tbeta-lactam\n"
    )
    hits = _parse_abricate_tsv(tsv, organism_names=["Staphylococcus aureus"])
    assert len(hits) == 1
    assert hits[0].gene == "mecA"
    assert hits[0].drug_class == "beta-lactam"
    assert hits[0].identity_pct == 99.5
    assert hits[0].coverage_pct == 100.0
    assert hits[0].organism_match == "Staphylococcus aureus"
    assert hits[0].database == "card"


def test_parse_abricate_tsv_organism_match_by_substring():
    tsv = (
        "#FILE\tSEQUENCE\tSTART\tEND\tGENE\tCOVERAGE\tGAPS\t%COVERAGE\t%IDENTITY\t"
        "DATABASE\tACCESSION\tPRODUCT\tRESISTANCE\n"
        "amr_reads.fa\tEscherichia_coli_contig_1\t1\t50\tblaCTX-M\t1-50/50\t0\t100\t98.0\t"
        "card\tNC_000002\tExtended-spectrum beta-lactamase\tbeta-lactam\n"
    )
    hits = _parse_abricate_tsv(tsv, organism_names=["Escherichia coli"])
    assert len(hits) == 1
    assert hits[0].organism_match == "Escherichia coli"


def test_parse_abricate_tsv_unknown_organism():
    tsv = (
        "#FILE\tSEQUENCE\tSTART\tEND\tGENE\tCOVERAGE\tGAPS\t%COVERAGE\t%IDENTITY\t"
        "DATABASE\tACCESSION\tPRODUCT\tRESISTANCE\n"
        "amr_reads.fa\tsome_random_contig\t1\t100\tvanA\t1-100/100\t0\t100\t98.1\t"
        "card\tNC_000003\tVancomycin resistance protein\tglycopeptide\n"
    )
    hits = _parse_abricate_tsv(tsv, organism_names=["Staphylococcus aureus"])
    meca = next((h for h in hits if h.gene == "vanA"), None)
    assert meca is not None
    assert meca.organism_match == "unknown"


def test_parse_abricate_tsv_attributes_contig_via_map():
    # contig-based input: SEQUENCE is a bare contig id with no organism name, so
    # substring matching fails -- the contig_to_org map is what resolves it.
    tsv = (
        "#FILE\tSEQUENCE\tSTART\tEND\tGENE\tCOVERAGE\tGAPS\t%COVERAGE\t%IDENTITY\t"
        "DATABASE\tACCESSION\tPRODUCT\tRESISTANCE\n"
        "contigs.fa\tk141_42\t1\t100\tblaCTX-M\t1-100/100\t0\t100\t99.0\t"
        "card\tNC\tESBL\tbeta-lactam\n"
    )
    # without the map: unknown (the demo bug)
    assert _parse_abricate_tsv(tsv, ["Escherichia coli"])[0].organism_match == "unknown"
    # with the map: attributed to the finding
    hits = _parse_abricate_tsv(tsv, ["Escherichia coli"], {"k141_42": "Escherichia coli"})
    assert hits[0].organism_match == "Escherichia coli"


def test_map_contigs_to_organisms_assigns_best_genome(tmp_path):
    cfg = _cfg(tmp_path)
    contigs = tmp_path / "contigs.fa"
    contigs.write_text(">k141_1\nACGT\n>k141_2\nTTTT\n")
    g0, g1 = tmp_path / "ecoli.fna", tmp_path / "saureus.fna"
    g0.write_text(">e\nA\n"); g1.write_text(">s\nA\n")
    hits = [SimpleNamespace(name="Escherichia coli", genome_path=g0),
            SimpleNamespace(name="Staphylococcus aureus", genome_path=g1)]
    # PAF: cols 0=query(contig) ... 9=residue matches. k141_1 aligns best to g0,
    # k141_2 best to g1 (later, higher-match line wins).
    paf_by_genome = {
        str(g0): "k141_1\t4\t0\t4\t+\te\t1\t0\t4\t40\t40\t60\nk141_2\t4\t0\t4\t+\te\t1\t0\t4\t5\t40\t60\n",
        str(g1): "k141_2\t4\t0\t4\t+\ts\t1\t0\t4\t50\t40\t60\n",
    }

    def fake_run(cmd, **kw):
        out = Path(cmd[cmd.index("-o") + 1])
        out.write_text(paf_by_genome[cmd[-2]])   # cmd[-2] = genome, cmd[-1] = contigs
        return SimpleNamespace(returncode=0)

    with patch("pathogeniq.amr.shutil.which", return_value="minimap2"), \
         patch("pathogeniq.amr.subprocess.run", side_effect=fake_run):
        mapping = map_contigs_to_organisms(cfg, contigs, hits)
    assert mapping == {"k141_1": "Escherichia coli", "k141_2": "Staphylococcus aureus"}


def test_map_contigs_to_organisms_non_blocking(tmp_path):
    # no minimap2 -> empty map, run survives
    with patch("pathogeniq.amr.shutil.which", return_value=None):
        assert map_contigs_to_organisms(_cfg(tmp_path), tmp_path / "c.fa", []) == {}


def test_run_amr_screen_skips_when_abricate_missing(tmp_path):
    cfg = _cfg(tmp_path)
    reads = tmp_path / "reads.fq"
    reads.write_text("@r1\nACGT\n+\nAAAA\n")
    with patch("pathogeniq.amr.shutil.which", return_value=None):
        hits = run_amr_screen(cfg, reads, organism_names=["E. coli"])
    assert hits == []


def test_run_amr_screen_with_mock_abricate(tmp_path):
    cfg = _cfg(tmp_path)
    reads = tmp_path / "reads.fq"
    reads.write_text("@r1\nACGT\n+\nAAAA\n")

    tsv_output = (
        "#FILE\tSEQUENCE\tSTART\tEND\tGENE\tCOVERAGE\tGAPS\t%COVERAGE\t%IDENTITY\t"
        "DATABASE\tACCESSION\tPRODUCT\tRESISTANCE\n"
        "amr_reads.fa\tStaphylococcus_aureus_read1\t1\t100\tmecA\t1-100/100\t0\t100\t99.5\t"
        "card\tNC_000001\tmethicillin resistance protein\tbeta-lactam\n"
    )

    with patch("pathogeniq.amr.shutil.which", return_value="/usr/bin/abricate"):
        with patch(
            "pathogeniq.amr.subprocess.run",
            return_value=type("obj", (object,), {"returncode": 0, "stdout": tsv_output})(),
        ):
            hits = run_amr_screen(cfg, reads, organism_names=["Staphylococcus aureus"])

    assert len(hits) == 1
    assert hits[0].gene == "mecA"


def test_run_virulence_screen_with_mock_vfdb(tmp_path):
    cfg = _cfg(tmp_path)
    reads = tmp_path / "reads.fq"
    reads.write_text("@r1\nACGT\n+\nAAAA\n")
    # VFDB hit: the factor lives in PRODUCT; RESISTANCE is empty for virulence
    tsv_output = (
        "#FILE\tSEQUENCE\tSTART\tEND\tGENE\tCOVERAGE\tGAPS\t%COVERAGE\t%IDENTITY\t"
        "DATABASE\tACCESSION\tPRODUCT\tRESISTANCE\n"
        "amr_reads.fa\tStaphylococcus_aureus_r1\t1\t100\thlb\t1-100/100\t0\t100\t97.0\t"
        "vfdb\tVFG001\tbeta-hemolysin\t\n"
    )
    with patch("pathogeniq.amr.shutil.which", return_value="/usr/bin/abricate"), \
         patch("pathogeniq.amr.subprocess.run",
               return_value=type("obj", (object,), {"returncode": 0, "stdout": tsv_output})()):
        hits = run_virulence_screen(cfg, reads, organism_names=["Staphylococcus aureus"])
    assert len(hits) == 1
    assert hits[0].gene == "hlb"
    assert hits[0].factor == "beta-hemolysin"      # PRODUCT, not RESISTANCE
    assert hits[0].organism_match == "Staphylococcus aureus"
    assert hits[0].database == "vfdb"


def test_run_virulence_screen_skips_when_abricate_missing(tmp_path):
    cfg = _cfg(tmp_path)
    reads = tmp_path / "reads.fq"
    reads.write_text("@r1\nACGT\n+\nAAAA\n")
    with patch("pathogeniq.amr.shutil.which", return_value=None):
        assert run_virulence_screen(cfg, reads, organism_names=["E. coli"]) == []
