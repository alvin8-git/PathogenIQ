from unittest.mock import patch

from pathogeniq.assembly import (
    _parse_checkm_table,
    _parse_gtdbtk_summary,
    _fasta_stats,
    run_megahit,
    run_assembly_stage,
)
from pathogeniq.config import PipelineConfig, ReadType, SpecimenType


def _cfg(tmp_path):
    return PipelineConfig(
        input_fastq=tmp_path / "r.fq.gz", read_type=ReadType.SHORT,
        specimen_type=SpecimenType.BLOOD, output_dir=tmp_path,
        db_tier1=tmp_path / "db", host_reference=tmp_path / "ref", threads=4,
    )


def test_parse_checkm_table():
    tsv = ("Bin Id\tMarker lineage\tCompleteness\tContamination\n"
           "bin.1\tk__Bacteria\t92.5\t1.3\n"
           "bin.2\tk__Bacteria\t40.0\t8.0\n")
    d = _parse_checkm_table(tsv)
    assert d["bin.1"] == (92.5, 1.3)
    assert d["bin.2"] == (40.0, 8.0)


def test_parse_gtdbtk_summary():
    tsv = ("user_genome\tclassification\tnote\n"
           "bin.1\td__Bacteria;p__Pseudomonadota;s__Escherichia coli\tok\n")
    d = _parse_gtdbtk_summary(tsv)
    assert d["bin.1"].startswith("d__Bacteria") and "Escherichia coli" in d["bin.1"]


def test_fasta_stats(tmp_path):
    fa = tmp_path / "bin.1.fa"
    fa.write_text(">c1\nACGTACGT\n>c2\nTTTT\n")
    assert _fasta_stats(fa) == (2, 12)   # 2 contigs, 8+4 bp


def test_run_megahit_missing_returns_none(tmp_path):
    with patch("pathogeniq.assembly.shutil.which", return_value=None):
        assert run_megahit(_cfg(tmp_path), tmp_path / "reads.fq.gz") is None


def test_run_assembly_stage_skips_without_megahit(tmp_path):
    with patch("pathogeniq.assembly.run_megahit", return_value=None):
        assert run_assembly_stage(_cfg(tmp_path), tmp_path / "reads.fq.gz") == []


def test_run_assembly_stage_filters_low_quality_and_builds_mags(tmp_path):
    bindir = tmp_path / "bins"
    bindir.mkdir()
    b1 = bindir / "bin.1.fa"
    b1.write_text(">c1\n" + "A" * 100 + "\n")
    b2 = bindir / "bin.2.fa"
    b2.write_text(">c1\n" + "A" * 50 + "\n")
    with patch("pathogeniq.assembly.run_megahit", return_value=tmp_path / "contigs.fa"), \
         patch("pathogeniq.assembly._contig_depths", return_value=tmp_path / "depth.txt"), \
         patch("pathogeniq.assembly.run_metabat2", return_value=[b1, b2]), \
         patch("pathogeniq.assembly.run_checkm",
               return_value={"bin.1": (80.0, 2.0), "bin.2": (40.0, 5.0)}), \
         patch("pathogeniq.assembly.run_gtdbtk",
               return_value={"bin.1": "d__Bacteria;p__Pseudomonadota;s__E coli"}):
        mags = run_assembly_stage(_cfg(tmp_path), tmp_path / "reads.fq.gz")
    assert [m.bin_id for m in mags] == ["bin.1"]   # bin.2 dropped: 40% < 50% completeness
    m = mags[0]
    assert m.completeness == 80.0 and m.contamination == 2.0
    assert m.taxonomy.startswith("d__Bacteria")
    assert m.total_bp == 100


def test_run_assembly_stage_keeps_bins_when_checkm_absent(tmp_path):
    bindir = tmp_path / "bins"
    bindir.mkdir()
    b1 = bindir / "bin.1.fa"
    b1.write_text(">c1\nACGT\n")
    with patch("pathogeniq.assembly.run_megahit", return_value=tmp_path / "contigs.fa"), \
         patch("pathogeniq.assembly._contig_depths", return_value=tmp_path / "depth.txt"), \
         patch("pathogeniq.assembly.run_metabat2", return_value=[b1]), \
         patch("pathogeniq.assembly.run_checkm", return_value={}), \
         patch("pathogeniq.assembly.run_gtdbtk", return_value={}):
        mags = run_assembly_stage(_cfg(tmp_path), tmp_path / "reads.fq.gz")
    assert len(mags) == 1 and mags[0].completeness is None   # no QC -> bin still kept
