from pathlib import Path
from unittest.mock import patch

import pytest

from pathogeniq.amr import (
    AMRHit,
    _fastq_to_fasta,
    _parse_abricate_tsv,
    run_amr_screen,
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


def test_run_amr_screen_skips_when_abricate_missing(tmp_path):
    cfg = _cfg(tmp_path)
    reads = tmp_path / "reads.fq"
    reads.write_text("@r1\nACGT\n+\nAAAA\n")
    with patch("pathogeniq.amr.shutil.which", return_value=None):
        hits = run_amr_screen(cfg, reads, organism_names=["E. coli"])
    assert hits == []


def test_fastq_to_fasta_conversion(tmp_path):
    fq = tmp_path / "reads.fq"
    fa = tmp_path / "reads.fa"
    fq.write_text("@read1\nACGTACGT\n+\nIIIIIIII\n@read2\nTGCATGCA\n+\nIIIIIIII\n")
    _fastq_to_fasta(fq, fa)
    content = fa.read_text()
    assert ">read1\nACGTACGT\n" in content
    assert ">read2\nTGCATGCA\n" in content


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
