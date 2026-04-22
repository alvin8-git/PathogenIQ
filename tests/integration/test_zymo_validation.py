"""
Integration test: run core pipeline on ZymoBIOMICS replicate 1 and verify
all 10 community members are detected at expected abundances.

Run with:
    HUMAN_REF=/path/to/GRCh38.fa TIER1_DB=/path/to/tier1.sbt.zip \
        pytest tests/integration/ -v -m integration

Requires:
  - fastp, minimap2, bwa-mem2, sourmash, samtools on PATH
  - HUMAN_REF env var: path to BWA-MEM2-indexed GRCh38
  - TIER1_DB env var: path to sourmash SBT database for Tier-1 pathogens
"""
import json
import os
import pytest
from pathlib import Path
from click.testing import CliRunner

from pathogeniq.cli import cli

ZYMO_FASTQ = Path(
    "/data/alvin/Metagenomics/ZymoStandardsFromMao/V350234554_L01_UDB-32.Zymo_Std.fq.gz"
)

EXPECTED_ORGANISMS = [
    "Pseudomonas aeruginosa",
    "Escherichia coli",
    "Salmonella enterica",
    "Lactobacillus fermentum",
    "Enterococcus faecalis",
    "Staphylococcus aureus",
    "Listeria monocytogenes",
    "Bacillus subtilis",
    "Saccharomyces cerevisiae",
    "Cryptococcus neoformans",
]

# Expected gDNA abundance: 12% bacteria, 2% yeasts (ZymoBIOMICS D6300)
EXPECTED_ABUNDANCE = {org: 0.12 for org in EXPECTED_ORGANISMS[:8]}
EXPECTED_ABUNDANCE.update({org: 0.02 for org in EXPECTED_ORGANISMS[8:]})


def _get_env():
    human_ref = os.environ.get("HUMAN_REF")
    tier1_db = os.environ.get("TIER1_DB")
    return human_ref, tier1_db


@pytest.mark.integration
@pytest.mark.skipif(not ZYMO_FASTQ.exists(), reason="ZymoBIOMICS FASTQ not found at expected path")
def test_zymo_all_organisms_detected(tmp_path):
    human_ref, tier1_db = _get_env()
    if not human_ref or not tier1_db:
        pytest.skip("Set HUMAN_REF and TIER1_DB env vars to run integration tests")

    runner = CliRunner()
    result = runner.invoke(cli, [
        "run",
        "--input", str(ZYMO_FASTQ),
        "--output", str(tmp_path / "out"),
        "--db", tier1_db,
        "--host-ref", human_ref,
        "--specimen", "blood",
        "--read-type", "short",
        "--threads", "16",
    ])

    assert result.exit_code == 0, f"Pipeline failed:\n{result.output}"

    report_json = tmp_path / "out" / "report" / "pathogeniq_report.json"
    assert report_json.exists(), "Report JSON not created"

    with open(report_json) as f:
        report = json.load(f)

    detected = {f["organism"] for f in report["findings"]}

    missing = [org for org in EXPECTED_ORGANISMS if org not in detected]
    assert not missing, f"Missing organisms: {missing}"


@pytest.mark.integration
@pytest.mark.skipif(not ZYMO_FASTQ.exists(), reason="ZymoBIOMICS FASTQ not found at expected path")
def test_zymo_abundance_accuracy(tmp_path):
    human_ref, tier1_db = _get_env()
    if not human_ref or not tier1_db:
        pytest.skip("Set HUMAN_REF and TIER1_DB env vars to run integration tests")

    runner = CliRunner()
    runner.invoke(cli, [
        "run",
        "--input", str(ZYMO_FASTQ),
        "--output", str(tmp_path / "out"),
        "--db", tier1_db,
        "--host-ref", human_ref,
        "--specimen", "blood",
        "--read-type", "short",
        "--threads", "16",
    ])

    report_json = tmp_path / "out" / "report" / "pathogeniq_report.json"
    if not report_json.exists():
        pytest.skip("Report not generated — check test_zymo_all_organisms_detected first")

    with open(report_json) as f:
        report = json.load(f)

    findings_by_org = {f["organism"]: f for f in report["findings"]}

    errors = []
    for org, expected_pct in EXPECTED_ABUNDANCE.items():
        if org not in findings_by_org:
            continue
        measured_pct = findings_by_org[org]["abundance_pct"] / 100
        rel_error = abs(measured_pct - expected_pct) / expected_pct
        if rel_error >= 0.20:
            errors.append(
                f"{org}: expected {expected_pct:.1%}, got {measured_pct:.1%} "
                f"(relative error {rel_error:.1%})"
            )

    assert not errors, "Abundance errors:\n" + "\n".join(errors)


@pytest.mark.integration
@pytest.mark.skipif(not ZYMO_FASTQ.exists(), reason="ZymoBIOMICS FASTQ not found at expected path")
def test_zymo_no_false_positives(tmp_path):
    """No organism outside the ZymoBIOMICS panel should receive Grade A."""
    human_ref, tier1_db = _get_env()
    if not human_ref or not tier1_db:
        pytest.skip("Set HUMAN_REF and TIER1_DB env vars to run integration tests")

    runner = CliRunner()
    runner.invoke(cli, [
        "run",
        "--input", str(ZYMO_FASTQ),
        "--output", str(tmp_path / "out"),
        "--db", tier1_db,
        "--host-ref", human_ref,
        "--specimen", "blood",
        "--read-type", "short",
        "--threads", "16",
    ])

    report_json = tmp_path / "out" / "report" / "pathogeniq_report.json"
    if not report_json.exists():
        pytest.skip("Report not generated")

    with open(report_json) as f:
        report = json.load(f)

    grade_a_fps = [
        f["organism"]
        for f in report["findings"]
        if f["grade"] == "A" and f["organism"] not in EXPECTED_ORGANISMS
    ]
    assert not grade_a_fps, f"Grade A false positives: {grade_a_fps}"
