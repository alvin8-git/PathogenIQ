import click
from dataclasses import replace
from pathlib import Path

from .amr import run_amr_screen, run_virulence_screen
from .background import (
    build_background,
    load_background_table,
    load_default_background,
)
from .config import PipelineConfig, ReadType, SpecimenType
from .em import bootstrap_ci, em_abundance
from .host_remove import run_host_removal, run_phix_removal
from .html_report import write_html_report
from .pdf_report import write_pdf_report
from .qc import run_qc
from .report import build_entries, write_report
from .sketch import run_sketch_screen
from .align import run_targeted_alignment


@click.group()
def cli():
    """PathogenIQ — clinical metagenomics pipeline."""


@cli.command()
@click.option("--input", "input_fastq", required=True, type=click.Path(exists=True, path_type=Path))
@click.option("--output", "output_dir", required=True, type=click.Path(path_type=Path))
@click.option("--db", "db_tier1", required=True, type=click.Path(exists=True, path_type=Path))
@click.option("--host-ref", "host_reference", required=True, type=click.Path(exists=True, path_type=Path))
@click.option("--specimen", type=click.Choice([s.value for s in SpecimenType]), required=True)
@click.option("--read-type", type=click.Choice([r.value for r in ReadType]), default="short", show_default=True)
@click.option("--threads", default=8, show_default=True)
@click.option("--sketch-threshold", default=0.003, show_default=True)
@click.option("--n-bootstrap", default=100, show_default=True)
@click.option("--amr-db", default="card", show_default=True, help="ABRicate database (card, resfinder, …)")
@click.option("--no-pdf", is_flag=True, default=False, help="Skip PDF report generation")
@click.option("--ntc", "ntc_fastq", type=click.Path(exists=True, path_type=Path), default=None,
              help="Batch-matched no-template control FASTQ (Tier 1)")
@click.option("--background", "background_table", type=click.Path(exists=True, path_type=Path), default=None,
              help="Precomputed pooled background table (Tier 2); takes precedence over --ntc")
@click.option("--no-background", is_flag=True, default=False,
              help="Disable NTC background correction (Tier 3, uncorrected)")
def run(input_fastq, output_dir, db_tier1, host_reference, specimen, read_type,
        threads, sketch_threshold, n_bootstrap, amr_db, no_pdf,
        ntc_fastq, background_table, no_background):
    """Run the full PathogenIQ pipeline."""
    output_dir.mkdir(parents=True, exist_ok=True)

    cfg = PipelineConfig(
        input_fastq=input_fastq,
        read_type=ReadType(read_type),
        specimen_type=SpecimenType(specimen),
        output_dir=output_dir,
        db_tier1=db_tier1,
        host_reference=host_reference,
        threads=threads,
        sketch_threshold=sketch_threshold,
        n_bootstrap=n_bootstrap,
        amr_db=amr_db,
    )

    click.echo("[1/6] QC & adapter trimming...")
    filtered, qc_metrics = run_qc(cfg)
    click.echo(f"      {qc_metrics.passing_reads:,} reads pass QC")

    click.echo("[2/6] Host removal...")
    nonhuman, hr_metrics = run_host_removal(cfg, filtered)
    click.echo(f"      Microbial fraction: {hr_metrics.microbial_fraction:.2%}")
    nonhuman, n_phix = run_phix_removal(cfg, nonhuman)
    if n_phix:
        click.echo(f"      {n_phix:,} PhiX spike-in reads removed")

    click.echo("[3/6] Sketch screening...")
    hits = run_sketch_screen(cfg, nonhuman)
    click.echo(f"      {len(hits)} candidate organisms shortlisted")

    if not hits:
        click.echo("No pathogens detected above threshold.")
        return

    click.echo("[4/6] Targeted alignment + EM abundance...")
    align_result = run_targeted_alignment(cfg, nonhuman, hits)
    em_result = em_abundance(align_result.alignment_matrix)
    ci_lower, ci_upper = bootstrap_ci(align_result.alignment_matrix, n_bootstrap=cfg.n_bootstrap)

    click.echo("[5/6] AMR & virulence screening...")
    amr_hits = run_amr_screen(cfg, nonhuman, organism_names=align_result.organism_names, db=cfg.amr_db)
    if amr_hits:
        click.echo(f"      {len(amr_hits)} AMR gene(s) detected")
    else:
        click.echo("      No AMR genes detected (or abricate not installed)")
    virulence_hits = run_virulence_screen(cfg, nonhuman, organism_names=align_result.organism_names)
    if virulence_hits:
        click.echo(f"      {len(virulence_hits)} virulence factor(s) detected (VFDB)")

    click.echo("[6/6] Background correction & report...")
    background = _resolve_background(cfg, ntc_fastq, background_table, no_background)
    tier = background.tier if background is not None else 3
    click.echo(f"      NTC background tier: {tier}"
               + ("" if background is not None else " (uncorrected)"))

    entries = build_entries(
        cfg, align_result.organism_names, em_result, ci_lower, ci_upper,
        align_result.taxon_ids, background=background,
    )
    removed = len(align_result.organism_names) - len(entries)
    if removed:
        click.echo(f"      {removed} taxon(s) removed as background")

    report_dir = write_report(cfg, entries, em_result, amr_hits=amr_hits, virulence_hits=virulence_hits)

    if not no_pdf:
        pdf_path = write_pdf_report(cfg, entries, amr_hits, virulence_hits=virulence_hits)
        click.echo(f"PDF report:        {pdf_path}")

    html_path = write_html_report(
        cfg, qc_metrics, hr_metrics, hits, entries, em_result,
        amr_hits=amr_hits, virulence_hits=virulence_hits,
    )
    click.echo(f"HTML report:       {html_path}")

    click.echo(f"Report written to: {report_dir}")


def _classify_taxon_counts(cfg: PipelineConfig, fastq: Path) -> tuple[dict[str, int], int]:
    """Classify an NTC FASTQ through the same QC->host->sketch->align stages the
    sample uses, returning (per-taxon read counts keyed by taxon_id, total
    classified reads). Returns ({}, 0) when nothing classifies."""
    ntc_cfg = replace(cfg, input_fastq=fastq, output_dir=cfg.output_dir / "ntc")
    ntc_cfg.output_dir.mkdir(parents=True, exist_ok=True)
    filtered, _ = run_qc(ntc_cfg)
    nonhuman, _ = run_host_removal(ntc_cfg, filtered)
    nonhuman, _ = run_phix_removal(ntc_cfg, nonhuman)
    hits = run_sketch_screen(ntc_cfg, nonhuman)
    if not hits:
        return {}, 0
    align = run_targeted_alignment(ntc_cfg, nonhuman, hits)
    total = len(align.read_ids)
    counts: dict[str, int] = {}
    for j, taxon_id in enumerate(align.taxon_ids):
        counts[taxon_id] = counts.get(taxon_id, 0) + int(align.alignment_matrix[:, j].sum())
    return counts, total


def _resolve_background(cfg, ntc_fastq, background_table, no_background):
    """Resolve the NTC background model. Precedence: --background > --ntc > none.
    Returns None for Tier 3 (no correction)."""
    if no_background:
        return None
    if background_table is not None:
        return load_background_table(background_table)
    if ntc_fastq is not None:
        counts, total = _classify_taxon_counts(cfg, ntc_fastq)
        return build_background([(counts, total)], tier=1)
    # Default: the packaged pooled background (Tier 2) if curated; else Tier 3.
    return load_default_background()
