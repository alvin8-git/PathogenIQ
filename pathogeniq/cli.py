import click
from pathlib import Path

from .amr import run_amr_screen
from .config import PipelineConfig, ReadType, SpecimenType
from .contaminants import flag_contaminants
from .em import bootstrap_ci, em_abundance
from .host_remove import run_host_removal
from .pdf_report import write_pdf_report
from .qc import run_qc
from .report import ReportEntry, write_report
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
def run(input_fastq, output_dir, db_tier1, host_reference, specimen, read_type,
        threads, sketch_threshold, n_bootstrap, amr_db, no_pdf):
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

    click.echo("[5/6] AMR screening...")
    amr_hits = run_amr_screen(cfg, nonhuman, organism_names=align_result.organism_names, db=cfg.amr_db)
    if amr_hits:
        click.echo(f"      {len(amr_hits)} AMR gene(s) detected")
    else:
        click.echo("      No AMR genes detected (or abricate not installed)")

    click.echo("[6/6] Generating report...")
    read_counts = (em_result.abundances * em_result.n_reads).astype(int)
    entries = [
        ReportEntry(
            organism=name,
            abundance=float(em_result.abundances[i]),
            ci_lower=float(ci_lower[i]),
            ci_upper=float(ci_upper[i]),
            read_count=int(read_counts[i]),
            specimen_type=cfg.specimen_type,
        )
        for i, name in enumerate(align_result.organism_names)
    ]
    entries = flag_contaminants(entries)

    report_dir = write_report(
        cfg, align_result.organism_names, em_result, ci_lower, ci_upper,
        amr_hits=amr_hits, entries_override=entries,
    )

    if not no_pdf:
        pdf_path = write_pdf_report(cfg, entries, amr_hits)
        click.echo(f"PDF report:        {pdf_path}")

    click.echo(f"Report written to: {report_dir}")
