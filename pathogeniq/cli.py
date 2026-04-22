import click
from pathlib import Path

from .config import PipelineConfig, ReadType, SpecimenType
from .qc import run_qc
from .host_remove import run_host_removal
from .sketch import run_sketch_screen
from .align import run_targeted_alignment
from .em import em_abundance, bootstrap_ci
from .report import write_report


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
@click.option("--sketch-threshold", default=0.01, show_default=True)
@click.option("--n-bootstrap", default=100, show_default=True)
def run(input_fastq, output_dir, db_tier1, host_reference, specimen, read_type, threads, sketch_threshold, n_bootstrap):
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
    )

    click.echo("[1/5] QC & adapter trimming...")
    filtered, qc_metrics = run_qc(cfg)
    click.echo(f"      {qc_metrics.passing_reads:,} reads pass QC")

    click.echo("[2/5] Host removal...")
    nonhuman, hr_metrics = run_host_removal(cfg, filtered)
    click.echo(f"      Microbial fraction: {hr_metrics.microbial_fraction:.2%}")

    click.echo("[3/5] Sketch screening...")
    hits = run_sketch_screen(cfg, nonhuman)
    click.echo(f"      {len(hits)} candidate organisms shortlisted")

    if not hits:
        click.echo("No pathogens detected above threshold.")
        return

    click.echo("[4/5] Targeted alignment + EM abundance...")
    align_result = run_targeted_alignment(cfg, nonhuman, hits)
    em_result = em_abundance(align_result.alignment_matrix)
    ci_lower, ci_upper = bootstrap_ci(align_result.alignment_matrix, n_bootstrap=cfg.n_bootstrap)

    click.echo("[5/5] Generating report...")
    report_dir = write_report(cfg, align_result.organism_names, em_result, ci_lower, ci_upper)
    click.echo(f"Report written to: {report_dir}")
