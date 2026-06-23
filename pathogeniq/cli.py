import click
import numpy as np
from dataclasses import replace
from pathlib import Path

from .amr import run_amr_screen, run_virulence_screen
from .background import (
    build_background,
    load_background_table,
    load_default_background,
)
from .config import PipelineConfig, ReadType, SpecimenType
from .em import EMResult, bootstrap_ci, em_abundance
from .host_remove import run_host_removal, run_phix_removal
from .quantify import quantify_entries
from .assembly import run_assembly_stage, run_megahit
from .novelty import assess_novelty
from .viral import run_viral_stage
from .pathogenicity import triage_mags
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
@click.option("--spike-taxon", default=None,
              help="Spike-in control taxon_id (must be in --db) for absolute quantification")
@click.option("--spike-copies", type=float, default=None,
              help="Known copies/cells of spike added to the sample (anchors absolute load)")
@click.option("--sample-volume", type=float, default=None,
              help="Volume the sample represents (e.g. mL air/blood) -> copies per volume")
@click.option("--assemble", is_flag=True, default=False,
              help="Run the de novo assembly + MAG arm (MEGAHIT/MetaBAT2/CheckM/GTDB-Tk) to "
                   "recover novel organisms not in the reference DB (slow; tools/DBs required)")
@click.option("--novelty", "novelty_flag", is_flag=True, default=False,
              help="Open-world novelty trigger: classify reads against a broad Kraken2 DB and "
                   "report the unclassified fraction (flags novel/uncatalogued content; "
                   "needs kraken2 + a Standard DB at $KRAKEN2_DB or databases/kraken2)")
@click.option("--viral", "viral_flag", is_flag=True, default=False,
              help="Run the viral identification arm (geNomad + CheckV) on the assembly to "
                   "detect airborne viruses not in the bacterial DBs (slow; tools/DBs required)")
def run(input_fastq, output_dir, db_tier1, host_reference, specimen, read_type,
        threads, sketch_threshold, n_bootstrap, amr_db, no_pdf,
        ntc_fastq, background_table, no_background,
        spike_taxon, spike_copies, sample_volume, assemble, novelty_flag, viral_flag):
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

    # Open-world novelty trigger runs on the non-host reads BEFORE the targeted
    # arm, so a sample with no DB hit but novel content is still flagged.
    novelty = None
    if novelty_flag:
        click.echo("      Open-world novelty screen (Kraken2, broad DB)...")
        novelty = assess_novelty(cfg, nonhuman)
        if novelty is None:
            click.echo("      Novelty screen skipped (kraken2 or broad DB missing)")
        else:
            click.echo(f"      {novelty.unclassified_fraction:.1%} of reads unclassified"
                       + (" — FLAGGED: novel/uncatalogued content, consider --assemble"
                          if novelty.flagged else ""))

    def _discovery_arms(contigs=None):
        """Open-world arms (assembly MAGs + viral) on the non-host reads. Run
        regardless of targeted hits so novel/viral content with no DB match is
        still recovered. ``contigs`` (the shared megahit assembly used for AMR) is
        reused for the viral arm to avoid re-assembling."""
        m = v = None
        if assemble:
            click.echo("      De novo assembly + MAG recovery (this is slow)...")
            m = run_assembly_stage(cfg, nonhuman)
            click.echo(f"      {len(m)} MAG(s) recovered" if m else
                       "      No MAGs recovered (assembly tools/DBs missing or no bins)")
        if viral_flag:
            click.echo("      Viral identification arm (geNomad + CheckV)...")
            if contigs is None:
                contigs = run_megahit(cfg, nonhuman)
            v = run_viral_stage(cfg, contigs) if contigs is not None else []
            click.echo(f"      {len(v)} viral contig(s) identified" if v else
                       "      No viral contigs (geNomad/DBs missing or none found)")
        # R4: triage recovered MAGs into pathogen vs environmental.
        t = None
        if m:
            t = triage_mags(cfg, m)
            n_path = sum(a.verdict != "ENVIRONMENTAL" for a in t)
            click.echo(f"      Pathogenicity triage: {n_path}/{len(t)} MAG(s) flagged "
                       f"pathogen-relevant (markers or pathogen-adjacent lineage)")
        return m, v, t

    click.echo("[3/6] Sketch screening...")
    hits = run_sketch_screen(cfg, nonhuman)
    click.echo(f"      {len(hits)} candidate organisms shortlisted")

    if not hits:
        click.echo("No targeted pathogens detected above threshold.")
        if novelty is not None and novelty.flagged:
            click.echo(f"NOTE: {novelty.unclassified_fraction:.1%} of reads are unclassified "
                       f"against the broad DB — possible novel organism.")
        # The open-world arms can still recover novel/viral genomes with no DB hit.
        mags, viral, triage = _discovery_arms()
        if any(x is not None for x in (novelty, mags, viral)):
            empty_em = EMResult(abundances=np.array([]), n_reads=0, n_organisms=0, iterations=0)
            report_dir = write_report(cfg, [], empty_em, novelty=novelty, mags=mags,
                                      viral=viral, pathogenicity=triage)
            click.echo(f"Report written to: {report_dir}")
        return

    click.echo("[4/6] Targeted alignment + EM abundance...")
    align_result = run_targeted_alignment(cfg, nonhuman, hits)
    em_result = em_abundance(align_result.alignment_matrix)
    ci_lower, ci_upper = bootstrap_ci(align_result.alignment_matrix, n_bootstrap=cfg.n_bootstrap)

    click.echo("[5/6] AMR & virulence screening...")
    # Assemble once: AMR/virulence (and the viral arm) screen the contigs, not the
    # millions of raw reads — abricate-on-reads was the pipeline's slowest stage.
    contigs = run_megahit(cfg, nonhuman)
    if contigs is None:
        click.echo("      Assembly unavailable (megahit missing) — AMR/virulence skipped")
    amr_hits = run_amr_screen(cfg, contigs, organism_names=align_result.organism_names, db=cfg.amr_db)
    if amr_hits:
        click.echo(f"      {len(amr_hits)} AMR gene(s) detected")
    else:
        click.echo("      No AMR genes detected (or abricate/assembly not available)")
    virulence_hits = run_virulence_screen(cfg, contigs, organism_names=align_result.organism_names)
    if virulence_hits:
        click.echo(f"      {len(virulence_hits)} virulence factor(s) detected (VFDB)")

    click.echo("[6/6] Background correction & report...")
    background = _resolve_background(cfg, ntc_fastq, background_table, no_background)
    tier = background.tier if background is not None else 3
    click.echo(f"      NTC background tier: {tier}"
               + ("" if background is not None else " (uncorrected)"))

    entries = build_entries(
        cfg, align_result.organism_names, em_result, ci_lower, ci_upper,
        align_result.taxon_ids, background=background, coverage=align_result.coverage,
    )
    removed = len(align_result.organism_names) - len(entries)
    if removed:
        click.echo(f"      {removed} taxon(s) removed as background")

    spike_info = None
    if spike_taxon and spike_copies:
        entries, spike_info = quantify_entries(
            entries, spike_taxon_id=spike_taxon, spike_copies=spike_copies,
            sample_volume=sample_volume,
        )
        if spike_info.found:
            click.echo(f"      Absolute quantification anchored on spike "
                       f"{spike_taxon} ({spike_info.spike_reads:,} reads = "
                       f"{spike_copies:g} copies)")
        else:
            click.echo(f"      WARNING: spike {spike_taxon} not detected — "
                       f"absolute quantification unavailable (check spike input)")

    mags, viral, triage = _discovery_arms(contigs)

    report_dir = write_report(cfg, entries, em_result, amr_hits=amr_hits,
                              virulence_hits=virulence_hits, spike_info=spike_info, mags=mags,
                              novelty=novelty, viral=viral, pathogenicity=triage)

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
