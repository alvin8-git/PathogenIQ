from __future__ import annotations

from datetime import datetime
from pathlib import Path

from reportlab.lib import colors
from reportlab.lib.pagesizes import letter
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import inch
from reportlab.platypus import (
    Paragraph,
    SimpleDocTemplate,
    Spacer,
    Table,
    TableStyle,
)

from .amr import AMRHit
from .config import PipelineConfig
from .report import EvidenceGrade, ReportEntry, grade_mag


_GRADE_COLORS: dict[EvidenceGrade, str] = {
    EvidenceGrade.A: "#2e7d32",
    EvidenceGrade.B: "#e65100",
    EvidenceGrade.C: "#757575",
    EvidenceGrade.X: "#c62828",
}


def write_pdf_report(
    cfg: PipelineConfig,
    entries: list[ReportEntry],
    amr_hits: list[AMRHit],
    virulence_hits: list | None = None,
    mags: list | None = None,
) -> Path:
    """Render a single-page clinical PDF report."""
    out = cfg.output_dir / "report"
    out.mkdir(parents=True, exist_ok=True)
    pdf_path = out / "pathogeniq_report.pdf"

    doc = SimpleDocTemplate(
        str(pdf_path),
        pagesize=letter,
        rightMargin=0.75 * inch,
        leftMargin=0.75 * inch,
        topMargin=0.75 * inch,
        bottomMargin=0.75 * inch,
    )

    styles = getSampleStyleSheet()
    title_style = ParagraphStyle(
        "Title",
        parent=styles["Heading1"],
        fontSize=18,
        spaceAfter=12,
    )
    subtitle_style = ParagraphStyle(
        "Subtitle",
        parent=styles["Normal"],
        fontSize=10,
        textColor=colors.grey,
        spaceAfter=18,
    )
    heading2_style = ParagraphStyle(
        "Heading2",
        parent=styles["Heading2"],
        fontSize=14,
        spaceAfter=8,
        spaceBefore=12,
    )

    story: list = []

    # Header
    story.append(Paragraph("PathogenIQ Clinical Report", title_style))
    story.append(
        Paragraph(
            f"Sample: <b>{cfg.input_fastq.name}</b> &nbsp;|&nbsp; "
            f"Specimen: <b>{cfg.specimen_type.value.upper()}</b> &nbsp;|&nbsp; "
            f"Date: {datetime.now().strftime('%Y-%m-%d %H:%M')}",
            subtitle_style,
        )
    )

    # Findings table
    story.append(Paragraph("Findings", heading2_style))
    if entries:
        # Absolute-copies column only when a spike-in anchored quantification (Plan 6 #2);
        # otherwise every row is empty, so drop it.
        show_copies = any(e.absolute_copies is not None for e in entries)
        header = ["Organism", "Abundance (%)", "CI Lower", "CI Upper", "Reads"]
        if show_copies:
            header.append("Abs. copies")
        header += ["Grade", "Contaminant"]
        findings_data = [header]
        for e in entries:
            grade_color = _GRADE_COLORS.get(e.grade, colors.black)
            row = [
                e.organism,
                f"{e.abundance * 100:.2f}",
                f"{e.ci_lower * 100:.2f}",
                f"{e.ci_upper * 100:.2f}",
                str(e.read_count),
            ]
            if show_copies:
                row.append("—" if e.absolute_copies is None else f"{e.absolute_copies:.3g}")
            row += [
                Paragraph(f'<font color="{grade_color}"><b>{e.grade.value}</b></font>', styles["Normal"]),
                "Yes" if e.contaminant_risk else "No",
            ]
            findings_data.append(row)
        findings_table = Table(findings_data, repeatRows=1)
        findings_table.setStyle(TableStyle([
            ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#f5f5f5")),
            ("TEXTCOLOR", (0, 0), (-1, 0), colors.black),
            ("ALIGN", (0, 0), (-1, 0), "CENTER"),
            ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
            ("FONTSIZE", (0, 0), (-1, 0), 10),
            ("BOTTOMPADDING", (0, 0), (-1, 0), 8),
            ("BACKGROUND", (0, 1), (-1, -1), colors.white),
            ("GRID", (0, 0), (-1, -1), 0.5, colors.grey),
            ("FONTSIZE", (0, 1), (-1, -1), 9),
            ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
        ]))
        story.append(findings_table)
    else:
        story.append(Paragraph("No organisms detected.", styles["Normal"]))

    # AMR table
    if amr_hits:
        story.append(Spacer(1, 12))
        story.append(Paragraph("Antimicrobial Resistance Genes", heading2_style))
        amr_data = [["Gene", "Drug Class", "Identity (%)", "Coverage (%)", "Organism Match"]]
        for h in amr_hits:
            amr_data.append([
                h.gene,
                h.drug_class,
                f"{h.identity_pct:.1f}",
                f"{h.coverage_pct:.1f}",
                h.organism_match,
            ])
        amr_table = Table(amr_data, repeatRows=1)
        amr_table.setStyle(TableStyle([
            ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#f5f5f5")),
            ("TEXTCOLOR", (0, 0), (-1, 0), colors.black),
            ("ALIGN", (0, 0), (-1, 0), "CENTER"),
            ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
            ("FONTSIZE", (0, 0), (-1, 0), 10),
            ("BOTTOMPADDING", (0, 0), (-1, 0), 8),
            ("BACKGROUND", (0, 1), (-1, -1), colors.white),
            ("GRID", (0, 0), (-1, -1), 0.5, colors.grey),
            ("FONTSIZE", (0, 1), (-1, -1), 9),
            ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
        ]))
        story.append(amr_table)

    # Virulence factor table (VFDB)
    if virulence_hits:
        story.append(Spacer(1, 12))
        story.append(Paragraph("Virulence Factors (VFDB)", heading2_style))
        vir_data = [["Gene", "Virulence Factor", "Identity (%)", "Coverage (%)", "Organism Match"]]
        for h in virulence_hits:
            vir_data.append([
                h.gene,
                h.factor,
                f"{h.identity_pct:.1f}",
                f"{h.coverage_pct:.1f}",
                h.organism_match,
            ])
        vir_table = Table(vir_data, repeatRows=1)
        vir_table.setStyle(TableStyle([
            ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#f5f5f5")),
            ("TEXTCOLOR", (0, 0), (-1, 0), colors.black),
            ("ALIGN", (0, 0), (-1, 0), "CENTER"),
            ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
            ("FONTSIZE", (0, 0), (-1, 0), 10),
            ("BOTTOMPADDING", (0, 0), (-1, 0), 8),
            ("BACKGROUND", (0, 1), (-1, -1), colors.white),
            ("GRID", (0, 0), (-1, -1), 0.5, colors.grey),
            ("FONTSIZE", (0, 1), (-1, -1), 9),
            ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
        ]))
        story.append(vir_table)

    # MAG table (open-world assembly arm, Plan 6 #3)
    if mags:
        story.append(Spacer(1, 12))
        story.append(Paragraph("Metagenome-Assembled Genomes (MAGs)", heading2_style))
        mag_data = [["Bin", "GTDB Taxonomy", "Complete (%)", "Contam (%)", "Contigs", "Size (Mb)", "Grade"]]
        for m in mags:
            # ponytail: renderer grades on CheckM QC only; the marker rescue
            # (completeness=None + pathogenicity markers -> C) is JSON-only.
            g = grade_mag(m)
            gc = _GRADE_COLORS.get(g, colors.black)
            mag_data.append([
                Paragraph(m.bin_id, styles["Normal"]),
                Paragraph(m.taxonomy or "unclassified", styles["Normal"]),
                "—" if m.completeness is None else f"{m.completeness:.1f}",
                "—" if m.contamination is None else f"{m.contamination:.1f}",
                str(m.n_contigs),
                f"{m.total_bp / 1e6:.2f}",
                Paragraph(f'<font color="{gc}"><b>{g.value}</b></font>', styles["Normal"]),
            ])
        mag_table = Table(mag_data, repeatRows=1)
        mag_table.setStyle(TableStyle([
            ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#f5f5f5")),
            ("TEXTCOLOR", (0, 0), (-1, 0), colors.black),
            ("ALIGN", (0, 0), (-1, 0), "CENTER"),
            ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
            ("FONTSIZE", (0, 0), (-1, 0), 10),
            ("BOTTOMPADDING", (0, 0), (-1, 0), 8),
            ("BACKGROUND", (0, 1), (-1, -1), colors.white),
            ("GRID", (0, 0), (-1, -1), 0.5, colors.grey),
            ("FONTSIZE", (0, 1), (-1, -1), 9),
            ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
        ]))
        story.append(mag_table)

    # Footer
    story.append(Spacer(1, 18))
    story.append(
        Paragraph(
            "<i>Grade definitions: A = high-confidence pathogen; B = probable pathogen; "
            "C = possible pathogen / contaminant; X = insufficient evidence. "
            "Contaminant flag indicates known specimen-specific environmental flora.</i>",
            ParagraphStyle(
                "Footer",
                parent=styles["Normal"],
                fontSize=8,
                textColor=colors.grey,
                leading=10,
            ),
        )
    )

    doc.build(story)
    return pdf_path
