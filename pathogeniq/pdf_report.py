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
    # reportlab only wraps text inside a Paragraph (bare strings overflow), so every
    # text cell is wrapped in _p() and each table gets explicit colWidths summing to
    # the usable page width — long gene names / GTDB lineages now wrap instead of
    # running off the page.
    cell_style = ParagraphStyle("Cell", parent=styles["Normal"], fontSize=8, leading=10)
    header_style = ParagraphStyle("CellHead", parent=styles["Normal"], fontSize=9, leading=11)
    usable_w = letter[0] - doc.leftMargin - doc.rightMargin

    def _p(text):
        return Paragraph(str(text), cell_style)

    _TABLE_STYLE = TableStyle([
        ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#f5f5f5")),
        ("ALIGN", (0, 0), (-1, 0), "CENTER"),
        ("BOTTOMPADDING", (0, 0), (-1, 0), 8),
        ("BACKGROUND", (0, 1), (-1, -1), colors.white),
        ("GRID", (0, 0), (-1, -1), 0.5, colors.grey),
        ("FONTSIZE", (0, 1), (-1, -1), 9),
        ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
    ])

    def _table(headers, rows, weights):
        """A width-safe table: headers/text cells wrap, columns sized by `weights`
        (proportional, normalized to the usable page width)."""
        head = [Paragraph(f"<b>{h}</b>", header_style) for h in headers]
        total = sum(weights)
        col_w = [usable_w * w / total for w in weights]
        t = Table([head] + rows, repeatRows=1, colWidths=col_w)
        t.setStyle(_TABLE_STYLE)
        return t

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
        rows = []
        for e in entries:
            grade_color = _GRADE_COLORS.get(e.grade, colors.black)
            row = [
                _p(f"<i>{e.organism}</i>"),
                f"{e.abundance * 100:.2f}",
                f"{e.ci_lower * 100:.2f}",
                f"{e.ci_upper * 100:.2f}",
                str(e.read_count),
            ]
            if show_copies:
                row.append("—" if e.absolute_copies is None else f"{e.absolute_copies:.3g}")
            row += [
                Paragraph(f'<font color="{grade_color}"><b>{e.grade.value}</b></font>', cell_style),
                "Yes" if e.contaminant_risk else "No",
            ]
            rows.append(row)
        weights = [3.0, 1.2, 1.0, 1.0, 0.9] + ([1.1] if show_copies else []) + [0.7, 1.1]
        story.append(_table(header, rows, weights))
    else:
        story.append(Paragraph("No organisms detected.", styles["Normal"]))

    # AMR table
    if amr_hits:
        story.append(Spacer(1, 12))
        story.append(Paragraph("Antimicrobial Resistance Genes", heading2_style))
        header = ["Gene", "Drug Class", "Identity (%)", "Coverage (%)", "Organism Match"]
        rows = [
            [_p(h.gene), _p(h.drug_class), f"{h.identity_pct:.1f}", f"{h.coverage_pct:.1f}",
             _p(f"<i>{h.organism_match}</i>")]
            for h in amr_hits
        ]
        story.append(_table(header, rows, [1.6, 2.2, 1.0, 1.0, 2.2]))

    # Virulence factor table (VFDB)
    if virulence_hits:
        story.append(Spacer(1, 12))
        story.append(Paragraph("Virulence Factors (VFDB)", heading2_style))
        header = ["Gene", "Virulence Factor", "Identity (%)", "Coverage (%)", "Organism Match"]
        rows = [
            [_p(h.gene), _p(h.factor), f"{h.identity_pct:.1f}", f"{h.coverage_pct:.1f}",
             _p(f"<i>{h.organism_match}</i>")]
            for h in virulence_hits
        ]
        story.append(_table(header, rows, [1.6, 2.2, 1.0, 1.0, 2.2]))

    # MAG table (open-world assembly arm, Plan 6 #3)
    if mags:
        story.append(Spacer(1, 12))
        story.append(Paragraph("Metagenome-Assembled Genomes (MAGs)", heading2_style))
        header = ["Bin", "GTDB Taxonomy", "Complete (%)", "Contam (%)", "Contigs", "Size (Mb)", "Grade"]
        rows = []
        for m in mags:
            # ponytail: renderer grades on CheckM QC only; the marker rescue
            # (completeness=None + pathogenicity markers -> C) is JSON-only.
            g = grade_mag(m)
            gc = _GRADE_COLORS.get(g, colors.black)
            rows.append([
                _p(m.bin_id),
                _p(f"<i>{m.taxonomy or 'unclassified'}</i>"),
                "—" if m.completeness is None else f"{m.completeness:.1f}",
                "—" if m.contamination is None else f"{m.contamination:.1f}",
                str(m.n_contigs),
                f"{m.total_bp / 1e6:.2f}",
                Paragraph(f'<font color="{gc}"><b>{g.value}</b></font>', cell_style),
            ])
        story.append(_table(header, rows, [0.9, 3.0, 1.0, 0.9, 0.8, 0.9, 0.7]))

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
