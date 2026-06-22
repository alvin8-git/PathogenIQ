from __future__ import annotations

from datetime import datetime
from pathlib import Path

from .amr import AMRHit
from .config import PipelineConfig
from .em import EMResult
from .host_remove import HostRemovalMetrics
from .qc import QCMetrics
from .report import EvidenceGrade, ReportEntry
from .sketch import SketchHit

_GRADE_BADGE: dict[EvidenceGrade, tuple[str, str]] = {
    EvidenceGrade.A: ("#2e7d32", "white"),
    EvidenceGrade.B: ("#e65100", "white"),
    EvidenceGrade.C: ("#757575", "white"),
    EvidenceGrade.X: ("#c62828", "white"),
}

_CSS = """
* { box-sizing: border-box; margin: 0; padding: 0; }
body { font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif;
       background: #f4f6f9; color: #1a1a2e; font-size: 14px; }
.page { max-width: 1100px; margin: 0 auto; padding: 24px; }
header { background: #1a237e; color: white; padding: 24px 32px;
         border-radius: 10px; margin-bottom: 24px; }
header h1 { font-size: 26px; font-weight: 700; letter-spacing: 1px; }
header .subtitle { opacity: 0.8; margin-top: 4px; font-size: 13px; }
.meta-row { display: flex; gap: 24px; margin-top: 14px; flex-wrap: wrap; }
.meta-item { background: rgba(255,255,255,0.15); border-radius: 6px;
             padding: 8px 14px; }
.meta-item .label { font-size: 11px; opacity: 0.7; text-transform: uppercase; letter-spacing: .5px; }
.meta-item .value { font-size: 14px; font-weight: 600; margin-top: 2px; }
.card { background: white; border-radius: 10px; padding: 22px 26px;
        margin-bottom: 20px; box-shadow: 0 1px 4px rgba(0,0,0,0.08); }
.card-title { font-size: 16px; font-weight: 700; color: #1a237e;
              border-bottom: 2px solid #e8eaf6; padding-bottom: 10px;
              margin-bottom: 16px; }
.stat-grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(160px, 1fr)); gap: 14px; }
.stat-box { background: #f8f9ff; border: 1px solid #e8eaf6; border-radius: 8px;
            padding: 14px 16px; }
.stat-box .num { font-size: 22px; font-weight: 700; color: #1a237e; }
.stat-box .desc { font-size: 12px; color: #666; margin-top: 4px; }
.stat-box.highlight { background: #e8f5e9; border-color: #a5d6a7; }
.stat-box.highlight .num { color: #2e7d32; }
.stat-box.warn { background: #fff3e0; border-color: #ffcc02; }
.stat-box.warn .num { color: #e65100; }
.bar-wrap { display: flex; align-items: center; gap: 10px; margin: 8px 0; }
.bar-bg { flex: 1; background: #e8eaf6; border-radius: 4px; height: 12px; overflow: hidden; }
.bar-fill { height: 100%; border-radius: 4px; transition: width 0.3s; }
.bar-label { font-size: 12px; width: 52px; text-align: right; color: #555; }
table { width: 100%; border-collapse: collapse; }
th { background: #e8eaf6; color: #1a237e; font-size: 12px; text-transform: uppercase;
     letter-spacing: .4px; padding: 10px 12px; text-align: left; }
td { padding: 10px 12px; border-bottom: 1px solid #f0f0f0; vertical-align: middle; }
tr:last-child td { border-bottom: none; }
tr:hover td { background: #fafbff; }
.grade-badge { display: inline-block; padding: 3px 10px; border-radius: 12px;
               font-size: 12px; font-weight: 700; }
.contam-flag { display: inline-block; padding: 2px 8px; border-radius: 10px;
               font-size: 11px; background: #fff3e0; color: #e65100;
               border: 1px solid #ffcc80; }
.abund-bar-wrap { display: flex; align-items: center; gap: 8px; }
.abund-bar-bg { width: 90px; background: #e8eaf6; border-radius: 3px; height: 8px; overflow: hidden; }
.abund-bar-fill { height: 100%; border-radius: 3px; }
em { font-style: italic; }
.sketch-chip { display: inline-block; background: #e8eaf6; border-radius: 12px;
               padding: 3px 10px; font-size: 12px; margin: 2px; }
.no-data { color: #999; font-style: italic; font-size: 13px; padding: 12px 0; }
footer { text-align: center; color: #999; font-size: 12px; margin-top: 32px; padding-top: 16px;
         border-top: 1px solid #e0e0e0; }
"""


def _badge(grade: EvidenceGrade) -> str:
    bg, fg = _GRADE_BADGE[grade]
    return f'<span class="grade-badge" style="background:{bg};color:{fg}">Grade {grade.value}</span>'


def _abund_bar(pct: float, grade: EvidenceGrade) -> str:
    bg, _ = _GRADE_BADGE[grade]
    width = min(100, pct)
    return (
        f'<div class="abund-bar-wrap">'
        f'<div class="abund-bar-bg"><div class="abund-bar-fill" style="width:{width:.1f}%;background:{bg}"></div></div>'
        f'<span>{pct:.2f}%</span></div>'
    )


def write_html_report(
    cfg: PipelineConfig,
    qc_metrics: QCMetrics,
    hr_metrics: HostRemovalMetrics,
    sketch_hits: list[SketchHit],
    entries: list[ReportEntry],
    em_result: EMResult,
    amr_hits: list[AMRHit] | None = None,
    virulence_hits: list | None = None,
) -> Path:
    out = cfg.output_dir / "report"
    out.mkdir(parents=True, exist_ok=True)

    # entries come pre-built (tier, contaminants, background filter applied) from
    # report.build_entries — the same list the JSON/TSV/PDF renderers use.
    entries = sorted(entries, key=lambda e: e.abundance, reverse=True)

    now = datetime.now().strftime("%Y-%m-%d %H:%M")
    sample_name = cfg.input_fastq.name

    # ── Header ──────────────────────────────────────────────────────────
    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>PathogenIQ Report — {sample_name}</title>
<style>{_CSS}</style>
</head>
<body>
<div class="page">
<header>
  <h1>PathogenIQ</h1>
  <div class="subtitle">Clinical Metagenomics Pipeline Report</div>
  <div class="meta-row">
    <div class="meta-item"><div class="label">Sample</div><div class="value">{sample_name}</div></div>
    <div class="meta-item"><div class="label">Specimen</div><div class="value">{cfg.specimen_type.value.upper()}</div></div>
    <div class="meta-item"><div class="label">Read Type</div><div class="value">{cfg.read_type.value.capitalize()}</div></div>
    <div class="meta-item"><div class="label">Generated</div><div class="value">{now}</div></div>
  </div>
</header>
"""

    # ── QC Card ──────────────────────────────────────────────────────────
    pass_pct = qc_metrics.pass_rate * 100
    qc_class = "highlight" if pass_pct >= 90 else "warn"
    html += f"""
<div class="card">
  <div class="card-title">Quality Control</div>
  <div class="stat-grid">
    <div class="stat-box"><div class="num">{qc_metrics.total_reads:,}</div><div class="desc">Input reads</div></div>
    <div class="stat-box {qc_class}"><div class="num">{qc_metrics.passing_reads:,}</div><div class="desc">Reads passing QC</div></div>
    <div class="stat-box {qc_class}"><div class="num">{pass_pct:.1f}%</div><div class="desc">Pass rate</div></div>
  </div>
  <div style="margin-top:14px">
    <div class="bar-wrap">
      <span style="font-size:12px;width:120px;color:#555">Passing reads</span>
      <div class="bar-bg"><div class="bar-fill" style="width:{pass_pct:.1f}%;background:#2e7d32"></div></div>
      <span class="bar-label">{pass_pct:.1f}%</span>
    </div>
  </div>
  <div style="margin-top:10px;font-size:12px;color:#888">
    Filters: Q≥20, length≥50bp, low-complexity removal
  </div>
</div>
"""

    # ── Host Removal Card ─────────────────────────────────────────────────
    mic_pct = hr_metrics.microbial_fraction * 100
    hum_pct = hr_metrics.human_fraction * 100
    html += f"""
<div class="card">
  <div class="card-title">Host Removal</div>
  <div class="stat-grid">
    <div class="stat-box"><div class="num">{hr_metrics.total_reads:,}</div><div class="desc">Total reads</div></div>
    <div class="stat-box warn"><div class="num">{hr_metrics.human_reads:,}</div><div class="desc">Human reads removed</div></div>
    <div class="stat-box highlight"><div class="num">{hr_metrics.nonhuman_reads:,}</div><div class="desc">Microbial reads</div></div>
    <div class="stat-box highlight"><div class="num">{mic_pct:.2f}%</div><div class="desc">Microbial fraction</div></div>
  </div>
  <div style="margin-top:14px">
    <div class="bar-wrap">
      <span style="font-size:12px;width:120px;color:#555">Microbial</span>
      <div class="bar-bg"><div class="bar-fill" style="width:{mic_pct:.2f}%;background:#2e7d32"></div></div>
      <span class="bar-label">{mic_pct:.1f}%</span>
    </div>
    <div class="bar-wrap">
      <span style="font-size:12px;width:120px;color:#555">Human</span>
      <div class="bar-bg"><div class="bar-fill" style="width:{hum_pct:.2f}%;background:#e65100"></div></div>
      <span class="bar-label">{hum_pct:.1f}%</span>
    </div>
  </div>
</div>
"""

    # ── Sketch Screening Card ─────────────────────────────────────────────
    hits_html = "".join(
        f'<span class="sketch-chip"><em>{h.name}</em> ({h.containment:.3f})</span>'
        for h in sorted(sketch_hits, key=lambda x: x.containment, reverse=True)
    ) if sketch_hits else '<span class="no-data">No candidates above threshold</span>'

    html += f"""
<div class="card">
  <div class="card-title">Sketch Screening</div>
  <div class="stat-grid" style="margin-bottom:14px">
    <div class="stat-box"><div class="num">{len(sketch_hits)}</div><div class="desc">Candidate organisms</div></div>
    <div class="stat-box"><div class="num">{cfg.sketch_threshold:.3f}</div><div class="desc">Containment threshold</div></div>
  </div>
  <div style="font-size:12px;color:#555;margin-bottom:6px">Candidates (containment score):</div>
  {hits_html}
</div>
"""

    # ── Microbial Findings Card ───────────────────────────────────────────
    if entries:
        rows = ""
        for e in entries:
            badge = _badge(e.grade)
            bar = _abund_bar(e.abundance * 100, e.grade)
            contam = '<span class="contam-flag">⚠ contaminant</span>' if e.contaminant_risk else ""
            ci_str = f"{e.ci_lower*100:.1f}–{e.ci_upper*100:.1f}%"
            rows += f"""<tr>
  <td><em>{e.organism}</em> {contam}</td>
  <td>{bar}</td>
  <td style="font-size:12px;color:#777">{ci_str}</td>
  <td>{e.read_count:,}</td>
  <td>{badge}</td>
</tr>"""
        findings_html = f"""<table>
  <thead><tr>
    <th>Organism</th><th>Abundance</th><th>95% CI</th><th>Reads</th><th>Grade</th>
  </tr></thead>
  <tbody>{rows}</tbody>
</table>"""
    else:
        findings_html = '<p class="no-data">No organisms detected above threshold.</p>'

    html += f"""
<div class="card">
  <div class="card-title">Microbial Findings</div>
  <div style="margin-bottom:14px">
    <div class="stat-grid">
      <div class="stat-box"><div class="num">{len(entries)}</div><div class="desc">Organisms detected</div></div>
      <div class="stat-box"><div class="num">{em_result.n_reads:,}</div><div class="desc">Classified reads</div></div>
      <div class="stat-box"><div class="num">{sum(1 for e in entries if e.grade == EvidenceGrade.A)}</div><div class="desc">Grade A detections</div></div>
    </div>
  </div>
  {findings_html}
  <div style="margin-top:12px;font-size:12px;color:#888">
    <strong>Grade key:</strong>
    <span class="grade-badge" style="background:#2e7d32;color:white;margin:0 4px">A</span> High confidence &nbsp;
    <span class="grade-badge" style="background:#e65100;color:white;margin:0 4px">B</span> Moderate confidence &nbsp;
    <span class="grade-badge" style="background:#757575;color:white;margin:0 4px">C</span> Low confidence / contaminant &nbsp;
    <span class="grade-badge" style="background:#c62828;color:white;margin:0 4px">X</span> Insufficient reads
  </div>
</div>
"""

    # ── AMR Card ──────────────────────────────────────────────────────────
    if amr_hits:
        amr_rows = "".join(
            f"<tr><td>{h.gene}</td><td>{h.drug_class}</td>"
            f"<td>{h.identity_pct:.1f}%</td><td>{h.coverage_pct:.1f}%</td>"
            f"<td><em>{h.organism_match}</em></td><td>{h.database}</td></tr>"
            for h in amr_hits
        )
        amr_html = f"""<table>
  <thead><tr>
    <th>Gene</th><th>Drug Class</th><th>Identity</th><th>Coverage</th><th>Organism</th><th>Database</th>
  </tr></thead>
  <tbody>{amr_rows}</tbody>
</table>"""
    else:
        amr_html = '<p class="no-data">No AMR genes detected (ABRicate not run or no hits above threshold).</p>'

    html += f"""
<div class="card">
  <div class="card-title">Antimicrobial Resistance Genes</div>
  {amr_html}
</div>
"""

    # ── Virulence Card (VFDB) ─────────────────────────────────────────────
    if virulence_hits:
        vir_rows = "".join(
            f"<tr><td>{h.gene}</td><td>{h.factor}</td>"
            f"<td>{h.identity_pct:.1f}%</td><td>{h.coverage_pct:.1f}%</td>"
            f"<td><em>{h.organism_match}</em></td><td>{h.database}</td></tr>"
            for h in virulence_hits
        )
        vir_html = f"""<table>
  <thead><tr>
    <th>Gene</th><th>Virulence Factor</th><th>Identity</th><th>Coverage</th><th>Organism</th><th>Database</th>
  </tr></thead>
  <tbody>{vir_rows}</tbody>
</table>"""
    else:
        vir_html = '<p class="no-data">No virulence factors detected (VFDB not run or no hits above threshold).</p>'

    html += f"""
<div class="card">
  <div class="card-title">Virulence Factors (VFDB)</div>
  {vir_html}
</div>
"""

    # ── Footer ────────────────────────────────────────────────────────────
    html += f"""
<footer>PathogenIQ v0.1 &nbsp;|&nbsp; Generated {now} &nbsp;|&nbsp; For research use only</footer>
</div>
</body>
</html>
"""

    html_path = out / "pathogeniq_report.html"
    html_path.write_text(html, encoding="utf-8")
    return html_path
