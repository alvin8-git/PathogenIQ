# PathogenIQ Plan 2 — Clinical Interpretation Engine

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Extend PathogenIQ with specimen-aware contaminant suppression, AMR gene detection via ABRicate+CARD, and a PDF clinical report.

**Architecture:** Three independent modules plug into the existing pipeline: (1) `contaminants.py` provides a registry of known per-specimen contaminants and downgrades Grade B findings to Grade C when the organism is a known environmental contaminant; (2) `amr.py` runs ABRicate against the CARD database on non-human reads and links AMR hits to detected organisms by genus/species substring match; (3) `pdf_report.py` uses ReportLab Platypus to render a formatted clinical PDF from the existing JSON report payload. All three are wired into `cli.py` and `report.py`.

**Tech Stack:** Python 3.11, ReportLab 4.x, ABRicate ≥ 1.0 (CARD DB), existing pathogeniq modules (config, report, em, cli).

---

## File Map

| File | Action | Responsibility |
|---|---|---|
| `pathogeniq/contaminants.py` | Create | Contaminant prior registry; `flag_contaminants()` |
| `pathogeniq/amr.py` | Create | ABRicate wrapper; `AMRHit`; `run_amr_screen()` |
| `pathogeniq/pdf_report.py` | Create | ReportLab PDF renderer; `write_pdf_report()` |
| `pathogeniq/report.py` | Modify | Add `contaminant_risk` field; add `amr_genes` to JSON; accept AMR hits |
| `pathogeniq/cli.py` | Modify | Wire contaminant flag, AMR step, PDF step; add `--no-pdf` flag |
| `pathogeniq/config.py` | Modify | Add `amr_db: str` field (default `"card"`) |
| `pyproject.toml` | Modify | Add `reportlab>=4.0` dependency |
| `tests/test_contaminants.py` | Create | Unit tests for contaminant flagging logic |
| `tests/test_amr.py` | Create | Unit tests for AMRHit parsing and run_amr_screen |
| `tests/test_pdf_report.py` | Create | Smoke-test PDF is written and non-empty |
| `tests/test_report.py` | Modify | Add tests for contaminant_risk field and amr_genes in JSON |

---

## Task 1: Contaminant Prior Registry

**Files:**
- Create: `pathogeniq/contaminants.py`
- Create: `tests/test_contaminants.py`

Background: Some organisms are clinically relevant pathogens in one specimen type but routine contaminants in another. *Cutibacterium acnes* is always a skin contaminant in blood cultures; *Streptococcus salivarius* is oral flora that frequently contaminates BAL. This module provides a registry of such organisms per specimen type and a function to flag report entries.

- [ ] **Step 1: Write the failing tests**

```python
# tests/test_contaminants.py
import pytest
from pathogeniq.contaminants import flag_contaminants, CONTAMINANT_PRIORS
from pathogeniq.config import SpecimenType
from pathogeniq.report import ReportEntry, EvidenceGrade


def _entry(organism: str, specimen: SpecimenType, read_count: int = 5,
           abundance: float = 0.05, ci_lower: float = 0.01, ci_upper: float = 0.10) -> ReportEntry:
    return ReportEntry(
        organism=organism,
        abundance=abundance,
        ci_lower=ci_lower,
        ci_upper=ci_upper,
        read_count=read_count,
        specimen_type=specimen,
        contaminant_risk=False,
    )


def test_cutibacterium_flagged_in_blood():
    entries = [_entry("Cutibacterium acnes", SpecimenType.BLOOD)]
    result = flag_contaminants(entries)
    assert result[0].contaminant_risk is True


def test_cutibacterium_not_flagged_in_tissue():
    entries = [_entry("Cutibacterium acnes", SpecimenType.TISSUE)]
    result = flag_contaminants(entries)
    assert result[0].contaminant_risk is False


def test_streptococcus_salivarius_flagged_in_bal():
    entries = [_entry("Streptococcus salivarius", SpecimenType.BAL)]
    result = flag_contaminants(entries)
    assert result[0].contaminant_risk is True


def test_staphylococcus_aureus_not_flagged_blood():
    # S. aureus is a true pathogen even in blood — must NOT be suppressed
    entries = [_entry("Staphylococcus aureus", SpecimenType.BLOOD)]
    result = flag_contaminants(entries)
    assert result[0].contaminant_risk is False


def test_staph_epidermidis_flagged_in_csf():
    entries = [_entry("Staphylococcus epidermidis", SpecimenType.CSF)]
    result = flag_contaminants(entries)
    assert result[0].contaminant_risk is True


def test_contaminant_registry_has_all_specimen_types():
    for stype in SpecimenType:
        assert stype in CONTAMINANT_PRIORS


def test_grade_b_contaminant_downgrades_to_c():
    # Grade B entry (enough reads, wide CI) that is a contaminant → reported as C
    entry = _entry("Cutibacterium acnes", SpecimenType.BLOOD,
                   read_count=10, ci_lower=0.00, ci_upper=0.20)
    flagged = flag_contaminants([entry])
    # Grade should fall through to C because contaminant_risk=True makes grade() return C
    assert flagged[0].grade == EvidenceGrade.C
```

- [ ] **Step 2: Run tests to confirm they fail**

```bash
pytest tests/test_contaminants.py -v 2>&1 | head -30
```

Expected: `ImportError` or `ModuleNotFoundError` for `pathogeniq.contaminants`.

- [ ] **Step 3: Create `pathogeniq/contaminants.py`**

```python
from __future__ import annotations

from pathogeniq.config import SpecimenType
from pathogeniq.report import ReportEntry


# Organisms that are routine contaminants in specific specimen types.
# Matching is substring (case-insensitive) against ReportEntry.organism.
CONTAMINANT_PRIORS: dict[SpecimenType, list[str]] = {
    SpecimenType.BLOOD: [
        "Cutibacterium acnes",
        "Staphylococcus epidermidis",
        "Staphylococcus capitis",
        "Staphylococcus hominis",
        "Bacillus cereus",
        "Bacillus subtilis",
        "Micrococcus luteus",
        "Corynebacterium striatum",
        "Corynebacterium jeikeium",
    ],
    SpecimenType.CSF: [
        "Cutibacterium acnes",
        "Staphylococcus epidermidis",
        "Staphylococcus capitis",
        "Corynebacterium striatum",
    ],
    SpecimenType.BAL: [
        "Streptococcus salivarius",
        "Streptococcus mitis",
        "Streptococcus oralis",
        "Prevotella melaninogenica",
        "Veillonella parvula",
        "Rothia mucilaginosa",
        "Neisseria sicca",
        "Fusobacterium nucleatum",
    ],
    SpecimenType.TISSUE: [],
}


def flag_contaminants(entries: list[ReportEntry]) -> list[ReportEntry]:
    """Set contaminant_risk=True on entries matching known contaminants for their specimen type."""
    for entry in entries:
        known = CONTAMINANT_PRIORS.get(entry.specimen_type, [])
        entry.contaminant_risk = any(
            c.lower() in entry.organism.lower() for c in known
        )
    return entries
```

- [ ] **Step 4: Update `ReportEntry` in `pathogeniq/report.py` to add `contaminant_risk` field and downgrade logic**

Replace the `ReportEntry` dataclass and `grade` property (lines 31–47):

```python
@dataclass
class ReportEntry:
    organism: str
    abundance: float
    ci_lower: float
    ci_upper: float
    read_count: int
    specimen_type: SpecimenType
    contaminant_risk: bool = False

    @property
    def grade(self) -> EvidenceGrade:
        min_reads = _MIN_READS.get(self.specimen_type, 5)
        ci_width = self.ci_upper - self.ci_lower
        if self.read_count >= min_reads and ci_width <= _MAX_CI_WIDTH_A and not self.contaminant_risk:
            return EvidenceGrade.A
        if self.read_count >= min_reads and not self.contaminant_risk:
            return EvidenceGrade.B
        if self.read_count >= min_reads:
            return EvidenceGrade.C
        return EvidenceGrade.X
```

- [ ] **Step 5: Run tests to confirm they pass**

```bash
pytest tests/test_contaminants.py -v
```

Expected: all 7 tests PASS.

- [ ] **Step 6: Run existing report tests to confirm no regression**

```bash
pytest tests/test_report.py -v
```

Expected: all PASS (contaminant_risk defaults to False).

- [ ] **Step 7: Commit**

```bash
git add pathogeniq/contaminants.py pathogeniq/report.py tests/test_contaminants.py
git commit -m "feat: add specimen-aware contaminant prior registry and grade downgrade"
```

---

## Task 2: Contaminant Field in Report Output

**Files:**
- Modify: `pathogeniq/report.py`
- Modify: `tests/test_report.py`

The JSON and TSV reports need to expose `contaminant_risk` so downstream consumers (dashboard, EHR) can filter or annotate.

- [ ] **Step 1: Add failing test for contaminant_risk in JSON output**

Open `tests/test_report.py` and add:

```python
def test_report_json_includes_contaminant_risk(tmp_path):
    from pathogeniq.contaminants import flag_contaminants

    cfg = PipelineConfig(
        input_fastq=tmp_path / "s.fq.gz",
        read_type=ReadType.SHORT,
        specimen_type=SpecimenType.BLOOD,
        output_dir=tmp_path,
        db_tier1=tmp_path / "db.sbt.zip",
        host_reference=tmp_path / "ref.fa",
    )
    n_orgs = 2
    em = EMResult(
        abundances=np.array([0.8, 0.2]),
        n_reads=100,
        converged=True,
        iterations=5,
    )
    entries_raw = [
        ReportEntry("Staphylococcus aureus", 0.8, 0.7, 0.9, 80, SpecimenType.BLOOD),
        ReportEntry("Cutibacterium acnes",   0.2, 0.1, 0.3, 20, SpecimenType.BLOOD),
    ]
    entries = flag_contaminants(entries_raw)

    write_report(cfg, ["Staphylococcus aureus", "Cutibacterium acnes"], em,
                 np.array([0.7, 0.1]), np.array([0.9, 0.3]),
                 amr_hits=[], entries_override=entries)

    with open(tmp_path / "report" / "pathogeniq_report.json") as f:
        payload = json.load(f)

    findings = {f["organism"]: f for f in payload["findings"]}
    assert findings["Staphylococcus aureus"]["contaminant_risk"] is False
    assert findings["Cutibacterium acnes"]["contaminant_risk"] is True
```

- [ ] **Step 2: Run to confirm failure**

```bash
pytest tests/test_report.py::test_report_json_includes_contaminant_risk -v
```

Expected: FAIL — `write_report` does not accept `amr_hits` or `entries_override`.

- [ ] **Step 3: Update `write_report` signature and output in `pathogeniq/report.py`**

Replace the full `write_report` function:

```python
def write_report(
    cfg: PipelineConfig,
    organism_names: list[str],
    em_result: EMResult,
    ci_lower: np.ndarray,
    ci_upper: np.ndarray,
    amr_hits: list | None = None,
    entries_override: list[ReportEntry] | None = None,
) -> Path:
    out = cfg.output_dir / "report"
    out.mkdir(parents=True, exist_ok=True)

    if entries_override is not None:
        entries = entries_override
    else:
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
            for i, name in enumerate(organism_names)
        ]
    entries.sort(key=lambda e: e.abundance, reverse=True)

    amr_by_org: dict[str, list[dict]] = {}
    for hit in (amr_hits or []):
        amr_by_org.setdefault(hit.organism_match, []).append({
            "gene": hit.gene,
            "drug_class": hit.drug_class,
            "identity_pct": hit.identity_pct,
            "coverage_pct": hit.coverage_pct,
        })

    tsv_path = out / "pathogeniq_report.tsv"
    with open(tsv_path, "w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=["organism", "abundance_pct", "ci_lower_pct", "ci_upper_pct",
                        "read_count", "grade", "contaminant_risk"],
            delimiter="\t",
        )
        writer.writeheader()
        for e in entries:
            writer.writerow({
                "organism": e.organism,
                "abundance_pct": f"{e.abundance * 100:.2f}",
                "ci_lower_pct": f"{e.ci_lower * 100:.2f}",
                "ci_upper_pct": f"{e.ci_upper * 100:.2f}",
                "read_count": e.read_count,
                "grade": e.grade.value,
                "contaminant_risk": e.contaminant_risk,
            })

    json_path = out / "pathogeniq_report.json"
    payload = {
        "sample": str(cfg.input_fastq.name),
        "specimen_type": cfg.specimen_type.value,
        "read_type": cfg.read_type.value,
        "total_classified_reads": em_result.n_reads,
        "findings": [
            {
                "organism": e.organism,
                "abundance_pct": round(e.abundance * 100, 2),
                "ci_lower_pct": round(e.ci_lower * 100, 2),
                "ci_upper_pct": round(e.ci_upper * 100, 2),
                "read_count": e.read_count,
                "grade": e.grade.value,
                "contaminant_risk": e.contaminant_risk,
                "amr_genes": amr_by_org.get(e.organism, []),
            }
            for e in entries
        ],
    }
    with open(json_path, "w") as f:
        json.dump(payload, f, indent=2)

    return out
```

- [ ] **Step 4: Run all report tests**

```bash
pytest tests/test_report.py -v
```

Expected: all PASS.

- [ ] **Step 5: Commit**

```bash
git add pathogeniq/report.py tests/test_report.py
git commit -m "feat: add contaminant_risk and amr_genes fields to report output"
```

---

## Task 3: AMR Detection Module

**Files:**
- Create: `pathogeniq/amr.py`
- Create: `tests/test_amr.py`

Background: ABRicate scans a FASTA file against curated AMR databases (CARD, ResFinder, VFDB). It outputs a TSV with columns: `#FILE`, `SEQUENCE`, `START`, `END`, `STRAND`, `GENE`, `COVERAGE`, `COVERAGE_MAP`, `GAPS`, `%COVERAGE`, `%IDENTITY`, `DATABASE`, `ACCESSION`, `PRODUCT`, `RESISTANCE`. We convert non-human reads from FASTQ to FASTA and run ABRicate on it. We then match AMR hits to detected organisms by substring matching the SEQUENCE header against organism names.

- [ ] **Step 1: Write failing tests**

```python
# tests/test_amr.py
import textwrap
from pathlib import Path
from unittest.mock import patch, MagicMock
import subprocess

import pytest

from pathogeniq.amr import AMRHit, _parse_abricate_tsv, run_amr_screen
from pathogeniq.config import PipelineConfig, ReadType, SpecimenType


FAKE_TSV = textwrap.dedent("""\
    #FILE\tSEQUENCE\tSTART\tEND\tSTRAND\tGENE\tCOVERAGE\tCOVERAGE_MAP\tGAPS\t%COVERAGE\t%IDENTITY\tDATABASE\tACCESSION\tPRODUCT\tRESISTANCE
    reads.fa\tStaphylococcus_aureus_read1\t1\t801\t+\tmecA\t1-801/801\t===============\t0\t100.00\t99.75\tcard\tAR3\tPBP2a\tbeta-lactam
    reads.fa\tKlebsiella_pneumoniae_read5\t1\t588\t+\tblaCTX-M-15\t1-588/588\t===============\t0\t100.00\t98.10\tcard\tAR4\tCTX-M-15\tbeta-lactam;cephalosporin
""")


def test_parse_abricate_tsv_returns_hits():
    hits = _parse_abricate_tsv(FAKE_TSV, organism_names=["Staphylococcus aureus", "Klebsiella pneumoniae"])
    assert len(hits) == 2


def test_parse_abricate_tsv_meca_hit():
    hits = _parse_abricate_tsv(FAKE_TSV, organism_names=["Staphylococcus aureus", "Klebsiella pneumoniae"])
    meca = next(h for h in hits if h.gene == "mecA")
    assert meca.drug_class == "beta-lactam"
    assert meca.identity_pct == pytest.approx(99.75)
    assert meca.coverage_pct == pytest.approx(100.0)
    assert meca.organism_match == "Staphylococcus aureus"


def test_parse_abricate_tsv_blactxm_hit():
    hits = _parse_abricate_tsv(FAKE_TSV, organism_names=["Staphylococcus aureus", "Klebsiella pneumoniae"])
    ctx = next(h for h in hits if h.gene == "blaCTX-M-15")
    assert ctx.organism_match == "Klebsiella pneumoniae"


def test_parse_abricate_no_match_organism():
    # SEQUENCE header doesn't match any organism name → organism_match is "unknown"
    tsv = FAKE_TSV.replace("Staphylococcus_aureus_read1", "some_random_contig")
    hits = _parse_abricate_tsv(tsv, organism_names=["Staphylococcus aureus"])
    meca = next((h for h in hits if h.gene == "mecA"), None)
    assert meca is not None
    assert meca.organism_match == "unknown"


def test_run_amr_screen_skips_when_abricate_missing(tmp_path):
    cfg = PipelineConfig(
        input_fastq=tmp_path / "s.fq.gz",
        read_type=ReadType.SHORT,
        specimen_type=SpecimenType.BLOOD,
        output_dir=tmp_path,
        db_tier1=tmp_path / "db.sbt.zip",
        host_reference=tmp_path / "ref.fa",
    )
    reads = tmp_path / "nonhuman.fq.gz"
    reads.touch()
    with patch("shutil.which", return_value=None):
        hits = run_amr_screen(cfg, reads, organism_names=["Staphylococcus aureus"])
    assert hits == []


def test_run_amr_screen_returns_hits(tmp_path):
    cfg = PipelineConfig(
        input_fastq=tmp_path / "s.fq.gz",
        read_type=ReadType.SHORT,
        specimen_type=SpecimenType.BLOOD,
        output_dir=tmp_path,
        db_tier1=tmp_path / "db.sbt.zip",
        host_reference=tmp_path / "ref.fa",
    )
    reads = tmp_path / "nonhuman.fq.gz"
    reads.touch()

    mock_result = MagicMock()
    mock_result.returncode = 0
    mock_result.stdout = FAKE_TSV

    with patch("shutil.which", return_value="/usr/bin/abricate"), \
         patch("subprocess.run", return_value=mock_result):
        hits = run_amr_screen(cfg, reads, organism_names=["Staphylococcus aureus", "Klebsiella pneumoniae"])

    assert len(hits) == 2
```

- [ ] **Step 2: Run to confirm failure**

```bash
pytest tests/test_amr.py -v 2>&1 | head -20
```

Expected: `ImportError` — `pathogeniq.amr` does not exist.

- [ ] **Step 3: Create `pathogeniq/amr.py`**

```python
from __future__ import annotations

import csv
import io
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path

from .config import PipelineConfig


@dataclass
class AMRHit:
    gene: str
    drug_class: str
    identity_pct: float
    coverage_pct: float
    organism_match: str  # matched organism name or "unknown"
    database: str


def _parse_abricate_tsv(tsv_text: str, organism_names: list[str]) -> list[AMRHit]:
    hits: list[AMRHit] = []
    reader = csv.DictReader(io.StringIO(tsv_text), delimiter="\t")
    for row in reader:
        if row.get("#FILE", "").startswith("#"):
            continue
        sequence = row.get("SEQUENCE", "")
        # Match sequence header to a known organism by substring (underscore-normalised)
        seq_norm = sequence.replace("_", " ").lower()
        matched = "unknown"
        for org in organism_names:
            if org.lower() in seq_norm:
                matched = org
                break
        hits.append(AMRHit(
            gene=row.get("GENE", ""),
            drug_class=row.get("RESISTANCE", "unknown"),
            identity_pct=float(row.get("%IDENTITY", 0)),
            coverage_pct=float(row.get("%COVERAGE", 0)),
            organism_match=matched,
            database=row.get("DATABASE", ""),
        ))
    return hits


def run_amr_screen(
    cfg: PipelineConfig,
    reads_path: Path,
    organism_names: list[str],
    db: str = "card",
    min_identity: float = 90.0,
    min_coverage: float = 80.0,
) -> list[AMRHit]:
    """Run ABRicate on non-human reads. Returns empty list if ABRicate not on PATH."""
    if not shutil.which("abricate"):
        return []

    # Convert FASTQ → FASTA in a temp file for abricate
    fasta_path = cfg.output_dir / "amr_reads.fa"
    _fastq_to_fasta(reads_path, fasta_path)

    result = subprocess.run(
        [
            "abricate",
            "--db", db,
            "--minid", str(min_identity),
            "--mincov", str(min_coverage),
            "--quiet",
            str(fasta_path),
        ],
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        return []

    return _parse_abricate_tsv(result.stdout, organism_names)


def _fastq_to_fasta(fastq_path: Path, fasta_path: Path) -> None:
    """Convert (optionally gzipped) FASTQ to FASTA."""
    import gzip

    open_fn = gzip.open if str(fastq_path).endswith(".gz") else open
    with open_fn(fastq_path, "rt") as fq, open(fasta_path, "w") as fa:
        while True:
            header = fq.readline()
            if not header:
                break
            seq = fq.readline()
            fq.readline()  # +
            fq.readline()  # quality
            fa.write(f">{header[1:]}{seq}")
```

- [ ] **Step 4: Run tests**

```bash
pytest tests/test_amr.py -v
```

Expected: all 6 tests PASS.

- [ ] **Step 5: Commit**

```bash
git add pathogeniq/amr.py tests/test_amr.py
git commit -m "feat: add AMR detection module wrapping ABRicate+CARD"
```

---

## Task 4: PDF Report Generator

**Files:**
- Create: `pathogeniq/pdf_report.py`
- Create: `tests/test_pdf_report.py`
- Modify: `pyproject.toml`

The PDF is a single-page clinical summary: header (sample, specimen, date), findings table (colour-coded by grade), AMR summary table, and a footer note. We use ReportLab's Platypus (high-level layout engine) — no external system dependencies.

Grade colour coding:
- Grade A → green (`#2e7d32`)
- Grade B → orange (`#e65100`)
- Grade C → grey (`#757575`)
- Grade X → red (`#c62828`)

- [ ] **Step 1: Add reportlab to pyproject.toml**

Open `pyproject.toml` and update `dependencies`:

```toml
dependencies = [
    "click>=8.1",
    "numpy>=1.26",
    "scipy>=1.12",
    "pysam>=0.22",
    "pandas>=2.2",
    "reportlab>=4.0",
]
```

- [ ] **Step 2: Install the new dependency**

```bash
pip install -e .
```

Expected: `Successfully installed reportlab-4.x.x` (or "already satisfied").

- [ ] **Step 3: Write the failing tests**

```python
# tests/test_pdf_report.py
from pathlib import Path
import pytest

from pathogeniq.config import PipelineConfig, ReadType, SpecimenType
from pathogeniq.report import ReportEntry
from pathogeniq.amr import AMRHit
from pathogeniq.pdf_report import write_pdf_report


def _cfg(tmp_path: Path) -> PipelineConfig:
    return PipelineConfig(
        input_fastq=tmp_path / "sample.fq.gz",
        read_type=ReadType.SHORT,
        specimen_type=SpecimenType.BLOOD,
        output_dir=tmp_path,
        db_tier1=tmp_path / "db.sbt.zip",
        host_reference=tmp_path / "ref.fa",
    )


def _entries() -> list[ReportEntry]:
    return [
        ReportEntry("Staphylococcus aureus", 0.80, 0.72, 0.88, 80, SpecimenType.BLOOD),
        ReportEntry("Cutibacterium acnes",   0.15, 0.08, 0.22, 15, SpecimenType.BLOOD, contaminant_risk=True),
        ReportEntry("Escherichia coli",      0.05, 0.01, 0.09,  5, SpecimenType.BLOOD),
    ]


def test_pdf_report_creates_file(tmp_path):
    cfg = _cfg(tmp_path)
    pdf_path = write_pdf_report(cfg, _entries(), amr_hits=[])
    assert pdf_path.exists()
    assert pdf_path.suffix == ".pdf"
    assert pdf_path.stat().st_size > 1000  # non-trivial PDF


def test_pdf_report_with_amr_hits(tmp_path):
    cfg = _cfg(tmp_path)
    hits = [
        AMRHit("mecA", "beta-lactam", 99.5, 100.0, "Staphylococcus aureus", "card"),
        AMRHit("vanA", "glycopeptide", 98.1,  95.0, "unknown", "card"),
    ]
    pdf_path = write_pdf_report(cfg, _entries(), amr_hits=hits)
    assert pdf_path.exists()
    assert pdf_path.stat().st_size > 1000


def test_pdf_report_output_path(tmp_path):
    cfg = _cfg(tmp_path)
    pdf_path = write_pdf_report(cfg, _entries(), amr_hits=[])
    assert pdf_path == tmp_path / "report" / "pathogeniq_report.pdf"
```

- [ ] **Step 4: Run to confirm failure**

```bash
pytest tests/test_pdf_report.py -v 2>&1 | head -20
```

Expected: `ImportError` — `pathogeniq.pdf_report` does not exist.

- [ ] **Step 5: Create `pathogeniq/pdf_report.py`**

```python
from __future__ import annotations

from datetime import datetime
from pathlib import Path

from reportlab.lib import colors
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.lib.units import cm
from reportlab.platypus import (
    Paragraph, SimpleDocTemplate, Spacer, Table, TableStyle
)

from .amr import AMRHit
from .config import PipelineConfig
from .report import EvidenceGrade, ReportEntry

_GRADE_COLORS = {
    EvidenceGrade.A: colors.HexColor("#2e7d32"),
    EvidenceGrade.B: colors.HexColor("#e65100"),
    EvidenceGrade.C: colors.HexColor("#757575"),
    EvidenceGrade.X: colors.HexColor("#c62828"),
}


def write_pdf_report(
    cfg: PipelineConfig,
    entries: list[ReportEntry],
    amr_hits: list[AMRHit],
) -> Path:
    out_dir = cfg.output_dir / "report"
    out_dir.mkdir(parents=True, exist_ok=True)
    pdf_path = out_dir / "pathogeniq_report.pdf"

    styles = getSampleStyleSheet()
    doc = SimpleDocTemplate(
        str(pdf_path),
        pagesize=A4,
        leftMargin=2 * cm,
        rightMargin=2 * cm,
        topMargin=2 * cm,
        bottomMargin=2 * cm,
    )
    story = []

    # ── Header ────────────────────────────────────────────────────────────────
    story.append(Paragraph("<b>PathogenIQ Clinical Metagenomics Report</b>", styles["Title"]))
    story.append(Spacer(1, 0.3 * cm))
    story.append(Paragraph(
        f"Sample: <b>{cfg.input_fastq.name}</b> &nbsp;&nbsp; "
        f"Specimen: <b>{cfg.specimen_type.value.upper()}</b> &nbsp;&nbsp; "
        f"Read type: <b>{cfg.read_type.value.upper()}</b> &nbsp;&nbsp; "
        f"Generated: <b>{datetime.now().strftime('%Y-%m-%d %H:%M')}</b>",
        styles["Normal"],
    ))
    story.append(Spacer(1, 0.5 * cm))

    # ── Findings table ────────────────────────────────────────────────────────
    story.append(Paragraph("<b>Findings</b>", styles["Heading2"]))
    story.append(Spacer(1, 0.2 * cm))

    header_row = ["Organism", "Abundance %", "95% CI", "Reads", "Grade", "Contaminant?"]
    data = [header_row]
    for e in entries:
        ci_str = f"{e.ci_lower * 100:.1f}–{e.ci_upper * 100:.1f}%"
        data.append([
            e.organism,
            f"{e.abundance * 100:.2f}%",
            ci_str,
            str(e.read_count),
            e.grade.value,
            "Yes" if e.contaminant_risk else "No",
        ])

    col_widths = [6.5 * cm, 2.5 * cm, 3 * cm, 1.5 * cm, 1.5 * cm, 2.5 * cm]
    table = Table(data, colWidths=col_widths, repeatRows=1)
    ts = TableStyle([
        ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#37474f")),
        ("TEXTCOLOR",  (0, 0), (-1, 0), colors.white),
        ("FONTNAME",   (0, 0), (-1, 0), "Helvetica-Bold"),
        ("FONTSIZE",   (0, 0), (-1, -1), 9),
        ("ALIGN",      (1, 0), (-1, -1), "CENTER"),
        ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, colors.HexColor("#f5f5f5")]),
        ("GRID", (0, 0), (-1, -1), 0.25, colors.HexColor("#cfd8dc")),
        ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
        ("TOPPADDING", (0, 0), (-1, -1), 4),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 4),
    ])
    # Colour-code Grade column
    for i, e in enumerate(entries, start=1):
        grade_color = _GRADE_COLORS.get(e.grade, colors.black)
        ts.add("TEXTCOLOR", (4, i), (4, i), grade_color)
        ts.add("FONTNAME",  (4, i), (4, i), "Helvetica-Bold")
        if e.contaminant_risk:
            ts.add("TEXTCOLOR", (5, i), (5, i), colors.HexColor("#e65100"))
    table.setStyle(ts)
    story.append(table)
    story.append(Spacer(1, 0.5 * cm))

    # ── AMR table ─────────────────────────────────────────────────────────────
    if amr_hits:
        story.append(Paragraph("<b>Antimicrobial Resistance Genes Detected</b>", styles["Heading2"]))
        story.append(Spacer(1, 0.2 * cm))
        amr_header = ["Gene", "Drug Class", "Identity %", "Coverage %", "Organism"]
        amr_data = [amr_header] + [
            [h.gene, h.drug_class, f"{h.identity_pct:.1f}", f"{h.coverage_pct:.1f}", h.organism_match]
            for h in amr_hits
        ]
        amr_col_widths = [3 * cm, 4 * cm, 2.5 * cm, 2.5 * cm, 5.5 * cm]
        amr_table = Table(amr_data, colWidths=amr_col_widths, repeatRows=1)
        amr_table.setStyle(TableStyle([
            ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#b71c1c")),
            ("TEXTCOLOR",  (0, 0), (-1, 0), colors.white),
            ("FONTNAME",   (0, 0), (-1, 0), "Helvetica-Bold"),
            ("FONTSIZE",   (0, 0), (-1, -1), 9),
            ("ALIGN",      (2, 0), (3, -1), "CENTER"),
            ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, colors.HexColor("#ffebee")]),
            ("GRID", (0, 0), (-1, -1), 0.25, colors.HexColor("#cfd8dc")),
            ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
            ("TOPPADDING", (0, 0), (-1, -1), 4),
            ("BOTTOMPADDING", (0, 0), (-1, -1), 4),
        ]))
        story.append(amr_table)
        story.append(Spacer(1, 0.5 * cm))

    # ── Grade legend ──────────────────────────────────────────────────────────
    legend_text = (
        "<b>Grade key:</b> "
        "<font color='#2e7d32'>A = High confidence</font> · "
        "<font color='#e65100'>B = Moderate confidence</font> · "
        "<font color='#757575'>C = Low confidence / contaminant risk</font> · "
        "<font color='#c62828'>X = Insufficient reads</font>"
    )
    story.append(Paragraph(legend_text, styles["Normal"]))
    story.append(Spacer(1, 0.3 * cm))
    story.append(Paragraph(
        "<i>This report is for research use only. Clinical decisions should not be based solely on this output.</i>",
        styles["Normal"],
    ))

    doc.build(story)
    return pdf_path
```

- [ ] **Step 6: Run PDF tests**

```bash
pytest tests/test_pdf_report.py -v
```

Expected: all 3 tests PASS.

- [ ] **Step 7: Commit**

```bash
git add pathogeniq/pdf_report.py tests/test_pdf_report.py pyproject.toml
git commit -m "feat: add ReportLab PDF clinical report generator"
```

---

## Task 5: Wire Everything into the CLI

**Files:**
- Modify: `pathogeniq/cli.py`
- Modify: `pathogeniq/config.py`
- Modify: `tests/test_cli.py`

The CLI needs to: (1) run contaminant flagging after EM, (2) run AMR screening in parallel, (3) pass AMR hits to `write_report`, (4) call `write_pdf_report`, (5) support `--no-pdf` and `--amr-db` flags.

- [ ] **Step 1: Add `amr_db` to config**

Open `pathogeniq/config.py` and add `amr_db` field:

```python
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path


class ReadType(str, Enum):
    SHORT = "short"
    LONG = "long"


class SpecimenType(str, Enum):
    BLOOD = "blood"
    CSF = "csf"
    BAL = "bal"
    TISSUE = "tissue"


@dataclass
class PipelineConfig:
    input_fastq: Path
    read_type: ReadType
    specimen_type: SpecimenType
    output_dir: Path
    db_tier1: Path
    host_reference: Path
    threads: int = 8
    sketch_threshold: float = 0.01
    n_bootstrap: int = 100
    amr_db: str = "card"
```

- [ ] **Step 2: Write failing CLI test for contaminant + AMR + PDF**

Open `tests/test_cli.py` and add:

```python
def test_cli_run_produces_pdf(tmp_path):
    """End-to-end CLI smoke test: PDF is produced when pipeline succeeds."""
    from unittest.mock import patch, MagicMock
    import numpy as np
    from click.testing import CliRunner
    from pathogeniq.cli import cli
    from pathogeniq.qc import QCMetrics
    from pathogeniq.host_remove import HostRemovalMetrics
    from pathogeniq.sketch import SketchHit
    from pathogeniq.align import AlignmentResult
    from pathogeniq.em import EMResult

    fq = tmp_path / "s.fq.gz"
    fq.touch()
    db = tmp_path / "db.sbt.zip"
    db.touch()
    ref = tmp_path / "ref.fa"
    ref.touch()

    fake_em = EMResult(abundances=np.array([1.0]), n_reads=10, converged=True, iterations=3)
    fake_align = AlignmentResult(
        alignment_matrix=np.ones((10, 1)),
        organism_names=["Staphylococcus aureus"],
        accessions=["GCF_000013425.1"],
    )

    runner = CliRunner()
    with patch("pathogeniq.cli.run_qc", return_value=(tmp_path / "filt.fq.gz", QCMetrics(100, 90, 0.1))), \
         patch("pathogeniq.cli.run_host_removal", return_value=(tmp_path / "nonhuman.fq.gz", HostRemovalMetrics(90, 10))), \
         patch("pathogeniq.cli.run_sketch_screen", return_value=[SketchHit("Staphylococcus aureus", "GCF_000013425.1", 0.05)]), \
         patch("pathogeniq.cli.run_targeted_alignment", return_value=fake_align), \
         patch("pathogeniq.cli.em_abundance", return_value=fake_em), \
         patch("pathogeniq.cli.bootstrap_ci", return_value=(np.array([0.9]), np.array([1.0]))), \
         patch("pathogeniq.cli.run_amr_screen", return_value=[]), \
         patch("pathogeniq.cli.write_pdf_report") as mock_pdf:
        mock_pdf.return_value = tmp_path / "report" / "pathogeniq_report.pdf"
        result = runner.invoke(cli, [
            "run",
            "--input", str(fq),
            "--output", str(tmp_path / "out"),
            "--db", str(db),
            "--host-ref", str(ref),
            "--specimen", "blood",
            "--read-type", "short",
        ])

    assert result.exit_code == 0, result.output
    mock_pdf.assert_called_once()


def test_cli_run_no_pdf_flag(tmp_path):
    """--no-pdf suppresses PDF generation."""
    from unittest.mock import patch, MagicMock
    import numpy as np
    from click.testing import CliRunner
    from pathogeniq.cli import cli
    from pathogeniq.qc import QCMetrics
    from pathogeniq.host_remove import HostRemovalMetrics
    from pathogeniq.sketch import SketchHit
    from pathogeniq.align import AlignmentResult
    from pathogeniq.em import EMResult

    fq = tmp_path / "s.fq.gz"
    fq.touch()
    db = tmp_path / "db.sbt.zip"
    db.touch()
    ref = tmp_path / "ref.fa"
    ref.touch()

    fake_em = EMResult(abundances=np.array([1.0]), n_reads=10, converged=True, iterations=3)
    fake_align = AlignmentResult(
        alignment_matrix=np.ones((10, 1)),
        organism_names=["Staphylococcus aureus"],
        accessions=["GCF_000013425.1"],
    )

    runner = CliRunner()
    with patch("pathogeniq.cli.run_qc", return_value=(tmp_path / "filt.fq.gz", QCMetrics(100, 90, 0.1))), \
         patch("pathogeniq.cli.run_host_removal", return_value=(tmp_path / "nonhuman.fq.gz", HostRemovalMetrics(90, 10))), \
         patch("pathogeniq.cli.run_sketch_screen", return_value=[SketchHit("Staphylococcus aureus", "GCF_000013425.1", 0.05)]), \
         patch("pathogeniq.cli.run_targeted_alignment", return_value=fake_align), \
         patch("pathogeniq.cli.em_abundance", return_value=fake_em), \
         patch("pathogeniq.cli.bootstrap_ci", return_value=(np.array([0.9]), np.array([1.0]))), \
         patch("pathogeniq.cli.run_amr_screen", return_value=[]), \
         patch("pathogeniq.cli.write_pdf_report") as mock_pdf:
        result = runner.invoke(cli, [
            "run",
            "--input", str(fq),
            "--output", str(tmp_path / "out"),
            "--db", str(db),
            "--host-ref", str(ref),
            "--specimen", "blood",
            "--no-pdf",
        ])

    assert result.exit_code == 0, result.output
    mock_pdf.assert_not_called()
```

- [ ] **Step 3: Run to confirm failure**

```bash
pytest tests/test_cli.py::test_cli_run_produces_pdf tests/test_cli.py::test_cli_run_no_pdf_flag -v 2>&1 | head -30
```

Expected: FAIL — `--no-pdf` not recognised, `write_pdf_report` not imported.

- [ ] **Step 4: Rewrite `pathogeniq/cli.py`**

```python
import click
from pathlib import Path

from .amr import run_amr_screen
from .config import PipelineConfig, ReadType, SpecimenType
from .contaminants import flag_contaminants
from .em import bootstrap_ci, em_abundance
from .host_remove import run_host_removal
from .pdf_report import write_pdf_report
from .qc import run_qc
from .report import write_report
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
@click.option("--sketch-threshold", default=0.01, show_default=True)
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
    from .report import ReportEntry
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
        click.echo(f"PDF report:    {pdf_path}")

    click.echo(f"Report written to: {report_dir}")
```

- [ ] **Step 5: Run all CLI tests**

```bash
pytest tests/test_cli.py -v
```

Expected: all PASS (including the two new tests).

- [ ] **Step 6: Run the full test suite to confirm no regressions**

```bash
pytest tests/ -v --ignore=tests/integration
```

Expected: all PASS.

- [ ] **Step 7: Commit**

```bash
git add pathogeniq/cli.py pathogeniq/config.py tests/test_cli.py
git commit -m "feat: wire contaminant flags, AMR overlay, and PDF into CLI"
```

---

## Task 6: Final Integration & Push

**Files:**
- Modify: `README.md` (add AMR + PDF to pipeline diagram)

- [ ] **Step 1: Run full unit test suite**

```bash
pytest tests/ -v --ignore=tests/integration
```

Expected: all PASS, zero failures.

- [ ] **Step 2: Quick smoke test of the PDF output (visual check)**

```bash
python -c "
from pathlib import Path
from pathogeniq.config import PipelineConfig, ReadType, SpecimenType
from pathogeniq.report import ReportEntry
from pathogeniq.pdf_report import write_pdf_report
import tempfile, os

with tempfile.TemporaryDirectory() as d:
    cfg = PipelineConfig(
        input_fastq=Path('sample.fq.gz'),
        read_type=ReadType.SHORT,
        specimen_type=SpecimenType.BLOOD,
        output_dir=Path(d),
        db_tier1=Path('db.sbt.zip'),
        host_reference=Path('ref.fa'),
    )
    entries = [
        ReportEntry('Staphylococcus aureus', 0.85, 0.78, 0.92, 85, SpecimenType.BLOOD),
        ReportEntry('Cutibacterium acnes',   0.10, 0.04, 0.16, 10, SpecimenType.BLOOD, contaminant_risk=True),
    ]
    pdf = write_pdf_report(cfg, entries, amr_hits=[])
    print(f'PDF: {pdf} ({pdf.stat().st_size:,} bytes)')
"
```

Expected: `PDF: /tmp/.../report/pathogeniq_report.pdf (XXXXX bytes)`.

- [ ] **Step 3: Push to GitHub**

```bash
git push origin master
```

- [ ] **Step 4: Verify CI passes**

Check: `https://github.com/alvin8-git/PathogenIQ/actions`

---

## Self-Review Checklist

- [x] Contaminant registry covers all 4 specimen types
- [x] Grade downgrade for known contaminants tested
- [x] S. aureus correctly NOT suppressed in blood (true pathogen exception via non-match)
- [x] AMR parser handles unknown organism header → "unknown"
- [x] AMR skips gracefully when abricate not on PATH
- [x] PDF created for findings with and without AMR hits
- [x] `--no-pdf` flag suppresses PDF generation
- [x] `write_report` remains backwards-compatible (amr_hits and entries_override default to None)
- [x] `contaminant_risk` field in JSON and TSV
- [x] No placeholder text anywhere in plan
- [x] Type consistency: `AMRHit`, `ReportEntry`, `flag_contaminants` all used consistently across tasks
