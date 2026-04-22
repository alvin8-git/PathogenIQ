# PathogenIQ Core Pipeline — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a working end-to-end metagenomics pipeline (QC → host removal → sketch screening → EM abundance → TSV/JSON report) that runs on the ZymoBIOMICS FASTQ files and correctly identifies all 10 community members.

**Architecture:** Five Python modules wrap external bioinformatics tools (fastp, minimap2, sourmash) and implement a pure-NumPy EM algorithm for abundance estimation. A Click CLI ties them together into a single `pathogeniq run` command.

**Tech Stack:** Python 3.11+, NumPy, SciPy, Click, pysam, pytest. External binaries: fastp, Chopper, minimap2, BWA-MEM2, sourmash.

---

## File Map

```
/data/alvin/Metagenomics/
├── pathogeniq/
│   ├── __init__.py          # version string only
│   ├── config.py            # PipelineConfig dataclass, ReadType/SpecimenType enums
│   ├── qc.py                # Stage 1: fastp (SR) and Chopper (LR) wrappers
│   ├── host_remove.py       # Stage 2: minimap2 wrapper, returns non-human FASTQ
│   ├── sketch.py            # Stage 3: sourmash sketch + containment search
│   ├── em.py                # Stage 4: EM abundance estimation + bootstrap CI
│   ├── align.py             # Stage 4: minimap2 targeted alignment → alignment matrix
│   ├── report.py            # Stage 5: TSV + JSON report writer
│   └── cli.py               # Click CLI: `pathogeniq run`
├── tests/
│   ├── conftest.py          # Fixtures: synthetic FASTQ, tmp paths, mock subprocess
│   ├── test_config.py
│   ├── test_qc.py
│   ├── test_host_remove.py
│   ├── test_sketch.py
│   ├── test_em.py
│   ├── test_align.py
│   ├── test_report.py
│   ├── test_cli.py
│   └── integration/
│       └── test_zymo_validation.py
└── pyproject.toml
```

---

## Task 1: Project Scaffold

**Files:**
- Create: `pyproject.toml`
- Create: `pathogeniq/__init__.py`
- Create: `tests/conftest.py`

- [ ] **Step 1: Write `pyproject.toml`**

```toml
[build-system]
requires = ["setuptools>=68"]
build-backend = "setuptools.backends.legacy:BuildBackend"

[project]
name = "pathogeniq"
version = "0.1.0"
requires-python = ">=3.11"
dependencies = [
    "click>=8.1",
    "numpy>=1.26",
    "scipy>=1.12",
    "pysam>=0.22",
    "pandas>=2.2",
]

[project.scripts]
pathogeniq = "pathogeniq.cli:cli"

[project.optional-dependencies]
dev = [
    "pytest>=8.0",
    "pytest-cov>=5.0",
]

[tool.pytest.ini_options]
testpaths = ["tests"]
markers = [
    "integration: requires external tools and large data files (deselect with '-m not integration')",
]
```

- [ ] **Step 2: Create `pathogeniq/__init__.py`**

```python
__version__ = "0.1.0"
```

- [ ] **Step 3: Write `tests/conftest.py`**

```python
import gzip
import pytest
from pathlib import Path


@pytest.fixture
def synthetic_sr_fastq(tmp_path) -> Path:
    """10 synthetic 150bp Illumina reads, all high quality."""
    reads = []
    for i in range(10):
        reads.append(f"@read_{i}")
        reads.append("ACGTACGTACGTACGT" * 9 + "ACGTACGTAC")  # 150 bp
        reads.append("+")
        reads.append("I" * 150)  # Q=40 throughout
    fastq = tmp_path / "synthetic_sr.fastq.gz"
    with gzip.open(fastq, "wt") as f:
        f.write("\n".join(reads) + "\n")
    return fastq


@pytest.fixture
def synthetic_lr_fastq(tmp_path) -> Path:
    """5 synthetic 1000bp Nanopore reads."""
    reads = []
    for i in range(5):
        reads.append(f"@read_lr_{i}")
        reads.append("ACGTACGTACGT" * 83 + "ACGT")  # 1000 bp
        reads.append("+")
        reads.append("I" * 1000)
    fastq = tmp_path / "synthetic_lr.fastq.gz"
    with gzip.open(fastq, "wt") as f:
        f.write("\n".join(reads) + "\n")
    return fastq
```

- [ ] **Step 4: Install in editable mode**

```bash
cd /data/alvin/Metagenomics
pip install -e ".[dev]"
```

Expected: `Successfully installed pathogeniq-0.1.0`

- [ ] **Step 5: Verify pytest collects zero tests (nothing to fail yet)**

```bash
pytest --co -q
```

Expected: `no tests ran` or `0 tests collected`

- [ ] **Step 6: Commit**

```bash
git add pyproject.toml pathogeniq/__init__.py tests/conftest.py
git commit -m "feat: project scaffold with fixtures"
```

---

## Task 2: Config Module

**Files:**
- Create: `pathogeniq/config.py`
- Create: `tests/test_config.py`

- [ ] **Step 1: Write the failing test**

```python
# tests/test_config.py
from pathlib import Path
import pytest
from pathogeniq.config import PipelineConfig, ReadType, SpecimenType


def test_pipeline_config_defaults(tmp_path):
    cfg = PipelineConfig(
        input_fastq=tmp_path / "sample.fastq.gz",
        read_type=ReadType.SHORT,
        specimen_type=SpecimenType.BLOOD,
        output_dir=tmp_path / "out",
        db_tier1=tmp_path / "db",
        host_reference=tmp_path / "human.fa",
    )
    assert cfg.threads == 8
    assert cfg.sketch_threshold == 0.01


def test_read_type_enum_values():
    assert ReadType.SHORT == "short"
    assert ReadType.LONG == "long"


def test_specimen_type_enum_values():
    assert SpecimenType.BLOOD == "blood"
    assert SpecimenType.CSF == "csf"
    assert SpecimenType.BAL == "bal"
    assert SpecimenType.TISSUE == "tissue"


def test_pipeline_config_custom_threads(tmp_path):
    cfg = PipelineConfig(
        input_fastq=tmp_path / "s.fastq.gz",
        read_type=ReadType.LONG,
        specimen_type=SpecimenType.CSF,
        output_dir=tmp_path,
        db_tier1=tmp_path,
        host_reference=tmp_path / "h.fa",
        threads=32,
    )
    assert cfg.threads == 32
```

- [ ] **Step 2: Run test to verify it fails**

```bash
pytest tests/test_config.py -v
```

Expected: `ImportError: cannot import name 'PipelineConfig'`

- [ ] **Step 3: Write `pathogeniq/config.py`**

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
```

- [ ] **Step 4: Run tests to verify they pass**

```bash
pytest tests/test_config.py -v
```

Expected: `4 passed`

- [ ] **Step 5: Commit**

```bash
git add pathogeniq/config.py tests/test_config.py
git commit -m "feat: PipelineConfig dataclass with ReadType/SpecimenType enums"
```

---

## Task 3: QC Module

**Files:**
- Create: `pathogeniq/qc.py`
- Create: `tests/test_qc.py`

- [ ] **Step 1: Write the failing tests**

```python
# tests/test_qc.py
import json
from pathlib import Path
from unittest.mock import MagicMock, patch, mock_open
import pytest
from pathogeniq.config import PipelineConfig, ReadType, SpecimenType
from pathogeniq.qc import run_qc, QCMetrics


def _cfg(tmp_path, read_type=ReadType.SHORT):
    return PipelineConfig(
        input_fastq=tmp_path / "in.fastq.gz",
        read_type=read_type,
        specimen_type=SpecimenType.BLOOD,
        output_dir=tmp_path,
        db_tier1=tmp_path / "db",
        host_reference=tmp_path / "h.fa",
        threads=4,
    )


FASTP_JSON = {
    "summary": {
        "before_filtering": {"total_reads": 1000},
        "after_filtering": {"total_reads": 950},
    }
}


def test_run_qc_sr_calls_fastp(tmp_path):
    cfg = _cfg(tmp_path, ReadType.SHORT)
    with patch("subprocess.run") as mock_run, \
         patch("builtins.open", mock_open(read_data=json.dumps(FASTP_JSON))):
        mock_run.return_value = MagicMock(returncode=0)
        filtered, metrics = run_qc(cfg)
    cmd = mock_run.call_args[0][0]
    assert cmd[0] == "fastp"
    assert "--qualified_quality_phred" in cmd
    assert "--low_complexity_filter" in cmd


def test_run_qc_sr_returns_metrics(tmp_path):
    cfg = _cfg(tmp_path, ReadType.SHORT)
    with patch("subprocess.run"), \
         patch("builtins.open", mock_open(read_data=json.dumps(FASTP_JSON))):
        filtered, metrics = run_qc(cfg)
    assert metrics.total_reads == 1000
    assert metrics.passing_reads == 950


def test_run_qc_lr_calls_chopper(tmp_path):
    cfg = _cfg(tmp_path, ReadType.LONG)
    with patch("subprocess.run") as mock_run:
        mock_run.return_value = MagicMock(returncode=0)
        filtered, metrics = run_qc(cfg)
    cmd = mock_run.call_args[0][0]
    assert "chopper" in " ".join(cmd).lower()


def test_run_qc_creates_output_dir(tmp_path):
    cfg = _cfg(tmp_path, ReadType.SHORT)
    with patch("subprocess.run"), \
         patch("builtins.open", mock_open(read_data=json.dumps(FASTP_JSON))):
        filtered, _ = run_qc(cfg)
    assert (tmp_path / "qc").is_dir()
```

- [ ] **Step 2: Run to verify failures**

```bash
pytest tests/test_qc.py -v
```

Expected: `ImportError: cannot import name 'run_qc'`

- [ ] **Step 3: Write `pathogeniq/qc.py`**

```python
import json
import subprocess
from dataclasses import dataclass
from pathlib import Path

from .config import PipelineConfig, ReadType


@dataclass
class QCMetrics:
    total_reads: int
    passing_reads: int

    @property
    def pass_rate(self) -> float:
        return self.passing_reads / self.total_reads if self.total_reads else 0.0


def run_qc(cfg: PipelineConfig) -> tuple[Path, QCMetrics]:
    out = cfg.output_dir / "qc"
    out.mkdir(parents=True, exist_ok=True)
    filtered = out / "filtered.fastq.gz"
    json_report = out / "fastp.json"

    if cfg.read_type == ReadType.SHORT:
        cmd = [
            "fastp",
            "--in1", str(cfg.input_fastq),
            "--out1", str(filtered),
            "--json", str(json_report),
            "--qualified_quality_phred", "20",
            "--length_required", "50",
            "--low_complexity_filter",
            "--thread", str(cfg.threads),
            "--disable_adapter_trimming",
        ]
        subprocess.run(cmd, capture_output=True, text=True, check=True)
        with open(json_report) as f:
            data = json.load(f)
        metrics = QCMetrics(
            total_reads=data["summary"]["before_filtering"]["total_reads"],
            passing_reads=data["summary"]["after_filtering"]["total_reads"],
        )
    else:
        cmd = [
            "bash", "-c",
            f"gunzip -c {cfg.input_fastq} "
            f"| chopper -q 8 -l 200 --threads {cfg.threads} "
            f"| gzip > {filtered}",
        ]
        subprocess.run(cmd, capture_output=True, text=True, check=True)
        metrics = QCMetrics(total_reads=0, passing_reads=0)

    return filtered, metrics
```

- [ ] **Step 4: Run tests to verify they pass**

```bash
pytest tests/test_qc.py -v
```

Expected: `4 passed`

- [ ] **Step 5: Commit**

```bash
git add pathogeniq/qc.py tests/test_qc.py
git commit -m "feat: QC module with fastp (SR) and Chopper (LR) wrappers"
```

---

## Task 4: Host Removal Module

**Files:**
- Create: `pathogeniq/host_remove.py`
- Create: `tests/test_host_remove.py`

- [ ] **Step 1: Write the failing tests**

```python
# tests/test_host_remove.py
from pathlib import Path
from unittest.mock import MagicMock, patch
import pytest
from pathogeniq.config import PipelineConfig, ReadType, SpecimenType
from pathogeniq.host_remove import run_host_removal, HostRemovalMetrics


def _cfg(tmp_path, read_type=ReadType.SHORT):
    return PipelineConfig(
        input_fastq=tmp_path / "filtered.fastq.gz",
        read_type=read_type,
        specimen_type=SpecimenType.BLOOD,
        output_dir=tmp_path,
        db_tier1=tmp_path / "db",
        host_reference=tmp_path / "human.fa",
        threads=4,
    )


def test_host_removal_sr_uses_bwa(tmp_path):
    cfg = _cfg(tmp_path, ReadType.SHORT)
    filtered = tmp_path / "filtered.fastq.gz"
    with patch("subprocess.run") as mock_run:
        mock_run.return_value = MagicMock(returncode=0)
        nonhuman, metrics = run_host_removal(cfg, filtered)
    calls = [" ".join(c[0][0]) for c in mock_run.call_args_list]
    assert any("bwa-mem2" in c or "minimap2" in c for c in calls)


def test_host_removal_lr_uses_minimap2(tmp_path):
    cfg = _cfg(tmp_path, ReadType.LONG)
    filtered = tmp_path / "filtered.fastq.gz"
    with patch("subprocess.run") as mock_run:
        mock_run.return_value = MagicMock(returncode=0)
        nonhuman, metrics = run_host_removal(cfg, filtered)
    calls = [" ".join(c[0][0]) for c in mock_run.call_args_list]
    assert any("minimap2" in c for c in calls)


def test_host_removal_creates_output_dir(tmp_path):
    cfg = _cfg(tmp_path)
    filtered = tmp_path / "filtered.fastq.gz"
    with patch("subprocess.run"):
        nonhuman, _ = run_host_removal(cfg, filtered)
    assert (tmp_path / "host_removal").is_dir()


def test_host_removal_metrics_fields():
    m = HostRemovalMetrics(total_reads=1000, human_reads=800, nonhuman_reads=200)
    assert m.microbial_fraction == pytest.approx(0.2)
    assert m.human_fraction == pytest.approx(0.8)
```

- [ ] **Step 2: Run to verify failures**

```bash
pytest tests/test_host_remove.py -v
```

Expected: `ImportError: cannot import name 'run_host_removal'`

- [ ] **Step 3: Write `pathogeniq/host_remove.py`**

```python
import subprocess
from dataclasses import dataclass
from pathlib import Path

from .config import PipelineConfig, ReadType


@dataclass
class HostRemovalMetrics:
    total_reads: int
    human_reads: int
    nonhuman_reads: int

    @property
    def microbial_fraction(self) -> float:
        return self.nonhuman_reads / self.total_reads if self.total_reads else 0.0

    @property
    def human_fraction(self) -> float:
        return self.human_reads / self.total_reads if self.total_reads else 0.0


def run_host_removal(cfg: PipelineConfig, filtered_fastq: Path) -> tuple[Path, HostRemovalMetrics]:
    out = cfg.output_dir / "host_removal"
    out.mkdir(parents=True, exist_ok=True)
    nonhuman_fastq = out / "nonhuman.fastq.gz"
    sam_file = out / "host_aligned.sam"

    if cfg.read_type == ReadType.SHORT:
        align_cmd = [
            "bwa-mem2", "mem",
            "-t", str(cfg.threads),
            str(cfg.host_reference),
            str(filtered_fastq),
            "-o", str(sam_file),
        ]
    else:
        align_cmd = [
            "minimap2",
            "-ax", "map-ont",
            "-t", str(cfg.threads),
            str(cfg.host_reference),
            str(filtered_fastq),
            "-o", str(sam_file),
        ]

    subprocess.run(align_cmd, capture_output=True, text=True, check=True)

    # Extract unmapped reads with samtools
    extract_cmd = [
        "bash", "-c",
        f"samtools view -f 4 -b {sam_file} "
        f"| samtools fastq - "
        f"| gzip > {nonhuman_fastq}",
    ]
    subprocess.run(extract_cmd, capture_output=True, text=True, check=True)

    # Placeholder metrics — real implementation parses samtools flagstat
    metrics = HostRemovalMetrics(total_reads=0, human_reads=0, nonhuman_reads=0)
    return nonhuman_fastq, metrics
```

- [ ] **Step 4: Run tests to verify they pass**

```bash
pytest tests/test_host_remove.py -v
```

Expected: `4 passed`

- [ ] **Step 5: Commit**

```bash
git add pathogeniq/host_remove.py tests/test_host_remove.py
git commit -m "feat: host removal module with BWA-MEM2 (SR) and minimap2 (LR)"
```

---

## Task 5: Sketch Screening Module

**Files:**
- Create: `pathogeniq/sketch.py`
- Create: `tests/test_sketch.py`

- [ ] **Step 1: Write the failing tests**

```python
# tests/test_sketch.py
from pathlib import Path
from unittest.mock import MagicMock, patch, mock_open
import pytest
from pathogeniq.config import PipelineConfig, ReadType, SpecimenType
from pathogeniq.sketch import run_sketch_screen, SketchHit


def _cfg(tmp_path):
    return PipelineConfig(
        input_fastq=tmp_path / "nonhuman.fastq.gz",
        read_type=ReadType.SHORT,
        specimen_type=SpecimenType.BLOOD,
        output_dir=tmp_path,
        db_tier1=tmp_path / "db.sbt.zip",
        host_reference=tmp_path / "h.fa",
        threads=4,
        sketch_threshold=0.01,
    )


SOURMASH_CSV = """name,containment,similarity,filename
GCF_000001405.40_GRCh38_ecoli.fna,0.35,0.30,ecoli.fna
GCF_000007545.1_staphaureus.fna,0.12,0.10,staph.fna
GCF_000001285.2_weak.fna,0.005,0.004,weak.fna
"""


def test_sketch_screen_returns_hits_above_threshold(tmp_path):
    cfg = _cfg(tmp_path)
    nonhuman = tmp_path / "nonhuman.fastq.gz"
    with patch("subprocess.run") as mock_run, \
         patch("builtins.open", mock_open(read_data=SOURMASH_CSV)):
        mock_run.return_value = MagicMock(returncode=0)
        hits = run_sketch_screen(cfg, nonhuman)
    # weak hit (0.005) is below threshold (0.01)
    assert len(hits) == 2
    names = [h.name for h in hits]
    assert any("ecoli" in n for n in names)
    assert all(h.containment >= 0.01 for h in hits)


def test_sketch_screen_hit_fields(tmp_path):
    cfg = _cfg(tmp_path)
    nonhuman = tmp_path / "nonhuman.fastq.gz"
    with patch("subprocess.run"), \
         patch("builtins.open", mock_open(read_data=SOURMASH_CSV)):
        hits = run_sketch_screen(cfg, nonhuman)
    hit = hits[0]
    assert hasattr(hit, "name")
    assert hasattr(hit, "containment")
    assert hasattr(hit, "genome_path")


def test_sketch_screen_calls_sourmash(tmp_path):
    cfg = _cfg(tmp_path)
    nonhuman = tmp_path / "nonhuman.fastq.gz"
    with patch("subprocess.run") as mock_run, \
         patch("builtins.open", mock_open(read_data=SOURMASH_CSV)):
        mock_run.return_value = MagicMock(returncode=0)
        run_sketch_screen(cfg, nonhuman)
    calls = [c[0][0] for c in mock_run.call_args_list]
    assert any(c[0] == "sourmash" for c in calls)


def test_sketch_screen_creates_output_dir(tmp_path):
    cfg = _cfg(tmp_path)
    nonhuman = tmp_path / "nonhuman.fastq.gz"
    with patch("subprocess.run"), \
         patch("builtins.open", mock_open(read_data=SOURMASH_CSV)):
        run_sketch_screen(cfg, nonhuman)
    assert (tmp_path / "sketch").is_dir()
```

- [ ] **Step 2: Run to verify failures**

```bash
pytest tests/test_sketch.py -v
```

Expected: `ImportError: cannot import name 'run_sketch_screen'`

- [ ] **Step 3: Write `pathogeniq/sketch.py`**

```python
import csv
import subprocess
from dataclasses import dataclass
from pathlib import Path

from .config import PipelineConfig


@dataclass
class SketchHit:
    name: str
    containment: float
    genome_path: Path


def run_sketch_screen(cfg: PipelineConfig, nonhuman_fastq: Path) -> list[SketchHit]:
    out = cfg.output_dir / "sketch"
    out.mkdir(parents=True, exist_ok=True)
    sample_sig = out / "sample.sig"
    results_csv = out / "results.csv"

    # Sketch the sample reads
    sketch_cmd = [
        "sourmash", "sketch", "dna",
        "-p", "k=31,scaled=1000",
        str(nonhuman_fastq),
        "-o", str(sample_sig),
    ]
    subprocess.run(sketch_cmd, capture_output=True, text=True, check=True)

    # Search against Tier-1 database
    search_cmd = [
        "sourmash", "search",
        str(sample_sig),
        str(cfg.db_tier1),
        "--containment",
        "--threshold", str(cfg.sketch_threshold),
        "-o", str(results_csv),
    ]
    subprocess.run(search_cmd, capture_output=True, text=True, check=True)

    hits: list[SketchHit] = []
    with open(results_csv) as f:
        reader = csv.DictReader(f)
        for row in reader:
            containment = float(row["containment"])
            if containment >= cfg.sketch_threshold:
                hits.append(SketchHit(
                    name=row["name"],
                    containment=containment,
                    genome_path=Path(row["filename"]),
                ))
    return hits
```

- [ ] **Step 4: Run tests to verify they pass**

```bash
pytest tests/test_sketch.py -v
```

Expected: `4 passed`

- [ ] **Step 5: Commit**

```bash
git add pathogeniq/sketch.py tests/test_sketch.py
git commit -m "feat: sketch screening module with sourmash containment search"
```

---

## Task 6: EM Abundance Estimation (Novel Core)

**Files:**
- Create: `pathogeniq/em.py`
- Create: `tests/test_em.py`

This module is pure NumPy — no subprocess mocking needed. Tests are deterministic.

- [ ] **Step 1: Write the failing tests**

```python
# tests/test_em.py
import numpy as np
import pytest
from pathogeniq.em import em_abundance, bootstrap_ci, EMResult


def test_em_all_reads_one_organism():
    """All reads uniquely map to organism 0."""
    matrix = np.array([[1, 0], [1, 0], [1, 0]], dtype=float)
    result = em_abundance(matrix)
    assert result.abundances[0] == pytest.approx(1.0, abs=1e-6)
    assert result.abundances[1] == pytest.approx(0.0, abs=1e-6)


def test_em_equal_split():
    """3 reads each to two organisms."""
    matrix = np.array([[1, 0], [1, 0], [1, 0], [0, 1], [0, 1], [0, 1]], dtype=float)
    result = em_abundance(matrix)
    assert result.abundances[0] == pytest.approx(0.5, abs=0.01)
    assert result.abundances[1] == pytest.approx(0.5, abs=0.01)


def test_em_multimapper_resolves_by_unique_evidence():
    """Multi-mappers are distributed proportional to unique read evidence."""
    # 6 unique to org 0, 2 unique to org 1, 2 multi-mappers
    matrix = np.array([
        [1, 0], [1, 0], [1, 0], [1, 0], [1, 0], [1, 0],
        [0, 1], [0, 1],
        [1, 1], [1, 1],
    ], dtype=float)
    result = em_abundance(matrix)
    assert result.abundances[0] > result.abundances[1]


def test_em_sums_to_one():
    """Abundances always sum to 1."""
    rng = np.random.default_rng(42)
    matrix = (rng.random((100, 10)) > 0.8).astype(float)
    result = em_abundance(matrix)
    assert result.abundances.sum() == pytest.approx(1.0, abs=1e-6)


def test_em_zymo_composition():
    """Simulate ZymoBIOMICS: 8 bacteria at 12%, 2 yeasts at 2%."""
    true_abundances = np.array([0.12] * 8 + [0.02] * 2)
    reads_per_org = (true_abundances * 10000).astype(int)
    rows = []
    for org_idx, count in enumerate(reads_per_org):
        row = np.zeros(10)
        row[org_idx] = 1.0
        rows.extend([row] * count)
    matrix = np.array(rows)
    result = em_abundance(matrix)
    expected = true_abundances / true_abundances.sum()
    assert np.allclose(result.abundances, expected, atol=0.02)


def test_em_returns_organism_count():
    matrix = np.ones((5, 3), dtype=float)
    result = em_abundance(matrix)
    assert len(result.abundances) == 3


def test_bootstrap_ci_shape():
    matrix = np.array([[1, 0], [1, 0], [0, 1], [0, 1]], dtype=float)
    lower, upper = bootstrap_ci(matrix, n_bootstrap=50)
    assert lower.shape == (2,)
    assert upper.shape == (2,)
    assert np.all(lower <= upper)


def test_bootstrap_ci_coverage():
    """95% CI lower bound should be <= point estimate <= upper bound."""
    matrix = np.array([[1, 0]] * 7 + [[0, 1]] * 3, dtype=float)
    result = em_abundance(matrix)
    lower, upper = bootstrap_ci(matrix, n_bootstrap=200)
    for i in range(2):
        assert lower[i] <= result.abundances[i] + 1e-6
        assert upper[i] >= result.abundances[i] - 1e-6
```

- [ ] **Step 2: Run to verify failures**

```bash
pytest tests/test_em.py -v
```

Expected: `ImportError: cannot import name 'em_abundance'`

- [ ] **Step 3: Write `pathogeniq/em.py`**

```python
from dataclasses import dataclass
import numpy as np


@dataclass
class EMResult:
    abundances: np.ndarray  # shape (n_organisms,), sums to 1
    n_reads: int
    n_organisms: int
    iterations: int


def em_abundance(
    alignment_matrix: np.ndarray,
    max_iter: int = 200,
    tol: float = 1e-8,
) -> EMResult:
    """
    EM algorithm for abundance estimation from multi-mapping alignments.

    alignment_matrix: shape (n_reads, n_organisms)
        Entry [i, j] = 1 if read i maps to organism j, else 0.
        Rows with all zeros are silently ignored.
    """
    n_reads, n_orgs = alignment_matrix.shape
    # Remove reads that map nowhere
    row_mask = alignment_matrix.sum(axis=1) > 0
    matrix = alignment_matrix[row_mask]
    n_effective = matrix.shape[0]

    theta = np.ones(n_orgs) / n_orgs
    iterations = 0

    for iterations in range(1, max_iter + 1):
        # E-step: weighted responsibilities
        weighted = matrix * theta[np.newaxis, :]
        row_sums = weighted.sum(axis=1, keepdims=True)
        row_sums = np.where(row_sums == 0, 1.0, row_sums)
        responsibilities = weighted / row_sums

        # M-step: update abundances
        theta_new = responsibilities.sum(axis=0) / n_effective

        if np.max(np.abs(theta_new - theta)) < tol:
            theta = theta_new
            break
        theta = theta_new

    return EMResult(
        abundances=theta,
        n_reads=n_effective,
        n_organisms=n_orgs,
        iterations=iterations,
    )


def bootstrap_ci(
    alignment_matrix: np.ndarray,
    n_bootstrap: int = 100,
    alpha: float = 0.05,
    seed: int = 42,
) -> tuple[np.ndarray, np.ndarray]:
    """Bootstrap (1-alpha) confidence interval on EM abundance estimates."""
    rng = np.random.default_rng(seed)
    n_reads = alignment_matrix.shape[0]
    estimates = np.zeros((n_bootstrap, alignment_matrix.shape[1]))

    for b in range(n_bootstrap):
        idx = rng.integers(0, n_reads, size=n_reads)
        boot = alignment_matrix[idx]
        estimates[b] = em_abundance(boot).abundances

    lower = np.percentile(estimates, 100 * alpha / 2, axis=0)
    upper = np.percentile(estimates, 100 * (1 - alpha / 2), axis=0)
    return lower, upper
```

- [ ] **Step 4: Run tests to verify they pass**

```bash
pytest tests/test_em.py -v
```

Expected: `8 passed`

- [ ] **Step 5: Commit**

```bash
git add pathogeniq/em.py tests/test_em.py
git commit -m "feat: EM abundance estimation with bootstrap CI"
```

---

## Task 7: Targeted Alignment Module

**Files:**
- Create: `pathogeniq/align.py`
- Create: `tests/test_align.py`

- [ ] **Step 1: Write the failing tests**

```python
# tests/test_align.py
import numpy as np
from pathlib import Path
from unittest.mock import MagicMock, patch
import pytest
from pathogeniq.config import PipelineConfig, ReadType, SpecimenType
from pathogeniq.sketch import SketchHit
from pathogeniq.align import run_targeted_alignment, AlignmentResult


def _cfg(tmp_path):
    return PipelineConfig(
        input_fastq=tmp_path / "nonhuman.fastq.gz",
        read_type=ReadType.SHORT,
        specimen_type=SpecimenType.BLOOD,
        output_dir=tmp_path,
        db_tier1=tmp_path / "db",
        host_reference=tmp_path / "h.fa",
        threads=4,
    )


def _hits(tmp_path):
    return [
        SketchHit("ecoli", 0.35, tmp_path / "ecoli.fa"),
        SketchHit("staph", 0.12, tmp_path / "staph.fa"),
    ]


def test_targeted_alignment_calls_minimap2(tmp_path):
    cfg = _cfg(tmp_path)
    nonhuman = tmp_path / "nonhuman.fastq.gz"
    hits = _hits(tmp_path)
    with patch("subprocess.run") as mock_run:
        mock_run.return_value = MagicMock(returncode=0)
        with patch("pathogeniq.align._parse_paf", return_value=np.zeros((0, 2))):
            run_targeted_alignment(cfg, nonhuman, hits)
    calls = [c[0][0] for c in mock_run.call_args_list]
    assert any(c[0] == "minimap2" for c in calls)


def test_targeted_alignment_result_shape(tmp_path):
    cfg = _cfg(tmp_path)
    nonhuman = tmp_path / "nonhuman.fastq.gz"
    hits = _hits(tmp_path)
    fake_matrix = np.array([[1, 0], [0, 1], [1, 1]], dtype=float)
    with patch("subprocess.run"), \
         patch("pathogeniq.align._parse_paf", return_value=fake_matrix):
        result = run_targeted_alignment(cfg, nonhuman, hits)
    assert result.alignment_matrix.shape[1] == len(hits)
    assert result.organism_names == ["ecoli", "staph"]


def test_alignment_result_fields(tmp_path):
    cfg = _cfg(tmp_path)
    nonhuman = tmp_path / "nonhuman.fastq.gz"
    hits = _hits(tmp_path)
    fake_matrix = np.eye(2, dtype=float)
    with patch("subprocess.run"), \
         patch("pathogeniq.align._parse_paf", return_value=fake_matrix):
        result = run_targeted_alignment(cfg, nonhuman, hits)
    assert hasattr(result, "alignment_matrix")
    assert hasattr(result, "organism_names")
    assert hasattr(result, "read_ids")
```

- [ ] **Step 2: Run to verify failures**

```bash
pytest tests/test_align.py -v
```

Expected: `ImportError: cannot import name 'run_targeted_alignment'`

- [ ] **Step 3: Write `pathogeniq/align.py`**

```python
import subprocess
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from .config import PipelineConfig, ReadType
from .sketch import SketchHit


@dataclass
class AlignmentResult:
    alignment_matrix: np.ndarray   # shape (n_reads, n_organisms)
    organism_names: list[str]
    read_ids: list[str]


def run_targeted_alignment(
    cfg: PipelineConfig,
    nonhuman_fastq: Path,
    hits: list[SketchHit],
) -> AlignmentResult:
    out = cfg.output_dir / "align"
    out.mkdir(parents=True, exist_ok=True)

    all_read_ids: list[str] = []
    # organism_idx -> set of read_ids
    org_reads: dict[int, set[str]] = {i: set() for i in range(len(hits))}

    preset = "sr" if cfg.read_type == ReadType.SHORT else "map-ont"

    for org_idx, hit in enumerate(hits):
        paf_file = out / f"org_{org_idx}.paf"
        cmd = [
            "minimap2",
            "-c",                     # output CIGAR in PAF
            "-x", preset,
            "-t", str(cfg.threads),
            str(hit.genome_path),
            str(nonhuman_fastq),
            "-o", str(paf_file),
        ]
        subprocess.run(cmd, capture_output=True, text=True, check=True)
        read_set = _parse_paf(paf_file)
        org_reads[org_idx] = read_set
        all_read_ids.extend(r for r in read_set if r not in all_read_ids)

    n_reads = len(all_read_ids)
    n_orgs = len(hits)
    matrix = np.zeros((n_reads, n_orgs), dtype=float)
    read_idx = {rid: i for i, rid in enumerate(all_read_ids)}

    for org_idx, read_set in org_reads.items():
        for rid in read_set:
            if rid in read_idx:
                matrix[read_idx[rid], org_idx] = 1.0

    return AlignmentResult(
        alignment_matrix=matrix,
        organism_names=[h.name for h in hits],
        read_ids=all_read_ids,
    )


def _parse_paf(paf_file: Path) -> set[str]:
    """Return set of read IDs that aligned (mapq > 0)."""
    read_ids: set[str] = set()
    if not paf_file.exists():
        return read_ids
    with open(paf_file) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 12:
                continue
            mapq = int(parts[11])
            if mapq > 0:
                read_ids.add(parts[0])
    return read_ids
```

- [ ] **Step 4: Run tests to verify they pass**

```bash
pytest tests/test_align.py -v
```

Expected: `3 passed`

- [ ] **Step 5: Commit**

```bash
git add pathogeniq/align.py tests/test_align.py
git commit -m "feat: targeted alignment module with PAF parsing"
```

---

## Task 8: Report Module

**Files:**
- Create: `pathogeniq/report.py`
- Create: `tests/test_report.py`

- [ ] **Step 1: Write the failing tests**

```python
# tests/test_report.py
import json
import csv
from pathlib import Path
import numpy as np
import pytest
from pathogeniq.config import PipelineConfig, ReadType, SpecimenType
from pathogeniq.em import EMResult
from pathogeniq.report import write_report, ReportEntry, EvidenceGrade


def _cfg(tmp_path):
    return PipelineConfig(
        input_fastq=tmp_path / "in.fastq.gz",
        read_type=ReadType.SHORT,
        specimen_type=SpecimenType.BLOOD,
        output_dir=tmp_path,
        db_tier1=tmp_path / "db",
        host_reference=tmp_path / "h.fa",
    )


def _em_result():
    return EMResult(
        abundances=np.array([0.70, 0.28, 0.02]),
        n_reads=1000,
        n_organisms=3,
        iterations=45,
    )


def test_write_report_creates_tsv(tmp_path):
    cfg = _cfg(tmp_path)
    em = _em_result()
    names = ["Escherichia coli", "Klebsiella pneumoniae", "Cutibacterium acnes"]
    lower = np.array([0.65, 0.23, 0.005])
    upper = np.array([0.75, 0.33, 0.04])
    write_report(cfg, names, em, lower, upper)
    tsv = tmp_path / "report" / "pathogeniq_report.tsv"
    assert tsv.exists()


def test_write_report_creates_json(tmp_path):
    cfg = _cfg(tmp_path)
    em = _em_result()
    names = ["Escherichia coli", "Klebsiella pneumoniae", "Cutibacterium acnes"]
    lower = np.array([0.65, 0.23, 0.005])
    upper = np.array([0.75, 0.33, 0.04])
    write_report(cfg, names, em, lower, upper)
    jsn = tmp_path / "report" / "pathogeniq_report.json"
    assert jsn.exists()
    with open(jsn) as f:
        data = json.load(f)
    assert "findings" in data
    assert len(data["findings"]) == 3


def test_evidence_grade_blood_high_confidence():
    entry = ReportEntry(
        organism="Escherichia coli",
        abundance=0.70,
        ci_lower=0.65,
        ci_upper=0.75,
        read_count=700,
        specimen_type=SpecimenType.BLOOD,
    )
    assert entry.grade == EvidenceGrade.A


def test_evidence_grade_blood_low_reads():
    entry = ReportEntry(
        organism="Cutibacterium acnes",
        abundance=0.001,
        ci_lower=0.0,
        ci_upper=0.01,
        read_count=1,
        specimen_type=SpecimenType.BLOOD,
    )
    assert entry.grade == EvidenceGrade.C


def test_tsv_columns(tmp_path):
    cfg = _cfg(tmp_path)
    em = _em_result()
    names = ["Organism A", "Organism B", "Organism C"]
    lower = np.array([0.60, 0.20, 0.00])
    upper = np.array([0.80, 0.36, 0.04])
    write_report(cfg, names, em, lower, upper)
    tsv = tmp_path / "report" / "pathogeniq_report.tsv"
    with open(tsv) as f:
        reader = csv.DictReader(f, delimiter="\t")
        rows = list(reader)
    assert "organism" in rows[0]
    assert "abundance_pct" in rows[0]
    assert "ci_lower_pct" in rows[0]
    assert "ci_upper_pct" in rows[0]
    assert "grade" in rows[0]
```

- [ ] **Step 2: Run to verify failures**

```bash
pytest tests/test_report.py -v
```

Expected: `ImportError: cannot import name 'write_report'`

- [ ] **Step 3: Write `pathogeniq/report.py`**

```python
import csv
import json
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path

import numpy as np

from .config import PipelineConfig, SpecimenType
from .em import EMResult


class EvidenceGrade(str, Enum):
    A = "A"
    B = "B"
    C = "C"
    X = "X"


# Min reads for Grade A by specimen type
_MIN_READS: dict[SpecimenType, int] = {
    SpecimenType.BLOOD: 3,
    SpecimenType.CSF: 2,
    SpecimenType.BAL: 10,
    SpecimenType.TISSUE: 10,
}

# CI width threshold for Grade A (tight CI = high confidence)
_MAX_CI_WIDTH_A = 0.15


@dataclass
class ReportEntry:
    organism: str
    abundance: float
    ci_lower: float
    ci_upper: float
    read_count: int
    specimen_type: SpecimenType

    @property
    def grade(self) -> EvidenceGrade:
        min_reads = _MIN_READS.get(self.specimen_type, 5)
        ci_width = self.ci_upper - self.ci_lower
        if self.read_count >= min_reads and ci_width <= _MAX_CI_WIDTH_A:
            return EvidenceGrade.A
        if self.read_count >= min_reads:
            return EvidenceGrade.B
        return EvidenceGrade.C


def write_report(
    cfg: PipelineConfig,
    organism_names: list[str],
    em_result: EMResult,
    ci_lower: np.ndarray,
    ci_upper: np.ndarray,
) -> Path:
    out = cfg.output_dir / "report"
    out.mkdir(parents=True, exist_ok=True)

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

    tsv_path = out / "pathogeniq_report.tsv"
    with open(tsv_path, "w", newline="") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=["organism", "abundance_pct", "ci_lower_pct", "ci_upper_pct", "read_count", "grade"],
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
            }
            for e in entries
        ],
    }
    with open(json_path, "w") as f:
        json.dump(payload, f, indent=2)

    return out
```

- [ ] **Step 4: Run tests to verify they pass**

```bash
pytest tests/test_report.py -v
```

Expected: `5 passed`

- [ ] **Step 5: Commit**

```bash
git add pathogeniq/report.py tests/test_report.py
git commit -m "feat: TSV + JSON report with evidence grading"
```

---

## Task 9: CLI

**Files:**
- Create: `pathogeniq/cli.py`
- Create: `tests/test_cli.py`

- [ ] **Step 1: Write the failing tests**

```python
# tests/test_cli.py
from click.testing import CliRunner
from unittest.mock import patch, MagicMock
import numpy as np
from pathogeniq.cli import cli
from pathogeniq.em import EMResult


def test_cli_run_requires_input():
    runner = CliRunner()
    result = runner.invoke(cli, ["run"])
    assert result.exit_code != 0
    assert "Missing" in result.output or "Error" in result.output


def test_cli_run_help():
    runner = CliRunner()
    result = runner.invoke(cli, ["run", "--help"])
    assert result.exit_code == 0
    assert "--input" in result.output
    assert "--specimen" in result.output
    assert "--read-type" in result.output


def test_cli_run_invokes_pipeline(tmp_path):
    runner = CliRunner()
    fake_fastq = tmp_path / "sample.fastq.gz"
    fake_fastq.touch()
    fake_db = tmp_path / "db.sbt.zip"
    fake_db.touch()
    fake_ref = tmp_path / "human.fa"
    fake_ref.touch()

    em_result = EMResult(
        abundances=np.array([1.0]),
        n_reads=100,
        n_organisms=1,
        iterations=10,
    )

    with patch("pathogeniq.cli.run_qc") as mock_qc, \
         patch("pathogeniq.cli.run_host_removal") as mock_hr, \
         patch("pathogeniq.cli.run_sketch_screen") as mock_sketch, \
         patch("pathogeniq.cli.run_targeted_alignment") as mock_align, \
         patch("pathogeniq.cli.em_abundance") as mock_em, \
         patch("pathogeniq.cli.bootstrap_ci") as mock_ci, \
         patch("pathogeniq.cli.write_report") as mock_report:

        from pathogeniq.qc import QCMetrics
        from pathogeniq.host_remove import HostRemovalMetrics
        from pathogeniq.sketch import SketchHit
        from pathogeniq.align import AlignmentResult
        import numpy as np

        mock_qc.return_value = (tmp_path / "filtered.fastq.gz", QCMetrics(100, 95))
        mock_hr.return_value = (tmp_path / "nonhuman.fastq.gz", HostRemovalMetrics(95, 80, 15))
        mock_sketch.return_value = [SketchHit("org1", 0.5, tmp_path / "org1.fa")]
        mock_align.return_value = AlignmentResult(
            alignment_matrix=np.array([[1.0]]),
            organism_names=["org1"],
            read_ids=["r1"],
        )
        mock_em.return_value = em_result
        mock_ci.return_value = (np.array([0.9]), np.array([1.0]))
        mock_report.return_value = tmp_path / "report"

        result = runner.invoke(cli, [
            "run",
            "--input", str(fake_fastq),
            "--output", str(tmp_path / "out"),
            "--db", str(fake_db),
            "--host-ref", str(fake_ref),
            "--specimen", "blood",
            "--read-type", "short",
        ])

    assert result.exit_code == 0, result.output
```

- [ ] **Step 2: Run to verify failures**

```bash
pytest tests/test_cli.py -v
```

Expected: `ImportError: cannot import name 'cli'`

- [ ] **Step 3: Write `pathogeniq/cli.py`**

```python
import click
import numpy as np
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
        click.echo("No candidates above threshold. No pathogens detected.")
        return

    click.echo("[4/5] Targeted alignment + EM abundance...")
    align_result = run_targeted_alignment(cfg, nonhuman, hits)
    em_result = em_abundance(align_result.alignment_matrix)
    ci_lower, ci_upper = bootstrap_ci(align_result.alignment_matrix, n_bootstrap=cfg.n_bootstrap)

    click.echo("[5/5] Generating report...")
    report_dir = write_report(cfg, align_result.organism_names, em_result, ci_lower, ci_upper)
    click.echo(f"Report written to: {report_dir}")
```

- [ ] **Step 4: Run tests to verify they pass**

```bash
pytest tests/test_cli.py -v
```

Expected: `3 passed`

- [ ] **Step 5: Run full test suite**

```bash
pytest -v --tb=short
```

Expected: all tests pass

- [ ] **Step 6: Commit**

```bash
git add pathogeniq/cli.py tests/test_cli.py
git commit -m "feat: Click CLI tying all pipeline stages together"
```

---

## Task 10: ZymoBIOMICS Integration Validation

**Files:**
- Create: `tests/integration/test_zymo_validation.py`
- Create: `tests/integration/__init__.py`

This task runs the real pipeline on `ZymoStandardsFromMao/V350234554_L01_UDB-32.Zymo_Std.fq.gz` and validates that all 10 ZymoBIOMICS organisms are detected at approximately correct abundances.

Prerequisites: fastp, minimap2, BWA-MEM2, sourmash installed; human reference indexed; Tier-1 DB built.

- [ ] **Step 1: Create integration test**

```python
# tests/integration/__init__.py
# (empty)
```

```python
# tests/integration/test_zymo_validation.py
"""
Integration test: run core pipeline on ZymoBIOMICS replicate 1 and verify
all 10 community members are detected at expected abundances.

Run with: pytest tests/integration/ -v -m integration

Requires:
  - fastp, minimap2, bwa-mem2, sourmash on PATH
  - HUMAN_REF env var: path to BWA-MEM2-indexed GRCh38
  - TIER1_DB env var: path to sourmash SBT database
"""
import os
import json
import pytest
from pathlib import Path
from click.testing import CliRunner
from pathogeniq.cli import cli

ZYMO_FASTQ = Path("/data/alvin/Metagenomics/ZymoStandardsFromMao/V350234554_L01_UDB-32.Zymo_Std.fq.gz")

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

# Expected gDNA abundance: 12% for bacteria, 2% for yeasts
EXPECTED_ABUNDANCE = {org: 0.12 for org in EXPECTED_ORGANISMS[:8]}
EXPECTED_ABUNDANCE.update({org: 0.02 for org in EXPECTED_ORGANISMS[8:]})


@pytest.mark.integration
@pytest.mark.skipif(not ZYMO_FASTQ.exists(), reason="ZymoBIOMICS FASTQ not found")
def test_zymo_all_organisms_detected(tmp_path):
    human_ref = os.environ.get("HUMAN_REF")
    tier1_db = os.environ.get("TIER1_DB")
    if not human_ref or not tier1_db:
        pytest.skip("HUMAN_REF and TIER1_DB env vars required")

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
    assert report_json.exists()

    with open(report_json) as f:
        report = json.load(f)

    detected = {f["organism"] for f in report["findings"]}

    for org in EXPECTED_ORGANISMS:
        assert org in detected, f"{org} not detected in report"


@pytest.mark.integration
@pytest.mark.skipif(not ZYMO_FASTQ.exists(), reason="ZymoBIOMICS FASTQ not found")
def test_zymo_abundance_accuracy(tmp_path):
    human_ref = os.environ.get("HUMAN_REF")
    tier1_db = os.environ.get("TIER1_DB")
    if not human_ref or not tier1_db:
        pytest.skip("HUMAN_REF and TIER1_DB env vars required")

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

    with open(tmp_path / "out" / "report" / "pathogeniq_report.json") as f:
        report = json.load(f)

    findings_by_org = {f["organism"]: f for f in report["findings"]}

    for org, expected_pct in EXPECTED_ABUNDANCE.items():
        if org not in findings_by_org:
            continue
        measured_pct = findings_by_org[org]["abundance_pct"] / 100
        # Allow 20% relative error (per validation spec)
        assert abs(measured_pct - expected_pct) / expected_pct < 0.20, (
            f"{org}: expected {expected_pct:.1%}, got {measured_pct:.1%}"
        )
```

- [ ] **Step 2: Verify integration tests are collected but skipped without env vars**

```bash
pytest tests/integration/ -v -m integration
```

Expected: tests collected but skipped with `HUMAN_REF and TIER1_DB env vars required`

- [ ] **Step 3: Run full unit test suite (no integration)**

```bash
pytest -v -m "not integration"
```

Expected: all unit tests pass

- [ ] **Step 4: Commit**

```bash
git add tests/integration/
git commit -m "test: ZymoBIOMICS integration validation harness"
```

- [ ] **Step 5: Push everything to GitHub**

```bash
git push origin master
```

---

## Running the Complete Test Suite

```bash
# Unit tests only (fast, no external tools required)
pytest -v -m "not integration"

# Integration tests (requires tools + env vars)
HUMAN_REF=/path/to/GRCh38.fa TIER1_DB=/path/to/tier1.sbt.zip \
  pytest tests/integration/ -v -m integration

# Coverage report
pytest -m "not integration" --cov=pathogeniq --cov-report=term-missing
```

---

## Next Steps (Plan 2)

Plan 2 covers:
- Specimen-aware contaminant priors (BAL normal flora model)
- AMR gene overlay via ABRicate + CARD
- PDF clinical report with full evidence grading
- Negative result interpretation and sensitivity statements
