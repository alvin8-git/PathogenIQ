from unittest.mock import MagicMock, patch, mock_open
from pathogeniq.config import PipelineConfig, ReadType, SpecimenType
from pathogeniq.sketch import run_sketch_screen


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


def test_sketch_screen_populates_taxon_id(tmp_path):
    # taxon_id (GCF accession) is parsed from the CSV name even on the fallback path
    cfg = _cfg(tmp_path)
    nonhuman = tmp_path / "nonhuman.fastq.gz"
    with patch("subprocess.run"), \
         patch("builtins.open", mock_open(read_data=SOURMASH_CSV)):
        hits = run_sketch_screen(cfg, nonhuman)
    taxon_ids = [h.taxon_id for h in hits]
    assert "GCF_000001405.40" in taxon_ids
    assert all(t.startswith("GCF_") for t in taxon_ids)


def test_sketch_screen_taxon_id_empty_without_accession(tmp_path):
    cfg = _cfg(tmp_path)
    csv_no_acc = "name,containment,similarity,filename\nsomeorg,0.5,0.4,someorg.fna\n"
    with patch("subprocess.run"), \
         patch("builtins.open", mock_open(read_data=csv_no_acc)):
        hits = run_sketch_screen(cfg, tmp_path / "n.fastq.gz")
    assert hits[0].taxon_id == ""


def test_sketch_screen_creates_output_dir(tmp_path):
    cfg = _cfg(tmp_path)
    nonhuman = tmp_path / "nonhuman.fastq.gz"
    with patch("subprocess.run"), \
         patch("builtins.open", mock_open(read_data=SOURMASH_CSV)):
        run_sketch_screen(cfg, nonhuman)
    assert (tmp_path / "sketch").is_dir()
