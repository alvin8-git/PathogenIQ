import json

import numpy as np

from pathogeniq.config import PipelineConfig, ReadType, SpecimenType
from pathogeniq.em import EMResult
from pathogeniq.quantify import absolute_copies, quantify_entries, SpikeInfo
from pathogeniq.report import ReportEntry, write_report


def _entry(name, taxon_id, reads):
    return ReportEntry(
        organism=name, abundance=0.0, ci_lower=0.0, ci_upper=0.0,
        read_count=reads, specimen_type=SpecimenType.BLOOD, taxon_id=taxon_id,
    )


def test_absolute_copies_formula():
    # 50 target reads / 100 spike reads * 1e6 spike copies = 500k copies
    assert absolute_copies(50, 100, 1_000_000) == 500_000


def test_absolute_copies_none_without_spike_reads():
    assert absolute_copies(50, 0, 1_000_000) is None


def test_quantify_removes_spike_and_anchors_load():
    entries = [_entry("Spike org", "GCF_spike", 100), _entry("E. coli", "GCF_ecoli", 50)]
    kept, info = quantify_entries(entries, spike_taxon_id="GCF_spike", spike_copies=1_000_000)
    assert [e.organism for e in kept] == ["E. coli"]      # spike dropped from findings
    assert info.found and info.spike_reads == 100
    assert kept[0].absolute_copies == 500_000             # 50/100 * 1e6


def test_quantify_missing_spike_flags_unfound():
    entries = [_entry("E. coli", "GCF_ecoli", 50)]
    kept, info = quantify_entries(entries, spike_taxon_id="GCF_spike", spike_copies=1_000_000)
    assert not info.found and info.spike_reads == 0
    assert kept[0].absolute_copies is None                # no anchor -> no absolute load


def test_report_emits_absolute_load(tmp_path):
    cfg = PipelineConfig(
        input_fastq=tmp_path / "s.fq.gz", read_type=ReadType.SHORT,
        specimen_type=SpecimenType.BLOOD, output_dir=tmp_path,
        db_tier1=tmp_path / "db", host_reference=tmp_path / "ref",
    )
    e = _entry("E. coli", "GCF_ecoli", 50)
    e.absolute_copies = 500_000.0
    em = EMResult(abundances=np.array([1.0]), n_reads=100, n_organisms=1, iterations=1)
    info = SpikeInfo(taxon_id="GCF_spike", copies_added=1e6, spike_reads=100,
                     sample_volume=5.0, found=True)
    write_report(cfg, [e], em, spike_info=info)
    j = json.loads((tmp_path / "report" / "pathogeniq_report.json").read_text())
    assert j["spike_in"]["found"] is True and j["spike_in"]["taxon_id"] == "GCF_spike"
    f = j["findings"][0]
    assert f["absolute_copies"] == 500_000.0
    assert f["copies_per_volume"] == 100_000.0            # 500000 / 5 mL
