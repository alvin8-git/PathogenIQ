import importlib.util
from pathlib import Path

_spec = importlib.util.spec_from_file_location(
    "select_kitome",
    Path(__file__).resolve().parent.parent / "scripts" / "05_select_kitome_controls.py",
)
sk = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(sk)


def test_spike_fraction_basic():
    # 90 spike reads of 100 total -> 0.9 (a high-biomass run, excluded by threshold)
    counts = {"GCF_000006945.2": 90, "GCF_pseudomonas": 10}
    assert abs(sk.spike_fraction(counts, 100, ["GCF_000006945"]) - 0.9) < 1e-9


def test_spike_fraction_prefix_version_agnostic():
    # matches the accession regardless of the .N version suffix
    counts = {"GCF_000006945.3": 5}
    assert sk.spike_fraction(counts, 10, ["GCF_000006945"]) == 0.5


def test_spike_fraction_kitome_dominated():
    # a diluted run: spike is a small share -> low fraction -> kept as a control
    counts = {"GCF_000006945.2": 2, "GCF_pseudomonas": 198}
    assert sk.spike_fraction(counts, 200, ["GCF_000006945"]) == 0.01


def test_spike_fraction_zero_total():
    assert sk.spike_fraction({}, 0, ["GCF_000006945"]) == 0.0


def test_spike_fraction_multiple_spikes():
    counts = {"GCF_AAA.1": 30, "GCF_BBB.1": 20, "GCF_other": 50}
    assert sk.spike_fraction(counts, 100, ["GCF_AAA", "GCF_BBB"]) == 0.5
