from pathogeniq.background import (
    BackgroundModel,
    build_background,
    background_pvalue,
    is_background,
    load_background_table,
)


def test_load_background_table_parses_rates_and_tier(tmp_path):
    table = tmp_path / "bg.tsv"
    table.write_text("# tier=2\ntaxon_id\trpm\nGCF_000001405.40\t12.3\nGCF_000007545.1\t4.5\n")
    model = load_background_table(table)
    assert model.tier == 2
    assert abs(model.rates["GCF_000001405.40"] - 12.3) < 1e-6
    assert abs(model.rates["GCF_000007545.1"] - 4.5) < 1e-6


def test_load_background_table_defaults_tier_2(tmp_path):
    table = tmp_path / "bg.tsv"
    table.write_text("taxon_id\trpm\nGCF_x\t1.0\n")
    assert load_background_table(table).tier == 2


def test_load_background_table_tier_1_header(tmp_path):
    table = tmp_path / "bg.tsv"
    table.write_text("# tier=1\ntaxon_id\trpm\nGCF_x\t1.0\n")
    assert load_background_table(table).tier == 1


def test_build_background_rate_is_mean_rpm_plus_pseudocount():
    # one NTC: 100 reads of Cutibacterium out of 1,000,000 -> 100 RPM + 0.5 pseudocount
    model = build_background([({"GCF_cuti": 100}, 1_000_000)], tier=1)
    assert model is not None
    assert abs(model.rates["GCF_cuti"] - 100.5) < 1e-6
    assert model.n_controls == 1
    assert model.tier == 1


def test_build_background_averages_across_controls():
    # two NTCs at the same depth: 100 and 300 RPM -> mean 200 + 0.5
    model = build_background(
        [({"GCF_cuti": 100}, 1_000_000), ({"GCF_cuti": 300}, 1_000_000)],
        tier=1,
    )
    assert abs(model.rates["GCF_cuti"] - 200.5) < 1e-6
    assert model.n_controls == 2


def test_build_background_normalizes_by_depth():
    # same count, different depth -> different RPM (depth-invariance via RPM)
    shallow = build_background([({"GCF_x": 10}, 100_000)], tier=2)   # 100 RPM
    deep = build_background([({"GCF_x": 10}, 10_000_000)], tier=2)   # 1 RPM
    assert shallow.rates["GCF_x"] > deep.rates["GCF_x"]


def test_build_background_none_when_all_controls_empty():
    # mandatory edge case: NTC with zero classified reads -> no model, caller flags
    assert build_background([({}, 0), ({"GCF_x": 5}, 0)], tier=3) is None


def test_real_signal_over_low_background_is_not_background():
    # taxon absent from NTC, strong sample signal -> tiny p, passes the filter
    model = build_background([({"GCF_bg": 50}, 1_000_000)], tier=1)
    p = background_pvalue("GCF_pathogen", sample_count=5000, sample_total=1_000_000, model=model)
    assert p < 0.01
    assert is_background("GCF_pathogen", 5000, 1_000_000, model) is False


def test_count_matching_background_is_background():
    # sample count ~ NTC background level -> high p, filtered out
    model = build_background([({"GCF_bg": 500}, 1_000_000)], tier=1)
    p = background_pvalue("GCF_bg", sample_count=505, sample_total=1_000_000, model=model)
    assert p >= 0.01
    assert is_background("GCF_bg", 505, 1_000_000, model) is True


def test_absent_taxon_uses_pseudocount_floor():
    # a taxon not in the NTC still gets a finite (pseudocount) background, so a
    # modest real signal clears it rather than dividing by zero
    model = build_background([({"GCF_bg": 50}, 1_000_000)], tier=2)
    assert "GCF_novel" not in model.rates
    assert background_pvalue("GCF_novel", 200, 1_000_000, model) < 0.01


def test_pvalue_monotonic_in_sample_count():
    # more reads of the same taxon -> stronger evidence -> non-increasing p
    model = build_background([({"GCF_bg": 100}, 1_000_000)], tier=1)
    p_low = background_pvalue("GCF_bg", 120, 1_000_000, model)
    p_high = background_pvalue("GCF_bg", 5000, 1_000_000, model)
    assert p_high <= p_low


def test_zero_sample_count_is_background():
    model = build_background([({"GCF_bg": 100}, 1_000_000)], tier=1)
    assert background_pvalue("GCF_bg", 0, 1_000_000, model) == 1.0
    assert is_background("GCF_bg", 0, 1_000_000, model) is True


def test_alpha_threshold_controls_decision():
    model = build_background([({"GCF_bg": 300}, 1_000_000)], tier=1)
    # a borderline count: stricter alpha keeps it, looser alpha drops it
    count = 700
    p = background_pvalue("GCF_bg", count, 1_000_000, model)
    # is_background is (p >= alpha): a smaller alpha is satisfied, a larger one is not
    assert is_background("GCF_bg", count, 1_000_000, model, alpha=p - 1e-9) is True
    assert is_background("GCF_bg", count, 1_000_000, model, alpha=p + 1e-9) is False
