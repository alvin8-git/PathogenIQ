from pathogeniq.coverage import coverage_from_paf, _merged_covered


def _paf_line(qn, ts, te, tlen=1_000_000, tname="ref", mapq=60):
    # qname qlen qs qe strand tname tlen tstart tend nmatch alen mapq
    alen = te - ts
    return f"{qn}\t{alen}\t0\t{alen}\t+\t{tname}\t{tlen}\t{ts}\t{te}\t{alen}\t{alen}\t{mapq}"


def test_merged_covered_unions_overlaps():
    assert _merged_covered([(0, 100), (50, 150), (200, 250)]) == 200  # 150 + 50


def test_clumped_reads_low_breadth_ratio():
    # 50 reads all piled at the same 150bp locus -> artifact
    paf = "\n".join(_paf_line(f"r{i}", 1000, 1150) for i in range(50))
    cov = coverage_from_paf(paf)
    assert cov.covered_bases == 150          # union of identical intervals
    assert cov.n_alignments == 50
    assert cov.breadth_ratio < 0.1           # far below uniform expectation


def test_spread_reads_high_breadth_ratio():
    # 50 reads at distinct, non-overlapping positions across the genome -> real
    paf = "\n".join(_paf_line(f"r{i}", i * 20000, i * 20000 + 150) for i in range(50))
    cov = coverage_from_paf(paf)
    assert cov.covered_bases == 50 * 150
    assert cov.breadth_ratio > 0.9           # matches uniform expectation


def test_same_count_opposite_verdict():
    # identical read COUNT, opposite breadth verdict -> the whole point
    clumped = coverage_from_paf("\n".join(_paf_line(f"c{i}", 5000, 5150) for i in range(50)))
    spread = coverage_from_paf("\n".join(_paf_line(f"s{i}", i * 20000, i * 20000 + 150) for i in range(50)))
    assert clumped.n_alignments == spread.n_alignments == 50
    assert clumped.breadth_ratio < 0.1 < spread.breadth_ratio


def test_skips_unmapped():
    paf = _paf_line("r1", 0, 150, mapq=0) + "\n" + _paf_line("r2", 200, 350, mapq=60)
    cov = coverage_from_paf(paf)
    assert cov.n_alignments == 1             # mapq=0 dropped


def test_multi_contig_aggregates():
    paf = _paf_line("r1", 0, 150, tlen=1000, tname="contig_a") + "\n" + \
          _paf_line("r2", 0, 150, tlen=2000, tname="contig_b")
    cov = coverage_from_paf(paf)
    assert cov.genome_length == 3000         # contigs summed
    assert cov.covered_bases == 300
