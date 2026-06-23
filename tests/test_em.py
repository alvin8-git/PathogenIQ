import numpy as np
import pytest
from pathogeniq.em import em_abundance, bootstrap_ci


def test_em_all_reads_one_organism():
    matrix = np.array([[1, 0], [1, 0], [1, 0]], dtype=float)
    result = em_abundance(matrix)
    assert result.abundances[0] == pytest.approx(1.0, abs=1e-6)
    assert result.abundances[1] == pytest.approx(0.0, abs=1e-6)


def test_em_equal_split():
    matrix = np.array([[1, 0], [1, 0], [1, 0], [0, 1], [0, 1], [0, 1]], dtype=float)
    result = em_abundance(matrix)
    assert result.abundances[0] == pytest.approx(0.5, abs=0.01)
    assert result.abundances[1] == pytest.approx(0.5, abs=0.01)


def test_em_multimapper_resolves_by_unique_evidence():
    # 6 unique to org 0, 2 unique to org 1, 2 multi-mappers
    matrix = np.array([
        [1, 0], [1, 0], [1, 0], [1, 0], [1, 0], [1, 0],
        [0, 1], [0, 1],
        [1, 1], [1, 1],
    ], dtype=float)
    result = em_abundance(matrix)
    assert result.abundances[0] > result.abundances[1]


def test_em_sums_to_one():
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


def test_bootstrap_ci_parallel_matches_serial():
    # per-iteration spawned seeds make the result independent of n_jobs / order,
    # so the multiprocessing path is bit-identical to the serial one.
    rng = np.random.default_rng(0)
    # >= 20k rows so the parallel path (Pool) is actually taken, not the serial guard
    matrix = (rng.random((25_000, 3)) > 0.4).astype(float)
    matrix[matrix.sum(axis=1) == 0, 0] = 1.0   # avoid all-zero rows
    s_lo, s_hi = bootstrap_ci(matrix, n_bootstrap=8, n_jobs=1)
    p_lo, p_hi = bootstrap_ci(matrix, n_bootstrap=8, n_jobs=4)   # exercises the Pool
    assert np.array_equal(s_lo, p_lo)
    assert np.array_equal(s_hi, p_hi)


def test_bootstrap_ci_empty_matrix():
    lower, upper = bootstrap_ci(np.zeros((0, 3)), n_bootstrap=10, n_jobs=2)
    assert lower.shape == (3,) and upper.shape == (3,)


def test_bootstrap_ci_coverage():
    matrix = np.array([[1, 0]] * 7 + [[0, 1]] * 3, dtype=float)
    result = em_abundance(matrix)
    lower, upper = bootstrap_ci(matrix, n_bootstrap=200)
    for i in range(2):
        assert lower[i] <= result.abundances[i] + 1e-6
        assert upper[i] >= result.abundances[i] - 1e-6
