from dataclasses import dataclass
from multiprocessing import Pool

import numpy as np


@dataclass
class EMResult:
    abundances: np.ndarray   # shape (n_organisms,), sums to 1
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


# Bootstrap workers (module-level so multiprocessing can pickle them). The
# alignment matrix is set once per worker via the initializer rather than pickled
# per task — it is large and read-only.
_BOOT_MATRIX: np.ndarray | None = None


def _boot_init(matrix: np.ndarray) -> None:
    global _BOOT_MATRIX
    _BOOT_MATRIX = matrix


def _boot_one(seed_seq: np.random.SeedSequence) -> np.ndarray:
    m = _BOOT_MATRIX
    rng = np.random.default_rng(seed_seq)
    idx = rng.integers(0, m.shape[0], size=m.shape[0])
    return em_abundance(m[idx]).abundances


def bootstrap_ci(
    alignment_matrix: np.ndarray,
    n_bootstrap: int = 100,
    alpha: float = 0.05,
    seed: int = 42,
    n_jobs: int = 1,
) -> tuple[np.ndarray, np.ndarray]:
    """Bootstrap (1-alpha) confidence interval on EM abundance estimates.

    The ``n_bootstrap`` iterations are independent (each resamples reads + runs EM),
    so with ``n_jobs > 1`` they run across processes. Each iteration draws from its
    own spawned RNG stream, so the result is **identical and reproducible**
    regardless of ``n_jobs`` or execution order — parallelism never changes the CI.
    """
    n_orgs = alignment_matrix.shape[1]
    if alignment_matrix.shape[0] == 0 or n_bootstrap <= 0:
        return np.zeros(n_orgs), np.zeros(n_orgs)
    seeds = np.random.SeedSequence(seed).spawn(n_bootstrap)

    # Only parallelize when the per-iteration EM work is large enough to amortize
    # process-spawn overhead — for a small (e.g. low-biomass clinical) matrix the
    # serial loop is faster. ponytail: 20k-read cutoff is a heuristic, tune if needed.
    use_parallel = (
        n_jobs and n_jobs > 1 and n_bootstrap > 1 and alignment_matrix.shape[0] >= 20_000
    )
    if use_parallel:
        with Pool(processes=min(n_jobs, n_bootstrap),
                  initializer=_boot_init, initargs=(alignment_matrix,)) as pool:
            estimates = np.array(pool.map(_boot_one, seeds))
    else:
        _boot_init(alignment_matrix)
        estimates = np.array([_boot_one(s) for s in seeds])

    lower = np.percentile(estimates, 100 * alpha / 2, axis=0)
    upper = np.percentile(estimates, 100 * (1 - alpha / 2), axis=0)
    return lower, upper
