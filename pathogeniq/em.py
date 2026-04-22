from dataclasses import dataclass

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
