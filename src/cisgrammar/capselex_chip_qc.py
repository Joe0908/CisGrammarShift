from __future__ import annotations

import numpy as np
from scipy.stats import spearmanr


def replicate_qc(first: np.ndarray, second: np.ndarray) -> dict[str, float]:
    first = np.asarray(first, dtype=float)
    second = np.asarray(second, dtype=float)
    if first.shape != second.shape:
        raise ValueError("replicates must have the same shape")
    finite = np.isfinite(first) & np.isfinite(second)
    if finite.sum() < 3:
        raise ValueError("at least three finite paired values are required")
    log_first = np.log1p(np.clip(first[finite], 0, None))
    log_second = np.log1p(np.clip(second[finite], 0, None))
    return {
        "paired_loci": int(finite.sum()),
        "spearman_log1p": float(spearmanr(log_first, log_second).statistic),
        "first_zero_fraction": float(np.mean(first[finite] == 0)),
        "second_zero_fraction": float(np.mean(second[finite] == 0)),
        "median_absolute_log1p_difference": float(np.median(np.abs(log_first - log_second))),
    }
