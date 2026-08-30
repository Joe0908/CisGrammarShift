from __future__ import annotations

import numpy as np
from sklearn.metrics import (
    average_precision_score,
    balanced_accuracy_score,
    brier_score_loss,
    matthews_corrcoef,
    roc_auc_score,
)


def select_threshold(y_true: np.ndarray, probabilities: np.ndarray) -> float:
    candidates = np.unique(np.concatenate([[0.0], probabilities, [1.0]]))
    scores = np.array(
        [balanced_accuracy_score(y_true, probabilities >= threshold) for threshold in candidates]
    )
    best = np.flatnonzero(scores == scores.max())
    tied = candidates[best]
    return float(tied[np.argmin(np.abs(tied - 0.5))])


def expected_calibration_error(
    y_true: np.ndarray,
    probabilities: np.ndarray,
    n_bins: int = 10,
) -> float:
    edges = np.linspace(0.0, 1.0, n_bins + 1)
    total = len(y_true)
    error = 0.0
    for index in range(n_bins):
        if index == n_bins - 1:
            selected = (probabilities >= edges[index]) & (probabilities <= edges[index + 1])
        else:
            selected = (probabilities >= edges[index]) & (probabilities < edges[index + 1])
        if not np.any(selected):
            continue
        observed = float(np.mean(y_true[selected]))
        predicted = float(np.mean(probabilities[selected]))
        error += float(np.sum(selected)) / total * abs(observed - predicted)
    return error


def binary_metrics(
    y_true: np.ndarray,
    probabilities: np.ndarray,
    threshold: float,
) -> dict[str, float]:
    y_true = np.asarray(y_true, dtype=int)
    probabilities = np.asarray(probabilities, dtype=float)
    if set(np.unique(y_true)) != {0, 1}:
        raise ValueError("binary metrics require both classes")
    predictions = probabilities >= threshold
    return {
        "auroc": float(roc_auc_score(y_true, probabilities)),
        "auprc": float(average_precision_score(y_true, probabilities)),
        "balanced_accuracy": float(balanced_accuracy_score(y_true, predictions)),
        "mcc": float(matthews_corrcoef(y_true, predictions)),
        "brier": float(brier_score_loss(y_true, probabilities)),
        "ece": expected_calibration_error(y_true, probabilities),
    }


def paired_metrics(
    y_true: np.ndarray,
    probabilities: np.ndarray,
    pair_ids: np.ndarray,
) -> dict[str, float]:
    deltas: list[float] = []
    correct: list[float] = []
    for pair_id in np.unique(pair_ids):
        selected = pair_ids == pair_id
        pair_y = y_true[selected]
        pair_p = probabilities[selected]
        if len(pair_y) != 2 or set(pair_y.astype(int)) != {0, 1}:
            raise ValueError(f"pair {pair_id!r} is not one-positive/one-negative")
        positive = float(pair_p[pair_y == 1][0])
        negative = float(pair_p[pair_y == 0][0])
        delta = positive - negative
        deltas.append(delta)
        correct.append(float(delta > 0) + 0.5 * float(delta == 0))
    return {
        "pairwise_accuracy": float(np.mean(correct)),
        "counterfactual_delta": float(np.mean(deltas)),
        "counterfactual_delta_std": float(np.std(deltas, ddof=1)) if len(deltas) > 1 else 0.0,
    }
