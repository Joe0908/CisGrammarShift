from __future__ import annotations

import numpy as np

from cisgrammar.motifs import Motif, reverse_complement_pwm


def _scan_best_log_odds(x: np.ndarray, pwm: np.ndarray, background: np.ndarray) -> np.ndarray:
    log_odds = np.log(np.clip(pwm, 1e-8, 1.0) / background[None, :])
    reverse = np.log(np.clip(reverse_complement_pwm(pwm), 1e-8, 1.0) / background[None, :])
    n_positions = x.shape[1] - pwm.shape[0] + 1
    if n_positions <= 0:
        raise ValueError("motif is longer than the input sequence")
    best = np.full(x.shape[0], -np.inf, dtype=np.float64)
    for start in range(n_positions):
        window = x[:, start : start + pwm.shape[0], :]
        forward_score = np.einsum("nlc,lc->n", window, log_odds)
        reverse_score = np.einsum("nlc,lc->n", window, reverse)
        best = np.maximum(best, np.maximum(forward_score, reverse_score))
    return best


class PWMPresenceBaseline:
    """Scores motif co-presence without representing their relative spacing."""

    def __init__(self, motif_a: Motif, motif_b: Motif, background_gc: float = 0.5):
        self.motif_a = motif_a
        self.motif_b = motif_b
        self.background = np.array(
            [(1 - background_gc) / 2, background_gc / 2, background_gc / 2, (1 - background_gc) / 2]
        )
        self.location_: float | None = None
        self.scale_: float | None = None

    def decision_function(self, x: np.ndarray) -> np.ndarray:
        score_a = _scan_best_log_odds(x, self.motif_a.pwm, self.background)
        score_b = _scan_best_log_odds(x, self.motif_b.pwm, self.background)
        return np.minimum(score_a, score_b)

    def fit(self, x: np.ndarray) -> PWMPresenceBaseline:
        scores = self.decision_function(x)
        self.location_ = float(np.median(scores))
        mad = float(np.median(np.abs(scores - self.location_)))
        self.scale_ = max(1.4826 * mad, 1e-6)
        return self

    def predict_proba(self, x: np.ndarray) -> np.ndarray:
        if self.location_ is None or self.scale_ is None:
            raise RuntimeError("fit the baseline before prediction")
        z = np.clip((self.decision_function(x) - self.location_) / self.scale_, -30, 30)
        return 1.0 / (1.0 + np.exp(-z))
