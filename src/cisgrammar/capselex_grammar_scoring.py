from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from cisgrammar.capselex import reverse_complement
from cisgrammar.capselex_monomer_motifs import PWM, scan_pwm


@dataclass(frozen=True)
class GrammarProfile:
    profile_id: str
    focal_pwm: PWM
    partner_pwm: PWM
    gaps: tuple[int, ...]
    orientations: tuple[str, ...] = ("FF", "FR", "RF", "RR")
    mechanism: str = "composite_motif"

    def __post_init__(self) -> None:
        if not self.gaps:
            raise ValueError("grammar profile requires at least one gap")
        if any(gap < 0 for gap in self.gaps):
            raise ValueError("grammar gaps cannot be negative")


def _oriented_pwm(pwm: PWM, reverse: bool) -> PWM:
    if not reverse:
        return pwm
    probabilities = pwm.probabilities[::-1, ::-1]
    return PWM(pwm.name + "_RC", probabilities)


def scan_grammar(sequence: str, profile: GrammarProfile) -> float:
    best = float("-inf")
    sequence = sequence.upper()
    for orientation in profile.orientations:
        focal = _oriented_pwm(profile.focal_pwm, orientation[0] == "R")
        partner = _oriented_pwm(profile.partner_pwm, orientation[1] == "R")
        focal_scores = scan_pwm(sequence, focal, double_stranded=False)
        partner_scores = scan_pwm(sequence, partner, double_stranded=False)
        for gap in profile.gaps:
            offset = focal.length + gap
            length = min(focal_scores.size, max(0, partner_scores.size - offset))
            if length:
                paired_scores = focal_scores[:length] + partner_scores[offset : offset + length]
                best = max(best, float(np.max(paired_scores)))
    return best


def raw_grammar_matrix(sequences: list[str], profiles: list[GrammarProfile]) -> np.ndarray:
    return np.array([[scan_grammar(sequence, profile) for profile in profiles] for sequence in sequences])


def binned_monomer_residuals(
    grammar: np.ndarray,
    focal_scores: np.ndarray,
    partner_scores: np.ndarray,
    train_mask: np.ndarray,
    bins: int = 20,
) -> np.ndarray:
    """Residualize grammar in training-derived focal/partner monomer quantile bins."""
    grammar = np.asarray(grammar, dtype=float)
    focal = np.asarray(focal_scores, dtype=float)
    partner = np.asarray(partner_scores, dtype=float)
    train = np.asarray(train_mask, dtype=bool)
    if not (grammar.shape == focal.shape == partner.shape == train.shape):
        raise ValueError("residualization inputs must have identical shape")
    focal_edges = np.unique(np.quantile(focal[train], np.linspace(0, 1, bins + 1)))
    partner_edges = np.unique(np.quantile(partner[train], np.linspace(0, 1, bins + 1)))
    focal_bin = np.clip(np.digitize(focal, focal_edges[1:-1]), 0, bins - 1)
    partner_bin = np.clip(np.digitize(partner, partner_edges[1:-1]), 0, bins - 1)
    residual = np.empty_like(grammar)
    global_mean = float(np.mean(grammar[train]))
    means: dict[tuple[int, int], float] = {}
    for first, second in zip(focal_bin[train], partner_bin[train], strict=True):
        key = (int(first), int(second))
        if key not in means:
            mask = train & (focal_bin == first) & (partner_bin == second)
            means[key] = float(np.mean(grammar[mask]))
    for index, (first, second) in enumerate(zip(focal_bin, partner_bin, strict=True)):
        residual[index] = grammar[index] - means.get((int(first), int(second)), global_mean)
    return residual


def strand_invariant_sequence(sequence: str) -> str:
    rc = reverse_complement(sequence)
    return min(sequence.upper(), rc)
