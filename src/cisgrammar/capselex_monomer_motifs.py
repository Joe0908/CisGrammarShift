from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from cisgrammar.capselex import reverse_complement

DNA = "ACGT"
DNA_INDEX = {base: index for index, base in enumerate(DNA)}


@dataclass(frozen=True)
class PWM:
    name: str
    probabilities: np.ndarray

    def __post_init__(self) -> None:
        matrix = np.asarray(self.probabilities, dtype=float)
        if matrix.ndim != 2 or matrix.shape[1] != 4:
            raise ValueError("PWM must have shape length x 4 in A,C,G,T order")
        if np.any(matrix < 0):
            raise ValueError("PWM probabilities cannot be negative")
        normalized = (matrix + 1e-6) / (matrix + 1e-6).sum(axis=1, keepdims=True)
        object.__setattr__(self, "probabilities", normalized)

    @property
    def length(self) -> int:
        return self.probabilities.shape[0]

    @property
    def log_odds(self) -> np.ndarray:
        return np.log2(self.probabilities / 0.25)


def read_jaspar(path: str | Path, name: str | None = None) -> PWM:
    rows: dict[str, list[float]] = {}
    motif_name = name or Path(path).stem
    with Path(path).open(encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if line.startswith(">") and name is None:
                motif_name = line[1:].strip().split(maxsplit=1)[-1]
            elif line and line[0] in DNA:
                values = line.split("[", 1)[-1].split("]", 1)[0].replace(",", " ")
                rows[line[0]] = [float(value) for value in values.split()]
    if set(rows) != set(DNA):
        raise ValueError(f"JASPAR motif must contain rows {DNA}")
    return PWM(motif_name, np.array([rows[base] for base in DNA]).T)


def score_window(sequence: str, pwm: PWM) -> float:
    if len(sequence) != pwm.length or any(base not in DNA_INDEX for base in sequence):
        return float("-inf")
    return float(sum(pwm.log_odds[offset, DNA_INDEX[base]] for offset, base in enumerate(sequence)))


def scan_pwm(sequence: str, pwm: PWM, double_stranded: bool = True) -> np.ndarray:
    sequence = sequence.upper()
    if len(sequence) < pwm.length:
        return np.empty(0, dtype=float)
    offsets = range(len(sequence) - pwm.length + 1)
    scores = np.array([score_window(sequence[offset : offset + pwm.length], pwm) for offset in offsets])
    if not double_stranded:
        return scores
    reverse_scores = np.array(
        [
            score_window(reverse_complement(sequence[offset : offset + pwm.length]), pwm)
            for offset in range(len(sequence) - pwm.length + 1)
        ]
    )
    return np.maximum(scores, reverse_scores)


def maximum_pwm_score(sequence: str, pwm: PWM) -> float:
    scores = scan_pwm(sequence, pwm)
    return float(scores.max()) if scores.size else float("-inf")
