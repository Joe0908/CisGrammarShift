from __future__ import annotations

from collections import defaultdict
from collections.abc import Sequence

import numpy as np

from cisgrammar.capselex_monomer_motifs import PWM

DNA_INDEX = np.full(256, 4, dtype=np.uint8)
for _index, _base in enumerate(b"ACGT"):
    DNA_INDEX[_base] = _index


def encode_dna(sequences: Sequence[str]) -> np.ndarray:
    """Encode equal-length DNA strings as A=0, C=1, G=2, T=3, other=4."""
    if not sequences:
        return np.empty((0, 0), dtype=np.uint8)
    width = len(sequences[0])
    if any(len(sequence) != width for sequence in sequences):
        raise ValueError("all sequences must have equal length")
    encoded = np.empty((len(sequences), width), dtype=np.uint8)
    for row, sequence in enumerate(sequences):
        raw = np.frombuffer(sequence.upper().encode("ascii"), dtype=np.uint8)
        encoded[row] = DNA_INDEX[raw]
    return encoded


def _reverse_complement_log_odds(pwm: PWM) -> np.ndarray:
    return pwm.log_odds[::-1, ::-1]


def _score_equal_length_pwms(encoded: np.ndarray, pwms: Sequence[PWM]) -> np.ndarray:
    if not pwms:
        return np.empty((encoded.shape[0], 0), dtype=float)
    width = pwms[0].length
    if any(pwm.length != width for pwm in pwms):
        raise ValueError("PWM group must have one common length")
    windows = encoded.shape[1] - width + 1
    if windows <= 0:
        return np.full((encoded.shape[0], len(pwms)), -np.inf, dtype=float)

    forward = np.stack([pwm.log_odds for pwm in pwms]).astype(np.float32, copy=False)
    reverse = np.stack([_reverse_complement_log_odds(pwm) for pwm in pwms]).astype(
        np.float32, copy=False
    )
    batch = encoded.shape[0]
    scores_forward = np.zeros((batch, len(pwms), windows), dtype=np.float32)
    scores_reverse = np.zeros_like(scores_forward)
    invalid_counts = np.zeros((batch, windows), dtype=np.uint16)

    for offset in range(width):
        codes = encoded[:, offset : offset + windows]
        invalid_counts += codes >= 4
        safe_codes = np.minimum(codes, 3)
        scores_forward += np.moveaxis(forward[:, offset, safe_codes], 0, 1)
        scores_reverse += np.moveaxis(reverse[:, offset, safe_codes], 0, 1)

    invalid = invalid_counts > 0
    scores_forward = np.maximum(scores_forward, scores_reverse)
    scores_forward[np.broadcast_to(invalid[:, None, :], scores_forward.shape)] = -np.inf
    return scores_forward.max(axis=2).astype(float)


def maximum_pwm_scores(
    sequences: Sequence[str],
    pwms: Sequence[PWM],
    batch_size: int = 512,
) -> np.ndarray:
    """Return exact double-stranded maximum log-odds scores for many PWMs.

    The implementation groups motifs by length and uses NumPy batches. Windows that
    contain a non-ACGT base are excluded, matching ``maximum_pwm_score`` exactly.
    Output columns preserve the input PWM order.
    """
    if batch_size < 1:
        raise ValueError("batch_size must be positive")
    if not pwms:
        return np.empty((len(sequences), 0), dtype=float)
    if not sequences:
        return np.empty((0, len(pwms)), dtype=float)

    grouped: dict[int, list[tuple[int, PWM]]] = defaultdict(list)
    for index, pwm in enumerate(pwms):
        grouped[pwm.length].append((index, pwm))
    output = np.full((len(sequences), len(pwms)), -np.inf, dtype=float)
    for start in range(0, len(sequences), batch_size):
        stop = min(start + batch_size, len(sequences))
        encoded = encode_dna(sequences[start:stop])
        for records in grouped.values():
            indices = [index for index, _ in records]
            values = _score_equal_length_pwms(encoded, [pwm for _, pwm in records])
            output[start:stop, indices] = values
    return output
