from __future__ import annotations

import numpy as np

from cisgrammar.capselex_monomer_motifs import PWM, maximum_pwm_score
from cisgrammar.capselex_pwm_batch import encode_dna, maximum_pwm_scores


def deterministic_pwm(name: str, motif: str) -> PWM:
    matrix = np.full((len(motif), 4), 0.01)
    index = {base: position for position, base in enumerate("ACGT")}
    for offset, base in enumerate(motif):
        matrix[offset, index[base]] = 0.97
    return PWM(name, matrix)


def test_encode_dna_marks_ambiguous_bases() -> None:
    encoded = encode_dna(["ACGTN", "acgt-"])
    assert encoded.tolist() == [[0, 1, 2, 3, 4], [0, 1, 2, 3, 4]]


def test_batch_scores_match_scalar_scanner() -> None:
    sequences = ["TTTGTTAAA", "NNNAACNNN", "CCCCCCCCC"]
    pwms = [deterministic_pwm("aac", "AAC"), deterministic_pwm("cccc", "CCCC")]
    observed = maximum_pwm_scores(sequences, pwms, batch_size=2)
    expected = np.array(
        [[maximum_pwm_score(sequence, pwm) for pwm in pwms] for sequence in sequences]
    )
    np.testing.assert_allclose(observed, expected)


def test_batch_scores_support_mixed_lengths_and_all_n() -> None:
    sequences = ["AACCGGTT", "NNNNNNNN"]
    pwms = [deterministic_pwm("aac", "AAC"), deterministic_pwm("ccgg", "CCGG")]
    scores = maximum_pwm_scores(sequences, pwms)
    assert np.isfinite(scores[0]).all()
    assert np.isneginf(scores[1]).all()
