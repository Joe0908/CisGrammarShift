from __future__ import annotations

import numpy as np
import pandas as pd

from cisgrammar.capselex_monomer_motifs import PWM, MonomerPWMRecord
from cisgrammar.capselex_nature_supplements import CapPWMRecord
from cisgrammar.capselex_primary_features import (
    FrozenPWMPanel,
    PairMechanismProfile,
    add_expression_supported_aggregates,
    score_frozen_pwm_panel,
    sequence_covariates,
)


def deterministic_pwm(name: str, motif: str) -> PWM:
    matrix = np.full((len(motif), 4), 0.01)
    index = {base: position for position, base in enumerate("ACGT")}
    for offset, base in enumerate(motif):
        matrix[offset, index[base]] = 0.97
    return PWM(name, matrix)


def cap_record(pair: str, motif: str, mechanism: str) -> CapPWMRecord:
    pwm = deterministic_pwm(pair, motif)
    return CapPWMRecord(pair, pair, "NA", "B", motif, 1, "C3", mechanism, True, pwm.probabilities, 1)


def test_sequence_covariates_are_explicit_about_n() -> None:
    result = sequence_covariates(["ACGTCGNN"])
    assert result.loc[0, "gc_fraction"] == 4 / 6
    assert result.loc[0, "cpg_per_100bp"] == 200 / 7
    assert result.loc[0, "n_fraction"] == 0.25


def test_frozen_panel_scoring_aggregates_pair_mechanisms() -> None:
    focal = MonomerPWMRecord("focal", "FOCAL", deterministic_pwm("focal", "AAA"))
    partner = MonomerPWMRecord("partner", "PARTNER", deterministic_pwm("partner", "CCC"))
    profile = PairMechanismProfile(
        "FOCAL_PARTNER",
        "PARTNER",
        "composite",
        (cap_record("FOCAL_PARTNER", "AAACCC", "composite"),),
        (partner,),
        "test",
    )
    panel = FrozenPWMPanel("FOCAL", focal, (profile,), ())
    sequences = ["AAACCC" + "T" * 14] * 25 + ["AAA" + "T" * 17] * 25
    features, manifest = score_frozen_pwm_panel(
        sequences,
        panel,
        chromosome_folds=np.tile(np.arange(5), 10),
        residual_bins=2,
    )
    assert "cap_composite_excess_mean" in features
    assert "cap_all_excess_mean" in features
    assert np.isfinite(features.to_numpy()).all()
    assert manifest["pair_mechanism_units"] == 1


def test_all_mechanism_score_does_not_double_weight_pair_with_two_mechanisms() -> None:
    focal = MonomerPWMRecord("focal", "FOCAL", deterministic_pwm("focal", "AAA"))
    first_partner = MonomerPWMRecord(
        "first_partner", "FIRST", deterministic_pwm("first_partner", "CCC")
    )
    second_partner = MonomerPWMRecord(
        "second_partner", "SECOND", deterministic_pwm("second_partner", "GGG")
    )
    profiles = (
        PairMechanismProfile(
            "FOCAL_FIRST",
            "FIRST",
            "composite",
            (cap_record("FOCAL_FIRST", "AAACCC", "composite"),),
            (first_partner,),
            "test",
        ),
        PairMechanismProfile(
            "FOCAL_FIRST",
            "FIRST",
            "spacing",
            (cap_record("FOCAL_FIRST", "CCCAAA", "spacing"),),
            (first_partner,),
            "test",
        ),
        PairMechanismProfile(
            "FOCAL_SECOND",
            "SECOND",
            "composite",
            (cap_record("FOCAL_SECOND", "AAAGGG", "composite"),),
            (second_partner,),
            "test",
        ),
    )
    panel = FrozenPWMPanel("FOCAL", focal, profiles, ())
    sequences = (["AAACCC" + "T" * 14] * 25) + (["AAAGGG" + "T" * 14] * 25)
    features, _ = score_frozen_pwm_panel(
        sequences,
        panel,
        chromosome_folds=np.tile(np.arange(5), 10),
        residual_bins=2,
    )
    first_pair = features[
        ["cap_excess::composite::FOCAL_FIRST", "cap_excess::spacing::FOCAL_FIRST"]
    ].mean(axis=1)
    second_pair = features["cap_excess::composite::FOCAL_SECOND"]
    np.testing.assert_allclose(features["cap_all_excess_mean"], (first_pair + second_pair) / 2)


def test_expression_supported_aggregation_excludes_unsupported_pairs() -> None:
    frame = pd.DataFrame(
        {
            "cap_excess::composite::A_B": [2.0, 4.0],
            "cap_excess::spacing::A_B": [4.0, 8.0],
            "cap_excess::composite::A_C": [100.0, 100.0],
            "partner_raw::A_B": [1.0, 2.0],
            "partner_raw::A_C": [3.0, 4.0],
        }
    )
    expression = pd.DataFrame(
        {
            "focal_tf": ["A", "A"],
            "pair": ["A_B", "A_C"],
            "expression_supported": ["True", "False"],
        }
    )
    result, summary = add_expression_supported_aggregates(frame, "A", expression)
    np.testing.assert_allclose(result["cap_expression_supported_excess_mean"], [3.0, 6.0])
    assert result["partner_monomer_expressed_max"].tolist() == [1.0, 2.0]
    assert summary["expression_supported_pairs"] == ["A_B"]
