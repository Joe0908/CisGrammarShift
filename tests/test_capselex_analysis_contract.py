from __future__ import annotations

from dataclasses import replace

import pytest

from cisgrammar.capselex_analysis_contract import (
    AnalysisContract,
    panel_result_conclusion,
    permutation_analysis_status,
)


def test_primary_contract_uses_frozen_gphn_pipeline() -> None:
    contract = AnalysisContract()
    contract.validate()
    assert contract.primary_chip_pipeline == "McGill_GPHN_only"

    with pytest.raises(ValueError, match="McGill GPHN"):
        replace(contract, primary_chip_pipeline="Toronto_GPZN").validate()


def test_permutation_status_requires_the_frozen_final_count() -> None:
    contract = AnalysisContract()
    assert permutation_analysis_status(contract.screening_grammar_permutations, contract) == "screening"
    assert permutation_analysis_status(contract.grammar_permutations, contract) == "final"
    assert permutation_analysis_status(500, contract) == "exploratory"

    with pytest.raises(ValueError, match="at least one permutation"):
        permutation_analysis_status(0, contract)


def test_failed_panel_is_reported_as_tf_specific() -> None:
    conclusion = panel_result_conclusion(1, 4)
    assert conclusion == {
        "scope": "TF-specific",
        "panel_gate": "not passed",
        "summary": "TF-specific result; panel gate not passed.",
    }


def test_passing_panel_is_not_labeled_tf_specific() -> None:
    conclusion = panel_result_conclusion(4, 4)
    assert conclusion["scope"] == "panel-supported"
    assert conclusion["panel_gate"] == "passed"
