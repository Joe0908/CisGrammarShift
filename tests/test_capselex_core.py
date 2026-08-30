from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from cisgrammar.capselex import chromosome_fold, normalize_tf_symbol, reverse_complement
from cisgrammar.capselex_analysis_contract import AnalysisContract, assert_primary_asset_allowed
from cisgrammar.capselex_dataset import assign_split, build_directed_pair_dataset, summarize_dataset
from cisgrammar.capselex_pair_baselines import protein_kmer_vector


def test_aliases_and_reverse_complement() -> None:
    assert normalize_tf_symbol("RFXDC2") == "RFX7"
    assert normalize_tf_symbol("TP73L") == "TP63"
    assert reverse_complement("ACGTN") == "NACGT"


def test_chromosome_fold_is_deterministic() -> None:
    assert chromosome_fold("chr1") == 0
    assert chromosome_fold("chr6") == 0
    with pytest.raises(ValueError):
        chromosome_fold("chrX")


def test_analysis_contract_blocks_circular_assets() -> None:
    AnalysisContract().validate()
    with pytest.raises(ValueError, match="circular"):
        assert_primary_asset_allowed("TOPs/example.bed")


def test_directed_dataset_preserves_asymmetry() -> None:
    interactions = pd.DataFrame(
        {
            "bait": ["A", "B"],
            "prey": ["B", "A"],
            "cooperative_signal": [1, 0],
            "composite_motif": [1, 0],
            "spacing_or_orientation": [0, 1],
        }
    )
    constructs = pd.DataFrame(
        {"tf": ["A", "B"], "family": ["x", "y"], "protein_sequence": ["AAAA", "CCCC"]}
    )
    dataset = build_directed_pair_dataset(interactions, constructs)
    assert dataset["directed_pair"].tolist() == ["A->B", "B->A"]
    assert summarize_dataset(dataset).directed_asymmetries == 1


def test_pair_split_keeps_duplicate_pair_together() -> None:
    frame = pd.DataFrame({"directed_pair": ["A->B", "A->B", "B->A"], "family_pair": ["x", "x", "y"]})
    split = assign_split(frame, "pair", seed=3, heldout_fraction=0.5)
    assert split.iloc[0] == split.iloc[1]


def test_protein_kmer_is_normalized() -> None:
    vector = protein_kmer_vector("ACAC", k=2)
    assert np.isclose(vector.sum(), 1.0)
    assert np.count_nonzero(vector) == 2
