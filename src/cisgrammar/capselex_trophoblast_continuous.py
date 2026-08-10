from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from cisgrammar.capselex import asset_manifest, chromosome_fold, write_json
from cisgrammar.capselex_genomic_model import NestedResult, chromosome_nested_continuous
from cisgrammar.capselex_trophoblast_deseq2 import read_deseq2_results, validate_deseq2_session_info

DEFAULT_BASELINE_PREFIXES = (
    "monomer_",
    "gcm1_chip",
    "ght_",
    "gc",
    "tss_distance",
    "wt_count",
    "direct",
    "loop_count",
    "loop_score",
    "linked_locus_count",
)


def infer_baseline_features(frame: pd.DataFrame) -> list[str]:
    features = [
        column
        for column in frame.columns
        if any(column == prefix or column.startswith(prefix) for prefix in DEFAULT_BASELINE_PREFIXES)
    ]
    if not features:
        raise ValueError("no baseline features were recognized")
    return features


def prepare_continuous_frame(
    linked_gene_features: pd.DataFrame,
    deseq2_results: pd.DataFrame,
    grammar_feature: str = "ets_grammar",
) -> pd.DataFrame:
    if grammar_feature not in linked_gene_features:
        raise ValueError(f"linked features need {grammar_feature}")
    frame = linked_gene_features.merge(
        deseq2_results[["gene", "gcm1_dependence", "baseMean", "padj"]], on="gene", how="inner"
    )
    frame = frame[np.isfinite(frame["gcm1_dependence"])].copy()
    frame["chromosome_fold"] = frame["chrom"].map(chromosome_fold)
    return frame


def run_continuous_state(
    linked_gene_features: pd.DataFrame,
    deseq2_results: pd.DataFrame,
    grammar_feature: str = "ets_grammar",
) -> NestedResult:
    frame = prepare_continuous_frame(linked_gene_features, deseq2_results, grammar_feature)
    baseline = infer_baseline_features(frame)
    return chromosome_nested_continuous(frame, "gcm1_dependence", baseline, [grammar_feature])


def run_trophoblast_deseq2_hichip_benchmark(
    evt_features_path: str | Path,
    st_features_path: str | Path,
    deseq2_directory: str | Path,
    output: str | Path,
    predictions: str | Path,
) -> dict[str, object]:
    directory = Path(deseq2_directory)
    software = validate_deseq2_session_info(directory / "sessionInfo.txt")
    state_paths = {"EVT": Path(evt_features_path), "ST": Path(st_features_path)}
    payload: dict[str, object] = {
        "schema_version": "gcm1_deseq2_hichip_continuous_benchmark_v1",
        "software": software,
        "design": {
            "outcome": "negative unshrunk KO-versus-WT DESeq2 log2FoldChange",
            "outcome_filter": "finite LFC only; no padj threshold",
            "outer_cv": "test chromosome fold i; validation i+1; remaining folds train",
            "selection_after_outcome_access": False,
        },
        "results": {},
    }
    prediction_frames = []
    assets: list[Path] = [directory / "sessionInfo.txt"]
    for state, feature_path in state_paths.items():
        features = pd.read_csv(feature_path, sep="\t")
        deseq_path = directory / f"{state}.deseq2.tsv"
        deseq = read_deseq2_results(deseq_path)
        prepared = prepare_continuous_frame(features, deseq)
        nested = chromosome_nested_continuous(
            prepared,
            "gcm1_dependence",
            infer_baseline_features(prepared),
            ["ets_grammar"],
        )
        state_predictions = prepared[["gene", "chrom", "locus_id"]].reset_index(drop=True)
        state_predictions = pd.concat([state_predictions, nested.predictions.reset_index(drop=True)], axis=1)
        state_predictions.insert(0, "state", state)
        prediction_frames.append(state_predictions)
        coefficient_folds = [
            fold["standardized_addition_coefficients"][0] for fold in nested.summary["fold_results"]
        ]
        payload["results"][state] = {
            **nested.summary,
            "finite_lfc_genes": int(np.isfinite(deseq["gcm1_dependence"]).sum()),
            "positive_coefficient_folds": int(np.sum(np.asarray(coefficient_folds) > 0)),
            "median_standardized_grammar_coefficient": float(np.median(coefficient_folds)),
        }
        assets.extend([feature_path, deseq_path])
    prediction_path = Path(predictions)
    prediction_path.parent.mkdir(parents=True, exist_ok=True)
    pd.concat(prediction_frames, ignore_index=True).to_csv(prediction_path, sep="\t", index=False)
    payload["assets"] = asset_manifest(assets)
    write_json(payload, output)
    return payload
