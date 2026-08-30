from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from cisgrammar.capselex import asset_manifest, sha256_file, stable_seed, write_json
from cisgrammar.capselex_genomic_assets import add_cross_assembly_bigwig_signal
from cisgrammar.capselex_genomic_model import (
    chromosome_bootstrap_partial_r2,
    chromosome_nested_continuous,
)
from cisgrammar.capselex_primary_features import add_expression_supported_aggregates

BASELINE = [
    "ght_score",
    "focal_monomer_score",
    "partner_monomer_expressed_max",
    "gc_fraction",
    "cpg_per_100bp",
    "hek293_dnase_log1p",
]
OUTCOMES = ["chip_mean_log1p", "chip_rep1_log1p", "chip_rep2_log1p"]


def _validate_dnase(directory: Path, manifest_path: Path) -> dict[str, Path]:
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    records = {record["filename"]: record for record in manifest["assets"]}
    required = {"GSM2902639_HEK293_DNase.bw", "hg38ToHg19.over.chain.gz", "liftOver"}
    paths = {name: directory / name for name in required}
    for name, path in paths.items():
        record = records[name]
        if path.stat().st_size != record["size_bytes"] or sha256_file(path) != record["sha256"]:
            raise RuntimeError(f"DNase asset integrity check failed: {path}")
    return paths


def _fit(frame: pd.DataFrame, addition: str, outcome: str) -> dict[str, object]:
    result = chromosome_nested_continuous(frame, outcome, BASELINE, [addition])
    coefficients = [
        record["standardized_addition_coefficients"][0]
        for record in result.summary["fold_results"]
    ]
    bootstrap = chromosome_bootstrap_partial_r2(
        frame,
        result.predictions,
        replicates=2000,
        seed=stable_seed(
            "pair_decomposition", frame["focal_tf"].iloc[0], addition, outcome
        ),
    )
    return {
        "partial_r2": result.summary["partial_r2"],
        "m0": result.summary["m0"],
        "m1": result.summary["m1"],
        "standardized_addition_coefficients": coefficients,
        "median_standardized_addition_coefficient": float(np.median(coefficients)),
        "positive_coefficient_in_all_folds": bool(np.all(np.asarray(coefficients) > 0)),
        "partial_r2_chromosome_bootstrap": bootstrap,
    }


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Exploratory pair-level decomposition after the primary panel gate"
    )
    parser.add_argument("--feature-directory", type=Path, required=True)
    parser.add_argument("--focal-tfs", nargs="+", required=True)
    parser.add_argument("--partner-expression", type=Path, required=True)
    parser.add_argument("--dnase-directory", type=Path, required=True)
    parser.add_argument("--dnase-resolved-manifest", type=Path, required=True)
    parser.add_argument("--primary-screening", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    screening = json.loads(args.primary_screening.read_text(encoding="utf-8"))
    if screening["panel_go_gate_passed"]:
        raise ValueError("pair decomposition script is reserved for diagnosing a failed panel gate")
    expression = pd.read_csv(args.partner_expression, sep="\t")
    dnase = _validate_dnase(args.dnase_directory, args.dnase_resolved_manifest)
    feature_paths = []
    tf_results = []
    for focal_tf in args.focal_tfs:
        feature_path = args.feature_directory / f"{focal_tf}.features.tsv.gz"
        feature_paths.append(feature_path)
        frame = pd.read_csv(feature_path, sep="\t")
        frame, expression_summary = add_expression_supported_aggregates(
            frame, focal_tf, expression
        )
        frame = add_cross_assembly_bigwig_signal(
            frame,
            dnase["GSM2902639_HEK293_DNase.bw"],
            dnase["hg38ToHg19.over.chain.gz"],
            dnase["liftOver"],
            "hek293_dnase",
            "hg19",
        )
        frame = frame.loc[frame["hek293_dnase_mapped"]].reset_index(drop=True)
        frame["hek293_dnase_log1p"] = np.log1p(frame["hek293_dnase"])
        pairs = expression_summary["expression_supported_pairs"]
        pair_columns: dict[str, list[str]] = {}
        for pair in pairs:
            columns = [
                column
                for column in frame
                if column.startswith("cap_excess::") and column.endswith(f"::{pair}")
            ]
            pair_columns[pair] = columns
            frame[f"pair_excess::{pair}"] = frame[columns].mean(axis=1)

        pair_results = []
        for pair in pairs:
            addition = f"pair_excess::{pair}"
            alone = {outcome: _fit(frame, addition, outcome) for outcome in OUTCOMES}
            other_pairs = [f"pair_excess::{other}" for other in pairs if other != pair]
            if other_pairs:
                leave_column = f"aggregate_without::{pair}"
                frame[leave_column] = frame[other_pairs].mean(axis=1)
                leave_out = _fit(frame, leave_column, "chip_mean_log1p")
            else:
                leave_out = None
            mechanisms = sorted(
                {column.split("::", 2)[1] for column in pair_columns[pair]}
            )
            pair_results.append(
                {
                    "pair": pair,
                    "mechanisms": mechanisms,
                    "pair_alone": alone,
                    "aggregate_without_pair": leave_out,
                }
            )
        pair_results.sort(
            key=lambda record: record["pair_alone"]["chip_mean_log1p"]["partial_r2"],
            reverse=True,
        )
        tf_results.append(
            {
                "focal_tf": focal_tf,
                "rows": len(frame),
                "expression_supported_pairs": len(pairs),
                "pair_results_ranked_by_mean_outcome_partial_r2": pair_results,
            }
        )
        top = pair_results[0]
        print(
            f"{focal_tf}\ttop={top['pair']}\t"
            f"partial_r2={top['pair_alone']['chip_mean_log1p']['partial_r2']:.6g}",
            flush=True,
        )

    write_json(
        {
            "schema_version": "capselex_pair_decomposition_v1",
            "analysis_role": "exploratory diagnosis after failure of the frozen primary panel gate",
            "primary_screening_gate_passed": False,
            "multiple_testing_claim": (
                "none; ranks and chromosome-bootstrap intervals are descriptive and do not define "
                "partner discoveries"
            ),
            "baseline_features": BASELINE,
            "outcomes": OUTCOMES,
            "tf_results": tf_results,
            "source_assets": asset_manifest(
                [
                    *feature_paths,
                    args.partner_expression,
                    args.dnase_resolved_manifest,
                    args.primary_screening,
                ]
            ),
            "claim_boundary": (
                "A high pair-level predictive increment does not establish that the named partner "
                "occupies the locus or physically cooperates with the focal TF."
            ),
        },
        args.output,
    )


if __name__ == "__main__":
    main()
