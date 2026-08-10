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
    chromosome_shift_null,
)
from cisgrammar.capselex_primary_features import add_expression_supported_aggregates

BASELINE_WITHOUT_ACCESSIBILITY = [
    "ght_score",
    "focal_monomer_score",
    "partner_monomer_expressed_max",
    "gc_fraction",
    "cpg_per_100bp",
]
OUTCOMES = ["chip_mean_log1p", "chip_rep1_log1p", "chip_rep2_log1p"]


def _coefficient_summary(model_summary: dict[str, object]) -> dict[str, object]:
    coefficients = np.array(
        [record["standardized_addition_coefficients"] for record in model_summary["fold_results"]],
        dtype=float,
    )
    return {
        "fold_coefficients": coefficients.tolist(),
        "median": np.median(coefficients, axis=0).tolist(),
        "positive_in_all_folds": np.all(coefficients > 0, axis=0).tolist(),
    }


def _fit_outcomes(
    frame: pd.DataFrame, baseline: list[str], additions: list[str]
) -> dict[str, dict[str, object]]:
    payload = {}
    for outcome in OUTCOMES:
        result = chromosome_nested_continuous(frame, outcome, baseline, additions)
        bootstrap = chromosome_bootstrap_partial_r2(
            frame,
            result.predictions,
            replicates=2000,
            seed=stable_seed(
                "chromosome_bootstrap",
                frame["focal_tf"].iloc[0],
                outcome,
                *additions,
            ),
        )
        payload[outcome] = {
            **result.summary,
            "addition_coefficients": _coefficient_summary(result.summary),
            "partial_r2_chromosome_bootstrap": bootstrap,
        }
    return payload


def _permutation_test(
    frame: pd.DataFrame,
    baseline: list[str],
    addition: str,
    permutations: int,
    seed: int,
) -> dict[str, object]:
    observed = chromosome_nested_continuous(frame, "chip_mean_log1p", baseline, [addition])
    observed_partial_r2 = float(observed.summary["partial_r2"])
    null = []
    for permutation in range(permutations):
        shifted = chromosome_shift_null(
            frame,
            [addition],
            seed=stable_seed(seed, permutation, addition),
        )
        result = chromosome_nested_continuous(shifted, "chip_mean_log1p", baseline, [addition])
        null.append(float(result.summary["partial_r2"]))
    null_array = np.asarray(null)
    exceedances = int(np.sum(null_array >= observed_partial_r2))
    return {
        "permutations": permutations,
        "null": "independent circular shift within each focal-TF chromosome",
        "observed_partial_r2": observed_partial_r2,
        "exceedances": exceedances,
        "one_sided_pvalue": (exceedances + 1) / (permutations + 1),
        "null_quantiles": {
            "minimum": float(np.min(null_array)),
            "median": float(np.median(null_array)),
            "p95": float(np.quantile(null_array, 0.95)),
            "maximum": float(np.max(null_array)),
        },
    }


def _validate_dnase_assets(
    directory: Path, manifest_path: Path
) -> tuple[dict[str, object], dict[str, Path]]:
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    records = {record["filename"]: record for record in manifest["assets"]}
    required = {"GSM2902639_HEK293_DNase.bw", "hg38ToHg19.over.chain.gz", "liftOver"}
    if not required.issubset(records):
        raise ValueError(f"DNase manifest lacks: {sorted(required - set(records))}")
    paths = {name: directory / name for name in required}
    for name, path in paths.items():
        record = records[name]
        if path.stat().st_size != record["size_bytes"] or sha256_file(path) != record["sha256"]:
            raise RuntimeError(f"HEK293 DNase asset integrity check failed: {path}")
    return manifest, paths


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Run chromosome-held-out CAP grammar addition models for the focal panel"
    )
    parser.add_argument("--feature-directory", type=Path, required=True)
    parser.add_argument("--focal-tfs", nargs="+", required=True)
    parser.add_argument("--sensitivity-tfs", nargs="*", default=[])
    parser.add_argument("--availability-negative-controls", nargs="*", default=[])
    parser.add_argument("--partner-expression", type=Path, required=True)
    parser.add_argument("--dnase-directory", type=Path, required=True)
    parser.add_argument("--dnase-resolved-manifest", type=Path, required=True)
    parser.add_argument("--permutations", type=int, default=100)
    parser.add_argument("--partial-r2-threshold", type=float, default=0.005)
    parser.add_argument("--minimum-positive-focal-tfs", type=int, default=4)
    parser.add_argument("--seed", type=int, default=20260809)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.permutations < 1:
        raise ValueError("at least one permutation is required")

    dnase_manifest, dnase_paths = _validate_dnase_assets(
        args.dnase_directory, args.dnase_resolved_manifest
    )
    focal_set = set(args.focal_tfs)
    sensitivity_set = set(args.sensitivity_tfs)
    negative_set = set(args.availability_negative_controls)
    panel_sets = [focal_set, sensitivity_set, negative_set]
    if any(panel_sets[first] & panel_sets[second] for first in range(3) for second in range(first + 1, 3)):
        raise ValueError("primary, sensitivity and availability-negative panels must not overlap")
    expression = pd.read_csv(args.partner_expression, sep="\t")

    records = []
    feature_paths = []
    for focal_tf in [
        *args.focal_tfs,
        *args.sensitivity_tfs,
        *args.availability_negative_controls,
    ]:
        feature_path = args.feature_directory / f"{focal_tf}.features.tsv.gz"
        feature_paths.append(feature_path)
        frame = pd.read_csv(feature_path, sep="\t")
        if frame.empty or set(frame["focal_tf"]) != {focal_tf}:
            raise ValueError(f"invalid feature table for {focal_tf}")
        frame, expression_summary = add_expression_supported_aggregates(
            frame, focal_tf, expression
        )
        frame = add_cross_assembly_bigwig_signal(
            frame,
            dnase_paths["GSM2902639_HEK293_DNase.bw"],
            dnase_paths["hg38ToHg19.over.chain.gz"],
            dnase_paths["liftOver"],
            "hek293_dnase",
            "hg19",
        )
        mapped_rows = int(frame["hek293_dnase_mapped"].sum())
        model_frame = frame.loc[frame["hek293_dnase_mapped"]].reset_index(drop=True)
        if model_frame.empty:
            raise RuntimeError(f"no {focal_tf} loci mapped from hg38 to the hg19 DNase assembly")
        model_frame["hek293_dnase_log1p"] = np.log1p(model_frame["hek293_dnase"])
        baseline = [*BASELINE_WITHOUT_ACCESSIBILITY, "hek293_dnase_log1p"]
        required = [
            *OUTCOMES,
            "cap_all_excess_mean",
            "partner_monomer_max",
            "ght_score",
            "focal_monomer_score",
            "gc_fraction",
            "cpg_per_100bp",
        ]
        if expression_summary["expression_supported_pair_count"]:
            required.extend([*baseline, "cap_expression_supported_excess_mean"])
        missing = sorted(set(required) - set(model_frame))
        if missing:
            raise ValueError(f"{focal_tf} feature table is missing: {', '.join(missing)}")

        all_pair_baseline_without_accessibility = [
            "ght_score",
            "focal_monomer_score",
            "partner_monomer_max",
            "gc_fraction",
            "cpg_per_100bp",
        ]
        all_pair_baseline = [*all_pair_baseline_without_accessibility, "hek293_dnase_log1p"]
        all_pair_sensitivity = _fit_outcomes(
            model_frame, all_pair_baseline, ["cap_all_excess_mean"]
        )
        overall = None
        no_accessibility = None
        no_accessibility_full_locus = None
        branches = {}
        permutation = None
        pass_threshold = False
        criterion = None
        if expression_summary["expression_supported_pair_count"]:
            overall = _fit_outcomes(
                model_frame, baseline, ["cap_expression_supported_excess_mean"]
            )
            no_accessibility = _fit_outcomes(
                model_frame,
                BASELINE_WITHOUT_ACCESSIBILITY,
                ["cap_expression_supported_excess_mean"],
            )
            no_accessibility_full_locus = _fit_outcomes(
                frame,
                BASELINE_WITHOUT_ACCESSIBILITY,
                ["cap_expression_supported_excess_mean"],
            )
            for mechanism in ("composite", "spacing"):
                column = f"cap_expression_supported_{mechanism}_excess_mean"
                if column in frame:
                    branches[mechanism] = _fit_outcomes(model_frame, baseline, [column])
            permutation = _permutation_test(
                model_frame,
                baseline,
                "cap_expression_supported_excess_mean",
                args.permutations,
                stable_seed(args.seed, focal_tf),
            )
            mean_model = overall["chip_mean_log1p"]
            rep1_median = overall["chip_rep1_log1p"]["addition_coefficients"]["median"][0]
            rep2_median = overall["chip_rep2_log1p"]["addition_coefficients"]["median"][0]
            mean_median = mean_model["addition_coefficients"]["median"][0]
            pass_threshold = bool(
                mean_model["partial_r2"] >= args.partial_r2_threshold
                and permutation["one_sided_pvalue"] <= 0.05
                and mean_median > 0
                and rep1_median > 0
                and rep2_median > 0
            )
            criterion = {
                "partial_r2_at_least": args.partial_r2_threshold,
                "permutation_p_at_most": 0.05,
                "positive_median_standardized_coefficient_mean_outcome": mean_median > 0,
                "positive_median_standardized_coefficient_rep1": rep1_median > 0,
                "positive_median_standardized_coefficient_rep2": rep2_median > 0,
                "passed": pass_threshold,
            }
        if focal_tf in focal_set:
            panel_role = "primary_expression_supported"
        elif focal_tf in sensitivity_set:
            panel_role = "low_CAP_coverage_sensitivity"
        else:
            panel_role = "partner_availability_negative_control"
        record = {
            "focal_tf": focal_tf,
            "panel_role": panel_role,
            "rows": len(frame),
            "rows_with_hg38_to_hg19_dnase_mapping": mapped_rows,
            "dnase_mapping_fraction": mapped_rows / len(frame),
            "partner_expression": expression_summary,
            "baseline_features": baseline,
            "addition_feature": "cap_expression_supported_excess_mean",
            "primary_outcome_models": overall,
            "mechanism_branch_models": branches,
            "no_accessibility_sensitivity": no_accessibility,
            "no_accessibility_full_locus_sensitivity": (
                no_accessibility_full_locus
                if expression_summary["expression_supported_pair_count"]
                else None
            ),
            "all_pair_sensitivity": all_pair_sensitivity,
            "screening_permutation": permutation,
            "positive_focal_tf_criterion": criterion,
        }
        records.append(record)
        if overall is None:
            print(f"{focal_tf}\tno expression-supported CAP partner", flush=True)
        else:
            print(
                f"{focal_tf}\tpartial_r2={overall['chip_mean_log1p']['partial_r2']:.6g}\t"
                f"p={permutation['one_sided_pvalue']:.6g}\tpassed={pass_threshold}",
                flush=True,
            )

    primary_positive = sum(
        record["positive_focal_tf_criterion"]["passed"]
        for record in records
        if record["panel_role"] == "primary_expression_supported"
    )
    write_json(
        {
            "schema_version": "capselex_primary_genomic_model_v1",
            "analysis_status": "screening" if args.permutations < 1000 else "final",
            "primary_panel": args.focal_tfs,
            "sensitivity_panel": args.sensitivity_tfs,
            "partner_availability_negative_controls": args.availability_negative_controls,
            "outcome": "mean_of_replicate_log1p_200bp_exact_mean_GPZN_signal",
            "split": "five chromosome folds with nested ridge-alpha selection",
            "accessibility_covariate": {
                "name": "independent HEK293 DNase-seq proxy",
                "source_manifest": dnase_manifest,
                "coordinate_projection": (
                    "each included hg38 200-bp locus mapped to hg19 with the frozen UCSC chain; "
                    "exact source-BigWig mean queried over the mapped interval"
                ),
                "required_sensitivity": (
                    "models are repeated without accessibility on both mapped and full locus sets"
                ),
            },
            "primary_partial_r2_threshold": args.partial_r2_threshold,
            "minimum_positive_focal_tfs": args.minimum_positive_focal_tfs,
            "positive_primary_focal_tfs": primary_positive,
            "panel_go_gate_passed": primary_positive >= args.minimum_positive_focal_tfs,
            "tf_results": records,
            "source_assets": asset_manifest(
                [*feature_paths, args.dnase_resolved_manifest, args.partner_expression]
            ),
            "claim_boundary": (
                "A positive computational increment is evidence of predictive information beyond the "
                "specified baseline, not proof of partner occupancy or causal cooperativity. The HEK293 "
                "DNase assay is an independent cell-line proxy rather than matched chromatin profiling."
            ),
        },
        args.output,
    )


if __name__ == "__main__":
    main()
