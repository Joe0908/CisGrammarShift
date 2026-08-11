from __future__ import annotations

import argparse
import json
from pathlib import Path

from cisgrammar.capselex import write_json

FROZEN_FIELDS = (
    "primary_panel",
    "sensitivity_panel",
    "partner_availability_negative_controls",
    "primary_partial_r2_threshold",
    "minimum_positive_focal_tfs",
    "split",
)


def _evaluable_records(report: dict) -> dict[str, dict]:
    return {
        record["focal_tf"]: record
        for record in report["tf_results"]
        if record.get("primary_outcome_models") is not None
    }


def build_comparison(gphn: dict, gpzn: dict) -> dict:
    for field in FROZEN_FIELDS:
        if gphn[field] != gpzn[field]:
            raise ValueError(f"pipeline reports differ in frozen field: {field}")

    gphn_records = _evaluable_records(gphn)
    gpzn_records = _evaluable_records(gpzn)
    if set(gphn_records) != set(gpzn_records):
        raise ValueError("pipeline reports contain different evaluable TF panels")

    records = []
    for tf in sorted(gphn_records):
        gphn_record = gphn_records[tf]
        gpzn_record = gpzn_records[tf]
        gphn_model = gphn_record["primary_outcome_models"]["chip_mean_log1p"]
        gpzn_model = gpzn_record["primary_outcome_models"]["chip_mean_log1p"]
        gphn_ci = gphn_model["partial_r2_chromosome_bootstrap"]["ci95_percentile"]
        gpzn_ci = gpzn_model["partial_r2_chromosome_bootstrap"]["ci95_percentile"]
        gphn_passed = gphn_record["positive_focal_tf_criterion"]["passed"]
        gpzn_passed = gpzn_record["positive_focal_tf_criterion"]["passed"]
        records.append(
            {
                "tf": tf,
                "panel_role": gphn_record["panel_role"],
                "partial_r2_gphn": gphn_model["partial_r2"],
                "partial_r2_gpzn": gpzn_model["partial_r2"],
                "ci95_low_gphn": gphn_ci[0],
                "ci95_low_gpzn": gpzn_ci[0],
                "ci95_high_gphn": gphn_ci[1],
                "ci95_high_gpzn": gpzn_ci[1],
                "median_standardized_coefficient_gphn": gphn_model[
                    "addition_coefficients"
                ]["median"][0],
                "median_standardized_coefficient_gpzn": gpzn_model[
                    "addition_coefficients"
                ]["median"][0],
                "passed_prespecified_threshold_gphn": gphn_passed,
                "passed_prespecified_threshold_gpzn": gpzn_passed,
                "partial_r2_difference_gphn_minus_gpzn": (
                    gphn_model["partial_r2"] - gpzn_model["partial_r2"]
                ),
                "positive_effect_both_pipelines": (
                    gphn_model["addition_coefficients"]["median"][0] > 0
                    and gpzn_model["addition_coefficients"]["median"][0] > 0
                ),
            }
        )
    return {
        "schema_version": "gphn_gpzn_pipeline_comparison_v1",
        "analysis_role": "prespecified ChIP-processing-pipeline sensitivity",
        "primary_pipeline": "McGill_GPHN",
        "sensitivity_pipeline": "Toronto_GPZN",
        "retuning_performed": False,
        "comparison_is_descriptive": True,
        "records": records,
        "claim_boundary": (
            "Pipeline concordance supports processing robustness of association direction "
            "but does not establish causal TF cooperation."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Compare frozen GPHN-primary and GPZN-sensitivity panel reports"
    )
    parser.add_argument("--gphn-screening", type=Path, required=True)
    parser.add_argument("--gpzn-screening", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    gphn = json.loads(args.gphn_screening.read_text(encoding="utf-8"))
    gpzn = json.loads(args.gpzn_screening.read_text(encoding="utf-8"))
    write_json(build_comparison(gphn, gpzn), args.output)


if __name__ == "__main__":
    main()
