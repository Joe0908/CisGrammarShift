from __future__ import annotations

import argparse
from pathlib import Path

from cisgrammar.capselex import asset_manifest, write_json
from cisgrammar.capselex_nature_supplements import (
    audit_supplement_coverage,
    read_interaction_matrix,
    read_pwm_models,
    read_spacing_counts,
)


def main() -> None:
    parser = argparse.ArgumentParser(description="Audit official CAP-SELEX Nature supplements")
    parser.add_argument("--interaction-table", type=Path, required=True)
    parser.add_argument("--pwm-table", type=Path, required=True)
    parser.add_argument("--spacing-table", type=Path, required=True)
    parser.add_argument(
        "--focal-tfs",
        nargs="+",
        default=["FLI1", "GABPA", "GCM1", "MAX", "PAX7", "RFX5"],
    )
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    report = audit_supplement_coverage(
        read_interaction_matrix(args.interaction_table),
        read_pwm_models(args.pwm_table),
        read_spacing_counts(args.spacing_table),
        args.focal_tfs,
    )
    report["source_assets"] = asset_manifest(
        [args.interaction_table, args.pwm_table, args.spacing_table]
    )
    write_json(report, args.output)
    for record in report["panel"]:
        print(
            f"{record['tf']}\tpositive_pairs={record['positive_directed_pairs']}\t"
            f"representative_pwms={record['representative_pwm_models']}\t"
            f"spacing_count_pairs={record['spacing_count_pairs']}"
        )


if __name__ == "__main__":
    main()
