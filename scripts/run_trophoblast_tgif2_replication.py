from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from cisgrammar.capselex import asset_manifest, sha256_file, stable_seed, write_json
from cisgrammar.capselex_genomic_assets import add_bigwig_signal
from cisgrammar.capselex_genomic_model import (
    chromosome_bootstrap_partial_r2,
    chromosome_nested_continuous,
    chromosome_shift_null,
)

PAIR = "TGIF2_GCM1"
ADDITION = f"cap_excess::composite::{PAIR}"
BASELINE = [
    "ght_score",
    "focal_monomer_score",
    f"partner_raw::{PAIR}",
    "gc_fraction",
    "cpg_per_100bp",
]


def _verify_assets(directory: Path, manifest_path: Path) -> tuple[dict, dict[tuple[str, str], Path]]:
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    records = {
        (record["state"], record["clone"]): record
        for record in manifest["assets"]
        if record.get("state") in {"EVT", "ST"}
    }
    expected = {(state, clone) for state in ("EVT", "ST") for clone in ("B31", "CT27")}
    if set(records) != expected:
        raise ValueError("external GCM1 ChIP manifest must contain two clones in EVT and ST")
    paths = {}
    for key, record in records.items():
        path = directory / record["filename"]
        if path.stat().st_size != record["size_bytes"] or sha256_file(path) != record["sha256"]:
            raise RuntimeError(f"external GCM1 ChIP integrity check failed: {path}")
        paths[key] = path
    return manifest, paths


def _fit(frame: pd.DataFrame, outcome: str) -> dict[str, object]:
    result = chromosome_nested_continuous(frame, outcome, BASELINE, [ADDITION])
    coefficients = np.array(
        [record["standardized_addition_coefficients"][0] for record in result.summary["fold_results"]]
    )
    return {
        **result.summary,
        "standardized_addition_coefficients": coefficients.tolist(),
        "median_standardized_addition_coefficient": float(np.median(coefficients)),
        "positive_coefficient_in_all_folds": bool(np.all(coefficients > 0)),
        "partial_r2_chromosome_bootstrap": chromosome_bootstrap_partial_r2(
            frame,
            result.predictions,
            replicates=2000,
            seed=stable_seed("trophoblast_bootstrap", outcome),
        ),
    }


def _permutation(
    frame: pd.DataFrame,
    outcome: str,
    observed: float,
    permutations: int,
    seed: int,
) -> dict[str, object]:
    null = []
    for index in range(permutations):
        shifted = chromosome_shift_null(
            frame,
            [ADDITION],
            seed=stable_seed(seed, outcome, index),
        )
        result = chromosome_nested_continuous(shifted, outcome, BASELINE, [ADDITION])
        null.append(float(result.summary["partial_r2"]))
    null_array = np.asarray(null)
    exceedances = int(np.sum(null_array >= observed))
    return {
        "permutations": permutations,
        "null": "candidate grammar circularly shifted within chromosome",
        "exceedances": exceedances,
        "one_sided_pvalue": (exceedances + 1) / (permutations + 1),
        "null_quantiles": {
            "minimum": float(np.min(null_array)),
            "median": float(np.median(null_array)),
            "p95": float(np.quantile(null_array, 0.95)),
            "maximum": float(np.max(null_array)),
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Test the frozen TGIF2-GCM1 composite hypothesis in independent trophoblast ChIP"
    )
    parser.add_argument("--features", type=Path, required=True)
    parser.add_argument("--asset-directory", type=Path, required=True)
    parser.add_argument("--resolved-manifest", type=Path, required=True)
    parser.add_argument("--expression-gate", type=Path, required=True)
    parser.add_argument("--pair-decomposition", type=Path, required=True)
    parser.add_argument("--permutations", type=int, default=100)
    parser.add_argument("--partial-r2-threshold", type=float, default=0.005)
    parser.add_argument("--seed", type=int, default=20260809)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    expression_gate = json.loads(args.expression_gate.read_text(encoding="utf-8"))
    if expression_gate["candidate"] != PAIR or not expression_gate["expression_gate_passed"]:
        raise ValueError("the frozen TGIF2 trophoblast expression gate did not pass")
    decomposition = json.loads(args.pair_decomposition.read_text(encoding="utf-8"))
    gcm1 = next(record for record in decomposition["tf_results"] if record["focal_tf"] == "GCM1")
    if gcm1["pair_results_ranked_by_mean_outcome_partial_r2"][0]["pair"] != PAIR:
        raise ValueError("TGIF2_GCM1 is not the frozen leading Codebook GCM1 pair")
    manifest, bigwigs = _verify_assets(args.asset_directory, args.resolved_manifest)
    frame = pd.read_csv(args.features, sep="\t")
    missing = sorted(set([*BASELINE, ADDITION]) - set(frame))
    if missing:
        raise ValueError(f"GCM1 feature table is missing: {', '.join(missing)}")

    state_results = []
    for state in ("EVT", "ST"):
        state_frame = frame.copy()
        for clone in ("B31", "CT27"):
            state_frame = add_bigwig_signal(
                state_frame, bigwigs[(state, clone)], f"{state}_{clone}_signal"
            )
            state_frame[f"{state}_{clone}_log1p"] = np.log1p(
                state_frame[f"{state}_{clone}_signal"]
            )
        clone_outcomes = [f"{state}_{clone}_log1p" for clone in ("B31", "CT27")]
        mean_outcome = f"{state}_mean_log1p"
        state_frame[mean_outcome] = state_frame[clone_outcomes].mean(axis=1)
        outcomes = {outcome: _fit(state_frame, outcome) for outcome in [mean_outcome, *clone_outcomes]}
        permutation = _permutation(
            state_frame,
            mean_outcome,
            outcomes[mean_outcome]["partial_r2"],
            args.permutations,
            args.seed,
        )
        clone_positive = all(
            outcomes[outcome]["median_standardized_addition_coefficient"] > 0
            for outcome in clone_outcomes
        )
        state_passed = bool(
            outcomes[mean_outcome]["partial_r2"] >= args.partial_r2_threshold
            and outcomes[mean_outcome]["median_standardized_addition_coefficient"] > 0
            and clone_positive
            and permutation["one_sided_pvalue"] <= 0.05
        )
        state_results.append(
            {
                "state": state,
                "rows": len(state_frame),
                "clone_log1p_spearman": float(
                    spearmanr(
                        state_frame[clone_outcomes[0]], state_frame[clone_outcomes[1]]
                    ).statistic
                ),
                "outcome_models": outcomes,
                "screening_permutation": permutation,
                "state_replication_passed": state_passed,
            }
        )
        print(
            f"{state}\tpartial_r2={outcomes[mean_outcome]['partial_r2']:.6g}\t"
            f"p={permutation['one_sided_pvalue']:.6g}\tpassed={state_passed}",
            flush=True,
        )

    overall_passed = all(record["state_replication_passed"] for record in state_results)
    write_json(
        {
            "schema_version": "trophoblast_tgif2_gcm1_replication_v1",
            "analysis_status": "screening" if args.permutations < 1000 else "final",
            "candidate": PAIR,
            "mechanism": "composite",
            "discovery_context": "Codebook GCM1 ChIP in HEK293",
            "replication_context": "GCM1 ChIP in differentiated human trophoblast stem cells",
            "locus_universe": "official Codebook v2 GCM1 GHT-only fixed 200-bp tiles",
            "baseline_features": BASELINE,
            "addition_feature": ADDITION,
            "partial_r2_threshold": args.partial_r2_threshold,
            "state_results": state_results,
            "replication_passed_in_both_states": overall_passed,
            "source_manifest_contract": manifest["frozen_hypothesis"],
            "source_assets": asset_manifest(
                [
                    args.features,
                    args.resolved_manifest,
                    args.expression_gate,
                    args.pair_decomposition,
                ]
            ),
            "claim_boundary": (
                "Cross-context predictive replication supports a transferable sequence association; "
                "it does not establish TGIF2 co-occupancy or causal cooperation with GCM1."
            ),
        },
        args.output,
    )


if __name__ == "__main__":
    main()
