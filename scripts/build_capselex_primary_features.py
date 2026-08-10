from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from cisgrammar.capselex import asset_manifest, audit_asset, sha256_file, write_json
from cisgrammar.capselex_genomic_assets import add_bigwig_signal, read_magix
from cisgrammar.capselex_primary_features import (
    extract_twobit_contexts,
    load_frozen_pwm_panel,
    longest_canonical_run,
    score_frozen_pwm_panel,
    sequence_covariates,
)


def _verify_asset(path: Path, record: dict[str, object]) -> None:
    if path.stat().st_size != record["size_bytes"]:
        raise RuntimeError(f"asset size check failed: {path}")
    if sha256_file(path) != record["sha256"]:
        raise RuntimeError(f"asset SHA-256 check failed: {path}")


def _reference_assets(manifest_path: Path, directory: Path) -> dict[str, Path]:
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    records = {record["filename"]: record for record in manifest["assets"]}
    required = {"hg38.2bit", "hg38.chrom.sizes", "twoBitToFa"}
    if set(records) != required:
        raise ValueError(f"reference manifest must contain exactly {sorted(required)}")
    paths = {filename: directory / filename for filename in required}
    for filename, path in paths.items():
        _verify_asset(path, records[filename])
    return paths


def _chip_assets(manifest_path: Path, directory: Path, focal_tf: str) -> list[tuple[Path, dict]]:
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    records = sorted(
        (record for record in manifest["assets"] if record["tf"] == focal_tf),
        key=lambda record: record["replicate"],
    )
    if [record["replicate"] for record in records] != [1, 2]:
        raise ValueError(f"{focal_tf} requires exactly ChIP replicates 1 and 2")
    assets = []
    for record in records:
        path = directory / record["filename"]
        _verify_asset(path, record)
        assets.append((path, record))
    return assets


def _attach_ght_scores(loci: pd.DataFrame, magix_path: Path, width_bp: int) -> pd.DataFrame:
    magix = read_magix(magix_path)
    magix["tile"] = (magix["midpoint"] // width_bp).astype(int)
    best = (
        magix.sort_values(
            ["chrom", "tile", "ght_score", "fdr", "pvalue"],
            ascending=[True, True, False, True, True],
        )
        .drop_duplicates(["chrom", "tile"])
        .loc[
            :,
            [
                "chrom",
                "tile",
                "ght_score",
                "ght_score_before_refinement",
                "fdr",
                "pvalue",
            ],
        ]
        .rename(columns={"fdr": "ght_fdr", "pvalue": "ght_pvalue"})
    )
    result = loci.merge(best, on=["chrom", "tile"], how="left", validate="one_to_one")
    if result["ght_score"].isna().any():
        raise RuntimeError("one or more frozen GHT loci did not map back to a selected MAGIX row")
    return result


def _sequence_digest(locus_ids: pd.Series, sequences: list[str]) -> str:
    digest = hashlib.sha256()
    for locus_id, sequence in zip(locus_ids, sequences, strict=True):
        digest.update(f"{locus_id}\t{sequence}\n".encode())
    return digest.hexdigest()


def _quantiles(values: pd.Series) -> dict[str, float]:
    array = values.to_numpy(dtype=float)
    return {
        "min": float(np.min(array)),
        "median": float(np.median(array)),
        "p90": float(np.quantile(array, 0.9)),
        "p99": float(np.quantile(array, 0.99)),
        "max": float(np.max(array)),
    }


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Build the frozen CAP x GHT x ChIP primary feature table for one focal TF"
    )
    parser.add_argument("--focal-tf", required=True)
    parser.add_argument("--loci", type=Path, required=True)
    parser.add_argument("--magix", type=Path, required=True)
    parser.add_argument("--chip-resolved-manifest", type=Path, required=True)
    parser.add_argument("--chip-directory", type=Path, required=True)
    parser.add_argument("--reference-resolved-manifest", type=Path, required=True)
    parser.add_argument("--reference-directory", type=Path, required=True)
    parser.add_argument("--cap-pwm-table", type=Path, required=True)
    parser.add_argument("--mex-top1", type=Path, required=True)
    parser.add_argument("--jaspar", type=Path, required=True)
    parser.add_argument("--monomer-audit", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    parser.add_argument("--locus-width", type=int, default=200)
    parser.add_argument("--context-width", type=int, default=400)
    parser.add_argument("--batch-size", type=int, default=512)
    parser.add_argument("--residual-bins", type=int, default=20)
    args = parser.parse_args()

    loci = pd.read_csv(args.loci, sep="\t")
    if loci.empty or set(loci["focal_tf"]) != {args.focal_tf}:
        raise ValueError("locus table focal TF does not match --focal-tf")
    if loci["locus_universe"].nunique() != 1 or loci["locus_universe"].iloc[0] != "ght-only":
        raise ValueError("primary feature construction requires the frozen GHT-only universe")
    loci = _attach_ght_scores(loci, args.magix, args.locus_width)

    references = _reference_assets(args.reference_resolved_manifest, args.reference_directory)
    sequences = extract_twobit_contexts(
        loci,
        references["hg38.2bit"],
        references["twoBitToFa"],
        references["hg38.chrom.sizes"],
        args.context_width,
    )
    sequence_sha256 = _sequence_digest(loci["locus_id"], sequences)
    panel = load_frozen_pwm_panel(
        args.focal_tf,
        args.cap_pwm_table,
        args.mex_top1,
        args.jaspar,
        args.monomer_audit,
    )
    maximum_pwm_length = max(
        [panel.focal_record.pwm.length]
        + [record.pwm.length for profile in panel.profiles for record in profile.partner_records]
        + [record.to_pwm().length for profile in panel.profiles for record in profile.cap_records]
    )
    scannable = np.array(
        [longest_canonical_run(sequence) >= maximum_pwm_length for sequence in sequences], dtype=bool
    )
    excluded_locus_ids = loci.loc[~scannable, "locus_id"].tolist()
    loci = loci.loc[scannable].reset_index(drop=True)
    sequences = [sequence for sequence, keep in zip(sequences, scannable, strict=True) if keep]

    sequence_features = sequence_covariates(sequences)
    pwm_features, pwm_manifest = score_frozen_pwm_panel(
        sequences,
        panel,
        loci["chromosome_fold"].to_numpy(),
        batch_size=args.batch_size,
        residual_bins=args.residual_bins,
    )
    result = pd.concat([loci, sequence_features, pwm_features], axis=1)

    chip_assets = _chip_assets(args.chip_resolved_manifest, args.chip_directory, args.focal_tf)
    result = add_bigwig_signal(result, chip_assets[0][0], "chip_rep1")
    result = add_bigwig_signal(result, chip_assets[1][0], "chip_rep2")
    result["chip_rep1_log1p"] = np.log1p(result["chip_rep1"])
    result["chip_rep2_log1p"] = np.log1p(result["chip_rep2"])
    result["chip_mean_log1p"] = (result["chip_rep1_log1p"] + result["chip_rep2_log1p"]) / 2
    if not np.isfinite(result.select_dtypes(include="number").to_numpy()).all():
        raise RuntimeError("primary feature table contains non-finite numeric values")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(args.output, sep="\t", index=False, compression="gzip", float_format="%.10g")
    report = {
        "schema_version": "capselex_primary_feature_table_v1",
        "focal_tf": args.focal_tf,
        "rows_before_reference_qc": int(len(scannable)),
        "rows": int(len(result)),
        "excluded_non_scannable_contexts": len(excluded_locus_ids),
        "excluded_non_scannable_locus_ids": excluded_locus_ids,
        "reference_assembly": "UCSC hg38",
        "context_width_bp": args.context_width,
        "sequence_context_sha256_before_qc": sequence_sha256,
        "sequence_exclusion_uses_chip_outcome": False,
        "ght_score_aggregation": (
            "maximum positive coefficient.ar per fixed 200-bp tile; ties resolved by FDR then p-value"
        ),
        "chip_outcome": "mean of replicate-specific log1p exact 200-bp mean GPZN signal",
        "chip_replicates": [
            {
                "replicate": record["replicate"],
                "gsm": record["gsm"],
                "filename": path.name,
                "sha256": record["sha256"],
            }
            for path, record in chip_assets
        ],
        "pwm_panel": pwm_manifest,
        "feature_columns": list(result.columns),
        "outcome_quantiles": _quantiles(result["chip_mean_log1p"]),
        "source_assets": asset_manifest(
            [
                args.loci,
                args.magix,
                args.chip_resolved_manifest,
                args.reference_resolved_manifest,
                args.cap_pwm_table,
                args.mex_top1,
                args.jaspar,
                args.monomer_audit,
            ]
        ),
        "output": audit_asset(args.output).__dict__,
    }
    write_json(report, args.report)
    print(
        f"{args.focal_tf}\trows={len(result)}\tpair_mechanisms={pwm_manifest['pair_mechanism_units']}\t"
        f"excluded_contexts={len(excluded_locus_ids)}"
    )


if __name__ == "__main__":
    main()
