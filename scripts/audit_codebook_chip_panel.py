from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from cisgrammar.capselex import write_json
from cisgrammar.capselex_chip_qc import replicate_qc
from cisgrammar.capselex_genomic_assets import add_bigwig_signal


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def quantiles(values: np.ndarray) -> dict[str, float]:
    return {
        "min": float(np.min(values)),
        "median": float(np.median(values)),
        "p90": float(np.quantile(values, 0.90)),
        "p99": float(np.quantile(values, 0.99)),
        "max": float(np.max(values)),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Audit focal GPZN ChIP replicate signals")
    parser.add_argument("--resolved-manifest", type=Path, required=True)
    parser.add_argument("--bigwig-directory", type=Path, required=True)
    parser.add_argument("--locus-directory", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    manifest = json.loads(args.resolved_manifest.read_text(encoding="utf-8"))
    assets = manifest["assets"]
    grouped_assets = {
        tf: sorted(
            (record for record in assets if record["tf"] == tf),
            key=lambda record: record["replicate"],
        )
        for tf in sorted({record["tf"] for record in assets})
    }
    locus_inputs: dict[str, tuple[Path, str]] = {}
    for tf, tf_assets in grouped_assets.items():
        if [record["replicate"] for record in tf_assets] != [1, 2]:
            raise ValueError(f"{tf} must have exactly replicates 1 and 2")
        for asset in tf_assets:
            path = args.bigwig_directory / asset["filename"]
            if path.stat().st_size != asset["size_bytes"] or sha256(path) != asset["sha256"]:
                raise RuntimeError(f"bigWig integrity check failed: {path}")
        locus_path = args.locus_directory / f"{tf}.ght_only.loci.tsv.gz"
        locus_manifest_path = locus_path.with_name(locus_path.name + ".manifest.json")
        locus_manifest = json.loads(locus_manifest_path.read_text(encoding="utf-8"))
        observed_locus_sha256 = sha256(locus_path)
        if observed_locus_sha256 != locus_manifest["output"]["sha256"]:
            raise RuntimeError(f"locus integrity check failed: {locus_path}")
        locus_inputs[tf] = (locus_path, observed_locus_sha256)

    records = []
    for tf, tf_assets in grouped_assets.items():
        paths = [args.bigwig_directory / record["filename"] for record in tf_assets]
        locus_path, locus_sha256 = locus_inputs[tf]
        loci = pd.read_csv(locus_path, sep="\t")
        first = add_bigwig_signal(loci, paths[0], "chip_rep1")["chip_rep1"].to_numpy()
        second = add_bigwig_signal(loci, paths[1], "chip_rep2")["chip_rep2"].to_numpy()
        mean_log1p = (np.log1p(first) + np.log1p(second)) / 2
        qc = replicate_qc(first, second)
        records.append(
            {
                "tf": tf,
                "loci": int(len(loci)),
                "locus_sha256": locus_sha256,
                "replicate_1_gsm": tf_assets[0]["gsm"],
                "replicate_2_gsm": tf_assets[1]["gsm"],
                **qc,
                "both_nonzero_fraction": float(np.mean((first > 0) & (second > 0))),
                "replicate_1_signal": quantiles(first),
                "replicate_2_signal": quantiles(second),
                "mean_log1p_outcome": quantiles(mean_log1p),
            }
        )
        print(f"{tf}\t{len(loci)}\t{qc['spearman_log1p']:.6f}", flush=True)

    write_json(
        {
            "schema_version": "codebook_chip_gpzn_smoke_qc_v1",
            "accession": manifest["accession"],
            "chip_assets_primary_eligible": True,
            "processing_pipeline": "Toronto_GPZN_only",
            "outcome_definition": "mean_of_replicate_log1p_200bp_exact_mean_signal",
            "replicate_confirmation": "same_incremental_effect_direction_in_both_replicates",
            "qc_locus_universe": "GSE278858 pre-v2 GHT-only smoke loci",
            "qc_locus_universe_primary_eligible": False,
            "reason_not_primary": "Codebook v2 revised GHT MAGIX candidate peak generation.",
            "tf_qc": records,
        },
        args.output,
    )


if __name__ == "__main__":
    main()
