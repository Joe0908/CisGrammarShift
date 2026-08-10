from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from cisgrammar.capselex import write_json
from cisgrammar.capselex_analysis_contract import AnalysisContract
from cisgrammar.capselex_codebook_metadata import audit_three_way_overlap
from cisgrammar.capselex_dataset import (
    build_directed_pair_dataset,
    read_construct_table,
    read_interaction_table,
)
from cisgrammar.capselex_genomic_assets import build_assay_union_loci, read_narrowpeak
from cisgrammar.capselex_tf_panel import run_tf_panel
from cisgrammar.capselex_trophoblast_continuous import run_trophoblast_deseq2_hichip_benchmark
from cisgrammar.capselex_trophoblast_deseq2 import export_deseq2_inputs
from cisgrammar.experiment import run_experiment


def _path(parser: argparse.ArgumentParser, name: str, required: bool = True) -> None:
    parser.add_argument(name, type=Path, required=required)


def _build_capselex_parser(subparsers: argparse._SubParsersAction[argparse.ArgumentParser]) -> None:
    cap = subparsers.add_parser("capselex", help="CAP/GHT/ChIP mechanistic reanalysis")
    commands = cap.add_subparsers(dest="capselex_command", required=True)

    contract = commands.add_parser("contract", help="write the frozen primary analysis contract")
    _path(contract, "--output")

    dataset = commands.add_parser("dataset", help="build a directed CAP pair dataset")
    _path(dataset, "--interaction-table")
    _path(dataset, "--construct-table")
    _path(dataset, "--output")

    panel = commands.add_parser("tf-panel", help="run leakage-resistant leave-one-TF-out baselines")
    _path(panel, "--dataset")
    panel.add_argument(
        "--target",
        choices=("cooperative_signal", "composite_motif", "spacing_or_orientation"),
        required=True,
    )
    _path(panel, "--output")

    audit = commands.add_parser("codebook-audit", help="audit CAP/Codebook/GHT/ChIP TF overlap")
    _path(audit, "--cap-dataset")
    _path(audit, "--plasmid-metadata")
    _path(audit, "--chip-metadata")
    _path(audit, "--ght-metadata")
    _path(audit, "--output")

    loci = commands.add_parser("build-loci", help="build outcome-independent ChIP/GHT assay-union loci")
    _path(loci, "--chip-peaks")
    _path(loci, "--ght-peaks")
    loci.add_argument("--focal-tf", required=True)
    loci.add_argument("--width-bp", type=int, default=200)
    _path(loci, "--output")

    deseq_export = commands.add_parser("trophoblast-deseq2-export", help="export exact DESeq2 count matrices")
    _path(deseq_export, "--rna-dataset")
    _path(deseq_export, "--raw-rsem-filter-audit", required=False)
    _path(deseq_export, "--output-directory")

    continuous = commands.add_parser(
        "trophoblast-deseq2-hichip-benchmark",
        help="run the final chromosome-held-out continuous functional sensitivity",
    )
    _path(continuous, "--evt-gene-features")
    _path(continuous, "--st-gene-features")
    _path(continuous, "--deseq2-directory")
    _path(continuous, "--output")
    _path(continuous, "--predictions")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="cisgrammar",
        description="Intrinsic specificity and cooperative DNA grammar under biological context shift",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)
    run_parser = subparsers.add_parser("run", help="run the original controlled synthetic benchmark")
    run_parser.add_argument("--config", type=Path, required=True)
    run_parser.add_argument("--output", type=Path, required=True)
    run_parser.add_argument("--device", default="auto")
    _build_capselex_parser(subparsers)
    return parser


def _run_capselex(args: argparse.Namespace) -> None:
    if args.capselex_command == "contract":
        AnalysisContract().write(args.output)
    elif args.capselex_command == "dataset":
        interactions = read_interaction_table(args.interaction_table)
        constructs = read_construct_table(args.construct_table)
        frame = build_directed_pair_dataset(interactions, constructs)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        frame.to_csv(args.output, sep="\t", index=False)
    elif args.capselex_command == "tf-panel":
        frame = pd.read_csv(args.dataset, sep="\t")
        write_json(run_tf_panel(frame, args.target), args.output)
    elif args.capselex_command == "codebook-audit":
        frame = pd.read_csv(args.cap_dataset, sep="\t")
        tfs = sorted(set(frame["bait"]) | set(frame["prey"]))
        payload = audit_three_way_overlap(
            tfs,
            args.plasmid_metadata,
            args.chip_metadata,
            args.ght_metadata,
        )
        write_json(payload, args.output)
    elif args.capselex_command == "build-loci":
        chip = read_narrowpeak(args.chip_peaks)
        ght = read_narrowpeak(args.ght_peaks)
        frame = build_assay_union_loci(chip, ght, args.focal_tf, args.width_bp)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        frame.to_csv(args.output, sep="\t", index=False)
    elif args.capselex_command == "trophoblast-deseq2-export":
        export_deseq2_inputs(args.rna_dataset, args.output_directory, args.raw_rsem_filter_audit)
    elif args.capselex_command == "trophoblast-deseq2-hichip-benchmark":
        run_trophoblast_deseq2_hichip_benchmark(
            args.evt_gene_features,
            args.st_gene_features,
            args.deseq2_directory,
            args.output,
            args.predictions,
        )
    else:
        raise ValueError(f"unsupported CAP-SELEX command: {args.capselex_command}")


def main() -> None:
    args = build_parser().parse_args()
    if args.command == "run":
        run_experiment(args.config, args.output, requested_device=args.device)
    elif args.command == "capselex":
        _run_capselex(args)
    else:
        raise ValueError(json.dumps(vars(args), default=str))


if __name__ == "__main__":
    main()
