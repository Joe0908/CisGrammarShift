from __future__ import annotations

import argparse
from pathlib import Path

from cisgrammar.experiment import run_experiment


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="cisgrammar",
        description="Counterfactual benchmark for cis-regulatory motif grammar",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)
    run_parser = subparsers.add_parser("run", help="run a configured benchmark")
    run_parser.add_argument("--config", type=Path, required=True)
    run_parser.add_argument("--output", type=Path, required=True)
    run_parser.add_argument("--device", default="auto", help="auto, cpu, cuda, or a torch device string")
    return parser


def main() -> None:
    args = build_parser().parse_args()
    if args.command == "run":
        run_experiment(args.config, args.output, requested_device=args.device)


if __name__ == "__main__":
    main()
