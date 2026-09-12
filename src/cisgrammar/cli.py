from __future__ import annotations

import argparse
from pathlib import Path


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="cisgrammar",
        description="Counterfactual benchmark of cis-regulatory grammar learning",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)
    run_parser = subparsers.add_parser("run", help="run the matched-pair grammar benchmark")
    run_parser.add_argument("--config", type=Path, required=True)
    run_parser.add_argument("--output", type=Path, required=True)
    run_parser.add_argument("--device", default="auto")
    return parser


def main() -> None:
    args = build_parser().parse_args()
    if args.command != "run":
        raise ValueError(f"unsupported command: {args.command}")
    try:
        from cisgrammar.experiment import run_experiment
    except ModuleNotFoundError as error:
        if error.name == "torch":
            raise SystemExit(
                'model training requires the ML extra: pip install -e ".[ml]"'
            ) from error
        raise
    run_experiment(args.config, args.output, requested_device=args.device)


if __name__ == "__main__":
    main()
