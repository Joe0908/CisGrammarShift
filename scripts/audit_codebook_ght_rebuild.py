from __future__ import annotations

import argparse
import json
from pathlib import Path

from cisgrammar.capselex import asset_manifest, write_json
from cisgrammar.capselex_ght_rebuild import build_ght_rebuild_audit


def verify_frozen_asset(path: Path, record: dict[str, object]) -> None:
    observed = asset_manifest([path])[0]
    if path.name != record["filename"]:
        raise RuntimeError(f"unexpected source filename: {path.name}")
    if observed["bytes"] != record["size_bytes"]:
        raise RuntimeError(f"source size does not match frozen config: {path.name}")
    if observed["sha256"] != record["sha256"]:
        raise RuntimeError(f"source SHA-256 does not match frozen config: {path.name}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Freeze focal GHT target FASTQs and audit exact MAGIX v2 rebuild readiness"
    )
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--ena-report", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    config = json.loads(args.config.read_text(encoding="utf-8"))
    verify_frozen_asset(args.metadata, config["metadata"])
    verify_frozen_asset(args.ena_report, config["ena_filereport"])

    report = build_ght_rebuild_audit(
        args.metadata,
        args.ena_report,
        config["primary_tfs"],
        config["sensitivity_tfs"],
        config["audit_date"],
        config["magix"]["commit"],
    )
    report["config_asset"] = asset_manifest([args.config])[0]
    write_json(report, args.output)
    print(
        f"primary={report['primary_target_fastq_gib']} GiB\t"
        f"all={report['all_target_fastq_gib']} GiB\t"
        f"exact_v2_ready={report['exact_author_v2_rebuild_ready']}"
    )


if __name__ == "__main__":
    main()
