from __future__ import annotations

import argparse
import gzip
import io
import json
import tarfile
from pathlib import Path

import pandas as pd

from cisgrammar.capselex import asset_manifest, sha256_file, write_json

SAMPLES = {
    "EVT": {
        "B31": "GSM7809996_B31_EVT.genes.results.gz",
        "CT27": "GSM7810002_CT27_EVT.genes.results.gz",
    },
    "ST_day2": {
        "B31": "GSM7809998_B31_ST_DAY2.genes.results.gz",
        "CT27": "GSM7810004_CT27_ST_DAY2.genes.results.gz",
    },
    "ST_day4": {
        "B31": "GSM7809999_B31_ST_DAY4.genes.results.gz",
        "CT27": "GSM7810005_CT27_ST_DAY4.genes.results.gz",
    },
}


def _read_gene(archive: tarfile.TarFile, member: str, gene: str) -> dict[str, float]:
    extracted = archive.extractfile(member)
    if extracted is None:
        raise ValueError(f"missing RSEM member: {member}")
    with gzip.GzipFile(fileobj=extracted) as compressed:
        frame = pd.read_csv(io.BytesIO(compressed.read()), sep="\t")
    selected = frame[frame["gene_id"].eq(gene)]
    if len(selected) != 1:
        raise ValueError(f"expected one {gene} row in {member}; observed {len(selected)}")
    return {
        "TPM": float(selected.iloc[0]["TPM"]),
        "FPKM": float(selected.iloc[0]["FPKM"]),
        "expected_count": float(selected.iloc[0]["expected_count"]),
    }


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Apply the frozen TGIF2 availability gate in independent trophoblast RNA-seq"
    )
    parser.add_argument("--archive", type=Path, required=True)
    parser.add_argument("--resolved-manifest", type=Path, required=True)
    parser.add_argument("--minimum-tpm", type=float, default=1.0)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    manifest = json.loads(args.resolved_manifest.read_text(encoding="utf-8"))
    record = next(record for record in manifest["assets"] if record["filename"] == args.archive.name)
    if args.archive.stat().st_size != record["size_bytes"] or sha256_file(args.archive) != record["sha256"]:
        raise RuntimeError("trophoblast RNA archive integrity check failed")

    states = {}
    with tarfile.open(args.archive) as archive:
        members = {member.name for member in archive.getmembers()}
        required = {member for samples in SAMPLES.values() for member in samples.values()}
        if not required.issubset(members):
            raise ValueError(f"archive is missing: {sorted(required - members)}")
        for state, samples in SAMPLES.items():
            states[state] = {
                clone: _read_gene(archive, member, "TGIF2")
                for clone, member in samples.items()
            }
    all_tpm = [record["TPM"] for state in states.values() for record in state.values()]
    passed = all(value >= args.minimum_tpm for value in all_tpm)
    write_json(
        {
            "schema_version": "trophoblast_tgif2_expression_gate_v1",
            "candidate": "TGIF2_GCM1",
            "gene": "TGIF2",
            "source": "GSE244254 submitted RSEM gene results",
            "states": states,
            "minimum_tpm_in_every_profile": args.minimum_tpm,
            "minimum_observed_tpm": min(all_tpm),
            "expression_gate_passed": passed,
            "source_assets": asset_manifest([args.archive, args.resolved_manifest]),
            "claim_boundary": (
                "RNA abundance supports TGIF2 availability but does not establish protein abundance "
                "or co-occupancy with GCM1."
            ),
        },
        args.output,
    )
    print(f"TGIF2\tminimum_TPM={min(all_tpm):.6g}\tpassed={passed}")


if __name__ == "__main__":
    main()
