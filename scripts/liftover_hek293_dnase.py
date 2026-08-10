from __future__ import annotations

import argparse
import gzip
import json
import os
import shutil
import subprocess
import tempfile
from pathlib import Path

import pyBigWig

from cisgrammar.capselex import audit_asset, sha256_file, write_json


def _verify(path: Path, record: dict[str, object]) -> None:
    if path.stat().st_size != record["size_bytes"] or sha256_file(path) != record["sha256"]:
        raise RuntimeError(f"integrity check failed: {path}")


def _bedgraph_summary(path: Path) -> dict[str, int | float]:
    records = 0
    bases = 0
    signal_sum = 0.0
    with path.open(encoding="utf-8") as handle:
        for raw in handle:
            if not raw or raw.startswith("#"):
                continue
            chrom, start, end, value = raw.split()[:4]
            del chrom
            width = int(end) - int(start)
            if width <= 0:
                raise ValueError(f"invalid bedGraph interval in {path}")
            records += 1
            bases += width
            signal_sum += width * float(value)
    return {"records": records, "covered_bases": bases, "base_weighted_signal_sum": signal_sum}


def main() -> None:
    parser = argparse.ArgumentParser(description="Lift the exact HEK293 DNase signal from hg19 to hg38")
    parser.add_argument("--source-manifest", type=Path, required=True)
    parser.add_argument("--source-directory", type=Path, required=True)
    parser.add_argument("--hg38-chrom-sizes", type=Path, required=True)
    parser.add_argument("--hg38-reference-manifest", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()

    source_manifest = json.loads(args.source_manifest.read_text(encoding="utf-8"))
    source_records = {record["filename"]: record for record in source_manifest["assets"]}
    paths = {name: args.source_directory / name for name in source_records}
    for name, path in paths.items():
        _verify(path, source_records[name])
    for name in ("bigWigToBedGraph", "liftOver", "bedGraphToBigWig"):
        paths[name].chmod(paths[name].stat().st_mode | 0o111)
    reference_manifest = json.loads(args.hg38_reference_manifest.read_text(encoding="utf-8"))
    chrom_record = next(
        record for record in reference_manifest["assets"] if record["filename"] == "hg38.chrom.sizes"
    )
    _verify(args.hg38_chrom_sizes, chrom_record)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="cisgrammar_dnase_", dir=args.output.parent) as temporary:
        directory = Path(temporary)
        source_bedgraph = directory / "hek293.hg19.bedGraph"
        target_unsorted = directory / "hek293.hg38.unsorted.bedGraph"
        target_bedgraph = directory / "hek293.hg38.bedGraph"
        unmapped = directory / "hek293.hg19.unmapped.bed"
        chain = directory / "hg19ToHg38.over.chain"

        subprocess.run(
            [
                str(paths["bigWigToBedGraph"].resolve()),
                str(paths["GSM2902639_HEK293_DNase.bw"]),
                str(source_bedgraph),
            ],
            check=True,
        )
        with gzip.open(paths["hg19ToHg38.over.chain.gz"], "rb") as source, chain.open("wb") as target:
            shutil.copyfileobj(source, target)
        subprocess.run(
            [
                str(paths["liftOver"].resolve()),
                "-bedPlus=4",
                "-tab",
                str(source_bedgraph),
                str(chain),
                str(target_unsorted),
                str(unmapped),
            ],
            check=True,
        )
        sort_environment = {**os.environ, "LC_ALL": "C", "TMPDIR": str(directory)}
        with target_bedgraph.open("w", encoding="utf-8") as handle:
            subprocess.run(
                ["sort", "-T", str(directory), "-k1,1", "-k2,2n", str(target_unsorted)],
                check=True,
                stdout=handle,
                env=sort_environment,
            )
        source_summary = _bedgraph_summary(source_bedgraph)
        target_summary = _bedgraph_summary(target_bedgraph)
        subprocess.run(
            [
                str(paths["bedGraphToBigWig"].resolve()),
                str(target_bedgraph),
                str(args.hg38_chrom_sizes),
                str(args.output),
            ],
            check=True,
        )

    with pyBigWig.open(str(args.output)) as bigwig:
        chromosomes = bigwig.chroms()
        if not all(f"chr{chromosome}" in chromosomes for chromosome in range(1, 23)):
            raise RuntimeError("lifted DNase bigWig does not contain every autosome")
    output_asset = audit_asset(args.output)
    write_json(
        {
            "schema_version": "hek293_dnase_hg38_liftover_v1",
            "source": source_manifest["biological_source"],
            "coordinate_conversion": {
                "source_assembly": "hg19",
                "target_assembly": "hg38",
                "chain": source_records["hg19ToHg38.over.chain.gz"],
                "liftOver_parameters": ["-bedPlus=4", "-tab"],
                "source_bedgraph": source_summary,
                "mapped_bedgraph": target_summary,
                "covered_base_fraction_retained": (
                    target_summary["covered_bases"] / source_summary["covered_bases"]
                ),
            },
            "output": output_asset.__dict__,
            "claim_boundary": source_manifest["biological_source"]["claim_boundary"],
        },
        args.report,
    )
    print(f"{args.output}\t{output_asset.bytes}\t{output_asset.sha256}")


if __name__ == "__main__":
    main()
