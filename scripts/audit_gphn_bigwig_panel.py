from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import pyBigWig

from cisgrammar.capselex import write_json

HG38_SENTINEL_LENGTHS = {
    "chr1": 248956422,
    "chr2": 242193529,
    "chr22": 50818468,
    "chrX": 156040895,
    "chrY": 57227415,
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(16 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Audit the author-supplied McGill-GPHN six-TF bigWig panel"
    )
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--bigwig-directory", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    manifest = json.loads(args.manifest.read_text(encoding="utf-8"))
    expected = {record["filename"]: record for record in manifest["assets"]}
    observed = sorted(path.name for path in args.bigwig_directory.glob("*.bw"))
    missing = sorted(set(expected) - set(observed))
    unexpected = sorted(set(observed) - set(expected))
    records = []

    for filename, asset in expected.items():
        path = args.bigwig_directory / filename
        record = {
            "filename": filename,
            "tf": asset["tf"],
            "replicate": int(asset["replicate"]),
            "exists": path.is_file(),
        }
        if not path.is_file():
            record["status"] = "FAIL_MISSING"
            records.append(record)
            continue

        record["size_bytes"] = path.stat().st_size
        record["sha256"] = sha256(path)
        integrity_ok = (
            record["size_bytes"] == int(asset["size_bytes"])
            and record["sha256"] == asset["sha256"]
        )
        try:
            with pyBigWig.open(str(path)) as bigwig:
                chroms = bigwig.chroms()
                header = bigwig.header()
                assembly_ok = all(
                    chroms.get(chromosome) == length
                    for chromosome, length in HG38_SENTINEL_LENGTHS.items()
                )
                is_bigwig = bool(bigwig.isBigWig())
                record.update(
                    {
                        "is_bigwig": is_bigwig,
                        "bigwig_version": header.get("version"),
                        "chromosome_count": len(chroms),
                        "n_bases_covered": header.get("nBasesCovered"),
                        "min_value": header.get("minVal"),
                        "max_value": header.get("maxVal"),
                        "hg38_primary_lengths_match": assembly_ok,
                        "status": (
                            "PASS"
                            if integrity_ok and is_bigwig and assembly_ok
                            else "FAIL_INTEGRITY_OR_FORMAT"
                        ),
                    }
                )
        except Exception as error:  # pragma: no cover - exercised on malformed external files
            record["status"] = "FAIL_UNREADABLE"
            record["error"] = repr(error)
        records.append(record)

    all_passed = (
        len(records) == 12
        and all(record.get("status") == "PASS" for record in records)
        and not missing
        and not unexpected
    )
    write_json(
        {
            "schema_version": "cisgrammarshift_gphn_bigwig_audit_v1",
            "processing_pipeline": manifest["processing_pipeline"],
            "data_directory": str(args.bigwig_directory),
            "expected_file_count": len(expected),
            "observed_bw_count": len(observed),
            "all_passed": all_passed,
            "missing_files": missing,
            "unexpected_files": unexpected,
            "files": records,
        },
        args.output,
    )
    if not all_passed:
        raise SystemExit("GPHN bigWig audit failed; see the output report")


if __name__ == "__main__":
    main()
