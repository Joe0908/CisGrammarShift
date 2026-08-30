from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import tarfile
from pathlib import Path

import pandas as pd

from cisgrammar.capselex import AUTOSOMES, write_json
from cisgrammar.capselex_genomic_assets import MAGIX_COLUMNS, read_magix


def file_digest(path: Path, algorithm: str) -> str:
    digest = hashlib.new(algorithm)
    with path.open("rb") as handle:
        while chunk := handle.read(1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def expected_md5(metadata_path: Path) -> dict[str, str]:
    frame = pd.read_excel(
        metadata_path,
        sheet_name="MD5 Checksums",
        skiprows=8,
        header=None,
        usecols=[0, 1],
        names=["filename", "md5"],
    ).dropna()
    return dict(zip(frame["filename"].astype(str), frame["md5"].astype(str), strict=True))


def extract_exact_panel(archive: Path, filenames: set[str], output_directory: Path) -> None:
    output_directory.mkdir(parents=True, exist_ok=True)
    extracted: set[str] = set()
    # Streaming reads the compressed archive exactly once. Repeated random seeks through
    # a multi-gigabyte gzip member are slow and can leave truncated files if interrupted.
    with tarfile.open(archive, "r|gz") as handle:
        for member in handle:
            basename = Path(member.name).name
            if member.isfile() and basename in filenames:
                if basename in extracted:
                    raise RuntimeError(f"archive contains duplicate panel file: {basename}")
                source = handle.extractfile(member)
                if source is None:
                    raise RuntimeError(f"cannot extract archive member: {member.name}")
                destination = output_directory / basename
                temporary = destination.with_suffix(destination.suffix + ".part")
                with source, temporary.open("wb") as target:
                    shutil.copyfileobj(source, target)
                temporary.replace(destination)
                extracted.add(basename)
    missing = sorted(filenames - extracted)
    if missing:
        raise RuntimeError(f"archive is missing panel files: {', '.join(missing)}")


def audit_file(path: Path, expected: str, fdr_max: float, require_positive: bool) -> dict[str, object]:
    raw = pd.read_csv(path, sep="\t")
    missing = [column for column in MAGIX_COLUMNS if column not in raw.columns]
    if missing:
        raise ValueError(f"{path.name} is missing MAGIX columns: {', '.join(missing)}")
    widths = pd.to_numeric(raw["stop"], errors="raise") - pd.to_numeric(
        raw["start"], errors="raise"
    )
    selected = read_magix(path, fdr_max=fdr_max, require_positive=require_positive)
    observed = file_digest(path, "md5")
    return {
        "filename": path.name,
        "size_bytes": path.stat().st_size,
        "expected_md5": expected,
        "observed_md5": observed,
        "md5_matches_metadata": observed == expected,
        "rows_total": int(len(raw)),
        "rows_autosomal": int(raw["chr"].isin(AUTOSOMES).sum()),
        "rows_selected": int(len(selected)),
        "interval_width_min": int(widths.min()),
        "interval_width_max": int(widths.max()),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Extract and audit a frozen focal MAGIX panel")
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--metadata-directory", type=Path, required=True)
    parser.add_argument("--archive", type=Path, required=True)
    parser.add_argument("--extraction-directory", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    config = json.loads(args.config.read_text(encoding="utf-8"))
    source = config["source"]
    if args.archive.stat().st_size != source["archive_size_bytes"]:
        raise RuntimeError("GEO MAGIX archive size does not match the frozen manifest")
    archive_sha256 = file_digest(args.archive, "sha256")
    if archive_sha256 != source["archive_sha256"]:
        raise RuntimeError("GEO MAGIX archive SHA-256 does not match the frozen manifest")

    filenames = {record["filename"] for record in config["panel"]}
    extract_exact_panel(args.archive, filenames, args.extraction_directory)
    metadata_path = args.metadata_directory / source["checksum_metadata_filename"]
    checksums = expected_md5(metadata_path)
    call_rule = config["call_rule"]
    records = []
    for panel_record in config["panel"]:
        filename = panel_record["filename"]
        if filename not in checksums:
            raise RuntimeError(f"metadata has no checksum for {filename}")
        record = audit_file(
            args.extraction_directory / filename,
            checksums[filename],
            call_rule["fdr_max"],
            call_rule["require_positive_score"],
        )
        records.append({**panel_record, **record})

    mismatches = [record["tf"] for record in records if not record["md5_matches_metadata"]]
    write_json(
        {
            "schema_version": "codebook_magix_asset_audit_v1",
            "source": source,
            "call_rule": call_rule,
            "archive_observed_sha256": archive_sha256,
            "all_files_match_metadata": not mismatches,
            "checksum_mismatch_tfs": mismatches,
            "primary_analysis_eligible": bool(source["primary_analysis_eligible"] and not mismatches),
            "files": records,
        },
        args.output,
    )


if __name__ == "__main__":
    main()
