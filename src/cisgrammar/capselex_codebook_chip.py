from __future__ import annotations

import gzip
import re
from pathlib import Path

from cisgrammar.capselex import asset_manifest


def _sample_records(soft_path: str | Path) -> list[dict[str, str]]:
    records: list[dict[str, str]] = []
    current: dict[str, str] | None = None
    with gzip.open(soft_path, "rt", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if line.startswith("^SAMPLE = "):
                if current is not None:
                    records.append(current)
                current = {"gsm": line.split(" = ", 1)[1]}
            elif current is not None and line.startswith("!Sample_title = "):
                current["title"] = line.split(" = ", 1)[1]
            elif current is not None and line.startswith("!Sample_description = "):
                current["description"] = line.split(" = ", 1)[1]
            elif current is not None and line.startswith("!Sample_supplementary_file_1 = "):
                current["bigwig_url"] = line.split(" = ", 1)[1].replace("ftp://", "https://")
        if current is not None:
            records.append(current)
    return records


def _file_sizes(filelist_path: str | Path) -> dict[str, int]:
    sizes: dict[str, int] = {}
    with Path(filelist_path).open(encoding="utf-8") as handle:
        for raw in handle:
            fields = raw.rstrip("\n").split("\t")
            if len(fields) == 5 and fields[0] == "File":
                sizes[fields[1]] = int(fields[3])
    return sizes


def build_gpzn_manifest(
    soft_path: str | Path,
    filelist_path: str | Path,
    focal_tfs: list[str] | tuple[str, ...],
) -> dict[str, object]:
    """Resolve two Toronto GPZN ChIP bigWigs per focal TF from GEO metadata."""
    focal = set(focal_tfs)
    sizes = _file_sizes(filelist_path)
    assets: list[dict[str, object]] = []
    for sample in _sample_records(soft_path):
        match = re.search(r"analysis of ([^ ]+) genomic target", sample.get("title", ""))
        if match is None or match.group(1) not in focal:
            continue
        tf = match.group(1)
        url = sample.get("bigwig_url", "")
        filename = url.rsplit("/", 1)[-1]
        replicate_match = re.search(r"ChIP(\d+)", sample.get("description", ""), flags=re.IGNORECASE)
        if replicate_match is None:
            raise ValueError(f"cannot resolve replicate number for {sample['gsm']}")
        if "GPZN" not in filename or not filename.endswith("_hg38_PE.bw"):
            raise ValueError(f"focal sample is not a Toronto GPZN hg38 PE bigWig: {filename}")
        if filename not in sizes:
            raise ValueError(f"GEO file list has no size for {filename}")
        assets.append(
            {
                "tf": tf,
                "replicate": int(replicate_match.group(1)),
                "gsm": sample["gsm"],
                "sample_description": sample["description"],
                "filename": filename,
                "url": url,
                "size_bytes": sizes[filename],
                "processing_pipeline": "Toronto_GPZN",
                "assembly": "hg38",
                "layout": "paired_end",
            }
        )
    assets.sort(key=lambda record: (str(record["tf"]), int(record["replicate"])))

    observed_tfs = {str(record["tf"]) for record in assets}
    missing = sorted(focal - observed_tfs)
    unexpected = sorted(observed_tfs - focal)
    if missing or unexpected:
        raise ValueError(f"TF panel mismatch; missing={missing}, unexpected={unexpected}")
    for tf in sorted(focal):
        records = [record for record in assets if record["tf"] == tf]
        replicates = [record["replicate"] for record in records]
        if replicates != [1, 2]:
            raise ValueError(f"{tf} requires exactly ChIP replicates 1 and 2; found {replicates}")
    if len({record["gsm"] for record in assets}) != len(assets):
        raise ValueError("GEO sample accessions must be unique")
    if len({record["filename"] for record in assets}) != len(assets):
        raise ValueError("GEO bigWig filenames must be unique")

    return {
        "schema_version": "codebook_chip_gpzn_manifest_v1",
        "accession": "GSE280248",
        "selection_rule": "two_biological_replicates_per_tf",
        "processing_pipeline": "Toronto_GPZN_only",
        "total_size_bytes": sum(int(record["size_bytes"]) for record in assets),
        "source_metadata": asset_manifest([soft_path, filelist_path]),
        "assets": assets,
    }
