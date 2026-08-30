from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pandas as pd

from cisgrammar.capselex import AUTOSOMES
from cisgrammar.capselex_genomic_assets import MAGIX_COLUMNS, read_magix


def file_digest(path: Path, algorithm: str = "sha256") -> str:
    digest = hashlib.new(algorithm)
    with path.open("rb") as handle:
        while chunk := handle.read(1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def audit_codebook_v2_focal_panel(
    config_path: str | Path,
    data_directory: str | Path,
) -> dict[str, object]:
    config_path = Path(config_path)
    data_directory = Path(data_directory)
    config = json.loads(config_path.read_text(encoding="utf-8"))
    call_rule = config["call_rule"]
    records: list[dict[str, object]] = []

    for asset in config["panel"]:
        path = data_directory / asset["local_filename"]
        observed_sha256 = file_digest(path)
        if path.stat().st_size != asset["size_bytes"]:
            raise RuntimeError(f"Codebook v2 MAGIX size mismatch: {path}")
        if observed_sha256 != asset["sha256"]:
            raise RuntimeError(f"Codebook v2 MAGIX SHA-256 mismatch: {path}")

        raw = pd.read_csv(path, sep="\t")
        missing = [column for column in MAGIX_COLUMNS if column not in raw.columns]
        if missing:
            raise ValueError(f"{path.name} is missing MAGIX columns: {', '.join(missing)}")
        numeric_start = pd.to_numeric(raw["start"], errors="raise")
        numeric_stop = pd.to_numeric(raw["stop"], errors="raise")
        widths = numeric_stop - numeric_start
        coordinates = raw[["chr", "start", "stop"]]
        expected_name = (
            raw["chr"].astype(str)
            + ":"
            + numeric_start.astype(str)
            + "-"
            + numeric_stop.astype(str)
        )
        selected = read_magix(
            path,
            fdr_max=float(call_rule["fdr_max"]),
            require_positive=bool(call_rule["require_positive_score"]),
        )
        records.append(
            {
                **asset,
                "observed_sha256": observed_sha256,
                "sha256_matches_manifest": True,
                "rows_total": int(len(raw)),
                "rows_autosomal": int(raw["chr"].isin(AUTOSOMES).sum()),
                "rows_selected": int(len(selected)),
                "interval_width_min": int(widths.min()),
                "interval_width_max": int(widths.max()),
                "duplicate_coordinates": int(coordinates.duplicated().sum()),
                "coordinate_name_matches": bool(raw["name"].astype(str).eq(expected_name).all()),
                "byte_identical_to_geo_member": (
                    observed_sha256 == asset.get("geo_member_sha256")
                ),
            }
        )

    all_valid = all(
        record["sha256_matches_manifest"]
        and record["interval_width_min"] == 200
        and record["interval_width_max"] == 200
        and record["duplicate_coordinates"] == 0
        and record["coordinate_name_matches"]
        for record in records
    )
    return {
        "schema_version": "codebook_ght_v2_focal_magix_audit_v1",
        "source": config["source"],
        "call_rule": call_rule,
        "all_files_valid": bool(all_valid),
        "primary_analysis_eligible": bool(
            config["source"]["primary_analysis_eligible"] and all_valid
        ),
        "focal_members_byte_identical_to_geo": all(
            record["byte_identical_to_geo_member"] for record in records
        ),
        "interpretation": (
            "The current official Codebook v2 web archive is a different tar archive from the GEO "
            "snapshot, but its six focal BED members are byte-identical to the corresponding GEO "
            "members. The official v2 archive provenance therefore resolves release eligibility "
            "without changing the focal numerical data."
        ),
        "files": records,
    }
