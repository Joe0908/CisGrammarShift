from __future__ import annotations

from pathlib import Path

import pandas as pd

from cisgrammar.capselex import normalize_tf_symbol


def read_all_excel_sheets(path: str | Path) -> pd.DataFrame:
    sheets = pd.read_excel(path, sheet_name=None)
    frames = []
    for sheet, frame in sheets.items():
        copy = frame.copy()
        copy.insert(0, "source_sheet", sheet)
        frames.append(copy)
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def infer_tf_column(frame: pd.DataFrame) -> str:
    for column in frame.columns:
        lowered = str(column).lower().strip()
        if lowered in {"tf", "gene", "gene name", "gene_name", "symbol", "tf name"}:
            return str(column)
    raise ValueError("metadata table has no recognizable TF/gene column")


def normalized_tf_inventory(path: str | Path) -> pd.DataFrame:
    frame = read_all_excel_sheets(path)
    tf_column = infer_tf_column(frame)
    result = frame.copy()
    result["tf"] = result[tf_column].map(normalize_tf_symbol)
    return result[result["tf"].ne("") & result["tf"].ne("NAN")].reset_index(drop=True)


def audit_three_way_overlap(
    cap_tfs: list[str],
    plasmid_metadata: str | Path,
    chip_metadata: str | Path,
    ght_metadata: str | Path,
) -> dict[str, object]:
    cap = {normalize_tf_symbol(tf) for tf in cap_tfs}
    plasmid = set(normalized_tf_inventory(plasmid_metadata)["tf"])
    chip = set(normalized_tf_inventory(chip_metadata)["tf"])
    ght = set(normalized_tf_inventory(ght_metadata)["tf"])
    overlap = sorted(cap & plasmid & chip & ght)
    return {
        "cap_tfs": len(cap),
        "plasmid_tfs": len(plasmid),
        "chip_tfs": len(chip),
        "ght_tfs": len(ght),
        "four_way_tfs": overlap,
        "four_way_count": len(overlap),
    }
