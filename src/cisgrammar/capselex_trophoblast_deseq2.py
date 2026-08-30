from __future__ import annotations

import re
from pathlib import Path

import pandas as pd

from cisgrammar.capselex import audit_asset, write_json

EXPECTED_STATES = ("EVT", "ST")


def _find_gene_column(frame: pd.DataFrame) -> str:
    lowered = {str(column).lower().strip(): str(column) for column in frame.columns}
    for candidate in ("gene", "gene_name", "symbol", "gene symbol"):
        if candidate in lowered:
            return lowered[candidate]
    raise ValueError("RNA table needs a gene-symbol column")


def read_dataset_s4(path: str | Path, sheet_name: int | str = 0) -> pd.DataFrame:
    frame = pd.read_excel(path, sheet_name=sheet_name)
    gene_column = _find_gene_column(frame)
    frame = frame.rename(columns={gene_column: "gene"})
    frame["gene"] = frame["gene"].astype(str).str.strip()
    return frame


def infer_count_columns(frame: pd.DataFrame, state: str) -> tuple[list[str], list[str]]:
    state_upper = state.upper()
    candidates = [column for column in frame.columns if state_upper in str(column).upper()]
    wt = [column for column in candidates if re.search(r"(^|[_ -])WT([_ -]|$)", str(column), re.I)]
    ko = [column for column in candidates if re.search(r"GCM1.*KO|KO.*GCM1|KNOCK.?OUT", str(column), re.I)]
    if len(wt) != 4 or len(ko) != 2:
        message = (
            f"{state} requires exactly four WT and two GCM1-KO count columns; "
            f"got {len(wt)} and {len(ko)}"
        )
        raise ValueError(message)
    return wt, ko


def export_deseq2_inputs(
    dataset: str | Path,
    output_directory: str | Path,
    raw_rsem_filter_audit: str | Path | None = None,
) -> dict[str, object]:
    output = Path(output_directory)
    output.mkdir(parents=True, exist_ok=True)
    frame = read_dataset_s4(dataset)
    manifest: dict[str, object] = {"schema_version": "gcm1_deseq2_input_v1", "states": {}, "assets": []}
    for state in EXPECTED_STATES:
        wt, ko = infer_count_columns(frame, state)
        columns = wt + ko
        counts = frame[["gene", *columns]].copy()
        for column in columns:
            counts[column] = pd.to_numeric(counts[column], errors="raise").round().astype("int64")
        counts = counts.drop_duplicates("gene").set_index("gene")
        metadata = pd.DataFrame(
            {
                "sample": columns,
                "condition": ["WT"] * len(wt) + ["GCM1_KO"] * len(ko),
                "state": state,
            }
        ).set_index("sample")
        count_path = output / f"{state}.counts.tsv"
        metadata_path = output / f"{state}.metadata.tsv"
        counts.to_csv(count_path, sep="\t")
        metadata.to_csv(metadata_path, sep="\t")
        manifest["states"][state] = {"genes": len(counts), "wt": len(wt), "ko": len(ko)}
        manifest["assets"].extend([audit_asset(count_path).__dict__, audit_asset(metadata_path).__dict__])
    if raw_rsem_filter_audit is not None:
        raw = pd.read_csv(raw_rsem_filter_audit, sep="\t")
        length_column = next(column for column in raw.columns if column.lower() == "length")
        retained = int((pd.to_numeric(raw[length_column], errors="coerce") >= 300).sum())
        manifest["raw_filter_audit"] = {
            "raw_genes": len(raw),
            "removed_length_lt_300": len(raw) - retained,
            "retained": retained,
        }
    write_json(manifest, output / "manifest.json")
    return manifest


def validate_deseq2_session_info(path: str | Path) -> dict[str, str]:
    text = Path(path).read_text(encoding="utf-8")
    if "R version 4.2.2" not in text:
        raise ValueError("analysis requires R 4.2.2")
    match = re.search(r"DESeq2_([0-9.]+)", text)
    if match is None or match.group(1) != "1.36.0":
        raise ValueError("analysis requires DESeq2 1.36.0")
    return {"r": "4.2.2", "deseq2": "1.36.0"}


def read_deseq2_results(path: str | Path) -> pd.DataFrame:
    frame = pd.read_csv(path, sep="\t")
    required = {"gene", "baseMean", "log2FoldChange", "lfcSE", "pvalue", "padj"}
    missing = required - set(frame.columns)
    if missing:
        raise ValueError(f"DESeq2 result missing columns: {sorted(missing)}")
    frame["gcm1_dependence"] = -pd.to_numeric(frame["log2FoldChange"], errors="coerce")
    return frame
