from __future__ import annotations

from pathlib import Path

import pandas as pd

from cisgrammar.capselex import chromosome_fold
from cisgrammar.capselex_genomic_model import NestedResult, chromosome_nested_binary


def read_target_genes(path: str | Path, state: str, sheet_name: int | str = 0) -> set[str]:
    frame = pd.read_excel(path, sheet_name=sheet_name)
    lowered = {str(column).lower(): column for column in frame.columns}
    gene_column = next((lowered[key] for key in ("gene", "gene_name", "symbol") if key in lowered), None)
    if gene_column is None:
        raise ValueError("target table requires a gene column")
    state_column = next((column for column in frame.columns if "state" in str(column).lower()), None)
    if state_column is not None:
        frame = frame[frame[state_column].astype(str).str.upper() == state.upper()]
    return set(frame[gene_column].dropna().astype(str).str.strip())


def run_target_benchmark(
    gene_features: pd.DataFrame,
    targets: set[str],
    baseline_features: list[str],
    grammar_feature: str = "ets_grammar",
) -> NestedResult:
    frame = gene_features.copy()
    frame["target"] = frame["gene"].isin(targets).astype("int8")
    if "chromosome_fold" not in frame:
        frame["chromosome_fold"] = frame["chrom"].map(chromosome_fold)
    return chromosome_nested_binary(frame, "target", baseline_features, [grammar_feature])
