from __future__ import annotations

import numpy as np
import pandas as pd

from cisgrammar.capselex import require_columns
from cisgrammar.capselex_trophoblast_hichip import choose_strongest_linked_locus


def assemble_locus_features(
    loci: pd.DataFrame,
    grammar: pd.DataFrame,
    monomers: pd.DataFrame,
    ght: pd.DataFrame | None = None,
) -> pd.DataFrame:
    require_columns(loci, ("locus_id", "chrom", "start", "end"), "loci")
    result = loci.merge(grammar, on="locus_id", how="left", validate="one_to_one")
    result = result.merge(monomers, on="locus_id", how="left", validate="one_to_one")
    if ght is not None:
        result = result.merge(ght, on="locus_id", how="left", validate="one_to_one")
    numeric = result.select_dtypes(include=["number"]).columns
    result[numeric] = result[numeric].replace([np.inf, -np.inf], np.nan).fillna(0.0)
    return result


def assemble_linked_gene_features(
    links: pd.DataFrame,
    locus_features: pd.DataFrame,
    promoters: pd.DataFrame,
    wt_expression: pd.DataFrame,
) -> pd.DataFrame:
    selected = choose_strongest_linked_locus(links, locus_features)
    selected = selected.merge(promoters[["gene", "chrom", "tss"]], on=["gene", "chrom"], how="left")
    selected["tss_distance"] = np.where(
        selected["tss"].notna(),
        np.abs(selected["midpoint"].to_numpy(dtype=float) - selected["tss"].to_numpy(dtype=float)),
        0.0,
    )
    expression = wt_expression[["gene", "wt_count"]].drop_duplicates("gene")
    selected = selected.merge(expression, on="gene", how="left")
    selected["wt_count"] = selected["wt_count"].fillna(0.0)
    selected["gcm1_chip"] = selected["mean_log1p_chip"]
    return selected
