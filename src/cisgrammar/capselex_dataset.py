from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Literal

import numpy as np
import pandas as pd

from cisgrammar.capselex import normalize_tf_symbol, require_columns, stable_seed

TARGETS = ("cooperative_signal", "composite_motif", "spacing_or_orientation")
SplitMode = Literal["pair", "tf", "family_pair"]


@dataclass(frozen=True)
class DatasetSummary:
    rows: int
    positive_pairs: int
    unique_tfs: int
    directed_asymmetries: int


def _as_binary(values: pd.Series) -> pd.Series:
    truthy = {"1", "true", "yes", "y", "+", "positive", "observed"}
    return values.map(lambda value: int(str(value).strip().lower() in truthy)).astype("int8")


def read_construct_table(path: str | Path) -> pd.DataFrame:
    frame = (
        pd.read_excel(path)
        if str(path).endswith((".xlsx", ".xls"))
        else pd.read_csv(path, sep=None, engine="python")
    )
    lowered = {column.lower().strip(): column for column in frame.columns}
    symbol = next((lowered[key] for key in ("tf", "gene", "symbol", "gene_symbol") if key in lowered), None)
    if symbol is None:
        raise ValueError("construct table needs a TF/gene/symbol column")
    family = next(
        (lowered[key] for key in ("structural_family", "family", "structural_class") if key in lowered),
        None,
    )
    sequence = next(
        (lowered[key] for key in ("protein_sequence", "sequence", "edbd_sequence") if key in lowered),
        None,
    )
    result = pd.DataFrame({"tf": frame[symbol].map(normalize_tf_symbol)})
    result["family"] = frame[family].fillna("unknown").astype(str) if family else "unknown"
    result["protein_sequence"] = frame[sequence].fillna("").astype(str) if sequence else ""
    return result.drop_duplicates("tf", keep="first").reset_index(drop=True)


def read_interaction_table(path: str | Path) -> pd.DataFrame:
    """Read either a long directed-pair table or a wide bait-by-prey matrix."""
    frame = (
        pd.read_excel(path)
        if str(path).endswith((".xlsx", ".xls"))
        else pd.read_csv(path, sep=None, engine="python")
    )
    lowered = {column.lower().strip(): column for column in frame.columns}
    bait = next((lowered[key] for key in ("bait", "row_tf", "tf1") if key in lowered), None)
    prey = next((lowered[key] for key in ("prey", "column_tf", "tf2", "partner") if key in lowered), None)
    if bait and prey:
        result = pd.DataFrame(
            {"bait": frame[bait].map(normalize_tf_symbol), "prey": frame[prey].map(normalize_tf_symbol)}
        )
        aliases = {
            "cooperative_signal": ("cooperative_signal", "observed", "interaction", "positive"),
            "composite_motif": ("composite_motif", "composite", "composite_pwm"),
            "spacing_or_orientation": ("spacing_or_orientation", "spacing", "gapped_kmer"),
        }
        for target, candidates in aliases.items():
            source = next((lowered[key] for key in candidates if key in lowered), None)
            result[target] = _as_binary(frame[source]) if source else 0
        return result.drop_duplicates(["bait", "prey"], keep="last").reset_index(drop=True)

    if len(frame.columns) < 2:
        raise ValueError("interaction matrix requires a row TF column and at least one prey column")
    row_column = frame.columns[0]
    melted = frame.melt(id_vars=row_column, var_name="prey", value_name="cooperative_signal")
    melted = melted.rename(columns={row_column: "bait"})
    melted["bait"] = melted["bait"].map(normalize_tf_symbol)
    melted["prey"] = melted["prey"].map(normalize_tf_symbol)
    melted["cooperative_signal"] = _as_binary(melted["cooperative_signal"])
    melted["composite_motif"] = 0
    melted["spacing_or_orientation"] = 0
    return melted


def build_directed_pair_dataset(interactions: pd.DataFrame, constructs: pd.DataFrame) -> pd.DataFrame:
    require_columns(interactions, ("bait", "prey", *TARGETS), "interaction table")
    require_columns(constructs, ("tf", "family", "protein_sequence"), "construct table")
    annotations = constructs.set_index("tf")
    result = interactions.copy()
    result["bait"] = result["bait"].map(normalize_tf_symbol)
    result["prey"] = result["prey"].map(normalize_tf_symbol)
    result = result[result["bait"].isin(annotations.index) & result["prey"].isin(annotations.index)].copy()
    result["bait_family"] = result["bait"].map(annotations["family"])
    result["prey_family"] = result["prey"].map(annotations["family"])
    result["bait_sequence"] = result["bait"].map(annotations["protein_sequence"])
    result["prey_sequence"] = result["prey"].map(annotations["protein_sequence"])
    result["directed_pair"] = result["bait"] + "->" + result["prey"]
    result["family_pair"] = result["bait_family"] + "->" + result["prey_family"]
    for target in TARGETS:
        result[target] = pd.to_numeric(result[target], errors="coerce").fillna(0).clip(0, 1).astype("int8")
    return result.sort_values(["bait", "prey"], kind="mergesort").reset_index(drop=True)


def summarize_dataset(frame: pd.DataFrame) -> DatasetSummary:
    reverse = frame[["bait", "prey", "cooperative_signal"]].rename(
        columns={"bait": "prey", "prey": "bait", "cooperative_signal": "reverse_signal"}
    )
    paired = frame.merge(reverse, on=["bait", "prey"], how="inner")
    return DatasetSummary(
        rows=len(frame),
        positive_pairs=int(frame["cooperative_signal"].sum()),
        unique_tfs=len(set(frame["bait"]) | set(frame["prey"])),
        directed_asymmetries=int((paired["cooperative_signal"] != paired["reverse_signal"]).sum() // 2),
    )


def assign_split(
    frame: pd.DataFrame,
    mode: SplitMode,
    seed: int,
    heldout_fraction: float = 0.2,
) -> pd.Series:
    """Assign leakage-resistant train/test labels; never split sequencing reads."""
    if not 0 < heldout_fraction < 1:
        raise ValueError("heldout_fraction must lie between zero and one")
    if mode == "pair":
        groups = frame["directed_pair"]
    elif mode == "family_pair":
        groups = frame["family_pair"]
    elif mode == "tf":
        tfs = sorted(set(frame["bait"]) | set(frame["prey"]))
        rng = np.random.default_rng(stable_seed(seed, mode))
        n_test = max(1, round(len(tfs) * heldout_fraction))
        test_tfs = set(rng.choice(tfs, size=n_test, replace=False))
        return pd.Series(
            np.where(frame["bait"].isin(test_tfs) | frame["prey"].isin(test_tfs), "test", "train"),
            index=frame.index,
            name="split",
        )
    else:
        raise ValueError(f"unsupported split mode: {mode}")
    unique = sorted(groups.unique())
    rng = np.random.default_rng(stable_seed(seed, mode))
    n_test = max(1, round(len(unique) * heldout_fraction))
    test_groups = set(rng.choice(unique, size=n_test, replace=False))
    return groups.map(lambda value: "test" if value in test_groups else "train").rename("split")


def leave_one_tf_out(frame: pd.DataFrame, heldout_tf: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    tf = normalize_tf_symbol(heldout_tf)
    mask = (frame["bait"] == tf) | (frame["prey"] == tf)
    if not mask.any():
        raise ValueError(f"held-out TF {tf} is absent")
    return frame.loc[~mask].copy(), frame.loc[mask].copy()
