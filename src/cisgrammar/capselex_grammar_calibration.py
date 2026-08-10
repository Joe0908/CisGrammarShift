from __future__ import annotations

import numpy as np
import pandas as pd

from cisgrammar.capselex_grammar_scoring import binned_monomer_residuals


def cross_fitted_residualize(
    frame: pd.DataFrame,
    raw_feature: str,
    focal_monomer: str,
    partner_monomer: str,
    fold_column: str = "chromosome_fold",
    bins: int = 20,
) -> np.ndarray:
    output = np.full(len(frame), np.nan)
    folds = frame[fold_column].to_numpy(dtype=int)
    grammar = frame[raw_feature].to_numpy(dtype=float)
    focal = frame[focal_monomer].to_numpy(dtype=float)
    partner = frame[partner_monomer].to_numpy(dtype=float)
    for fold in sorted(np.unique(folds)):
        train = folds != fold
        residuals = binned_monomer_residuals(grammar, focal, partner, train, bins=bins)
        output[folds == fold] = residuals[folds == fold]
    return output


def build_ets_grammar_score(
    frame: pd.DataFrame,
    profile_columns: list[str],
    focal_monomer: str,
    partner_monomer_columns: list[str],
) -> pd.Series:
    if len(profile_columns) != len(partner_monomer_columns):
        raise ValueError("each grammar profile needs its partner monomer control")
    residuals = []
    for profile, partner in zip(profile_columns, partner_monomer_columns, strict=True):
        residuals.append(cross_fitted_residualize(frame, profile, focal_monomer, partner))
    return pd.Series(np.mean(np.column_stack(residuals), axis=1), index=frame.index, name="ets_grammar")
