from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr
from sklearn.linear_model import LogisticRegression, Ridge
from sklearn.metrics import (
    average_precision_score,
    log_loss,
    mean_absolute_error,
    mean_squared_error,
    roc_auc_score,
)
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

RIDGE_ALPHAS = (0.01, 0.1, 1.0, 10.0, 100.0)


@dataclass(frozen=True)
class NestedResult:
    predictions: pd.DataFrame
    summary: dict[str, object]


def _continuous_metrics(y: np.ndarray, prediction: np.ndarray) -> dict[str, float]:
    mse = float(mean_squared_error(y, prediction))
    return {
        "mse": mse,
        "rmse": float(np.sqrt(mse)),
        "mae": float(mean_absolute_error(y, prediction)),
        "r2": float(1 - np.sum((y - prediction) ** 2) / np.sum((y - y.mean()) ** 2)),
        "pearson": float(pearsonr(y, prediction).statistic),
        "spearman": float(spearmanr(y, prediction).statistic),
    }


def _binary_metrics(y: np.ndarray, prediction: np.ndarray) -> dict[str, float]:
    return {
        "log_loss": float(log_loss(y, prediction, labels=[0, 1])),
        "auprc": float(average_precision_score(y, prediction)),
        "auroc": float(roc_auc_score(y, prediction)),
    }


def _select_ridge_alpha(x: np.ndarray, y: np.ndarray, validation: np.ndarray, train: np.ndarray) -> float:
    losses = {}
    for alpha in RIDGE_ALPHAS:
        model = make_pipeline(StandardScaler(), Ridge(alpha=alpha))
        model.fit(x[train], y[train])
        losses[alpha] = mean_squared_error(y[validation], model.predict(x[validation]))
    return float(min(losses, key=losses.get))


def chromosome_nested_continuous(
    frame: pd.DataFrame,
    outcome: str,
    baseline_features: list[str],
    addition_features: list[str],
    fold_column: str = "chromosome_fold",
) -> NestedResult:
    y = frame[outcome].to_numpy(dtype=float)
    folds = frame[fold_column].to_numpy(dtype=int)
    x0 = frame[baseline_features].to_numpy(dtype=float)
    x1 = frame[baseline_features + addition_features].to_numpy(dtype=float)
    pred0 = np.full(len(frame), np.nan)
    pred1 = np.full(len(frame), np.nan)
    fold_results = []
    for fold in sorted(np.unique(folds)):
        test = folds == fold
        validation = folds == ((fold + 1) % (int(folds.max()) + 1))
        train = ~(test | validation)
        alpha0 = _select_ridge_alpha(x0, y, validation, train)
        alpha1 = _select_ridge_alpha(x1, y, validation, train)
        model0 = make_pipeline(StandardScaler(), Ridge(alpha=alpha0)).fit(x0[train], y[train])
        model1 = make_pipeline(StandardScaler(), Ridge(alpha=alpha1)).fit(x1[train], y[train])
        pred0[test] = model0.predict(x0[test])
        pred1[test] = model1.predict(x1[test])
        ridge = model1.named_steps["ridge"]
        fold_results.append(
            {
                "outer_fold": int(fold),
                "genes_or_loci": int(test.sum()),
                "m0_alpha": alpha0,
                "m1_alpha": alpha1,
                "standardized_addition_coefficients": ridge.coef_[-len(addition_features) :].tolist(),
                "m0": _continuous_metrics(y[test], pred0[test]),
                "m1": _continuous_metrics(y[test], pred1[test]),
            }
        )
    metrics0 = _continuous_metrics(y, pred0)
    metrics1 = _continuous_metrics(y, pred1)
    partial_r2 = 1 - metrics1["mse"] / metrics0["mse"]
    predictions = frame[[fold_column]].copy()
    predictions["outcome"] = y
    predictions["m0_prediction"] = pred0
    predictions["m1_prediction"] = pred1
    return NestedResult(
        predictions,
        {
            "rows": len(frame),
            "m0": metrics0,
            "m1": metrics1,
            "partial_r2": float(partial_r2),
            "fold_results": fold_results,
        },
    )


def chromosome_nested_binary(
    frame: pd.DataFrame,
    outcome: str,
    baseline_features: list[str],
    addition_features: list[str],
    fold_column: str = "chromosome_fold",
    regularization: float = 1.0,
) -> NestedResult:
    y = frame[outcome].to_numpy(dtype=int)
    folds = frame[fold_column].to_numpy(dtype=int)
    matrices = [
        frame[baseline_features].to_numpy(dtype=float),
        frame[baseline_features + addition_features].to_numpy(dtype=float),
    ]
    predictions = [np.full(len(frame), np.nan), np.full(len(frame), np.nan)]
    fold_results = []
    for fold in sorted(np.unique(folds)):
        test = folds == fold
        train = ~test
        fold_payload: dict[str, object] = {"outer_fold": int(fold), "genes_or_loci": int(test.sum())}
        for index, matrix in enumerate(matrices):
            model = make_pipeline(
                StandardScaler(),
                LogisticRegression(C=regularization, max_iter=2000, class_weight="balanced"),
            )
            model.fit(matrix[train], y[train])
            predictions[index][test] = model.predict_proba(matrix[test])[:, 1]
            fold_payload[f"m{index}"] = _binary_metrics(y[test], predictions[index][test])
        fold_results.append(fold_payload)
    metrics0 = _binary_metrics(y, predictions[0])
    metrics1 = _binary_metrics(y, predictions[1])
    output = frame[[fold_column]].copy()
    output["outcome"] = y
    output["m0_prediction"] = predictions[0]
    output["m1_prediction"] = predictions[1]
    return NestedResult(
        output,
        {
            "rows": len(frame),
            "positives": int(y.sum()),
            "m0": metrics0,
            "m1": metrics1,
            "m1_minus_m0": {key: metrics1[key] - metrics0[key] for key in metrics0},
            "fold_results": fold_results,
        },
    )


def chromosome_shift_null(
    frame: pd.DataFrame,
    columns: list[str],
    chromosome_column: str = "chrom",
    seed: int = 20260809,
) -> pd.DataFrame:
    """Circularly shift grammar features within chromosomes while keeping context fixed."""
    rng = np.random.default_rng(seed)
    result = frame.copy()
    for _, indices in result.groupby(chromosome_column, sort=False).groups.items():
        positions = np.asarray(list(indices), dtype=int)
        if positions.size < 2:
            continue
        shift = int(rng.integers(1, positions.size))
        result.loc[positions, columns] = np.roll(result.loc[positions, columns].to_numpy(), shift, axis=0)
    return result
