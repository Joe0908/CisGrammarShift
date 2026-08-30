from __future__ import annotations

from collections.abc import Callable

import numpy as np
import pandas as pd

from cisgrammar.capselex_dataset import TARGETS, leave_one_tf_out
from cisgrammar.capselex_pair_baselines import (
    ProteinKmerBaseline,
    TFIdentityBaseline,
    binary_metrics,
    family_rate_baseline,
)


def eligible_panel_tfs(
    frame: pd.DataFrame,
    min_observed_positives: int = 10,
    min_composite_positives: int = 5,
    min_spacing_positives: int = 5,
) -> list[str]:
    thresholds = dict(
        cooperative_signal=min_observed_positives,
        composite_motif=min_composite_positives,
        spacing_or_orientation=min_spacing_positives,
    )
    candidates = sorted(set(frame["bait"]) | set(frame["prey"]))
    eligible = []
    for tf in candidates:
        subset = frame[(frame["bait"] == tf) | (frame["prey"] == tf)]
        if all(int(subset[target].sum()) >= threshold for target, threshold in thresholds.items()):
            eligible.append(tf)
    return eligible


def _model_factories() -> dict[str, Callable[[], object]]:
    return {
        "family_rate": family_rate_baseline,
        "tf_identity": TFIdentityBaseline,
        "protein_kmer": ProteinKmerBaseline,
    }


def run_leave_one_tf_out(frame: pd.DataFrame, heldout_tf: str, target: str) -> dict[str, object]:
    if target not in TARGETS:
        raise ValueError(f"target must be one of {TARGETS}")
    train, test = leave_one_tf_out(frame, heldout_tf)
    labels = test[target].to_numpy(dtype=int)
    results: dict[str, object] = {
        "heldout_tf": heldout_tf,
        "target": target,
        "train_rows": len(train),
        "test_rows": len(test),
        "models": {},
    }
    for name, factory in _model_factories().items():
        if train[target].nunique() < 2:
            continue
        model = factory()
        model.fit(train, target)
        scores = model.predict_proba(test)
        results["models"][name] = binary_metrics(labels, scores)
    return results


def bootstrap_mean_interval(
    values: list[float],
    seed: int = 20260809,
    repetitions: int = 2000,
) -> dict[str, float]:
    clean = np.asarray([value for value in values if np.isfinite(value)], dtype=float)
    if clean.size == 0:
        return {"mean": float("nan"), "lower": float("nan"), "upper": float("nan")}
    rng = np.random.default_rng(seed)
    samples = rng.choice(clean, size=(repetitions, clean.size), replace=True).mean(axis=1)
    return {
        "mean": float(clean.mean()),
        "lower": float(np.quantile(samples, 0.025)),
        "upper": float(np.quantile(samples, 0.975)),
    }


def run_tf_panel(frame: pd.DataFrame, target: str) -> dict[str, object]:
    tfs = eligible_panel_tfs(frame)
    per_tf = [run_leave_one_tf_out(frame, tf, target) for tf in tfs]
    summary: dict[str, object] = {}
    for model_name in _model_factories():
        improvements = []
        for result in per_tf:
            model = result["models"].get(model_name)
            if model is not None:
                improvements.append(float(model["auprc_minus_prevalence"]))
        summary[model_name] = bootstrap_mean_interval(improvements)
    return {"target": target, "eligible_tfs": tfs, "per_tf": per_tf, "summary": summary}
