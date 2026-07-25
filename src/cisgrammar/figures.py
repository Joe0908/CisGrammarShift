from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np
import pandas as pd


def plot_metric(
    summary: pd.DataFrame,
    metric: str,
    output_path: Path,
) -> None:
    mean_column = f"{metric}_mean"
    standard_deviation_column = f"{metric}_std"
    if mean_column not in summary:
        raise KeyError(f"summary does not contain {mean_column}")

    models = list(summary["model"].drop_duplicates())
    conditions = list(summary["condition"].drop_duplicates())
    x = np.arange(len(conditions), dtype=float)
    width = min(0.18, 0.8 / max(len(models), 1))

    figure, axis = plt.subplots(figsize=(max(8, len(conditions) * 1.6), 5))
    for model_index, model_name in enumerate(models):
        model_rows = summary.set_index(["model", "condition"]).loc[model_name]
        means = np.array([model_rows.loc[condition, mean_column] for condition in conditions])
        errors = np.array([model_rows.loc[condition, standard_deviation_column] for condition in conditions])
        errors = np.nan_to_num(errors, nan=0.0)
        offset = (model_index - (len(models) - 1) / 2) * width
        axis.bar(
            x + offset,
            means,
            width=width,
            yerr=errors,
            capsize=3,
            label=model_name,
        )

    axis.axhline(0.5, color="black", linestyle="--", linewidth=1, alpha=0.6)
    axis.set_xticks(x, conditions, rotation=20, ha="right")
    axis.set_ylim(0.0, 1.0)
    axis.set_ylabel(metric.replace("_", " ").title())
    axis.set_title(f"{metric.replace('_', ' ').title()} across distribution shifts")
    axis.legend(frameon=False, ncols=min(4, len(models)))
    axis.spines[["top", "right"]].set_visible(False)
    figure.tight_layout()
    figure.savefig(output_path, dpi=200)
    plt.close(figure)


def write_standard_figures(summary: pd.DataFrame, output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    plot_metric(summary, "auroc", output_dir / "auroc_by_condition.png")
    plot_metric(summary, "pairwise_accuracy", output_dir / "pairwise_accuracy_by_condition.png")
