from __future__ import annotations

import hashlib
import json
import platform
import sys
from functools import partial
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import torch
import yaml

from cisgrammar import __version__
from cisgrammar.baselines import PWMPresenceBaseline
from cisgrammar.data import (
    Condition,
    GrammarRule,
    SequenceDataset,
    derive_seed,
    generate_matched_dataset,
)
from cisgrammar.figures import write_standard_figures
from cisgrammar.interpretation import attribution_localisation
from cisgrammar.metrics import binary_metrics, paired_metrics, select_threshold
from cisgrammar.models import build_model, count_parameters
from cisgrammar.motifs import Motif, get_motif
from cisgrammar.training import (
    fit_model,
    make_loader,
    predict_probabilities,
    seed_everything,
)


def _build_dataset(
    *,
    data_config: dict[str, Any],
    condition_mapping: dict[str, Any],
    n_pairs: int,
    seed: int,
    condition_name: str,
    motif_a: Motif,
    motif_b: Motif,
    rule: GrammarRule,
) -> SequenceDataset:
    return generate_matched_dataset(
        n_pairs=n_pairs,
        sequence_length=int(data_config["sequence_length"]),
        motif_a=motif_a,
        motif_b=motif_b,
        rule=rule,
        condition=Condition.from_mapping(condition_mapping),
        anchor_min=int(data_config["anchor_min"]),
        anchor_max=int(data_config["anchor_max"]),
        seed=seed,
        condition_name=condition_name,
    )


def _save_manifest(dataset: SequenceDataset, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(dataset.records).to_csv(path, index=False)


def _evaluate(
    *,
    y: np.ndarray,
    probabilities: np.ndarray,
    pair_ids: np.ndarray,
    threshold: float,
) -> dict[str, float]:
    return {
        **binary_metrics(y, probabilities, threshold),
        **paired_metrics(y, probabilities, pair_ids),
    }


def _resolve_device(requested: str) -> torch.device:
    if requested == "auto":
        return torch.device("cuda" if torch.cuda.is_available() else "cpu")
    device = torch.device(requested)
    if device.type == "cuda" and not torch.cuda.is_available():
        raise RuntimeError("CUDA was requested but is unavailable")
    return device


def _write_metadata(
    output_dir: Path,
    config_text: str,
    motif_a: Motif,
    motif_b: Motif,
    device: torch.device,
) -> None:
    metadata = {
        "cisgrammar_version": __version__,
        "config_sha256": hashlib.sha256(config_text.encode()).hexdigest(),
        "python": sys.version,
        "platform": platform.platform(),
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "torch": torch.__version__,
        "device": str(device),
        "motifs": [
            {
                "name": motif.name,
                "accession": motif.accession,
                "source": motif.source,
                "source_url": motif.source_url,
            }
            for motif in (motif_a, motif_b)
        ],
    }
    (output_dir / "run_metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")


def _summarise(metrics: pd.DataFrame) -> pd.DataFrame:
    identity_columns = {"model", "condition", "seed"}
    numeric_columns = [
        column
        for column in metrics.select_dtypes(include=[np.number]).columns
        if column not in identity_columns
    ]
    summary = metrics.groupby(["model", "condition"], sort=True)[numeric_columns].agg(["mean", "std"])
    summary.columns = [f"{name}_{statistic}" for name, statistic in summary.columns]
    return summary.reset_index()


def run_experiment(
    config_path: Path,
    output_dir: Path,
    *,
    requested_device: str = "auto",
) -> pd.DataFrame:
    config_text = config_path.read_text()
    config = yaml.safe_load(config_text)
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "config.snapshot.yaml").write_text(yaml.safe_dump(config, sort_keys=False))

    data_config = config["dataset"]
    experiment_config = config["experiment"]
    training_config = config["training"]
    model_section = config["models"]
    motif_a = get_motif(data_config["motif_a"])
    motif_b = get_motif(data_config["motif_b"])
    rule = GrammarRule(
        period=int(data_config["rule"]["period"]),
        allowed_phases=tuple(int(value) for value in data_config["rule"]["allowed_phases"]),
    )
    device = _resolve_device(requested_device)
    _write_metadata(output_dir, config_text, motif_a, motif_b, device)

    checkpoint_dir = output_dir / "checkpoints"
    checkpoint_dir.mkdir(exist_ok=True)
    rows: list[dict[str, Any]] = []

    for seed in experiment_config["seeds"]:
        seed = int(seed)
        print(f"\n=== Seed {seed} ===")
        train_data = _build_dataset(
            data_config=data_config,
            condition_mapping=data_config["train"],
            n_pairs=int(data_config["train_pairs"]),
            seed=derive_seed(seed, "train"),
            condition_name="train",
            motif_a=motif_a,
            motif_b=motif_b,
            rule=rule,
        )
        validation_data = _build_dataset(
            data_config=data_config,
            condition_mapping=data_config["train"],
            n_pairs=int(data_config["validation_pairs"]),
            seed=derive_seed(seed, "validation"),
            condition_name="validation",
            motif_a=motif_a,
            motif_b=motif_b,
            rule=rule,
        )
        evaluation_data = {
            condition_name: _build_dataset(
                data_config=data_config,
                condition_mapping=condition_mapping,
                n_pairs=int(data_config["test_pairs"]),
                seed=derive_seed(seed, "evaluation", condition_name),
                condition_name=condition_name,
                motif_a=motif_a,
                motif_b=motif_b,
                rule=rule,
            )
            for condition_name, condition_mapping in data_config["evaluations"].items()
        }

        _save_manifest(train_data, output_dir / "manifests" / f"seed_{seed}_train.csv")
        _save_manifest(validation_data, output_dir / "manifests" / f"seed_{seed}_validation.csv")
        for condition_name, dataset in evaluation_data.items():
            _save_manifest(dataset, output_dir / "manifests" / f"seed_{seed}_{condition_name}.csv")

        for model_name in model_section["names"]:
            print(f"{model_name}: fitting and evaluating")
            if model_name == "pwm_presence":
                baseline = PWMPresenceBaseline(
                    motif_a,
                    motif_b,
                    background_gc=float(data_config["train"]["background_gc"]),
                ).fit(train_data.x)
                validation_probabilities = baseline.predict_proba(validation_data.x)
                threshold = select_threshold(validation_data.y, validation_probabilities)
                predictor = baseline.predict_proba
                neural_model = None
                parameter_count = 0
                best_validation_loss = float("nan")
                epochs_completed = 0
            else:
                seed_everything(
                    derive_seed(seed, "model", model_name),
                    deterministic=bool(experiment_config.get("deterministic", True)),
                )
                neural_model = build_model(
                    model_name,
                    sequence_length=int(data_config["sequence_length"]),
                    model_config=dict(model_section.get(model_name, {})),
                )
                parameter_count = count_parameters(neural_model)
                train_loader = make_loader(
                    train_data.x,
                    train_data.y,
                    batch_size=int(training_config["batch_size"]),
                    shuffle=True,
                    seed=derive_seed(seed, "loader", model_name),
                )
                validation_loader = make_loader(
                    validation_data.x,
                    validation_data.y,
                    batch_size=int(training_config["batch_size"]),
                    shuffle=False,
                    seed=derive_seed(seed, "validation_loader", model_name),
                )
                trained = fit_model(
                    neural_model,
                    train_loader,
                    validation_loader,
                    epochs=int(training_config["epochs"]),
                    patience=int(training_config["patience"]),
                    learning_rate=float(training_config["learning_rate"]),
                    weight_decay=float(training_config["weight_decay"]),
                    device=device,
                )
                neural_model = trained.model
                best_validation_loss = trained.best_validation_loss
                epochs_completed = trained.epochs_completed
                validation_probabilities = predict_probabilities(
                    neural_model,
                    validation_data.x,
                    batch_size=int(training_config["batch_size"]),
                    device=device,
                )
                threshold = select_threshold(validation_data.y, validation_probabilities)

                predictor = partial(
                    predict_probabilities,
                    neural_model,
                    batch_size=int(training_config["batch_size"]),
                    device=device,
                )

                torch.save(
                    {
                        "state_dict": neural_model.state_dict(),
                        "model_name": model_name,
                        "model_config": model_section.get(model_name, {}),
                        "sequence_length": int(data_config["sequence_length"]),
                        "seed": seed,
                        "threshold": threshold,
                    },
                    checkpoint_dir / f"{model_name}_seed_{seed}.pt",
                )

            validation_scores = binary_metrics(
                validation_data.y,
                validation_probabilities,
                threshold,
            )
            for condition_name, dataset in evaluation_data.items():
                probabilities = predictor(dataset.x)
                prediction_path = (
                    output_dir / "predictions" / f"{model_name}_seed_{seed}_{condition_name}.csv"
                )
                prediction_path.parent.mkdir(parents=True, exist_ok=True)
                prediction_frame = pd.DataFrame(dataset.records)
                prediction_frame["probability"] = probabilities
                prediction_frame.to_csv(prediction_path, index=False)
                metrics = _evaluate(
                    y=dataset.y,
                    probabilities=probabilities,
                    pair_ids=dataset.pair_ids,
                    threshold=threshold,
                )
                if neural_model is None:
                    interpretation = {
                        "attribution_mass_mean": float("nan"),
                        "attribution_mass_std": float("nan"),
                        "sequence_coverage": float(np.mean(dataset.motif_mask)),
                    }
                else:
                    interpretation = attribution_localisation(
                        neural_model,
                        dataset.x,
                        dataset.motif_mask,
                        max_examples=int(training_config["attribution_examples"]),
                        device=device,
                    )
                rows.append(
                    {
                        "seed": seed,
                        "model": model_name,
                        "condition": condition_name,
                        "threshold": threshold,
                        "parameter_count": parameter_count,
                        "epochs_completed": epochs_completed,
                        "best_validation_loss": best_validation_loss,
                        "validation_auroc": validation_scores["auroc"],
                        **metrics,
                        **interpretation,
                    }
                )

    metrics_frame = pd.DataFrame(rows)
    metrics_frame.to_csv(output_dir / "metrics.csv", index=False)
    summary = _summarise(metrics_frame)
    summary.to_csv(output_dir / "summary.csv", index=False)
    write_standard_figures(summary, output_dir / "figures")
    print(f"\nSaved {len(metrics_frame)} evaluation rows to {output_dir / 'metrics.csv'}")
    return metrics_frame
