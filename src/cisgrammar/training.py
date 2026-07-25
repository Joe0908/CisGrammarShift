from __future__ import annotations

import copy
import random
from dataclasses import dataclass

import numpy as np
import torch
import torch.nn as nn
from torch.utils.data import DataLoader, TensorDataset


@dataclass
class TrainingResult:
    model: nn.Module
    best_validation_loss: float
    epochs_completed: int


def seed_everything(seed: int, deterministic: bool = True) -> None:
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)
    if deterministic:
        torch.use_deterministic_algorithms(True, warn_only=True)
        torch.backends.cudnn.benchmark = False


def make_loader(
    x: np.ndarray,
    y: np.ndarray,
    batch_size: int,
    shuffle: bool,
    seed: int,
) -> DataLoader:
    generator = torch.Generator()
    generator.manual_seed(seed)
    dataset = TensorDataset(torch.from_numpy(x).float(), torch.from_numpy(y).float())
    return DataLoader(dataset, batch_size=batch_size, shuffle=shuffle, generator=generator)


def _mean_loss(model: nn.Module, loader: DataLoader, loss_function: nn.Module, device: torch.device) -> float:
    model.eval()
    losses: list[float] = []
    with torch.no_grad():
        for x_batch, y_batch in loader:
            logits = model(x_batch.to(device))
            loss = loss_function(logits, y_batch.to(device))
            losses.append(float(loss.detach().cpu()))
    return float(np.mean(losses))


def fit_model(
    model: nn.Module,
    train_loader: DataLoader,
    validation_loader: DataLoader,
    *,
    epochs: int,
    patience: int,
    learning_rate: float,
    weight_decay: float,
    device: torch.device,
) -> TrainingResult:
    model = model.to(device)
    optimiser = torch.optim.AdamW(
        model.parameters(),
        lr=learning_rate,
        weight_decay=weight_decay,
    )
    loss_function = nn.BCEWithLogitsLoss()
    best_loss = float("inf")
    best_state: dict[str, torch.Tensor] | None = None
    stale_epochs = 0
    epochs_completed = 0

    for epoch in range(1, epochs + 1):
        model.train()
        for x_batch, y_batch in train_loader:
            optimiser.zero_grad(set_to_none=True)
            logits = model(x_batch.to(device))
            loss = loss_function(logits, y_batch.to(device))
            if not torch.isfinite(loss):
                raise FloatingPointError(f"non-finite training loss at epoch {epoch}")
            loss.backward()
            optimiser.step()

        validation_loss = _mean_loss(model, validation_loader, loss_function, device)
        epochs_completed = epoch
        if validation_loss < best_loss - 1e-7:
            best_loss = validation_loss
            best_state = copy.deepcopy(model.state_dict())
            stale_epochs = 0
        else:
            stale_epochs += 1
            if stale_epochs >= patience:
                break

    if best_state is None:
        raise RuntimeError("training completed without a finite validation state")
    model.load_state_dict(best_state)
    return TrainingResult(model=model, best_validation_loss=best_loss, epochs_completed=epochs_completed)


def predict_probabilities(
    model: nn.Module,
    x: np.ndarray,
    batch_size: int,
    device: torch.device,
) -> np.ndarray:
    model.eval()
    probabilities: list[np.ndarray] = []
    loader = DataLoader(torch.from_numpy(x).float(), batch_size=batch_size, shuffle=False)
    with torch.no_grad():
        for x_batch in loader:
            logits = model(x_batch.to(device))
            probabilities.append(torch.sigmoid(logits).cpu().numpy())
    return np.concatenate(probabilities)
