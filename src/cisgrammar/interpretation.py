from __future__ import annotations

import numpy as np
import torch
import torch.nn as nn


def attribution_localisation(
    model: nn.Module,
    x: np.ndarray,
    motif_mask: np.ndarray,
    *,
    max_examples: int,
    device: torch.device,
    batch_size: int = 64,
) -> dict[str, float]:
    n_examples = min(max_examples, len(x))
    if n_examples <= 0:
        return {
            "attribution_mass_mean": float("nan"),
            "attribution_mass_std": float("nan"),
            "sequence_coverage": float("nan"),
        }

    model.eval()
    fractions: list[np.ndarray] = []
    for start in range(0, n_examples, batch_size):
        stop = min(start + batch_size, n_examples)
        inputs = torch.from_numpy(x[start:stop]).float().to(device)
        inputs.requires_grad_(True)
        model.zero_grad(set_to_none=True)
        logits = model(inputs)
        gradients = torch.autograd.grad(logits.sum(), inputs)[0]
        attribution = torch.abs(gradients * inputs).sum(dim=-1)
        masks = torch.from_numpy(motif_mask[start:stop]).to(device)
        inside = (attribution * masks).sum(dim=1)
        total = attribution.sum(dim=1).clamp_min(1e-12)
        fractions.append((inside / total).detach().cpu().numpy())

    values = np.concatenate(fractions)
    return {
        "attribution_mass_mean": float(np.mean(values)),
        "attribution_mass_std": float(np.std(values, ddof=1)) if len(values) > 1 else 0.0,
        "sequence_coverage": float(np.mean(motif_mask[:n_examples])),
    }
