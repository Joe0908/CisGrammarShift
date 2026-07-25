from __future__ import annotations

from typing import Any

import torch
import torch.nn as nn


class LocalCNN(nn.Module):
    """Motif detector with a deliberately limited interaction receptive field."""

    def __init__(self, channels: int = 64, kernel_size: int = 19):
        super().__init__()
        if kernel_size % 2 == 0:
            raise ValueError("kernel_size must be odd")
        self.features = nn.Sequential(
            nn.Conv1d(4, channels, kernel_size=kernel_size, padding=kernel_size // 2),
            nn.BatchNorm1d(channels),
            nn.GELU(),
        )
        self.pool = nn.AdaptiveMaxPool1d(1)
        self.classifier = nn.Linear(channels, 1)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        hidden = self.features(x.transpose(1, 2))
        pooled = self.pool(hidden).squeeze(-1)
        return self.classifier(pooled).squeeze(-1)


class ResidualDilatedBlock(nn.Module):
    def __init__(self, channels: int, dilation: int, dropout: float = 0.1):
        super().__init__()
        self.block = nn.Sequential(
            nn.Conv1d(channels, channels, kernel_size=3, padding=dilation, dilation=dilation),
            nn.BatchNorm1d(channels),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Conv1d(channels, channels, kernel_size=1),
        )
        self.activation = nn.GELU()

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return self.activation(x + self.block(x))


class DilatedCNN(nn.Module):
    """Long-range convolutional model with an explicit sequence-wide receptive field."""

    def __init__(
        self,
        channels: int = 64,
        kernel_size: int = 19,
        dilations: list[int] | tuple[int, ...] = (1, 2, 4, 8, 16, 32),
        dropout: float = 0.1,
    ):
        super().__init__()
        if kernel_size % 2 == 0:
            raise ValueError("kernel_size must be odd")
        self.stem = nn.Sequential(
            nn.Conv1d(4, channels, kernel_size=kernel_size, padding=kernel_size // 2),
            nn.BatchNorm1d(channels),
            nn.GELU(),
        )
        self.blocks = nn.Sequential(
            *(ResidualDilatedBlock(channels, int(dilation), dropout) for dilation in dilations)
        )
        self.pool = nn.AdaptiveMaxPool1d(1)
        self.classifier = nn.Sequential(nn.Dropout(dropout), nn.Linear(channels, 1))

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        hidden = self.stem(x.transpose(1, 2))
        hidden = self.blocks(hidden)
        pooled = self.pool(hidden).squeeze(-1)
        return self.classifier(pooled).squeeze(-1)


class SequenceTransformer(nn.Module):
    def __init__(
        self,
        sequence_length: int,
        embedding_dim: int = 64,
        attention_heads: int = 4,
        layers: int = 3,
        feedforward_dim: int = 128,
        dropout: float = 0.1,
    ):
        super().__init__()
        self.sequence_length = sequence_length
        self.embedding = nn.Linear(4, embedding_dim)
        self.cls_token = nn.Parameter(torch.zeros(1, 1, embedding_dim))
        self.position = nn.Parameter(torch.empty(1, sequence_length + 1, embedding_dim))
        encoder_layer = nn.TransformerEncoderLayer(
            d_model=embedding_dim,
            nhead=attention_heads,
            dim_feedforward=feedforward_dim,
            dropout=dropout,
            activation="gelu",
            batch_first=True,
            norm_first=True,
        )
        self.encoder = nn.TransformerEncoder(encoder_layer, num_layers=layers)
        self.normalise = nn.LayerNorm(embedding_dim)
        self.classifier = nn.Linear(embedding_dim, 1)
        nn.init.trunc_normal_(self.position, std=0.02)
        nn.init.trunc_normal_(self.cls_token, std=0.02)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        if x.shape[1] != self.sequence_length:
            raise ValueError(f"expected sequence length {self.sequence_length}, received {x.shape[1]}")
        tokens = self.embedding(x)
        cls = self.cls_token.expand(x.shape[0], -1, -1)
        hidden = torch.cat([cls, tokens], dim=1) + self.position
        encoded = self.encoder(hidden)
        return self.classifier(self.normalise(encoded[:, 0])).squeeze(-1)


def build_model(name: str, sequence_length: int, model_config: dict[str, Any]) -> nn.Module:
    if name == "local_cnn":
        return LocalCNN(**model_config)
    if name == "dilated_cnn":
        return DilatedCNN(**model_config)
    if name == "transformer":
        return SequenceTransformer(sequence_length=sequence_length, **model_config)
    raise KeyError(f"unknown neural model {name!r}")


def count_parameters(model: nn.Module) -> int:
    return sum(parameter.numel() for parameter in model.parameters() if parameter.requires_grad)
