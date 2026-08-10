from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from itertools import product
from typing import Protocol

import numpy as np
import pandas as pd
from sklearn.feature_extraction import DictVectorizer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler


class PairModel(Protocol):
    def fit(self, frame: pd.DataFrame, target: str) -> PairModel: ...

    def predict_proba(self, frame: pd.DataFrame) -> np.ndarray: ...


@dataclass
class SmoothedRateBaseline:
    group_columns: tuple[str, ...]
    alpha: float = 1.0
    beta: float = 1.0

    def fit(self, frame: pd.DataFrame, target: str) -> SmoothedRateBaseline:
        self.target = target
        self.global_rate = (frame[target].sum() + self.alpha) / (len(frame) + self.alpha + self.beta)
        grouped = frame.groupby(list(self.group_columns), dropna=False)[target].agg(["sum", "count"])
        self.rates = ((grouped["sum"] + self.alpha) / (grouped["count"] + self.alpha + self.beta)).to_dict()
        return self

    def predict_proba(self, frame: pd.DataFrame) -> np.ndarray:
        keys = frame[list(self.group_columns)].itertuples(index=False, name=None)
        return np.array([self.rates.get(tuple(key), self.global_rate) for key in keys], dtype=float)


class TFIdentityBaseline:
    def __init__(self, regularization: float = 1.0) -> None:
        self.model = make_pipeline(
            DictVectorizer(sparse=True),
            LogisticRegression(C=regularization, max_iter=2000, class_weight="balanced"),
        )

    @staticmethod
    def _records(frame: pd.DataFrame) -> list[dict[str, str]]:
        return [
            {"bait": str(bait), "prey": str(prey), "bait_family": str(bf), "prey_family": str(pf)}
            for bait, prey, bf, pf in frame[["bait", "prey", "bait_family", "prey_family"]].itertuples(
                index=False, name=None
            )
        ]

    def fit(self, frame: pd.DataFrame, target: str) -> TFIdentityBaseline:
        self.model.fit(self._records(frame), frame[target].to_numpy(dtype=int))
        return self

    def predict_proba(self, frame: pd.DataFrame) -> np.ndarray:
        return self.model.predict_proba(self._records(frame))[:, 1]


def protein_kmer_vector(sequence: str, k: int = 2) -> np.ndarray:
    alphabet = "ACDEFGHIKLMNPQRSTVWY"
    vocabulary = {"".join(chars): index for index, chars in enumerate(product(alphabet, repeat=k))}
    counts = np.zeros(len(vocabulary), dtype=np.float32)
    cleaned = "".join(residue for residue in str(sequence).upper() if residue in alphabet)
    for offset in range(max(0, len(cleaned) - k + 1)):
        counts[vocabulary[cleaned[offset : offset + k]]] += 1
    total = counts.sum()
    return counts / total if total else counts


class ProteinKmerBaseline:
    def __init__(self, k: int = 2, regularization: float = 1.0) -> None:
        self.k = k
        self.model = make_pipeline(
            StandardScaler(),
            LogisticRegression(C=regularization, max_iter=2000, class_weight="balanced"),
        )

    def _matrix(self, frame: pd.DataFrame) -> np.ndarray:
        bait = np.stack([protein_kmer_vector(sequence, self.k) for sequence in frame["bait_sequence"]])
        prey = np.stack([protein_kmer_vector(sequence, self.k) for sequence in frame["prey_sequence"]])
        return np.concatenate([bait, prey, np.abs(bait - prey), bait * prey], axis=1)

    def fit(self, frame: pd.DataFrame, target: str) -> ProteinKmerBaseline:
        self.model.fit(self._matrix(frame), frame[target].to_numpy(dtype=int))
        return self

    def predict_proba(self, frame: pd.DataFrame) -> np.ndarray:
        return self.model.predict_proba(self._matrix(frame))[:, 1]


def binary_metrics(labels: np.ndarray, scores: np.ndarray) -> dict[str, float]:
    labels = np.asarray(labels, dtype=int)
    scores = np.asarray(scores, dtype=float)
    prevalence = float(labels.mean())
    auprc = float(average_precision_score(labels, scores)) if labels.sum() else float("nan")
    auroc = float(roc_auc_score(labels, scores)) if np.unique(labels).size == 2 else float("nan")
    return {
        "genes_or_pairs": int(labels.size),
        "positives": int(labels.sum()),
        "prevalence": prevalence,
        "auprc": auprc,
        "auroc": auroc,
        "auprc_minus_prevalence": auprc - prevalence,
        "auprc_over_prevalence": auprc / prevalence if prevalence else float("nan"),
    }


def family_rate_baseline() -> SmoothedRateBaseline:
    return SmoothedRateBaseline(("bait_family", "prey_family"))


def partner_frequency(frame: pd.DataFrame, target: str) -> Counter[str]:
    positives = frame.loc[frame[target] == 1]
    return Counter(positives["prey"])
