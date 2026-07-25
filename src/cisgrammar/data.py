from __future__ import annotations

import hashlib
from dataclasses import dataclass
from typing import Any

import numpy as np

from cisgrammar.motifs import Motif, reverse_complement

BASES = np.array(list("ACGT"))
BASE_TO_INDEX = {base: index for index, base in enumerate(BASES)}


@dataclass(frozen=True)
class GrammarRule:
    period: int = 10
    allowed_phases: tuple[int, ...] = (0, 1, 9)

    def __post_init__(self) -> None:
        if self.period <= 1:
            raise ValueError("period must be greater than one")
        phases = tuple(sorted({int(phase) % self.period for phase in self.allowed_phases}))
        if not phases or len(phases) == self.period:
            raise ValueError("allowed phases must define both positive and negative gaps")
        object.__setattr__(self, "allowed_phases", phases)

    def is_positive(self, gap: int) -> bool:
        return int(gap) % self.period in self.allowed_phases


@dataclass(frozen=True)
class Condition:
    gap_min: int
    gap_max: int
    background_gc: float = 0.5
    background_persistence: float = 0.15
    allow_reverse: bool = False
    motif_temperature: float = 1.0

    def __post_init__(self) -> None:
        if self.gap_min < 0 or self.gap_max < self.gap_min:
            raise ValueError("gap range must satisfy 0 <= gap_min <= gap_max")
        if not 0 < self.background_gc < 1:
            raise ValueError("background_gc must be between zero and one")
        if not 0 <= self.background_persistence < 1:
            raise ValueError("background_persistence must be in [0, 1)")
        if self.motif_temperature <= 0:
            raise ValueError("motif_temperature must be positive")

    @classmethod
    def from_mapping(cls, mapping: dict[str, Any]) -> Condition:
        return cls(**mapping)


@dataclass
class SequenceDataset:
    x: np.ndarray
    y: np.ndarray
    pair_ids: np.ndarray
    motif_mask: np.ndarray
    records: list[dict[str, Any]]

    def validate(self, rule: GrammarRule) -> None:
        n = len(self.y)
        if self.x.shape[0] != n or self.pair_ids.shape[0] != n or self.motif_mask.shape[0] != n:
            raise ValueError("dataset fields have inconsistent row counts")
        if self.x.ndim != 3 or self.x.shape[2] != 4:
            raise ValueError(f"x must have shape N×L×4, received {self.x.shape}")
        if self.motif_mask.shape != self.x.shape[:2]:
            raise ValueError("motif mask does not match sequence dimensions")
        labels, counts = np.unique(self.y, return_counts=True)
        if labels.tolist() != [0.0, 1.0] or counts[0] != counts[1]:
            raise ValueError("dataset must be exactly balanced")
        unique_pairs, pair_counts = np.unique(self.pair_ids, return_counts=True)
        if np.any(pair_counts != 2) or len(unique_pairs) * 2 != n:
            raise ValueError("each pair ID must occur exactly twice")
        by_pair: dict[str, list[dict[str, Any]]] = {}
        for record in self.records:
            by_pair.setdefault(record["pair_id"], []).append(record)
            expected = rule.is_positive(record["gap"])
            if bool(record["label"]) != expected:
                raise ValueError("record label violates the configured grammar rule")
        for pair_records in by_pair.values():
            if sorted(record["label"] for record in pair_records) != [0, 1]:
                raise ValueError("each pair must contain one positive and one negative")
            first, second = pair_records
            controlled_fields = (
                "motif_a_instance",
                "motif_b_instance",
                "orientation_a",
                "orientation_b",
                "background_sha256",
            )
            for field in controlled_fields:
                if first[field] != second[field]:
                    raise ValueError(f"matched pair differs in controlled field {field}")


def derive_seed(base_seed: int, *parts: str) -> int:
    payload = "::".join([str(base_seed), *parts]).encode()
    return int.from_bytes(hashlib.sha256(payload).digest()[:8], "little") % (2**32)


def one_hot_encode(sequences: list[str]) -> np.ndarray:
    if not sequences:
        raise ValueError("cannot encode an empty sequence list")
    length = len(sequences[0])
    encoded = np.zeros((len(sequences), length, 4), dtype=np.float32)
    for row, sequence in enumerate(sequences):
        if len(sequence) != length:
            raise ValueError("all sequences must have equal length")
        for column, base in enumerate(sequence):
            try:
                encoded[row, column, BASE_TO_INDEX[base]] = 1.0
            except KeyError as exc:
                raise ValueError(f"unsupported base {base!r}") from exc
    return encoded


def sample_background(
    rng: np.random.Generator,
    length: int,
    gc: float,
    persistence: float,
) -> str:
    probabilities = np.array([(1 - gc) / 2, gc / 2, gc / 2, (1 - gc) / 2], dtype=np.float64)
    sequence = [str(rng.choice(BASES, p=probabilities))]
    for _ in range(1, length):
        if rng.random() < persistence:
            sequence.append(sequence[-1])
        else:
            sequence.append(str(rng.choice(BASES, p=probabilities)))
    return "".join(sequence)


def _implant(
    background: str,
    motif_a: str,
    motif_b: str,
    start_a: int,
    start_b: int,
) -> tuple[str, np.ndarray]:
    sequence = list(background)
    mask = np.zeros(len(background), dtype=bool)
    end_a = start_a + len(motif_a)
    end_b = start_b + len(motif_b)
    if end_a > start_b or end_b > len(background):
        raise ValueError("motif intervals overlap or extend outside the sequence")
    sequence[start_a:end_a] = motif_a
    sequence[start_b:end_b] = motif_b
    mask[start_a:end_a] = True
    mask[start_b:end_b] = True
    return "".join(sequence), mask


def generate_matched_dataset(
    *,
    n_pairs: int,
    sequence_length: int,
    motif_a: Motif,
    motif_b: Motif,
    rule: GrammarRule,
    condition: Condition,
    anchor_min: int,
    anchor_max: int,
    seed: int,
    condition_name: str,
) -> SequenceDataset:
    if n_pairs <= 0:
        raise ValueError("n_pairs must be positive")
    if not 0 <= anchor_min <= anchor_max:
        raise ValueError("invalid anchor range")

    positive_gaps = [gap for gap in range(condition.gap_min, condition.gap_max + 1) if rule.is_positive(gap)]
    negative_gaps = [
        gap for gap in range(condition.gap_min, condition.gap_max + 1) if not rule.is_positive(gap)
    ]
    if not positive_gaps or not negative_gaps:
        raise ValueError("condition gap range must contain both positive and negative grammar phases")

    worst_end = anchor_max + motif_a.length + condition.gap_max + motif_b.length
    if worst_end > sequence_length:
        raise ValueError(
            f"sequence length {sequence_length} is too short; maximum motif end would be {worst_end}"
        )

    rng = np.random.default_rng(seed)
    sequences: list[str] = []
    labels: list[float] = []
    pair_ids: list[str] = []
    masks: list[np.ndarray] = []
    records: list[dict[str, Any]] = []

    for pair_index in range(n_pairs):
        background = sample_background(
            rng,
            sequence_length,
            condition.background_gc,
            condition.background_persistence,
        )
        background_hash = hashlib.sha256(background.encode()).hexdigest()
        instance_a = motif_a.sample(rng, condition.motif_temperature)
        instance_b = motif_b.sample(rng, condition.motif_temperature)
        orientation_a = "-" if condition.allow_reverse and rng.random() < 0.5 else "+"
        orientation_b = "-" if condition.allow_reverse and rng.random() < 0.5 else "+"
        inserted_a = reverse_complement(instance_a) if orientation_a == "-" else instance_a
        inserted_b = reverse_complement(instance_b) if orientation_b == "-" else instance_b
        start_a = int(rng.integers(anchor_min, anchor_max + 1))
        positive_gap = int(rng.choice(positive_gaps))
        negative_gap = int(rng.choice(negative_gaps))
        pair_id = f"{condition_name}-{seed}-{pair_index:07d}"

        for label, gap in ((1, positive_gap), (0, negative_gap)):
            start_b = start_a + motif_a.length + gap
            sequence, motif_mask = _implant(background, inserted_a, inserted_b, start_a, start_b)
            sequences.append(sequence)
            labels.append(float(label))
            pair_ids.append(pair_id)
            masks.append(motif_mask)
            records.append(
                {
                    "pair_id": pair_id,
                    "condition": condition_name,
                    "label": label,
                    "gap": gap,
                    "phase": gap % rule.period,
                    "start_a": start_a,
                    "end_a": start_a + motif_a.length,
                    "start_b": start_b,
                    "end_b": start_b + motif_b.length,
                    "orientation_a": orientation_a,
                    "orientation_b": orientation_b,
                    "motif_a_instance": instance_a,
                    "motif_b_instance": instance_b,
                    "background_sha256": background_hash,
                }
            )

    order = rng.permutation(len(sequences))
    dataset = SequenceDataset(
        x=one_hot_encode(sequences)[order],
        y=np.asarray(labels, dtype=np.float32)[order],
        pair_ids=np.asarray(pair_ids)[order],
        motif_mask=np.asarray(masks, dtype=bool)[order],
        records=[records[index] for index in order],
    )
    dataset.validate(rule)
    return dataset
