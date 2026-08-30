from __future__ import annotations

import re
import zipfile
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from cisgrammar.capselex import reverse_complement

DNA = "ACGT"
DNA_INDEX = {base: index for index, base in enumerate(DNA)}


@dataclass(frozen=True)
class PWM:
    name: str
    probabilities: np.ndarray

    def __post_init__(self) -> None:
        matrix = np.asarray(self.probabilities, dtype=float)
        if matrix.ndim != 2 or matrix.shape[1] != 4:
            raise ValueError("PWM must have shape length x 4 in A,C,G,T order")
        if np.any(matrix < 0):
            raise ValueError("PWM probabilities cannot be negative")
        normalized = (matrix + 1e-6) / (matrix + 1e-6).sum(axis=1, keepdims=True)
        object.__setattr__(self, "probabilities", normalized)

    @property
    def length(self) -> int:
        return self.probabilities.shape[0]

    @property
    def log_odds(self) -> np.ndarray:
        return np.log2(self.probabilities / 0.25)


@dataclass(frozen=True)
class MonomerPWMRecord:
    source_id: str
    symbol: str
    pwm: PWM


def read_jaspar(path: str | Path, name: str | None = None) -> PWM:
    rows: dict[str, list[float]] = {}
    motif_name = name or Path(path).stem
    with Path(path).open(encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if line.startswith(">") and name is None:
                motif_name = line[1:].strip().split(maxsplit=1)[-1]
            elif line and line[0] in DNA:
                values = line.split("[", 1)[-1].split("]", 1)[0].replace(",", " ")
                rows[line[0]] = [float(value) for value in values.split()]
    if set(rows) != set(DNA):
        raise ValueError(f"JASPAR motif must contain rows {DNA}")
    return PWM(motif_name, np.array([rows[base] for base in DNA]).T)


def read_jaspar_collection(path: str | Path) -> list[MonomerPWMRecord]:
    records: list[MonomerPWMRecord] = []
    matrix_id: str | None = None
    symbol: str | None = None
    rows: dict[str, list[float]] = {}

    def finish() -> None:
        if matrix_id is None:
            return
        if set(rows) != set(DNA) or symbol is None:
            raise ValueError(f"incomplete JASPAR matrix {matrix_id}")
        records.append(
            MonomerPWMRecord(
                source_id=matrix_id,
                symbol=symbol.upper(),
                pwm=PWM(f"{matrix_id} {symbol}", np.array([rows[base] for base in DNA]).T),
            )
        )

    with Path(path).open(encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if line.startswith(">"):
                finish()
                fields = line[1:].split(maxsplit=1)
                if len(fields) != 2:
                    raise ValueError(f"invalid JASPAR header: {line!r}")
                matrix_id, symbol = fields
                rows = {}
            elif line and line[0] in DNA:
                values = line.split("[", 1)[-1].split("]", 1)[0].replace(",", " ")
                rows[line[0]] = [float(value) for value in values.split()]
    finish()
    if not records:
        raise ValueError("JASPAR collection is empty")
    return records


def _read_mex_text(text: str, source_id: str) -> MonomerPWMRecord:
    lines = [line.strip() for line in text.splitlines() if line.strip()]
    if len(lines) < 2 or not lines[0].startswith(">"):
        raise ValueError(f"invalid MEX matrix {source_id}")
    name = lines[0][1:].strip()
    matrix = np.array([[float(value) for value in line.split()] for line in lines[1:]], dtype=float)
    if matrix.ndim != 2 or matrix.shape[1] != 4:
        raise ValueError(f"MEX matrix {source_id} must have four A,C,G,T columns")
    symbol = re.split(r"[.@]", name, maxsplit=1)[0].upper()
    return MonomerPWMRecord(source_id=source_id, symbol=symbol, pwm=PWM(name, matrix))


def read_mex_top1_archive(path: str | Path) -> list[MonomerPWMRecord]:
    records: list[MonomerPWMRecord] = []
    with zipfile.ZipFile(path) as archive:
        for info in sorted(archive.infolist(), key=lambda member: member.filename):
            if info.is_dir():
                continue
            if Path(info.filename).name != info.filename:
                raise ValueError(f"MEX archive member must be a flat filename: {info.filename}")
            text = archive.read(info).decode("utf-8")
            records.append(_read_mex_text(text, info.filename))
    if not records:
        raise ValueError("MEX top-1 archive is empty")
    symbols = [record.symbol for record in records]
    if len(set(symbols)) != len(symbols):
        raise ValueError("MEX top-1 archive contains duplicate TF symbols")
    return records


def pwm_consensus(pwm: PWM) -> str:
    return "".join(DNA[index] for index in np.argmax(pwm.probabilities, axis=1))


def pwm_information_bits(pwm: PWM) -> float:
    probabilities = pwm.probabilities
    entropy = -np.sum(probabilities * np.log2(probabilities), axis=1)
    return float(np.sum(2.0 - entropy))


def score_window(sequence: str, pwm: PWM) -> float:
    if len(sequence) != pwm.length or any(base not in DNA_INDEX for base in sequence):
        return float("-inf")
    return float(sum(pwm.log_odds[offset, DNA_INDEX[base]] for offset, base in enumerate(sequence)))


def scan_pwm(sequence: str, pwm: PWM, double_stranded: bool = True) -> np.ndarray:
    sequence = sequence.upper()
    if len(sequence) < pwm.length:
        return np.empty(0, dtype=float)
    offsets = range(len(sequence) - pwm.length + 1)
    scores = np.array([score_window(sequence[offset : offset + pwm.length], pwm) for offset in offsets])
    if not double_stranded:
        return scores
    reverse_scores = np.array(
        [
            score_window(reverse_complement(sequence[offset : offset + pwm.length]), pwm)
            for offset in range(len(sequence) - pwm.length + 1)
        ]
    )
    return np.maximum(scores, reverse_scores)


def maximum_pwm_score(sequence: str, pwm: PWM) -> float:
    scores = scan_pwm(sequence, pwm)
    return float(scores.max()) if scores.size else float("-inf")
