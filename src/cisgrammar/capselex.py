from __future__ import annotations

import hashlib
import json
import re
from collections.abc import Iterable
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

TF_ALIASES = {
    "TP73L": "TP63",
    "RFXDC2": "RFX7",
    "RXF5": "RFX5",
    "POU5F1B": "POU5F1",
}

AUTOSOMES = tuple(f"chr{i}" for i in range(1, 23))


@dataclass(frozen=True)
class Asset:
    path: str
    bytes: int
    sha256: str


def normalize_tf_symbol(value: object, aliases: dict[str, str] | None = None) -> str:
    """Return an uppercase HGNC-like symbol while preserving explicit aliases."""
    text = str(value).strip().upper().replace(" ", "")
    text = re.sub(r"^(HUMAN|HS)[_:.-]", "", text)
    text = re.sub(r"[-_](FULL|EDBD|DBD|ISOFORM\d+)$", "", text)
    mapping = TF_ALIASES if aliases is None else {**TF_ALIASES, **aliases}
    return mapping.get(text, text)


def sha256_file(path: str | Path, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        while chunk := handle.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def audit_asset(path: str | Path) -> Asset:
    resolved = Path(path)
    return Asset(str(resolved), resolved.stat().st_size, sha256_file(resolved))


def write_json(payload: Any, path: str | Path) -> None:
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def asset_manifest(paths: Iterable[str | Path]) -> list[dict[str, object]]:
    return [asdict(audit_asset(path)) for path in paths]


def chromosome_fold(chromosome: object, n_folds: int = 5) -> int:
    match = re.fullmatch(r"(?:chr)?(\d+)", str(chromosome))
    if match is None or not 1 <= int(match.group(1)) <= 22:
        raise ValueError(f"autosomal chromosome required, got {chromosome!r}")
    return (int(match.group(1)) - 1) % n_folds


def add_chromosome_folds(
    frame: pd.DataFrame,
    chromosome_column: str = "chrom",
    n_folds: int = 5,
) -> pd.DataFrame:
    result = frame.copy()
    result["chromosome_fold"] = [chromosome_fold(value, n_folds) for value in result[chromosome_column]]
    return result


def reverse_complement(sequence: str) -> str:
    return sequence.upper().translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def stable_seed(*parts: object) -> int:
    digest = hashlib.blake2b("|".join(map(str, parts)).encode(), digest_size=8).digest()
    return int.from_bytes(digest, "little") % (2**32 - 1)


def require_columns(frame: pd.DataFrame, columns: Iterable[str], label: str = "table") -> None:
    missing = sorted(set(columns) - set(frame.columns))
    if missing:
        raise ValueError(f"{label} missing required columns: {', '.join(missing)}")


def finite_numeric(series: pd.Series, fill: float = 0.0) -> np.ndarray:
    values = pd.to_numeric(series, errors="coerce").to_numpy(dtype=float)
    values[~np.isfinite(values)] = fill
    return values
