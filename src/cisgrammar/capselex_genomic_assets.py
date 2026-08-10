from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from cisgrammar.capselex import AUTOSOMES, add_chromosome_folds, require_columns
from cisgrammar.capselex_analysis_contract import assert_primary_asset_allowed

NARROWPEAK_COLUMNS = (
    "chrom",
    "start",
    "end",
    "name",
    "score",
    "strand",
    "signal",
    "pvalue",
    "qvalue",
    "summit_offset",
)


def read_narrowpeak(path: str | Path) -> pd.DataFrame:
    assert_primary_asset_allowed(path)
    frame = pd.read_csv(path, sep="\t", header=None, comment="#")
    if frame.shape[1] < 3:
        raise ValueError("narrowPeak requires at least three BED columns")
    frame = frame.iloc[:, : min(10, frame.shape[1])]
    frame.columns = list(NARROWPEAK_COLUMNS[: frame.shape[1]])
    for column in ("start", "end"):
        frame[column] = pd.to_numeric(frame[column], errors="raise").astype(int)
    if "summit_offset" in frame:
        summit = frame["start"] + pd.to_numeric(frame["summit_offset"], errors="coerce").fillna(0).astype(int)
    else:
        summit = ((frame["start"] + frame["end"]) // 2).astype(int)
    frame["midpoint"] = summit
    return frame[frame["chrom"].isin(AUTOSOMES)].reset_index(drop=True)


def fixed_loci(points: pd.DataFrame, width_bp: int = 200, source: str = "unknown") -> pd.DataFrame:
    require_columns(points, ("chrom", "midpoint"), "point table")
    if width_bp <= 0 or width_bp % 2:
        raise ValueError("locus width must be a positive even number")
    half = width_bp // 2
    result = points[["chrom", "midpoint"]].copy()
    result["start"] = (result["midpoint"].astype(int) - half).clip(lower=0)
    result["end"] = result["start"] + width_bp
    result["source"] = source
    return result[["chrom", "start", "end", "midpoint", "source"]]


def build_assay_union_loci(
    chip_peaks: pd.DataFrame,
    ght_peaks: pd.DataFrame,
    focal_tf: str,
    width_bp: int = 200,
) -> pd.DataFrame:
    chip = fixed_loci(chip_peaks, width_bp, "chip")
    ght = fixed_loci(ght_peaks, width_bp, "ght")
    combined = pd.concat([chip, ght], ignore_index=True)
    combined["tile"] = (combined["midpoint"] // width_bp).astype(int)
    membership = combined.pivot_table(
        index=["chrom", "tile"], columns="source", values="midpoint", aggfunc="size", fill_value=0
    )
    selected = (
        combined.sort_values(["chrom", "tile", "source", "midpoint"])
        .drop_duplicates(["chrom", "tile"])
        .drop(columns="source")
    )
    selected = selected.merge(membership.reset_index(), on=["chrom", "tile"], how="left")
    selected["chip_source"] = selected.get("chip", 0).gt(0)
    selected["ght_source"] = selected.get("ght", 0).gt(0)
    selected["focal_tf"] = focal_tf
    selected["locus_id"] = (
        selected["focal_tf"]
        + ":"
        + selected["chrom"]
        + ":"
        + selected["start"].astype(str)
        + "-"
        + selected["end"].astype(str)
    )
    return add_chromosome_folds(selected).sort_values(["chrom", "start"]).reset_index(drop=True)


def add_bigwig_signal(
    loci: pd.DataFrame,
    bigwig_path: str | Path,
    column: str,
    statistic: str = "mean",
) -> pd.DataFrame:
    import pyBigWig

    result = loci.copy()
    with pyBigWig.open(str(bigwig_path)) as bigwig:
        values = []
        for chrom, start, end in result[["chrom", "start", "end"]].itertuples(index=False, name=None):
            value = bigwig.stats(str(chrom), int(start), int(end), type=statistic, exact=True)[0]
            values.append(0.0 if value is None or not np.isfinite(value) else float(value))
    result[column] = values
    return result


def read_fasta(path: str | Path) -> dict[str, str]:
    sequences: dict[str, list[str]] = {}
    name: str | None = None
    with Path(path).open(encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if line.startswith(">"):
                name = line[1:].split()[0]
                sequences[name] = []
            elif name is not None:
                sequences[name].append(line.upper())
    return {name: "".join(parts) for name, parts in sequences.items()}


def extract_contexts(loci: pd.DataFrame, genome: dict[str, str], width_bp: int = 400) -> list[str]:
    half = width_bp // 2
    sequences = []
    for chrom, midpoint in loci[["chrom", "midpoint"]].itertuples(index=False, name=None):
        start = int(midpoint) - half
        end = start + width_bp
        chromosome = genome.get(str(chrom), "")
        sequence = chromosome[max(0, start) : max(0, end)]
        if start < 0:
            sequence = "N" * -start + sequence
        sequences.append(sequence.ljust(width_bp, "N")[:width_bp])
    return sequences
