from __future__ import annotations

import gzip
import shutil
import subprocess
import tempfile
from collections.abc import Iterator
from pathlib import Path
from typing import Literal

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

MAGIX_COLUMNS = (
    "chr",
    "start",
    "stop",
    "name",
    "coefficient.br",
    "coefficient.ar",
    "full_LL",
    "reduced_LL",
    "pvalue",
    "fdr",
)

LocusUniverse = Literal["ght-only", "assay-union", "fixed-genome", "legacy-assay-union"]
LOCUS_UNIVERSES: tuple[LocusUniverse, ...] = (
    "ght-only",
    "assay-union",
    "fixed-genome",
    "legacy-assay-union",
)
LOCUS_UNIVERSE_AUDIT: dict[LocusUniverse, dict[str, str | bool]] = {
    "ght-only": {
        "selection_rule": "ght_tile_only",
        "center_rule": "fixed_tile_center",
        "outcome_used_for_selection": False,
        "outcome_used_for_centering": False,
    },
    "assay-union": {
        "selection_rule": "chip_or_ght_tile",
        "center_rule": "fixed_tile_center",
        "outcome_used_for_selection": True,
        "outcome_used_for_centering": False,
    },
    "fixed-genome": {
        "selection_rule": "all_full_autosomal_tiles",
        "center_rule": "fixed_tile_center",
        "outcome_used_for_selection": False,
        "outcome_used_for_centering": False,
    },
    "legacy-assay-union": {
        "selection_rule": "chip_or_ght_tile_chip_priority",
        "center_rule": "legacy_selected_peak_midpoint",
        "outcome_used_for_selection": True,
        "outcome_used_for_centering": True,
    },
}


def _validate_width(width_bp: int) -> None:
    if width_bp <= 0 or width_bp % 2:
        raise ValueError("locus width must be a positive even number")


def _tile_counts(points: pd.DataFrame | None, width_bp: int, source: str) -> pd.DataFrame:
    if points is None or points.empty:
        return pd.DataFrame(columns=["chrom", "tile", f"{source}_count"])
    require_columns(points, ("chrom", "midpoint"), f"{source} point table")
    frame = points.loc[points["chrom"].isin(AUTOSOMES), ["chrom", "midpoint"]].copy()
    frame["midpoint"] = pd.to_numeric(frame["midpoint"], errors="raise").astype(int)
    frame = frame[frame["midpoint"] >= 0]
    frame["tile"] = (frame["midpoint"] // width_bp).astype(int)
    return (
        frame.groupby(["chrom", "tile"], as_index=False)
        .size()
        .rename(columns={"size": f"{source}_count"})
    )


def _membership_table(
    chip_peaks: pd.DataFrame | None,
    ght_peaks: pd.DataFrame | None,
    width_bp: int,
) -> pd.DataFrame:
    chip = _tile_counts(chip_peaks, width_bp, "chip")
    ght = _tile_counts(ght_peaks, width_bp, "ght")
    membership = chip.merge(ght, on=["chrom", "tile"], how="outer")
    for column in ("chip_count", "ght_count"):
        if column not in membership:
            membership[column] = 0
        membership[column] = membership[column].fillna(0).astype(int)
    return membership


def _finalize_fixed_tiles(
    tiles: pd.DataFrame,
    focal_tf: str,
    width_bp: int,
    locus_universe: LocusUniverse,
) -> pd.DataFrame:
    require_columns(tiles, ("chrom", "tile", "chip_count", "ght_count"), "tile table")
    result = tiles.copy()
    result["tile"] = pd.to_numeric(result["tile"], errors="raise").astype(int)
    result["start"] = result["tile"] * width_bp
    result["end"] = result["start"] + width_bp
    result["midpoint"] = result["start"] + width_bp // 2
    result["chip_source"] = result["chip_count"].gt(0)
    result["ght_source"] = result["ght_count"].gt(0)
    result["focal_tf"] = focal_tf
    result["locus_universe"] = locus_universe
    for key, value in LOCUS_UNIVERSE_AUDIT[locus_universe].items():
        result[key] = value
    result["locus_id"] = (
        result["focal_tf"]
        + ":"
        + result["chrom"]
        + ":"
        + result["start"].astype(str)
        + "-"
        + result["end"].astype(str)
    )
    columns = [
        "locus_id",
        "focal_tf",
        "chrom",
        "start",
        "end",
        "midpoint",
        "tile",
        "chip_source",
        "ght_source",
        "chip_count",
        "ght_count",
        "locus_universe",
        "selection_rule",
        "center_rule",
        "outcome_used_for_selection",
        "outcome_used_for_centering",
    ]
    result = add_chromosome_folds(result[columns])
    result["_chrom_order"] = result["chrom"].map({chrom: i for i, chrom in enumerate(AUTOSOMES)})
    return result.sort_values(["_chrom_order", "start", "locus_id"]).drop(columns="_chrom_order").reset_index(
        drop=True
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
    if frame["start"].lt(0).any() or frame["end"].le(frame["start"]).any():
        raise ValueError("narrowPeak intervals require 0 <= start < end")
    if "summit_offset" in frame:
        offset = pd.to_numeric(frame["summit_offset"], errors="coerce")
        invalid_offset = offset.ge(0) & offset.ge(frame["end"] - frame["start"])
        if invalid_offset.any():
            raise ValueError("narrowPeak summit offsets must fall inside their intervals")
        summit = frame["start"] + offset.fillna(-1).astype(int)
        summit = summit.where(offset.ge(0), (frame["start"] + frame["end"]) // 2)
    else:
        summit = ((frame["start"] + frame["end"]) // 2).astype(int)
    frame["midpoint"] = summit
    return frame[frame["chrom"].isin(AUTOSOMES)].reset_index(drop=True)


def read_magix(
    path: str | Path,
    fdr_max: float = 0.05,
    require_positive: bool = True,
) -> pd.DataFrame:
    """Read a MAGIX result table and retain outcome-independent GHT calls.

    ``coefficient.ar`` is the refined coefficient emitted after MAGIX's
    likelihood-ratio testing step. The default call requires a Benjamini-Hochberg
    FDR no greater than 0.05 and a positive refined coefficient. Both criteria use
    only GHT-SELEX data; ChIP is never used to select or center these intervals.
    """
    assert_primary_asset_allowed(path)
    if not 0 <= fdr_max <= 1:
        raise ValueError("MAGIX FDR threshold must fall between zero and one")
    frame = pd.read_csv(path, sep="\t")
    missing = [column for column in MAGIX_COLUMNS if column not in frame.columns]
    if missing:
        raise ValueError(f"MAGIX table is missing columns: {', '.join(missing)}")
    frame = frame.loc[:, MAGIX_COLUMNS].copy()
    frame = frame.rename(
        columns={
            "chr": "chrom",
            "stop": "end",
            "coefficient.br": "ght_score_before_refinement",
            "coefficient.ar": "ght_score",
        }
    )
    numeric = (
        "start",
        "end",
        "ght_score_before_refinement",
        "ght_score",
        "full_LL",
        "reduced_LL",
        "pvalue",
        "fdr",
    )
    for column in numeric:
        frame[column] = pd.to_numeric(frame[column], errors="raise")
    frame[["start", "end"]] = frame[["start", "end"]].astype(int)
    if frame["start"].lt(0).any() or frame["end"].le(frame["start"]).any():
        raise ValueError("MAGIX intervals require 0 <= start < end")
    if frame[["pvalue", "fdr"]].lt(0).any().any() or frame[["pvalue", "fdr"]].gt(1).any().any():
        raise ValueError("MAGIX p-values and FDR values must fall between zero and one")
    frame["midpoint"] = ((frame["start"] + frame["end"]) // 2).astype(int)
    selected = frame["fdr"].le(fdr_max)
    if require_positive:
        selected &= frame["ght_score"].gt(0)
    frame = frame.loc[selected & frame["chrom"].isin(AUTOSOMES)].copy()
    frame["ght_fdr_threshold"] = float(fdr_max)
    frame["ght_positive_required"] = bool(require_positive)
    return frame.reset_index(drop=True)


def fixed_loci(points: pd.DataFrame, width_bp: int = 200, source: str = "unknown") -> pd.DataFrame:
    require_columns(points, ("chrom", "midpoint"), "point table")
    _validate_width(width_bp)
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
    """Build the ChIP/GHT tile union without allowing either assay to set the centre.

    ChIP still contributes to inclusion, so this is an outcome-informed *selection*
    sensitivity. Every sequence window is nevertheless centred on a fixed genomic tile,
    removing the more serious summit-centering shortcut in the historical implementation.
    """
    _validate_width(width_bp)
    membership = _membership_table(chip_peaks, ght_peaks, width_bp)
    return _finalize_fixed_tiles(
        membership,
        focal_tf,
        width_bp,
        locus_universe="assay-union",
    )


def build_ght_only_loci(
    ght_peaks: pd.DataFrame,
    focal_tf: str,
    width_bp: int = 200,
    chip_peaks: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """Build the primary outcome-independent universe from GHT peak tiles only."""
    _validate_width(width_bp)
    membership = _membership_table(chip_peaks, ght_peaks, width_bp)
    selected = membership[membership["ght_count"].gt(0)].copy()
    return _finalize_fixed_tiles(
        selected,
        focal_tf,
        width_bp,
        locus_universe="ght-only",
    )


def build_legacy_assay_union_loci(
    chip_peaks: pd.DataFrame,
    ght_peaks: pd.DataFrame,
    focal_tf: str,
    width_bp: int = 200,
) -> pd.DataFrame:
    """Reproduce the historical ChIP-prioritized summit-centred universe explicitly.

    This function exists only for a documented legacy sensitivity. It must not be used as
    the primary analysis because ChIP affects both locus inclusion and sequence centering.
    """
    _validate_width(width_bp)
    chip = fixed_loci(chip_peaks, width_bp, "chip")
    ght = fixed_loci(ght_peaks, width_bp, "ght")
    combined = pd.concat([chip, ght], ignore_index=True)
    combined["tile"] = (combined["midpoint"] // width_bp).astype(int)
    membership = _membership_table(chip_peaks, ght_peaks, width_bp)
    selected = (
        combined.sort_values(["chrom", "tile", "source", "midpoint"])
        .drop_duplicates(["chrom", "tile"])
        .drop(columns="source")
        .merge(membership, on=["chrom", "tile"], how="left", validate="one_to_one")
    )
    selected["chip_source"] = selected["chip_count"].gt(0)
    selected["ght_source"] = selected["ght_count"].gt(0)
    selected["focal_tf"] = focal_tf
    selected["locus_universe"] = "legacy-assay-union"
    for key, value in LOCUS_UNIVERSE_AUDIT["legacy-assay-union"].items():
        selected[key] = value
    selected["locus_id"] = (
        selected["focal_tf"]
        + ":"
        + selected["chrom"]
        + ":"
        + selected["start"].astype(str)
        + "-"
        + selected["end"].astype(str)
    )
    columns = [
        "locus_id",
        "focal_tf",
        "chrom",
        "start",
        "end",
        "midpoint",
        "tile",
        "chip_source",
        "ght_source",
        "chip_count",
        "ght_count",
        "locus_universe",
        "selection_rule",
        "center_rule",
        "outcome_used_for_selection",
        "outcome_used_for_centering",
    ]
    selected = add_chromosome_folds(selected[columns])
    selected["_chrom_order"] = selected["chrom"].map(
        {chrom: i for i, chrom in enumerate(AUTOSOMES)}
    )
    return (
        selected.sort_values(["_chrom_order", "start", "locus_id"])
        .drop(columns="_chrom_order")
        .reset_index(drop=True)
    )


def read_chrom_sizes(path: str | Path) -> pd.DataFrame:
    """Read a UCSC-style two-column chromosome-size file and retain autosomes."""
    frame = pd.read_csv(
        path,
        sep=r"\s+",
        header=None,
        comment="#",
        usecols=[0, 1],
        names=["chrom", "size"],
    )
    frame["size"] = pd.to_numeric(frame["size"], errors="raise").astype(int)
    frame = frame[frame["chrom"].isin(AUTOSOMES)].copy()
    if frame["chrom"].duplicated().any():
        duplicated = ", ".join(sorted(frame.loc[frame["chrom"].duplicated(), "chrom"].unique()))
        raise ValueError(f"chromosome sizes contain duplicates: {duplicated}")
    if frame.empty or frame["size"].le(0).any():
        raise ValueError("positive autosomal chromosome sizes are required")
    order = {chrom: index for index, chrom in enumerate(AUTOSOMES)}
    frame["_order"] = frame["chrom"].map(order)
    return frame.sort_values("_order").drop(columns="_order").reset_index(drop=True)


def iter_fixed_genome_loci(
    chrom_sizes: pd.DataFrame,
    focal_tf: str,
    width_bp: int = 200,
    chip_peaks: pd.DataFrame | None = None,
    ght_peaks: pd.DataFrame | None = None,
) -> Iterator[pd.DataFrame]:
    """Yield full-width autosomal tiles in chromosome-sized chunks.

    Streaming avoids materializing roughly fifteen million hg38 200-bp tiles at once.
    The final chromosome remainder shorter than ``width_bp`` is excluded so every locus
    has exactly the frozen width.
    """
    _validate_width(width_bp)
    require_columns(chrom_sizes, ("chrom", "size"), "chromosome sizes")
    membership = _membership_table(chip_peaks, ght_peaks, width_bp)
    for chrom, size in chrom_sizes[["chrom", "size"]].itertuples(index=False, name=None):
        if chrom not in AUTOSOMES:
            continue
        n_tiles = int(size) // width_bp
        if n_tiles < 1:
            continue
        tiles = pd.DataFrame({"chrom": chrom, "tile": np.arange(n_tiles, dtype=np.int64)})
        chrom_membership = membership[membership["chrom"].eq(chrom)].drop(columns="chrom")
        tiles = tiles.merge(chrom_membership, on="tile", how="left", validate="one_to_one")
        for column in ("chip_count", "ght_count"):
            tiles[column] = tiles[column].fillna(0).astype(int)
        yield _finalize_fixed_tiles(
            tiles,
            focal_tf,
            width_bp,
            locus_universe="fixed-genome",
        )


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


def add_cross_assembly_bigwig_signal(
    loci: pd.DataFrame,
    source_bigwig: str | Path,
    chain_gz: str | Path,
    lift_over_binary: str | Path,
    column: str,
    source_assembly: str,
) -> pd.DataFrame:
    """Project locus intervals through a UCSC chain and query the source-assembly BigWig."""
    import pyBigWig

    require_columns(loci, ("chrom", "start", "end"), "locus table")
    result = loci.copy()
    binary = Path(lift_over_binary).resolve()
    if not binary.stat().st_mode & 0o111:
        binary.chmod(binary.stat().st_mode | 0o111)
    mapped_coordinates: dict[int, tuple[str, int, int]] = {}
    with tempfile.TemporaryDirectory(prefix="cisgrammar_liftover_") as temporary:
        directory = Path(temporary)
        input_bed = directory / "input.bed"
        mapped_bed = directory / "mapped.bed"
        unmapped_bed = directory / "unmapped.bed"
        chain = directory / "mapping.over.chain"
        with input_bed.open("w", encoding="utf-8") as handle:
            for index, (chrom, start, end) in enumerate(
                result[["chrom", "start", "end"]].itertuples(index=False, name=None)
            ):
                handle.write(f"{chrom}\t{int(start)}\t{int(end)}\t{index}\n")
        with gzip.open(chain_gz, "rb") as source, chain.open("wb") as target:
            shutil.copyfileobj(source, target)
        subprocess.run(
            [
                str(binary),
                "-bedPlus=4",
                "-tab",
                str(input_bed),
                str(chain),
                str(mapped_bed),
                str(unmapped_bed),
            ],
            check=True,
            capture_output=True,
            text=True,
        )
        with mapped_bed.open(encoding="utf-8") as handle:
            for raw in handle:
                chrom, start, end, raw_index = raw.rstrip("\n").split("\t")[:4]
                index = int(raw_index)
                if index in mapped_coordinates:
                    raise RuntimeError(f"locus mapped more than once: {index}")
                mapped_coordinates[index] = (chrom, int(start), int(end))

    values = np.full(len(result), np.nan, dtype=float)
    mapped_chrom = np.full(len(result), "", dtype=object)
    mapped_start = np.full(len(result), -1, dtype=np.int64)
    mapped_end = np.full(len(result), -1, dtype=np.int64)
    with pyBigWig.open(str(source_bigwig)) as bigwig:
        chromosomes = bigwig.chroms()
        for index, (chrom, start, end) in mapped_coordinates.items():
            if chrom not in chromosomes or start < 0 or end <= start or end > chromosomes[chrom]:
                continue
            value = bigwig.stats(chrom, start, end, type="mean", exact=True)[0]
            values[index] = 0.0 if value is None or not np.isfinite(value) else float(value)
            mapped_chrom[index] = chrom
            mapped_start[index] = start
            mapped_end[index] = end
    prefix = f"{column}_{source_assembly}"
    result[column] = values
    result[f"{column}_mapped"] = np.isfinite(values)
    result[f"{prefix}_chrom"] = mapped_chrom
    result[f"{prefix}_start"] = mapped_start
    result[f"{prefix}_end"] = mapped_end
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
