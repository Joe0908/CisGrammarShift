from __future__ import annotations

import gzip
import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from cisgrammar.capselex import AUTOSOMES, require_columns


@dataclass(frozen=True)
class HiChIPSummary:
    loops: int
    promoters: int
    direct_links: int
    distal_links: int


def _attributes(text: str) -> dict[str, str]:
    return dict(re.findall(r'(\w+) "([^"]+)"', text))


def read_gencode_genes(path: str | Path, promoter_half_width_bp: int = 2000) -> pd.DataFrame:
    opener = gzip.open if str(path).endswith(".gz") else open
    records = []
    with opener(path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            if len(fields) != 9 or fields[2] != "gene" or fields[0] not in AUTOSOMES:
                continue
            attrs = _attributes(fields[8])
            start, end = int(fields[3]) - 1, int(fields[4])
            tss = start if fields[6] == "+" else end - 1
            records.append(
                {
                    "chrom": fields[0],
                    "tss": tss,
                    "promoter_start": max(0, tss - promoter_half_width_bp),
                    "promoter_end": tss + promoter_half_width_bp,
                    "gene_id": attrs.get("gene_id", "").split(".")[0],
                    "gene": attrs.get("gene_name", attrs.get("gene_id", "")),
                }
            )
    result = pd.DataFrame(records)
    return result.drop_duplicates(["gene", "chrom", "tss"]).reset_index(drop=True)


def read_fithichip_bedpe(path: str | Path, anchor_width_bp: int = 5000) -> pd.DataFrame:
    frame = pd.read_csv(path, sep="\t", header=None, comment="#")
    if frame.shape[1] < 6:
        raise ValueError("FitHiChIP BEDPE requires six coordinate columns")
    result = frame.iloc[:, :6].copy()
    result.columns = ["chrom1", "start1", "end1", "chrom2", "start2", "end2"]
    for column in ("start1", "end1", "start2", "end2"):
        result[column] = pd.to_numeric(result[column], errors="raise").astype(int)
    for side in (1, 2):
        midpoint = (result[f"start{side}"] + result[f"end{side}"]) // 2
        result[f"start{side}"] = (midpoint - anchor_width_bp // 2).clip(lower=0)
        result[f"end{side}"] = result[f"start{side}"] + anchor_width_bp
    result["loop_score"] = (
        pd.to_numeric(frame.iloc[:, 6], errors="coerce").fillna(1.0) if frame.shape[1] > 6 else 1.0
    )
    result = result[result["chrom1"].isin(AUTOSOMES) & result["chrom2"].isin(AUTOSOMES)]
    return result.reset_index(drop=True)


def _overlaps(intervals: pd.DataFrame, starts: np.ndarray, ends: np.ndarray) -> list[np.ndarray]:
    matches = []
    interval_starts = intervals["start"].to_numpy(dtype=int)
    interval_ends = intervals["end"].to_numpy(dtype=int)
    for start, end in zip(starts, ends, strict=True):
        matches.append(np.flatnonzero((interval_starts < end) & (interval_ends > start)))
    return matches


def link_loci_to_promoters(
    loci: pd.DataFrame,
    promoters: pd.DataFrame,
    loops: pd.DataFrame,
) -> tuple[pd.DataFrame, HiChIPSummary]:
    require_columns(loci, ("locus_id", "chrom", "start", "end"), "locus table")
    links: list[dict[str, object]] = []
    direct = 0
    distal = 0
    for chrom in AUTOSOMES:
        chrom_loci = loci[loci["chrom"] == chrom].reset_index()
        chrom_genes = promoters[promoters["chrom"] == chrom]
        if chrom_loci.empty or chrom_genes.empty:
            continue
        locus_intervals = chrom_loci[["start", "end"]]
        direct_hits = _overlaps(
            locus_intervals,
            chrom_genes["promoter_start"].to_numpy(),
            chrom_genes["promoter_end"].to_numpy(),
        )
        for gene_row, hits in zip(chrom_genes.itertuples(index=False), direct_hits, strict=True):
            for hit in hits:
                links.append(
                    {
                        "gene": gene_row.gene,
                        "locus_id": chrom_loci.loc[hit, "locus_id"],
                        "direct": 1,
                        "loop_score": 0.0,
                    }
                )
                direct += 1
        chrom_loops = loops[(loops["chrom1"] == chrom) | (loops["chrom2"] == chrom)]
        for loop in chrom_loops.itertuples(index=False):
            for promoter_side, distal_side in ((1, 2), (2, 1)):
                promoter_chrom = getattr(loop, f"chrom{promoter_side}")
                distal_chrom = getattr(loop, f"chrom{distal_side}")
                genes = promoters[promoters["chrom"] == promoter_chrom]
                distal_loci = loci[loci["chrom"] == distal_chrom]
                gene_mask = (genes["promoter_start"] < getattr(loop, f"end{promoter_side}")) & (
                    genes["promoter_end"] > getattr(loop, f"start{promoter_side}")
                )
                locus_mask = (distal_loci["start"] < getattr(loop, f"end{distal_side}")) & (
                    distal_loci["end"] > getattr(loop, f"start{distal_side}")
                )
                for gene in genes.loc[gene_mask, "gene"]:
                    for locus_id in distal_loci.loc[locus_mask, "locus_id"]:
                        links.append(
                            {
                                "gene": gene,
                                "locus_id": locus_id,
                                "direct": 0,
                                "loop_score": float(loop.loop_score),
                            }
                        )
                        distal += 1
    link_frame = pd.DataFrame(links, columns=["gene", "locus_id", "direct", "loop_score"])
    if not link_frame.empty:
        link_frame = (
            link_frame.groupby(["gene", "locus_id"], as_index=False)
            .agg(direct=("direct", "max"), loop_count=("direct", "size"), loop_score=("loop_score", "max"))
        )
    return link_frame, HiChIPSummary(len(loops), len(promoters), direct, distal)


def choose_strongest_linked_locus(
    links: pd.DataFrame,
    locus_features: pd.DataFrame,
    chip_columns: tuple[str, str] = ("external_chip_rep1", "external_chip_rep2"),
) -> pd.DataFrame:
    require_columns(locus_features, ("locus_id", *chip_columns), "locus feature table")
    features = locus_features.copy()
    features["mean_log1p_chip"] = np.mean(
        np.log1p(features[list(chip_columns)].to_numpy(dtype=float)), axis=1
    )
    merged = links.merge(features, on="locus_id", how="inner")
    opportunity = merged.groupby("gene")["locus_id"].nunique().rename("linked_locus_count")
    selected = (
        merged.sort_values(["gene", "mean_log1p_chip", "locus_id"], ascending=[True, False, True])
        .drop_duplicates("gene")
        .merge(opportunity, on="gene", how="left")
    )
    return selected.reset_index(drop=True)
