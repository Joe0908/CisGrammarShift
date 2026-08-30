from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from cisgrammar.capselex_trophoblast_deseq2 import (
    export_deseq2_inputs,
    read_deseq2_results,
    validate_deseq2_session_info,
)
from cisgrammar.capselex_trophoblast_hichip import link_loci_to_promoters


def test_hichip_linking_keeps_direct_and_distal_channels() -> None:
    loci = pd.DataFrame(
        {
            "locus_id": ["direct", "distal"],
            "chrom": ["chr1", "chr1"],
            "start": [90, 5000],
            "end": [110, 5100],
        }
    )
    promoters = pd.DataFrame(
        {"gene": ["G"], "chrom": ["chr1"], "promoter_start": [80], "promoter_end": [120]}
    )
    loops = pd.DataFrame(
        {
            "chrom1": ["chr1"],
            "start1": [50],
            "end1": [150],
            "chrom2": ["chr1"],
            "start2": [4900],
            "end2": [5200],
            "loop_score": [4.0],
        }
    )
    links, summary = link_loci_to_promoters(loci, promoters, loops)
    assert set(links["locus_id"]) == {"direct", "distal"}
    assert summary.direct_links == 1
    assert summary.distal_links == 1


def test_validate_exact_deseq2_runtime(tmp_path: Path) -> None:
    session = tmp_path / "sessionInfo.txt"
    session.write_text("R version 4.2.2 Patched\nother attached packages:\nDESeq2_1.36.0\n", encoding="utf-8")
    assert validate_deseq2_session_info(session) == {"r": "4.2.2", "deseq2": "1.36.0"}
    session.write_text("R version 4.3.0\nDESeq2_1.36.0\n", encoding="utf-8")
    with pytest.raises(ValueError, match="R 4.2.2"):
        validate_deseq2_session_info(session)


def test_read_deseq2_normalizes_direction(tmp_path: Path) -> None:
    path = tmp_path / "EVT.deseq2.tsv"
    pd.DataFrame(
        {
            "gene": ["A"],
            "baseMean": [10],
            "log2FoldChange": [-2.0],
            "lfcSE": [0.2],
            "pvalue": [0.01],
            "padj": [0.02],
        }
    ).to_csv(path, sep="\t", index=False)
    assert read_deseq2_results(path).loc[0, "gcm1_dependence"] == 2.0


def test_export_deseq2_inputs(tmp_path: Path) -> None:
    columns = {"gene": ["A", "B"]}
    for state in ("EVT", "ST"):
        for replicate in range(1, 5):
            columns[f"{state}_WT_{replicate}"] = [replicate, replicate + 1]
        for replicate in range(1, 3):
            columns[f"{state}_GCM1_KO_{replicate}"] = [replicate + 1, replicate + 2]
    source = tmp_path / "dataset.xlsx"
    pd.DataFrame(columns).to_excel(source, index=False)
    manifest = export_deseq2_inputs(source, tmp_path / "out")
    assert manifest["states"]["EVT"] == {"genes": 2, "wt": 4, "ko": 2}
    counts = pd.read_csv(tmp_path / "out" / "EVT.counts.tsv", sep="\t")
    assert counts.shape == (2, 7)
