from __future__ import annotations

import csv

import pytest
from openpyxl import Workbook

from cisgrammar.capselex_ght_rebuild import build_ght_rebuild_audit, read_selected_ght_runs


def _write_metadata(path) -> None:
    workbook = Workbook()
    sheet = workbook.active
    sheet.title = "Metadata"
    sheet.append(["preamble"])
    sheet.append(
        [
            "*library name",
            "*title",
            "*SRA Experiment or Run",
            "*SRA Experiment or Run",
            "*processed data file ",
            "processed data file ",
            "processed data file ",
        ]
    )
    sheet.append(
        [
            None,
            "Genomic HT-SELEX analysis of GCM1 target sites, selection Cycle1 Exp. GHT00206_1",
            "ERR1",
            None,
            "all.bw",
            "cycle.bw",
            "GCM1_MAGIX_GHT00206.bed",
        ]
    )
    sheet.append(
        [
            None,
            "Genomic HT-SELEX analysis of GCM1 target sites, selection Cycle2 Exp. GHT00206_2",
            "ERR2",
            None,
            "all.bw",
            "cycle.bw",
            "GCM1_MAGIX_GHT00206.bed",
        ]
    )
    sheet.append(
        [
            None,
            "Genomic HT-SELEX analysis of GCM1 target sites, selection Cycle3 Exp. GHT00999_3",
            "ERR_UNUSED",
            None,
            "all.bw",
            "cycle.bw",
            None,
        ]
    )
    sheet.append(
        [
            None,
            "Genomic HT-SELEX analysis of MAX target sites, selection Cycle1 Exp. GHT00107_1",
            "ERR3",
            None,
            "all.bw",
            "cycle.bw",
            "MAX_MAGIX_GHT00107.bed",
        ]
    )
    workbook.save(path)


def _write_ena(path) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["run_accession", "fastq_ftp", "fastq_md5", "fastq_bytes"])
        writer.writerow(["ERR1", "host/ERR1_1.fastq.gz;host/ERR1_2.fastq.gz", "a;b", "10;11"])
        writer.writerow(["ERR2", "host/ERR2.fastq.gz", "c", "12"])
        writer.writerow(["ERR3", "host/ERR3_1.fastq.gz;host/ERR3_2.fastq.gz", "d;e", "13;14"])


def test_selected_runs_follow_final_magix_column(tmp_path) -> None:
    metadata = tmp_path / "metadata.xlsx"
    _write_metadata(metadata)
    records = read_selected_ght_runs(metadata, ["GCM1"])
    assert [record.run_accession for record in records] == ["ERR1", "ERR2"]
    assert [record.cycle for record in records] == [1, 2]


def test_rebuild_audit_freezes_sizes_and_marks_missing_production_design(tmp_path) -> None:
    metadata = tmp_path / "metadata.xlsx"
    ena = tmp_path / "ena.tsv"
    _write_metadata(metadata)
    _write_ena(ena)
    report = build_ght_rebuild_audit(
        metadata,
        ena,
        primary_tfs=["GCM1"],
        sensitivity_tfs=["MAX"],
        audit_date="2026-08-10",
        magix_commit="abc123",
    )
    assert report["primary_target_fastq_bytes"] == 33
    assert report["all_target_fastq_bytes"] == 60
    assert report["target_fastq_download_ready"] is True
    assert report["exact_author_v2_rebuild_ready"] is False
    assert report["panel"][0]["runs"][0]["read_layout"] == "paired"
    assert report["panel"][0]["runs"][1]["read_layout"] == "single"


def test_selected_runs_require_every_requested_tf(tmp_path) -> None:
    metadata = tmp_path / "metadata.xlsx"
    _write_metadata(metadata)
    with pytest.raises(ValueError, match="no author-selected GHT runs"):
        read_selected_ght_runs(metadata, ["RFX5"])
