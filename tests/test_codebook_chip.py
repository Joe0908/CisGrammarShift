from __future__ import annotations

import gzip

import pytest

from cisgrammar.capselex_codebook_chip import build_gpzn_manifest


def _write_soft(path, samples: list[tuple[str, str, str]]) -> None:
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        for gsm, tf, replicate in samples:
            handle.write(f"^SAMPLE = {gsm}\n")
            handle.write(f"!Sample_title = ChIP-seq analysis of {tf} genomic target locations\n")
            handle.write(f"!Sample_description = THC000{replicate}TANN_{tf}_ChIP{replicate}\n")
            handle.write(
                "!Sample_supplementary_file_1 = "
                f"ftp://ftp.ncbi.nlm.nih.gov/{gsm}_THC000{replicate}GPZN_"
                f"{tf}_ChIP{replicate}_hg38_PE.bw\n"
            )


def _write_filelist(path, samples: list[tuple[str, str, str]]) -> None:
    lines = ["#Archive/File\tName\tTime\tSize\tType"]
    for gsm, tf, replicate in samples:
        filename = f"{gsm}_THC000{replicate}GPZN_{tf}_ChIP{replicate}_hg38_PE.bw"
        lines.append(f"File\t{filename}\t01/01/2026 00:00:00\t{100 + int(replicate)}\tBW")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def test_build_gpzn_manifest_requires_and_resolves_two_replicates(tmp_path) -> None:
    samples = [
        ("GSM1", "GCM1", "1"),
        ("GSM2", "GCM1", "2"),
        ("GSM3", "MAX", "1"),
        ("GSM4", "MAX", "2"),
    ]
    soft = tmp_path / "family.soft.gz"
    filelist = tmp_path / "filelist.txt"
    _write_soft(soft, samples)
    _write_filelist(filelist, samples)
    manifest = build_gpzn_manifest(soft, filelist, ["GCM1", "MAX"])
    assert len(manifest["assets"]) == 4
    assert manifest["total_size_bytes"] == 406
    assert manifest["assets"][0]["url"].startswith("https://")
    assert [record["replicate"] for record in manifest["assets"][:2]] == [1, 2]


def test_build_gpzn_manifest_rejects_missing_replicate(tmp_path) -> None:
    samples = [("GSM1", "GCM1", "1")]
    soft = tmp_path / "family.soft.gz"
    filelist = tmp_path / "filelist.txt"
    _write_soft(soft, samples)
    _write_filelist(filelist, samples)
    with pytest.raises(ValueError, match="replicates 1 and 2"):
        build_gpzn_manifest(soft, filelist, ["GCM1"])
