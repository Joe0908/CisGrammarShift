from __future__ import annotations

import gzip
import json

import numpy as np
import pandas as pd
import pytest

from cisgrammar.capselex_chip_qc import replicate_qc
from cisgrammar.capselex_genomic_assets import (
    add_cross_assembly_bigwig_signal,
    build_assay_union_loci,
    build_ght_only_loci,
    build_legacy_assay_union_loci,
    iter_fixed_genome_loci,
    read_chrom_sizes,
    read_magix,
)
from cisgrammar.capselex_genomic_model import (
    chromosome_bootstrap_partial_r2,
    chromosome_nested_continuous,
    chromosome_shift_null,
)
from cisgrammar.cli import _run_capselex, build_parser


def test_assay_union_deduplicates_fixed_tiles() -> None:
    chip = pd.DataFrame({"chrom": ["chr1", "chr1"], "midpoint": [190, 210]})
    ght = pd.DataFrame({"chrom": ["chr1", "chr2"], "midpoint": [110, 300]})
    loci = build_assay_union_loci(chip, ght, "GCM1", width_bp=200)
    assert len(loci) == 3
    first = loci[loci["chrom"] == "chr1"].iloc[0]
    assert first["chip_source"] and first["ght_source"]
    assert (first["start"], first["midpoint"], first["end"]) == (0, 100, 200)
    assert first["center_rule"] == "fixed_tile_center"
    assert bool(first["outcome_used_for_selection"])
    assert not bool(first["outcome_used_for_centering"])


def test_narrowpeak_without_summit_uses_interval_midpoint(tmp_path) -> None:
    path = tmp_path / "no_summit.narrowPeak"
    path.write_text("chr1\t100\t180\tpeak\t10\t.\t2\t3\t4\t-1\n", encoding="utf-8")
    from cisgrammar.capselex_genomic_assets import read_narrowpeak

    assert read_narrowpeak(path).loc[0, "midpoint"] == 140


def test_magix_uses_refined_positive_fdr_calls(tmp_path) -> None:
    path = tmp_path / "GCM1.magix.bed"
    path.write_text(
        "chr\tstart\tstop\tname\tcoefficient.br\tcoefficient.ar\tfull_LL\t"
        "reduced_LL\tpvalue\tfdr\n"
        "chr1\t100\t300\tpositive\t0.2\t1.4\t-2\t-4\t0.001\t0.01\n"
        "chr1\t300\t500\tnegative\t0.2\t-1.4\t-2\t-4\t0.001\t0.01\n"
        "chr1\t500\t700\tnonsignificant\t0.2\t1.4\t-2\t-4\t0.2\t0.2\n"
        "chrX\t100\t300\tnonautosomal\t0.2\t1.4\t-2\t-4\t0.001\t0.01\n",
        encoding="utf-8",
    )
    frame = read_magix(path)
    assert frame["name"].tolist() == ["positive"]
    assert frame["midpoint"].tolist() == [200]
    assert frame["ght_score"].tolist() == [1.4]
    assert frame["ght_score_before_refinement"].tolist() == [0.2]


def test_magix_rejects_invalid_interval(tmp_path) -> None:
    path = tmp_path / "invalid.magix.bed"
    path.write_text(
        "chr\tstart\tstop\tname\tcoefficient.br\tcoefficient.ar\tfull_LL\t"
        "reduced_LL\tpvalue\tfdr\n"
        "chr1\t300\t100\tinvalid\t0.2\t1.4\t-2\t-4\t0.001\t0.01\n",
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="0 <= start < end"):
        read_magix(path)


def test_ght_only_excludes_chip_only_tiles_and_is_outcome_independent() -> None:
    chip = pd.DataFrame({"chrom": ["chr1", "chr1"], "midpoint": [190, 410]})
    ght = pd.DataFrame({"chrom": ["chr1"], "midpoint": [110]})
    loci = build_ght_only_loci(ght, "GCM1", width_bp=200, chip_peaks=chip)
    assert loci["tile"].tolist() == [0]
    assert loci["chip_count"].tolist() == [1]
    assert not loci["outcome_used_for_selection"].any()
    assert not loci["outcome_used_for_centering"].any()


def test_legacy_assay_union_preserves_chip_priority_for_auditing() -> None:
    chip = pd.DataFrame({"chrom": ["chr1"], "midpoint": [190]})
    ght = pd.DataFrame({"chrom": ["chr1"], "midpoint": [110]})
    loci = build_legacy_assay_union_loci(chip, ght, "GCM1", width_bp=200)
    assert len(loci) == 1
    assert (loci.loc[0, "start"], loci.loc[0, "midpoint"], loci.loc[0, "end"]) == (90, 190, 290)
    assert loci.loc[0, "center_rule"] == "legacy_selected_peak_midpoint"
    assert bool(loci.loc[0, "outcome_used_for_selection"])
    assert bool(loci.loc[0, "outcome_used_for_centering"])


def test_fixed_genome_streams_full_width_tiles(tmp_path) -> None:
    sizes_path = tmp_path / "tiny.chrom.sizes"
    sizes_path.write_text("chr1\t450\nchr2\t210\nchrX\t1000\n", encoding="utf-8")
    chrom_sizes = read_chrom_sizes(sizes_path)
    chip = pd.DataFrame({"chrom": ["chr1"], "midpoint": [210]})
    ght = pd.DataFrame({"chrom": ["chr2"], "midpoint": [100]})
    loci = pd.concat(
        iter_fixed_genome_loci(chrom_sizes, "GCM1", 200, chip_peaks=chip, ght_peaks=ght),
        ignore_index=True,
    )
    assert len(loci) == 3
    assert loci[["chrom", "start", "end"]].values.tolist() == [
        ["chr1", 0, 200],
        ["chr1", 200, 400],
        ["chr2", 0, 200],
    ]
    assert loci["chip_source"].sum() == 1
    assert loci["ght_source"].sum() == 1
    assert not loci["outcome_used_for_selection"].any()
    assert not loci["outcome_used_for_centering"].any()


def test_build_loci_cli_writes_auditable_manifest(tmp_path) -> None:
    chip_path = tmp_path / "GCM1.chip.narrowPeak"
    ght_path = tmp_path / "GCM1.ght.narrowPeak"
    output = tmp_path / "GCM1.ght_only.tsv.gz"
    chip_path.write_text("chr1\t150\t230\tchip\t10\t.\t2\t3\t4\t40\n", encoding="utf-8")
    ght_path.write_text("chr1\t80\t160\tght\t10\t.\t2\t3\t4\t30\n", encoding="utf-8")
    args = build_parser().parse_args(
        [
            "capselex",
            "build-loci",
            "--universe",
            "ght-only",
            "--chip-peaks",
            str(chip_path),
            "--ght-peaks",
            str(ght_path),
            "--focal-tf",
            "GCM1",
            "--output",
            str(output),
        ]
    )
    _run_capselex(args)
    loci = pd.read_csv(output, sep="\t")
    manifest = json.loads((tmp_path / "GCM1.ght_only.tsv.gz.manifest.json").read_text())
    assert loci["locus_universe"].tolist() == ["ght-only"]
    assert loci["midpoint"].tolist() == [100]
    assert manifest["counts"] == {"chip_tiles": 1, "ght_tiles": 1, "rows": 1}
    assert manifest["outcome_used_for_selection"] is False
    assert manifest["outcome_used_for_centering"] is False
    assert manifest["output"]["sha256"]


def test_build_loci_cli_records_magix_filter(tmp_path) -> None:
    ght_path = tmp_path / "GCM1.magix.bed"
    output = tmp_path / "GCM1.ght_only.tsv.gz"
    ght_path.write_text(
        "chr\tstart\tstop\tname\tcoefficient.br\tcoefficient.ar\tfull_LL\t"
        "reduced_LL\tpvalue\tfdr\n"
        "chr1\t80\t280\tght\t0.2\t1.4\t-2\t-4\t0.001\t0.01\n",
        encoding="utf-8",
    )
    args = build_parser().parse_args(
        [
            "capselex",
            "build-loci",
            "--universe",
            "ght-only",
            "--ght-magix",
            str(ght_path),
            "--focal-tf",
            "GCM1",
            "--output",
            str(output),
        ]
    )
    _run_capselex(args)
    manifest = json.loads((tmp_path / "GCM1.ght_only.tsv.gz.manifest.json").read_text())
    assert manifest["counts"]["rows"] == 1
    assert manifest["ght_processing"] == {
        "fdr_max": 0.05,
        "format": "magix",
        "require_positive_score": True,
        "score_column": "coefficient.ar",
        "uses_chip_outcome": False,
    }


def test_replicate_qc() -> None:
    qc = replicate_qc(np.arange(10), np.arange(10) * 2)
    assert qc["spearman_log1p"] == pytest.approx(1.0)


def test_cross_assembly_bigwig_signal_preserves_mapping_audit(tmp_path) -> None:
    import pyBigWig

    bigwig_path = tmp_path / "source.bw"
    with pyBigWig.open(str(bigwig_path), "w") as bigwig:
        bigwig.addHeader([("chr1", 1000)])
        bigwig.addEntries(["chr1"], [100], ends=[200], values=[3.0])
    chain_path = tmp_path / "test.chain.gz"
    with gzip.open(chain_path, "wt", encoding="utf-8") as handle:
        handle.write("unused by fake liftOver\n")
    binary = tmp_path / "liftOver"
    binary.write_text(
        "#!/usr/bin/env python3\n"
        "import pathlib, sys\n"
        "source, chain, mapped, unmapped = map(pathlib.Path, sys.argv[-4:])\n"
        "lines = source.read_text().splitlines()\n"
        "pathlib.Path(mapped).write_text(lines[0] + '\\n')\n"
        "pathlib.Path(unmapped).write_text('#unmapped\\n' + lines[1] + '\\n')\n",
        encoding="utf-8",
    )
    binary.chmod(0o755)
    loci = pd.DataFrame(
        {"chrom": ["chr1", "chr1"], "start": [100, 300], "end": [200, 400]}
    )
    result = add_cross_assembly_bigwig_signal(
        loci, bigwig_path, chain_path, binary, "dnase", "hg19"
    )
    assert result["dnase_mapped"].tolist() == [True, False]
    assert result.loc[0, "dnase"] == pytest.approx(3.0)
    assert np.isnan(result.loc[1, "dnase"])
    assert result.loc[0, "dnase_hg19_start"] == 100


def test_nested_continuous_recovers_increment() -> None:
    rng = np.random.default_rng(9)
    n = 500
    baseline = rng.normal(size=n)
    grammar = rng.normal(size=n)
    frame = pd.DataFrame(
        {
            "chromosome_fold": np.tile(np.arange(5), n // 5),
            "baseline": baseline,
            "grammar": grammar,
            "outcome": baseline + 0.8 * grammar + rng.normal(scale=0.2, size=n),
        }
    )
    result = chromosome_nested_continuous(frame, "outcome", ["baseline"], ["grammar"])
    assert result.summary["partial_r2"] > 0.3
    assert result.predictions.notna().all().all()
    frame["chrom"] = [f"chr{index % 10 + 1}" for index in range(n)]
    bootstrap = chromosome_bootstrap_partial_r2(frame, result.predictions, replicates=100, seed=4)
    assert bootstrap["ci95_percentile"][0] > 0.3


def test_chromosome_shift_preserves_within_chromosome_values() -> None:
    frame = pd.DataFrame({"chrom": ["chr1"] * 4 + ["chr2"] * 4, "grammar": np.arange(8)})
    shifted = chromosome_shift_null(frame, ["grammar"], seed=2)
    assert sorted(shifted.loc[:3, "grammar"]) == [0, 1, 2, 3]
    assert not shifted["grammar"].equals(frame["grammar"])


def test_chromosome_shift_respects_tf_within_chromosome_blocks() -> None:
    frame = pd.DataFrame(
        {
            "focal_tf": ["A"] * 3 + ["B"] * 3,
            "chrom": ["chr1"] * 6,
            "grammar": np.arange(6),
        }
    )
    shifted = chromosome_shift_null(frame, ["grammar"], group_columns=["focal_tf"], seed=2)
    assert sorted(shifted.loc[:2, "grammar"]) == [0, 1, 2]
    assert sorted(shifted.loc[3:, "grammar"]) == [3, 4, 5]


def test_chromosome_shift_supports_noninteger_index() -> None:
    frame = pd.DataFrame(
        {"chrom": ["chr1"] * 4, "grammar": np.arange(4)}, index=["a", "b", "c", "d"]
    )
    shifted = chromosome_shift_null(frame, ["grammar"], seed=2)
    assert sorted(shifted["grammar"]) == [0, 1, 2, 3]
