from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from cisgrammar.capselex_chip_qc import replicate_qc
from cisgrammar.capselex_genomic_assets import build_assay_union_loci
from cisgrammar.capselex_genomic_model import chromosome_nested_continuous, chromosome_shift_null


def test_assay_union_deduplicates_fixed_tiles() -> None:
    chip = pd.DataFrame({"chrom": ["chr1", "chr1"], "midpoint": [100, 110]})
    ght = pd.DataFrame({"chrom": ["chr1", "chr2"], "midpoint": [105, 300]})
    loci = build_assay_union_loci(chip, ght, "GCM1", width_bp=200)
    assert len(loci) == 2
    first = loci[loci["chrom"] == "chr1"].iloc[0]
    assert first["chip_source"] and first["ght_source"]


def test_replicate_qc() -> None:
    qc = replicate_qc(np.arange(10), np.arange(10) * 2)
    assert qc["spearman_log1p"] == pytest.approx(1.0)


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


def test_chromosome_shift_preserves_within_chromosome_values() -> None:
    frame = pd.DataFrame({"chrom": ["chr1"] * 4 + ["chr2"] * 4, "grammar": np.arange(8)})
    shifted = chromosome_shift_null(frame, ["grammar"], seed=2)
    assert sorted(shifted.loc[:3, "grammar"]) == [0, 1, 2, 3]
    assert not shifted["grammar"].equals(frame["grammar"])
