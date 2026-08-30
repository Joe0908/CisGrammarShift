from __future__ import annotations

import numpy as np
import pandas as pd

from cisgrammar.capselex_grammar_calibration import build_ets_grammar_score, cross_fitted_residualize
from cisgrammar.capselex_grammar_scoring import GrammarProfile, binned_monomer_residuals, scan_grammar
from cisgrammar.capselex_monomer_motifs import (
    PWM,
    maximum_pwm_score,
    pwm_consensus,
    pwm_information_bits,
    read_jaspar_collection,
    read_mex_top1_archive,
    scan_pwm,
)


def deterministic_pwm(name: str, motif: str) -> PWM:
    matrix = np.full((len(motif), 4), 0.01)
    index = {base: position for position, base in enumerate("ACGT")}
    for offset, base in enumerate(motif):
        matrix[offset, index[base]] = 0.97
    return PWM(name, matrix)


def test_pwm_scans_both_strands() -> None:
    pwm = deterministic_pwm("test", "AAC")
    assert maximum_pwm_score("TTTGTTAAA", pwm) == scan_pwm("TTTGTTAAA", pwm).max()
    assert np.isfinite(maximum_pwm_score("TTTGTTAAA", pwm))


def test_grammar_prefers_specified_gap() -> None:
    focal = deterministic_pwm("focal", "AAA")
    partner = deterministic_pwm("partner", "CCC")
    profile = GrammarProfile("A_C", focal, partner, gaps=(2,), orientations=("FF",))
    good = scan_grammar("GGGAAATTCCCGGG", profile)
    bad = scan_grammar("GGGAAATTTTCCCG", profile)
    assert good > bad


def test_read_monomer_pwm_archives(tmp_path) -> None:
    import zipfile

    mex = tmp_path / "top1.zip"
    with zipfile.ZipFile(mex, "w") as archive:
        archive.writestr(
            "GCM1.NA@PBM@test.pcm",
            ">GCM1.NA@PBM@test\n10 0 0 0\n0 10 0 0\n0 0 10 0\n0 0 0 10\n",
        )
    mex_records = read_mex_top1_archive(mex)
    assert mex_records[0].symbol == "GCM1"
    assert pwm_consensus(mex_records[0].pwm) == "ACGT"
    assert pwm_information_bits(mex_records[0].pwm) > 7.9

    jaspar = tmp_path / "jaspar.txt"
    jaspar.write_text(
        ">MA0001.1 GCM1\n"
        "A [10 0]\n"
        "C [0 10]\n"
        "G [0 0]\n"
        "T [0 0]\n",
        encoding="utf-8",
    )
    jaspar_records = read_jaspar_collection(jaspar)
    assert jaspar_records[0].source_id == "MA0001.1"
    assert jaspar_records[0].symbol == "GCM1"
    assert pwm_consensus(jaspar_records[0].pwm) == "AC"


def test_training_only_residualization_removes_bin_mean() -> None:
    grammar = np.arange(20, dtype=float)
    monomer = np.repeat(np.arange(4), 5).astype(float)
    train = np.arange(20) < 15
    residual = binned_monomer_residuals(grammar, monomer, monomer, train, bins=4)
    assert np.isfinite(residual).all()
    assert abs(residual[train].mean()) < 1e-9


def test_cross_fitted_ets_score_has_no_missing_values() -> None:
    frame = pd.DataFrame(
        {
            "chromosome_fold": np.tile(np.arange(5), 10),
            "raw_a": np.linspace(0, 2, 50),
            "raw_b": np.linspace(1, 3, 50),
            "focal": np.linspace(0, 1, 50),
            "partner_a": np.linspace(0.2, 1.2, 50),
            "partner_b": np.linspace(0.3, 1.3, 50),
        }
    )
    residual = cross_fitted_residualize(frame, "raw_a", "focal", "partner_a", bins=5)
    assert np.isfinite(residual).all()
    score = build_ets_grammar_score(
        frame,
        ["raw_a", "raw_b"],
        "focal",
        ["partner_a", "partner_b"],
    )
    assert score.shape == (50,)
