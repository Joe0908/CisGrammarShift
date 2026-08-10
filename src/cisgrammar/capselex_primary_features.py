from __future__ import annotations

import json
import os
import subprocess
import tempfile
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from cisgrammar.capselex import normalize_tf_symbol, require_columns
from cisgrammar.capselex_genomic_assets import read_chrom_sizes
from cisgrammar.capselex_grammar_calibration import cross_fitted_residualize
from cisgrammar.capselex_monomer_motifs import (
    PWM,
    MonomerPWMRecord,
    read_jaspar_collection,
    read_mex_top1_archive,
)
from cisgrammar.capselex_nature_supplements import CapPWMRecord, read_pwm_models
from cisgrammar.capselex_pwm_batch import maximum_pwm_scores


@dataclass(frozen=True)
class PairMechanismProfile:
    pair: str
    partner: str
    mechanism: str
    cap_records: tuple[CapPWMRecord, ...]
    partner_records: tuple[MonomerPWMRecord, ...]
    partner_source: str

    @property
    def feature_id(self) -> str:
        return f"{self.mechanism}::{self.pair}"


@dataclass(frozen=True)
class FrozenPWMPanel:
    focal_tf: str
    focal_record: MonomerPWMRecord
    profiles: tuple[PairMechanismProfile, ...]
    excluded_pairs: tuple[str, ...]


def _profile_lookup(records: list[MonomerPWMRecord]) -> dict[str, MonomerPWMRecord]:
    lookup = {record.source_id: record for record in records}
    if len(lookup) != len(records):
        raise ValueError("monomer source identifiers must be unique")
    return lookup


def load_frozen_pwm_panel(
    focal_tf: str,
    cap_pwm_table: str | Path,
    mex_top1: str | Path,
    jaspar: str | Path,
    monomer_audit: str | Path,
) -> FrozenPWMPanel:
    """Load exactly the author-representative CAP and pre-audited monomer profiles."""
    focal_tf = normalize_tf_symbol(focal_tf)
    audit = json.loads(Path(monomer_audit).read_text(encoding="utf-8"))
    try:
        focal_audit = next(record for record in audit["panel"] if record["tf"] == focal_tf)
    except StopIteration as error:
        raise ValueError(f"{focal_tf} is absent from the frozen monomer audit") from error

    mex_records = read_mex_top1_archive(mex_top1)
    jaspar_records = [
        record for record in read_jaspar_collection(jaspar) if "::" not in record.symbol
    ]
    cap_records = read_pwm_models(cap_pwm_table)
    cap_individual = [
        MonomerPWMRecord(record.model_id, normalize_tf_symbol(record.pair), record.to_pwm())
        for record in cap_records
        if record.members is None
    ]
    source_records = {
        "Codebook_MEX_top1": mex_records,
        "CAP_Supplementary_Table_3_individual": cap_individual,
        "JASPAR2024_CORE_vertebrates": jaspar_records,
    }
    source_lookups = {source: _profile_lookup(records) for source, records in source_records.items()}
    focal_id = focal_audit["focal_mex_profile"]
    if focal_id is None or focal_id not in source_lookups["Codebook_MEX_top1"]:
        raise ValueError(f"{focal_tf} lacks its frozen Codebook MEX focal profile")

    complete_pairs = set(focal_audit["representative_pairs_with_complete_monomers"])
    excluded_pairs = tuple(focal_audit["excluded_representative_pairs_without_partner_monomer"])
    partner_audit = {record["partner"]: record for record in focal_audit["partner_sources"]}
    representative: dict[tuple[str, str], list[CapPWMRecord]] = defaultdict(list)
    for record in cap_records:
        if record.representative and record.pair in complete_pairs:
            if record.interaction_type not in {"composite", "spacing"}:
                raise ValueError(f"representative CAP PWM lacks a mechanism: {record.model_id}")
            representative[(record.pair, record.interaction_type)].append(record)

    found_pairs = {pair for pair, _ in representative}
    if found_pairs != complete_pairs:
        missing = ", ".join(sorted(complete_pairs - found_pairs))
        unexpected = ", ".join(sorted(found_pairs - complete_pairs))
        raise ValueError(f"CAP representative pair mismatch; missing={missing}; unexpected={unexpected}")

    profiles = []
    for (pair, mechanism), models in sorted(representative.items()):
        first, second = pair.split("_", 1)
        partner = normalize_tf_symbol(second if normalize_tf_symbol(first) == focal_tf else first)
        partner_record = partner_audit[partner]
        source = partner_record["source"]
        if source not in source_lookups:
            raise ValueError(f"unsupported frozen monomer source for {partner}: {source}")
        ids = partner_record["profile_ids"]
        missing_ids = [source_id for source_id in ids if source_id not in source_lookups[source]]
        if missing_ids:
            raise ValueError(f"missing frozen {source} profiles for {partner}: {missing_ids}")
        profiles.append(
            PairMechanismProfile(
                pair=pair,
                partner=partner,
                mechanism=mechanism,
                cap_records=tuple(sorted(models, key=lambda record: record.model_id)),
                partner_records=tuple(source_lookups[source][source_id] for source_id in ids),
                partner_source=source,
            )
        )
    return FrozenPWMPanel(
        focal_tf=focal_tf,
        focal_record=source_lookups["Codebook_MEX_top1"][focal_id],
        profiles=tuple(profiles),
        excluded_pairs=excluded_pairs,
    )


def _read_ordered_fasta(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    name: str | None = None
    parts: list[str] = []
    with path.open(encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(parts).upper()))
                name = line[1:].split()[0]
                parts = []
            elif name is not None:
                parts.append(line)
    if name is not None:
        records.append((name, "".join(parts).upper()))
    return records


def extract_twobit_contexts(
    loci: pd.DataFrame,
    two_bit: str | Path,
    two_bit_to_fa: str | Path,
    chrom_sizes_path: str | Path,
    width_bp: int = 400,
) -> list[str]:
    """Extract midpoint-centred contexts from a UCSC twoBit with audited edge padding."""
    if width_bp <= 0 or width_bp % 2:
        raise ValueError("context width must be a positive even number")
    require_columns(loci, ("chrom", "midpoint"), "locus table")
    sizes = dict(read_chrom_sizes(chrom_sizes_path)[["chrom", "size"]].itertuples(index=False))
    half = width_bp // 2
    intervals: list[tuple[str, int, int, int, int]] = []
    for chrom, midpoint in loci[["chrom", "midpoint"]].itertuples(index=False, name=None):
        if chrom not in sizes:
            raise ValueError(f"chromosome absent from size table: {chrom}")
        requested_start = int(midpoint) - half
        requested_end = requested_start + width_bp
        start = max(0, requested_start)
        end = min(int(sizes[chrom]), requested_end)
        intervals.append((str(chrom), start, end, start - requested_start, requested_end - end))

    with tempfile.TemporaryDirectory(prefix="cisgrammar_twobit_") as temporary:
        directory = Path(temporary)
        bed = directory / "contexts.bed"
        fasta = directory / "contexts.fa"
        with bed.open("w", encoding="utf-8") as handle:
            for index, (chrom, start, end, _, _) in enumerate(intervals):
                handle.write(f"{chrom}\t{start}\t{end}\tcontext_{index}\n")
        two_bit_binary = Path(two_bit_to_fa).resolve()
        if not os.access(two_bit_binary, os.X_OK):
            two_bit_binary.chmod(two_bit_binary.stat().st_mode | 0o111)
        command = [
            str(two_bit_binary),
            f"-bed={bed}",
            "-bedPos",
            str(Path(two_bit).resolve()),
            str(fasta),
        ]
        subprocess.run(command, check=True, capture_output=True, text=True)
        records = _read_ordered_fasta(fasta)

    if len(records) != len(intervals):
        raise RuntimeError(f"twoBitToFa returned {len(records)} contexts for {len(intervals)} loci")
    sequences = []
    for (header, sequence), (chrom, start, end, left_pad, right_pad) in zip(
        records, intervals, strict=True
    ):
        expected_header = f"{chrom}:{start}-{end}"
        if header != expected_header or len(sequence) != end - start:
            raise RuntimeError(f"unexpected twoBitToFa record: {header}")
        sequence = "N" * left_pad + sequence + "N" * right_pad
        if len(sequence) != width_bp:
            raise RuntimeError("twoBit context padding produced the wrong width")
        sequences.append("".join(base if base in "ACGT" else "N" for base in sequence))
    return sequences


def sequence_covariates(sequences: list[str]) -> pd.DataFrame:
    rows = []
    for sequence in sequences:
        canonical = sum(base in "ACGT" for base in sequence)
        denominator = max(canonical, 1)
        rows.append(
            {
                "gc_fraction": (sequence.count("G") + sequence.count("C")) / denominator,
                "cpg_per_100bp": 100.0 * sequence.count("CG") / max(len(sequence) - 1, 1),
                "n_fraction": 1.0 - canonical / max(len(sequence), 1),
            }
        )
    return pd.DataFrame(rows)


def longest_canonical_run(sequence: str) -> int:
    longest = 0
    current = 0
    for base in sequence:
        if base in "ACGT":
            current += 1
            longest = max(longest, current)
        else:
            current = 0
    return longest


def add_expression_supported_aggregates(
    frame: pd.DataFrame,
    focal_tf: str,
    expression: pd.DataFrame,
) -> tuple[pd.DataFrame, dict[str, object]]:
    """Aggregate CAP residuals using only independently expression-supported partners."""
    focal = expression[expression["focal_tf"].eq(focal_tf)].copy()
    values = focal["expression_supported"]
    if values.dtype != bool:
        normalized = values.astype(str).str.lower()
        if not normalized.isin({"true", "false"}).all():
            raise ValueError("expression_supported must contain only true/false values")
        focal["expression_supported"] = normalized.eq("true")
    supported = set(focal.loc[focal["expression_supported"], "pair"])
    pair_residuals: dict[str, list[str]] = {}
    mechanism_residuals: dict[str, list[str]] = {"composite": [], "spacing": []}
    partner_columns = []
    for pair in sorted(supported):
        pair_columns = [
            column
            for column in frame
            if column.startswith("cap_excess::") and column.endswith(f"::{pair}")
        ]
        if not pair_columns:
            raise ValueError(f"expression-supported pair lacks CAP residual features: {pair}")
        pair_residuals[pair] = pair_columns
        for column in pair_columns:
            mechanism = column.split("::", 2)[1]
            mechanism_residuals[mechanism].append(column)
        partner_column = f"partner_raw::{pair}"
        if partner_column not in frame:
            raise ValueError(f"expression-supported pair lacks partner monomer feature: {pair}")
        partner_columns.append(partner_column)

    result = frame.copy()
    if supported:
        pair_values = [result[columns].mean(axis=1).to_numpy() for columns in pair_residuals.values()]
        result["cap_expression_supported_excess_mean"] = np.mean(
            np.column_stack(pair_values), axis=1
        )
        result["partner_monomer_expressed_max"] = result[partner_columns].max(axis=1)
        for mechanism, columns in mechanism_residuals.items():
            if columns:
                result[f"cap_expression_supported_{mechanism}_excess_mean"] = result[
                    columns
                ].mean(axis=1)
    return result, {
        "expression_supported_pairs": sorted(supported),
        "expression_supported_pair_count": len(supported),
        "unsupported_pairs": sorted(set(focal["pair"]) - supported),
        "aggregation": "mean across mechanisms within pair then equal-weight mean across pairs",
    }


def score_frozen_pwm_panel(
    sequences: list[str],
    panel: FrozenPWMPanel,
    chromosome_folds: np.ndarray,
    batch_size: int = 512,
    residual_bins: int = 20,
) -> tuple[pd.DataFrame, dict[str, object]]:
    """Score representative CAP PWMs and construct leakage-controlled excess scores."""
    pwm_by_id: dict[str, PWM] = {panel.focal_record.source_id: panel.focal_record.pwm}
    for profile in panel.profiles:
        for record in profile.cap_records:
            pwm_by_id[record.model_id] = record.to_pwm()
        for record in profile.partner_records:
            pwm_by_id[record.source_id] = record.pwm
    pwm_ids = list(pwm_by_id)
    score_matrix = maximum_pwm_scores(sequences, [pwm_by_id[source_id] for source_id in pwm_ids], batch_size)
    if not np.isfinite(score_matrix).all():
        raise ValueError("all contexts must contain a canonical window for every frozen PWM")
    score_index = {source_id: index for index, source_id in enumerate(pwm_ids)}

    frame = pd.DataFrame({"chromosome_fold": np.asarray(chromosome_folds, dtype=int)})
    focal_column = "focal_monomer_score"
    frame[focal_column] = score_matrix[:, score_index[panel.focal_record.source_id]]
    residual_columns: dict[str, list[str]] = defaultdict(list)
    pair_residual_columns: dict[str, list[str]] = defaultdict(list)
    partner_columns = []
    profile_manifest = []
    for profile in panel.profiles:
        partner_column = f"partner_raw::{profile.pair}"
        if partner_column not in frame:
            partner_indices = [score_index[record.source_id] for record in profile.partner_records]
            frame[partner_column] = np.max(score_matrix[:, partner_indices], axis=1)
            partner_columns.append(partner_column)
        raw_column = f"cap_raw::{profile.feature_id}"
        cap_indices = [score_index[record.model_id] for record in profile.cap_records]
        frame[raw_column] = np.max(score_matrix[:, cap_indices], axis=1)
        residual_column = f"cap_excess::{profile.feature_id}"
        frame[residual_column] = cross_fitted_residualize(
            frame,
            raw_column,
            focal_column,
            partner_column,
            bins=residual_bins,
        )
        residual_columns[profile.mechanism].append(residual_column)
        pair_residual_columns[profile.pair].append(residual_column)
        profile_manifest.append(
            {
                "pair": profile.pair,
                "partner": profile.partner,
                "mechanism": profile.mechanism,
                "cap_model_ids": [record.model_id for record in profile.cap_records],
                "partner_source": profile.partner_source,
                "partner_profile_ids": [record.source_id for record in profile.partner_records],
                "raw_aggregation": "maximum_across_author_representative_models_within_pair_mechanism",
            }
        )

    # Defragment before aggregate operations because focal panels can contain dozens of pairs.
    frame = frame.copy()
    frame["partner_monomer_max"] = frame[partner_columns].max(axis=1)
    for mechanism, columns in sorted(residual_columns.items()):
        aggregate = f"cap_{mechanism}_excess_mean"
        frame[aggregate] = frame[columns].mean(axis=1)
    pair_values = [frame[columns].mean(axis=1).to_numpy() for columns in pair_residual_columns.values()]
    frame["cap_all_excess_mean"] = np.mean(np.column_stack(pair_values), axis=1)
    manifest = {
        "focal_tf": panel.focal_tf,
        "focal_monomer_profile_id": panel.focal_record.source_id,
        "focal_monomer_source": "Codebook_MEX_top1",
        "excluded_pairs_without_partner_monomer": list(panel.excluded_pairs),
        "profile_units": profile_manifest,
        "pair_mechanism_units": len(profile_manifest),
        "all_mechanism_aggregation": (
            "mean_across_mechanisms_within_pair_then_equal_weight_mean_across_pairs"
        ),
        "mechanism_branch_aggregation": "equal_weight_mean_across_pairs",
        "residualization": {
            "method": "training_chromosome_fold_derived_2d_monomer_quantile_bins",
            "bins_per_monomer": residual_bins,
            "fold_column": "chromosome_fold",
        },
    }
    return frame.drop(columns="chromosome_fold"), manifest
