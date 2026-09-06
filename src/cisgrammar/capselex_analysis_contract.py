from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path

from cisgrammar.capselex import write_json


@dataclass(frozen=True)
class AnalysisContract:
    locus_width_bp: int = 200
    sequence_context_bp: int = 400
    chromosome_folds: int = 5
    primary_locus_universe: str = "ght-only"
    primary_ght_asset_version: str = "codebook-v2"
    primary_ght_format: str = "magix"
    primary_ght_fdr_max: float = 0.05
    primary_ght_require_positive_refined_coefficient: bool = True
    primary_chip_pipeline: str = "McGill_GPHN_only"
    chip_replicates_per_tf: int = 2
    primary_chip_outcome: str = "mean_of_replicate_log1p_200bp_exact_mean_signal"
    replicate_confirmation: str = "same_incremental_effect_direction_in_both_replicates"
    sensitivity_locus_universes: tuple[str, ...] = (
        "assay-union",
        "fixed-genome",
        "legacy-assay-union",
    )
    screening_grammar_permutations: int = 100
    grammar_permutations: int = 1000
    primary_partial_r2: float = 0.005
    minimum_positive_focal_tfs: int = 4
    minimum_training_spacing_hits: int = 20
    promoter_half_width_bp: int = 2000
    claim: str = "incremental occupancy association"

    def validate(self) -> None:
        if self.sequence_context_bp < self.locus_width_bp:
            raise ValueError("sequence context must contain the outcome locus")
        if self.chromosome_folds < 3:
            raise ValueError("at least three chromosome folds are required")
        if self.primary_locus_universe != "ght-only":
            raise ValueError("the primary universe must be selected without ChIP outcome information")
        if self.primary_ght_asset_version != "codebook-v2":
            raise ValueError("the primary GHT asset must use the corrected Codebook v2 release")
        if self.primary_ght_format != "magix":
            raise ValueError("the primary GHT asset must use MAGIX output")
        if not 0 <= self.primary_ght_fdr_max <= 1:
            raise ValueError("the primary GHT FDR threshold must fall between zero and one")
        if not self.primary_ght_require_positive_refined_coefficient:
            raise ValueError("the primary GHT call must require a positive refined coefficient")
        if self.primary_chip_pipeline != "McGill_GPHN_only":
            raise ValueError("primary ChIP outcomes must use the frozen McGill GPHN pipeline")
        if self.chip_replicates_per_tf != 2:
            raise ValueError("the primary ChIP outcome requires exactly two biological replicates")
        if self.replicate_confirmation != "same_incremental_effect_direction_in_both_replicates":
            raise ValueError("replicate-resolved effect-direction confirmation is required")
        if self.screening_grammar_permutations < 1:
            raise ValueError("at least one null permutation is required")
        if self.grammar_permutations < self.screening_grammar_permutations:
            raise ValueError("the final null must use at least as many permutations as screening")

    def write(self, path: str | Path) -> None:
        self.validate()
        write_json({"schema_version": "analysis_contract_v4", **asdict(self)}, path)


PRIMARY_FOCAL_TFS = ("FLI1", "GABPA", "GCM1", "RFX5")
PARTNER_AVAILABILITY_NEGATIVE_CONTROLS = ("PAX7",)
LOW_CAP_COVERAGE_SENSITIVITY_TFS = ("MAX",)
FORBIDDEN_PRIMARY_ASSETS = (
    "TOPs",
    "CTOPs",
    "MOODS_Triple_Optimized",
    "Peaks_Triple_Optimized",
    "Peaks_Triple_Overlap",
)


def permutation_analysis_status(permutations: int, contract: AnalysisContract | None = None) -> str:
    """Classify a run without allowing non-frozen counts to masquerade as final."""
    active_contract = contract or AnalysisContract()
    if permutations < 1:
        raise ValueError("at least one permutation is required")
    if permutations == active_contract.grammar_permutations:
        return "final"
    if permutations == active_contract.screening_grammar_permutations:
        return "screening"
    return "exploratory"


def panel_result_conclusion(
    positive_focal_tfs: int,
    minimum_positive_focal_tfs: int,
) -> dict[str, object]:
    """Return the frozen panel gate and its publication-facing interpretation."""
    if positive_focal_tfs < 0 or minimum_positive_focal_tfs < 1:
        raise ValueError("panel counts must be non-negative with a positive minimum")
    panel_gate_passed = positive_focal_tfs >= minimum_positive_focal_tfs
    if panel_gate_passed:
        return {
            "scope": "panel-supported",
            "panel_gate": "passed",
            "summary": "Panel-supported result; panel gate passed.",
        }
    return {
        "scope": "TF-specific",
        "panel_gate": "not passed",
        "summary": "TF-specific result; panel gate not passed.",
    }


def assert_primary_asset_allowed(path: str | Path) -> None:
    text = str(path).lower()
    forbidden = [name for name in FORBIDDEN_PRIMARY_ASSETS if name.lower() in text]
    if forbidden:
        raise ValueError(f"circular primary asset is forbidden: {', '.join(forbidden)}")
