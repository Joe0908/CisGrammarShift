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
        if self.screening_grammar_permutations < 1:
            raise ValueError("at least one null permutation is required")
        if self.grammar_permutations < self.screening_grammar_permutations:
            raise ValueError("the final null must use at least as many permutations as screening")

    def write(self, path: str | Path) -> None:
        self.validate()
        write_json({"schema_version": "analysis_contract_v2", **asdict(self)}, path)


PRIMARY_FOCAL_TFS = ("FLI1", "GABPA", "GCM1", "PAX7", "RFX5")
CAP_NULL_CONTROL = "MAX"
FORBIDDEN_PRIMARY_ASSETS = (
    "TOPs",
    "CTOPs",
    "MOODS_Triple_Optimized",
    "Peaks_Triple_Optimized",
    "Peaks_Triple_Overlap",
)


def assert_primary_asset_allowed(path: str | Path) -> None:
    text = str(path).lower()
    forbidden = [name for name in FORBIDDEN_PRIMARY_ASSETS if name.lower() in text]
    if forbidden:
        raise ValueError(f"circular primary asset is forbidden: {', '.join(forbidden)}")
