from __future__ import annotations

import argparse
import json
from collections import defaultdict
from pathlib import Path

from cisgrammar.capselex import asset_manifest, normalize_tf_symbol, write_json
from cisgrammar.capselex_monomer_motifs import (
    pwm_consensus,
    pwm_information_bits,
    read_jaspar_collection,
    read_mex_top1_archive,
)
from cisgrammar.capselex_nature_supplements import read_pwm_models


def main() -> None:
    parser = argparse.ArgumentParser(description="Audit monomer PWM coverage for the CAP focal panel")
    parser.add_argument("--mex-top1", type=Path, required=True)
    parser.add_argument("--jaspar", type=Path, required=True)
    parser.add_argument("--cap-pwm-table", type=Path, required=True)
    parser.add_argument("--cap-audit", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    mex_records = read_mex_top1_archive(args.mex_top1)
    mex = {normalize_tf_symbol(record.symbol): record for record in mex_records}
    jaspar: dict[str, list] = defaultdict(list)
    for record in read_jaspar_collection(args.jaspar):
        if "::" not in record.symbol:
            jaspar[normalize_tf_symbol(record.symbol)].append(record)
    cap_individual: dict[str, list] = defaultdict(list)
    for record in read_pwm_models(args.cap_pwm_table):
        if record.members is None:
            cap_individual[normalize_tf_symbol(record.pair)].append(record)

    cap = json.loads(args.cap_audit.read_text(encoding="utf-8"))
    panel = []
    for focal in cap["panel"]:
        tf = normalize_tf_symbol(focal["tf"])
        representative_partners = set()
        for pair in focal["representative_pwm_pairs"]:
            first, second = pair.split("_", 1)
            representative_partners.add(second if normalize_tf_symbol(first) == tf else first)
        partner_sources = []
        for raw_partner in sorted(representative_partners):
            partner = normalize_tf_symbol(raw_partner)
            if partner in mex:
                source = "Codebook_MEX_top1"
                profiles = [mex[partner].source_id]
            elif partner in cap_individual:
                source = "CAP_Supplementary_Table_3_individual"
                profiles = sorted(record.model_id for record in cap_individual[partner])
            elif partner in jaspar:
                source = "JASPAR2024_CORE_vertebrates"
                profiles = sorted(record.source_id for record in jaspar[partner])
            else:
                source = None
                profiles = []
            partner_sources.append(
                {"partner": partner, "source": source, "profile_ids": profiles}
            )
        missing_partners = {
            record["partner"] for record in partner_sources if record["source"] is None
        }
        complete_pairs = []
        excluded_pairs = []
        for pair in focal["representative_pwm_pairs"]:
            first, second = pair.split("_", 1)
            partner = normalize_tf_symbol(second if normalize_tf_symbol(first) == tf else first)
            (excluded_pairs if partner in missing_partners else complete_pairs).append(pair)
        focal_record = mex.get(tf)
        panel.append(
            {
                "tf": tf,
                "focal_mex_profile": focal_record.source_id if focal_record else None,
                "focal_pwm_length": focal_record.pwm.length if focal_record else None,
                "focal_pwm_consensus": pwm_consensus(focal_record.pwm) if focal_record else None,
                "focal_pwm_information_bits": (
                    pwm_information_bits(focal_record.pwm) if focal_record else None
                ),
                "representative_partners": len(representative_partners),
                "partners_with_monomer_pwm": sum(record["source"] is not None for record in partner_sources),
                "partners_missing_monomer_pwm": [
                    record["partner"] for record in partner_sources if record["source"] is None
                ],
                "representative_pairs_with_complete_monomers": sorted(complete_pairs),
                "excluded_representative_pairs_without_partner_monomer": sorted(excluded_pairs),
                "partner_sources": partner_sources,
            }
        )

    report = {
        "schema_version": "monomer_pwm_panel_audit_v1",
        "source_assets": asset_manifest(
            [args.mex_top1, args.jaspar, args.cap_pwm_table, args.cap_audit]
        ),
        "source_policy": (
            "Codebook MEX top-1 is preferred because motifs were ranked on held-out data; "
            "individual HT-SELEX PWMs in CAP Supplementary Table 3 are second; JASPAR 2024 "
            "CORE vertebrate monomer profiles are the frozen fallback."
        ),
        "mex_top1_profiles": len(mex_records),
        "cap_individual_pwm_models": sum(len(records) for records in cap_individual.values()),
        "jaspar_monomer_profiles": sum(len(records) for records in jaspar.values()),
        "panel": panel,
        "claim_boundary": (
            "Profiles lacking a frozen partner monomer are excluded before outcome analysis. Coverage "
            "does not guarantee that a partner is expressed or present in the ChIP cell state."
        ),
    }
    write_json(report, args.output)
    for record in panel:
        print(
            f"{record['tf']}\tfocal_mex={record['focal_mex_profile'] is not None}\t"
            f"partner_pwm={record['partners_with_monomer_pwm']}/{record['representative_partners']}"
        )


if __name__ == "__main__":
    main()
