from __future__ import annotations

import argparse
import gzip
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd

from cisgrammar.capselex import asset_manifest, normalize_tf_symbol, write_json

ATTRIBUTE = re.compile(r'(\S+) "([^"]+)";')


def read_transcript_to_gene(path: Path) -> dict[str, str]:
    mapping: dict[str, str] = {}
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for raw in handle:
            if raw.startswith("#"):
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "transcript":
                continue
            attributes = dict(ATTRIBUTE.findall(fields[8]))
            if "transcript_id" not in attributes or "gene_name" not in attributes:
                continue
            transcript = attributes["transcript_id"].split(".", 1)[0]
            gene = normalize_tf_symbol(attributes["gene_name"])
            previous = mapping.setdefault(transcript, gene)
            if previous != gene:
                raise ValueError(f"transcript maps to multiple genes: {transcript}")
    if not mapping:
        raise ValueError("GENCODE transcript-to-gene mapping is empty")
    return mapping


def read_gene_tpm(path: Path, transcript_to_gene: dict[str, str]) -> tuple[pd.Series, dict[str, float]]:
    frame = pd.read_csv(path, sep="\t", compression="gzip")
    required = {"Name", "TPM"}
    if not required.issubset(frame):
        raise ValueError(f"Salmon table lacks columns: {sorted(required - set(frame))}")
    frame["transcript"] = frame["Name"].astype(str).str.split(".").str[0]
    frame["gene"] = frame["transcript"].map(transcript_to_gene)
    frame["TPM"] = pd.to_numeric(frame["TPM"], errors="raise")
    if frame["TPM"].lt(0).any():
        raise ValueError("Salmon TPM values cannot be negative")
    mapped = frame["gene"].notna()
    summary = {
        "transcripts": int(len(frame)),
        "mapped_transcripts": int(mapped.sum()),
        "transcript_mapping_fraction": float(mapped.mean()),
        "total_tpm": float(frame["TPM"].sum()),
        "mapped_tpm": float(frame.loc[mapped, "TPM"].sum()),
        "tpm_mapping_fraction": float(frame.loc[mapped, "TPM"].sum() / frame["TPM"].sum()),
    }
    return frame.loc[mapped].groupby("gene")["TPM"].sum(), summary


def main() -> None:
    parser = argparse.ArgumentParser(description="Audit HEK293 expression of frozen CAP partners")
    parser.add_argument("--replicate-1", type=Path, required=True)
    parser.add_argument("--replicate-2", type=Path, required=True)
    parser.add_argument("--gencode", type=Path, required=True)
    parser.add_argument("--monomer-audit", type=Path, required=True)
    parser.add_argument("--minimum-tpm-each-replicate", type=float, default=1.0)
    parser.add_argument("--output-table", type=Path, required=True)
    parser.add_argument("--output-report", type=Path, required=True)
    args = parser.parse_args()
    if args.minimum_tpm_each_replicate < 0:
        raise ValueError("expression threshold cannot be negative")

    transcript_to_gene = read_transcript_to_gene(args.gencode)
    first, first_summary = read_gene_tpm(args.replicate_1, transcript_to_gene)
    second, second_summary = read_gene_tpm(args.replicate_2, transcript_to_gene)
    audit = json.loads(args.monomer_audit.read_text(encoding="utf-8"))
    rows = []
    for focal in audit["panel"]:
        complete_pairs = set(focal["representative_pairs_with_complete_monomers"])
        for pair in sorted(complete_pairs):
            first_tf, second_tf = pair.split("_", 1)
            focal_tf = normalize_tf_symbol(focal["tf"])
            partner = normalize_tf_symbol(
                second_tf if normalize_tf_symbol(first_tf) == focal_tf else first_tf
            )
            rep1 = float(first.get(partner, 0.0))
            rep2 = float(second.get(partner, 0.0))
            rows.append(
                {
                    "focal_tf": focal_tf,
                    "pair": pair,
                    "partner": partner,
                    "replicate_1_tpm": rep1,
                    "replicate_2_tpm": rep2,
                    "mean_tpm": (rep1 + rep2) / 2,
                    "minimum_replicate_tpm": min(rep1, rep2),
                    "expression_supported": (
                        min(rep1, rep2) >= args.minimum_tpm_each_replicate
                    ),
                }
            )
    table = pd.DataFrame(rows).drop_duplicates(["focal_tf", "pair"])
    if table.empty:
        raise ValueError("no frozen CAP pairs were found")
    args.output_table.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(args.output_table, sep="\t", index=False, float_format="%.10g")

    panel = []
    for focal_tf, group in table.groupby("focal_tf", sort=True):
        panel.append(
            {
                "focal_tf": focal_tf,
                "pairs": len(group),
                "expression_supported_pairs": int(group["expression_supported"].sum()),
                "unsupported_pairs": group.loc[~group["expression_supported"], "pair"].tolist(),
                "supported_pairs": group.loc[group["expression_supported"], "pair"].tolist(),
            }
        )
    replicate_log_correlation = float(
        np.corrcoef(np.log1p(table["replicate_1_tpm"]), np.log1p(table["replicate_2_tpm"]))[0, 1]
    )
    write_json(
        {
            "schema_version": "hek293_partner_expression_audit_v1",
            "biosample": "independent wild-type HEK293 RNA-seq",
            "geo_sample": "GSM3611199",
            "transcript_mapping": "GENCODE v29 transcript IDs with version suffix removed",
            "replicate_1_mapping": first_summary,
            "replicate_2_mapping": second_summary,
            "partner_replicate_log1p_pearson": replicate_log_correlation,
            "expression_supported_definition": (
                f"gene-level summed transcript TPM >= {args.minimum_tpm_each_replicate:g} "
                "in each of two independent HEK293 replicates"
            ),
            "threshold_role": (
                "outcome-independent audit; primary-versus-sensitivity placement is decided only "
                "after coverage is inspected and documented"
            ),
            "panel": panel,
            "source_assets": asset_manifest(
                [args.replicate_1, args.replicate_2, args.gencode, args.monomer_audit]
            ),
            "output_table": asset_manifest([args.output_table])[0],
            "claim_boundary": (
                "Expression in an independent HEK293 culture supports partner availability but does "
                "not establish protein abundance, simultaneous occupancy, or a physical interaction."
            ),
        },
        args.output_report,
    )
    for record in panel:
        print(
            f"{record['focal_tf']}\t{record['expression_supported_pairs']}/{record['pairs']} "
            "expression-supported pairs"
        )


if __name__ == "__main__":
    main()
