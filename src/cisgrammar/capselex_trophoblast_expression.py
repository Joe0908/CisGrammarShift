from __future__ import annotations

from pathlib import Path

import pandas as pd

from cisgrammar.capselex import normalize_tf_symbol


def read_rsem_gene_result(path: str | Path) -> pd.DataFrame:
    frame = pd.read_csv(path, sep="\t")
    lowered = {column.lower(): column for column in frame.columns}
    gene = lowered.get("gene_id") or lowered.get("gene")
    tpm = lowered.get("tpm")
    expected = lowered.get("expected_count")
    length = lowered.get("length")
    if gene is None or tpm is None:
        raise ValueError("RSEM gene result needs gene_id and TPM")
    return pd.DataFrame(
        {
            "gene": frame[gene].astype(str).str.split(".").str[0].map(normalize_tf_symbol),
            "tpm": pd.to_numeric(frame[tpm], errors="coerce").fillna(0.0),
            "expected_count": (
                pd.to_numeric(frame[expected], errors="coerce").fillna(0.0) if expected else 0.0
            ),
            "length": pd.to_numeric(frame[length], errors="coerce") if length else float("nan"),
        }
    )


def audit_partner_expression(
    files: list[str | Path],
    partners: list[str],
    threshold_tpm: float = 1.0,
) -> pd.DataFrame:
    rows = []
    for path in files:
        sample = Path(path).name.split(".genes.results")[0]
        frame = read_rsem_gene_result(path).set_index("gene")
        for partner in partners:
            symbol = normalize_tf_symbol(partner)
            tpm = float(frame.loc[symbol, "tpm"]) if symbol in frame.index else 0.0
            rows.append({"sample": sample, "partner": symbol, "tpm": tpm, "expressed": tpm >= threshold_tpm})
    return pd.DataFrame(rows)
