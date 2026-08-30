from __future__ import annotations

import hashlib
import json

from cisgrammar.capselex_codebook_v2 import audit_codebook_v2_focal_panel


def test_codebook_v2_panel_audit_validates_official_member(tmp_path) -> None:
    bed = tmp_path / "GCM1.bed"
    bed.write_text(
        "chr\tstart\tstop\tname\tcoefficient.br\tcoefficient.ar\tfull_LL\t"
        "reduced_LL\tpvalue\tfdr\n"
        "chr1\t100\t300\tchr1:100-300\t0.1\t1.2\t-1\t-2\t0.001\t0.01\n",
        encoding="utf-8",
    )
    digest = hashlib.sha256(bed.read_bytes()).hexdigest()
    config = {
        "source": {"primary_analysis_eligible": True},
        "call_rule": {"fdr_max": 0.05, "require_positive_score": True},
        "panel": [
            {
                "tf": "GCM1",
                "local_filename": "GCM1.bed",
                "size_bytes": bed.stat().st_size,
                "sha256": digest,
                "geo_member_sha256": digest,
            }
        ],
    }
    config_path = tmp_path / "config.json"
    config_path.write_text(json.dumps(config), encoding="utf-8")
    result = audit_codebook_v2_focal_panel(config_path, tmp_path)
    assert result["primary_analysis_eligible"] is True
    assert result["focal_members_byte_identical_to_geo"] is True
    assert result["files"][0]["rows_selected"] == 1

