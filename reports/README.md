# Reference outputs

This directory contains machine-readable acquisition, QC, feature, model, and
external-context records from the frozen analysis. The files are retained for
numerical verification and do not constitute a prose results section.

The preferred entry points are:

- `gphn_capselex_primary_genomic_model_screening.json` — primary panel summary;
- `trophoblast_tgif2_gcm1_replication_screening.json` — external-context screen;
- `capselex_nature_supplement_audit.json` and
  `codebook_ght_v2_rebuild_audit.json` — source/QC audits;
- `*_resolved.json` — downloaded-asset hashes and resolved manifests.

Large feature matrices and raw sequencing/track files are excluded from Git and
must be reconstructed from the manifests in `configs/`.
