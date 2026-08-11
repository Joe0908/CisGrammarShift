# Toronto-GPZN fallback and processing sensitivity

This directory documents the retained fallback used before the author-supplied McGill-GPHN panel became
available. It is provenance and a no-retuning processing-pipeline sensitivity, not the default analysis.

## What is retained

The fallback does **not** generate 12 bigWigs from the 175.39-GB merged Codebook archive. It resolves and
downloads 12 existing, per-sample Toronto-GPZN bigWigs from GEO `GSE280248` (six TFs by two biological
replicates):

- `scripts/build_codebook_chip_manifest.py` derives the whitelist from GEO metadata;
- `configs/codebook_chip_gpzn_panel_manifest.json` freezes the 12 public files;
- `reports/codebook_chip_gpzn_panel_resolved.json` records downloaded-file hashes;
- `reports/capselex_primary_genomic_model_screening.json` is the historical GPZN model report.

The large merged archive is therefore unnecessary for both the primary and fallback workflows.

## Current role

The author-supplied McGill-GPHN panel in `configs/codebook_chip_gphn_panel_manifest.json` is the primary
outcome source. The same frozen loci, features, expression gate, model specification, chromosome splits,
seed, permutation count, and `partial R2 >= 0.005` decision rule were then applied to Toronto-GPZN without
retuning. `reports/gphn_gpzn_pipeline_comparison.json` records that comparison.

Do not mix GPHN and GPZN replicates within a TF or use one pipeline to choose features for the other. The
comparison supports robustness to processing choice only; it does not establish co-occupancy or causality.

## Re-running the fallback

Follow the GEO metadata, manifest-building, and download commands in `docs/data_access.md`, then run the
feature and model scripts with the GPZN directories and:

```bash
python scripts/run_capselex_primary_models.py \
  --feature-directory results/capselex/primary_features_gpzn \
  --chip-processing-pipeline Toronto_GPZN_only \
  --focal-tfs FLI1 GABPA GCM1 RFX5 \
  --sensitivity-tfs MAX \
  --availability-negative-controls PAX7 \
  --partner-expression results/capselex/hek293_partner_expression.tsv \
  --dnase-directory data/reference/hek293_dnase \
  --dnase-resolved-manifest reports/hek293_dnase_resolved.json \
  --permutations 100 \
  --partial-r2-threshold 0.005 \
  --minimum-positive-focal-tfs 4 \
  --seed 20260809 \
  --output reports/capselex_primary_genomic_model_screening.json
```
