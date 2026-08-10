# Data access, storage and circularity rules

## Primary sources

| Assay | Public record | Use |
|---|---|---|
| CAP-SELEX | ENA `PRJEB66722` | Directed pair labels and grammar profiles |
| Codebook ChIP-seq | ENA `PRJEB78913`; GEO `GSE280248` | Focal-TF occupancy |
| GHT-SELEX | ENA `PRJEB76622`; Codebook v2 | Intrinsic binding |
| hTSC GCM1 ChIP | GEO `GSE244252` | Independent occupancy |
| hTSC H3K4me3 HiChIP | GEO `GSE244253` | Promoter–distal linkage |
| hTSC knockout RNA | GEO `GSE244254` | Functional outcome |

Codebook v2 metadata files are the canonical source for corrected labels. All joins use explicit normalized
TF symbols and retain original identifiers for audit.

## Avoiding the 175.39-GB ChIP archive

The monolithic Codebook merged-bigWig archive is unnecessary. GEO `GSE280248` provides Toronto-processed
`GPZN` bigWigs per sample, allowing only the required replicates to be downloaded. All focal samples must use
the same processing pipeline; Toronto and McGill absolute signals are never mixed in one comparison.

For the frozen six-TF panel, the 12 selected GPZN files total 4,616,661,075 bytes (4.30 GiB). Their GSM
accessions and GEO URLs are derived from the series SOFT rather than manually copied. The resolved manifest
records a locally computed SHA-256 for every downloaded bigWig. Both biological replicates are retained:
the mean of replicate `log1p` signals is the primary outcome, while replicate-resolved models must agree in
incremental-effect direction.

The checksum-aware `scripts/download_manifest.py` downloads a whitelisted subset. URLs, SHA-256 hashes and
expected byte sizes are frozen in `configs/codebook_geo_metadata_manifest.json`. Large assay files remain
under `data/`, which is ignored by Git.

## GHT MAGIX version boundary

The GEO `GSE278858_BED_files.tar.gz` archive is a reproducible pre-v2 snapshot. Its SHA-256 is frozen in
`configs/codebook_geo_v1_magix_panel.json`, and `scripts/audit_codebook_magix.py` extracts only the six-panel
files, checks each against the GEO workbook MD5, validates the MAGIX schema, and reports locus counts. It is
explicitly **smoke-test only** because Codebook v2 changed MAGIX preprocessing from fixed genomic bins to
peak-first candidate regions and consequently changed locations and scores.

The manuscript primary analysis therefore requires the corrected Codebook v2 MAGIX files. It uses the
outcome-independent rule:

- Benjamini-Hochberg `fdr <= 0.05`;
- refined `coefficient.ar > 0`;
- autosomes only;
- fixed 200-bp genomic tiles whose centers do not depend on ChIP.

This primary universe estimates the incremental CAP grammar association among loci with supported intrinsic
GHT binding. The `fixed-genome` sensitivity is required to address GHT-negative/ChIP-positive targeting,
because a GHT-only universe cannot by construction estimate that class.

## Primary whitelist

- corrected Codebook v2 TF/plasmid and assay metadata;
- Toronto or one consistently chosen independent peak pipeline;
- Toronto `GPZN` bigWigs for all focal replicates;
- revised GHT MAGIX peaks/scores;
- hg38 chromosome sizes from a versioned reference for the fixed-genome sensitivity;
- raw CAP composite PWMs and spacing annotations;
- external hTSC files from the GEO series above.

## Permanent primary blacklist

The following integrate motif, ChIP or GHT calls and would create circularity if used as primary predictors
or outcomes:

- TOPs and CTOPs;
- MOODS triple-optimized intersections;
- ChIP/GHT triple-optimized or triple-overlap peak sets;
- any peak or feature selected after reading the held-out outcome.

ChIP-derived loci are also forbidden in the primary universe. The `assay-union` implementation is retained
only as an explicitly outcome-informed sensitivity, while `legacy-assay-union` reproduces the former
ChIP-prioritized centering behavior for audit purposes.

They may be described as external descriptive resources only, never used to construct the primary universe
or CAP branch.

## Files in Git

Git contains code, tests, configs, data contracts and small result summaries. It does not contain raw FASTQ,
bigWig, BEDPE, large peak archives, count matrices, derived per-locus tables or model checkpoints. Generated
outputs record SHA-256 hashes so the analysis can be reproduced without publishing controlled or bulky data.
