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

The checksum-aware `scripts/download_manifest.py` downloads a whitelisted subset. URLs and expected hashes
belong in a local manifest under `data/`, which is ignored by Git.

## Primary whitelist

- corrected Codebook v2 TF/plasmid and assay metadata;
- Toronto or one consistently chosen independent peak pipeline;
- Toronto `GPZN` bigWigs for all focal replicates;
- revised GHT MAGIX peaks/scores;
- raw CAP composite PWMs and spacing annotations;
- external hTSC files from the GEO series above.

## Permanent primary blacklist

The following integrate motif, ChIP or GHT calls and would create circularity if used as primary predictors
or outcomes:

- TOPs and CTOPs;
- MOODS triple-optimized intersections;
- ChIP/GHT triple-optimized or triple-overlap peak sets;
- any peak or feature selected after reading the held-out outcome.

They may be described as external descriptive resources only, never used to construct the primary universe
or CAP branch.

## Files in Git

Git contains code, tests, configs, data contracts and small result summaries. It does not contain raw FASTQ,
bigWig, BEDPE, large peak archives, count matrices, derived per-locus tables or model checkpoints. Generated
outputs record SHA-256 hashes so the analysis can be reproduced without publishing controlled or bulky data.
