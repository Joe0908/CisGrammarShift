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

The focused trophoblast test downloads only four processed GCM1 BigWigs from `GSE244252` (B31 and CT27 in
EVT and ST) plus the 22.5-MiB `GSE244254_RAW.tar` archive used for the pre-outcome TGIF2 expression gate.
Together the four BigWigs are 1,181,654,675 bytes. Exact URLs, byte sizes and SHA-256 values are frozen in
`configs/trophoblast_gcm1_replication_manifest.json`; raw assay-wide archives are not required.

## CAP-SELEX supplementary assets

The official Nature article supplements are downloaded individually rather than reconstructed from figures:

- Supplementary Table 1 (`MOESM3`): DNA and protein sequence metadata;
- Supplementary Table 2 (`MOESM4`): directed interaction matrix;
- Supplementary Table 3 (`MOESM5`): 1,348 PWM models;
- Supplementary Table 7 (`MOESM9`): 3,914 oriented 6-mer count rows across gaps 0--28.

`configs/capselex_nature_supplement_manifest.json` freezes the publisher URLs and byte sizes, while
`reports/capselex_nature_supplement_resolved.json` records locally computed SHA-256 hashes. The parser
normalizes documented aliases such as `RFXDC2 -> RFX7`; after normalization, all 913 composite-positive
matrix pairs have a Table 3 PWM and all 1,336 spacing-positive pairs have Table 7 counts. The raw interaction
matrix contains 2,223 positive directed cells. That literal table-cell count is not silently substituted for
the paper's separately reported, deduplicated headline count.

Focal monomer PWMs come from the Codebook Motif Explorer top-1 set archived at Zenodo record `15667805`;
these models were ranked using held-out assay data. Partner controls use the same source when available,
then individual HT-SELEX PWMs in CAP Supplementary Table 3, then the frozen JASPAR 2024 CORE vertebrate
non-redundant release. The resulting representative-profile partner coverage is audited before sequence
scoring; no ChIP outcome is used to choose a monomer model.
The only uncovered representative partner is ZBTB20 in ZBTB20–FLI1. That CAP profile is excluded by the
frozen contract; it is not assigned an inferred or outcome-selected surrogate motif.

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

The GEO `GSE278858_BED_files.tar.gz` archive is a reproducible independently deposited snapshot. Its
SHA-256 is frozen in `configs/codebook_geo_v1_magix_panel.json`, and
`scripts/audit_codebook_magix.py` extracts only the six-panel files, checks each against the GEO workbook
MD5, validates the MAGIX schema, and reports locus counts.

The current official Codebook v2 `Peaks_MAGIX_McGill.tar.gz` archive was downloaded manually from the v2
data page on 2026-08-10. Its archive SHA-256 is
`e1f225936529f8f0b2b338fcdfd8913e653012da1380e9622f200cce1c269a9b`. Direct extraction from that tar,
without reusing an existing directory, established that the six focal v2 BED members are byte-identical to
the corresponding GEO files. The v2 tar differs in packaging and member names, but the focal numerical
content does not. This is recorded as an observed release-specific equivalence, not assumed from dates.
`scripts/audit_codebook_v2_magix_panel.py` verifies each member's SHA-256, schema, width, coordinate names,
duplicates and selected-locus count before primary use.

The manuscript primary analysis therefore requires the corrected Codebook v2 MAGIX files. It uses the
outcome-independent rule:

- Benjamini-Hochberg `fdr <= 0.05`;
- refined `coefficient.ar > 0`;
- autosomes only;
- fixed 200-bp genomic tiles whose centers do not depend on ChIP.

This primary universe estimates the incremental CAP grammar association among loci with supported intrinsic
GHT binding. The `fixed-genome` sensitivity is required to address GHT-negative/ChIP-positive targeting,
because a GHT-only universe cannot by construction estimate that class.

### Corrected-v2 recovery audit

MAGIX is now public under GPL-3.0 at `csglab/MAGIX`; the audit pins commit
`bc800e825d19686aa5b73f3090d1eb31dddabbd3` and the authors' Zenodo v1.0.1 DOI. The corrected GEO workbook
identifies the exact runs that contributed to each published merged MAGIX peak set. The formal paper's
Supplementary Table 3 independently supplies experiment approval and batch membership. Joining these tables
to the ENA `PRJEB76622` file report gives the following selective target-read download sizes:

| TF | Role | Runs | FASTQs | Compressed size (GiB) |
|---|---|---:|---:|---:|
| FLI1 | primary | 19 | 35 | 0.842 |
| GABPA | primary | 8 | 16 | 0.385 |
| GCM1 | primary | 12 | 24 | 0.589 |
| PAX7 | availability control | 8 | 16 | 0.280 |
| RFX5 | primary | 3 | 6 | 0.164 |
| MAX | sensitivity | 16 | 32 | 0.636 |

Thus, focal target FASTQs are not the storage bottleneck: the five-TF acquisition panel totals 2.260 GiB,
or 2.896 GiB with MAX. The four expression-evaluable primary TFs total 1.980 GiB. A from-FASTQ exact
reproduction is nevertheless not declared ready. The public production scripts fit
batch-aggregate covariates and a genome-wide library-size model over approximately 13 million 200-bp bins;
they request 100 GB RAM for that step. The nine batches containing the primary focal experiments contain
79.476 GiB of all GHT experiment FASTQs, or 26.870 GiB if restricted to approved experiments. The production
design matrices and file lists are not tracked in the MAGIX source repository, and the released methods do
not disambiguate these two aggregate membership rules. Focal-only refitting would change the published
model and remains forbidden. The acquired official processed v2 BED members are used instead, so this raw
reconstruction ambiguity is no longer a primary-analysis blocker.

`configs/codebook_ght_v2_rebuild.json` freezes the source versions, and
`reports/codebook_ght_v2_rebuild_audit.json` records all selected ERR accessions, per-file ENA URLs, MD5s,
byte sizes, and the unresolved exact-rebuild boundary. The processed-asset audit is separately recorded in
`reports/codebook_ght_v2_magix_focal_audit.json`.

## Primary whitelist

- corrected Codebook v2 TF/plasmid and assay metadata;
- Toronto or one consistently chosen independent peak pipeline;
- Toronto `GPZN` bigWigs for all focal replicates;
- revised GHT MAGIX peaks/scores;
- hg38 chromosome sizes from a versioned reference for the fixed-genome sensitivity;
- UCSC hg38.2bit and the checksummed `twoBitToFa` binary for exact 400-bp sequence contexts;
- HEK293 DNase-seq `GSM2902639` plus the UCSC hg38-to-hg19 chain for locus-specific accessibility queries;
- two wild-type HEK293 Salmon quantifications from `GSM3611199` and GENCODE v29 for the independent
  partner-expression audit;
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
