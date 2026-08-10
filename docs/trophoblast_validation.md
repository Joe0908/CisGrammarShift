# GCM1 trophoblast association replication

## Design

The historical 204,916-locus GCM1 universe was projected into human trophoblast data:

- `GSE244252`: two-clone GCM1 ChIP-seq in EVT and ST differentiation;
- `GSE244253`: H3K4me3 HiChIP interactions for promoter–distal linkage;
- `GSE244254`: WT and GCM1-knockout RNA-seq;
- GENCODE v48 basic annotation on GRCh38 for promoters.

External hTSC ChIP did not define a new peak-selected test set. Only the two outcome columns were replaced;
loci, sequences, CAP profiles, monomer scores and chromosome folds remained frozen. However, the frozen
universe had originally been constructed using Codebook ChIP, so this is an external-outcome association
replication, not a fully assay-independent locus-selection design. It will be rerun on GHT-only tiles.

## Historical occupancy association

| State | Full CAP partial R² | Binned-monomer partial R² | Empirical p |
|---|---:|---:|---:|
| EVT | 0.001532 | 0.001538 | 1/101 |
| ST | 0.001759 | 0.001493 | 1/101 |

Composite profiles retain most of the signal; spacing-only chromosome-block intervals cross zero.

The leading ETV1, GABPA and ELK3 CAP models are correlated ETS-like grammars. RNA expression argues against a
literal ETV1-partner interpretation in hTSCs: ETV1 is essentially absent, while GABPA, ELK3 and several
other ETS factors are expressed. The supported unit is a transferable GCM1–ETS-like sequence grammar.

## Functional test 1: published target lists

Dataset S5 defines GCM1-dependent target genes by combining nearby GCM1 ChIP and knockout expression.
Using a 100-kb outcome-blind nearest-locus design, adding frozen ETS grammar slightly worsens log loss,
AUPRC and AUROC in both lineages.

| State | Targets / genes | M0 log loss | M1 log loss | M0 AUPRC | M1 AUPRC |
|---|---:|---:|---:|---:|---:|
| EVT | 1,157 / 22,029 | 0.169814 | 0.169907 | 0.205245 | 0.204408 |
| ST | 1,846 / 22,029 | 0.227138 | 0.227380 | 0.314821 | 0.313469 |

## Functional test 2: HiChIP-linked targets

Promoters are TSS ±2 kb. Promoter-overlapping loci and state-matched significant FitHiChIP distal links are
kept as separate channels. Among eligible loci, external-ChIP strength selects the representative locus
without reading the target label.

| State | Loops | Targets / linked genes | M0 log loss | M1 log loss | M0 AUPRC | M1 AUPRC |
|---|---:|---:|---:|---:|---:|---:|
| EVT | 107,409 | 1,075 / 17,808 | 0.194598 | 0.194609 | 0.194101 | 0.194254 |
| ST | 68,685 | 1,700 / 17,274 | 0.265723 | 0.265788 | 0.303030 | 0.302787 |

The tiny EVT AUPRC increase is not accompanied by improved log loss or AUROC. All ST metrics worsen.

## Functional test 3: exact continuous DESeq2 effect

Dataset S4 provides integer expected counts for four WT and two GCM1-KO clones per state. The source filter
was reproduced exactly: 28,250 raw RSEM genes minus 2,662 genes shorter than 300 bp equals the deposited
25,588-gene matrix. Differential expression uses exact R 4.2.2 and DESeq2 1.36.0. The outcome is negative,
unshrunk KO-versus-WT log2 fold change, and no `padj` threshold is applied.

| State | Finite LFC genes | Linked genes | M0 MSE | M1 MSE | Partial R² | Positive coefficient folds |
|---|---:|---:|---:|---:|---:|---:|
| EVT | 18,089 | 14,705 | 1.977103 | 1.977575 | -0.000239 | 1/5 |
| ST | 18,410 | 14,631 | 3.243476 | 3.244969 | -0.000460 | 2/5 |

## Claim boundary

All three functional endpoints are negative. ChIP, HiChIP and knockout RNA were collected on different
differentiation days, so the null does not prove absence of every possible transcriptional effect.
Establishing causality requires a time-matched ETS-partner or composite-motif perturbation. Any occupancy
claim remains provisional until the GHT-only primary analysis is complete.

No alternative profile weighting, promoter window, contact threshold, expression filter or locus aggregator
will be selected after observing these endpoints.
