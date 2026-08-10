# GCM1 trophoblast association replication

## Design

The historical 204,916-locus GCM1 universe was projected into human trophoblast data:

- `GSE244252`: two-clone GCM1 ChIP-seq in EVT and ST differentiation;
- `GSE244253`: H3K4me3 HiChIP interactions for promoter–distal linkage;
- `GSE244254`: WT and GCM1-knockout RNA-seq;
- GENCODE v48 basic annotation on GRCh38 for promoters.

The original external hTSC analysis used a Codebook-ChIP-informed historical universe. It is retained below
for provenance. The focused TGIF2–GCM1 test instead uses the same official Codebook v2 GCM1 GHT-only tiles as
the primary analysis; external hTSC ChIP affects neither locus selection nor feature construction.

## Frozen TGIF2–GCM1 external-context test

TGIF2–GCM1 was the leading exploratory GCM1 pair in Codebook HEK293 ChIP. Before reading the hTSC ChIP
outcomes, the candidate, composite mechanism, two states, two-clone agreement rule, partial-R² threshold
(`0.005`) and two-stage permutation rule were frozen. TGIF2 then passed the outcome-independent expression
gate in all six relevant RNA profiles:

| State | B31 TGIF2 TPM | CT27 TGIF2 TPM |
|---|---:|---:|
| EVT | 8.61 | 4.44 |
| ST day 2 | 3.50 | 4.03 |
| ST day 4 | 1.19 | 1.31 |

RNA detection establishes availability only; it does not establish TGIF2 protein abundance or co-occupancy.
The frozen model tested the composite CAP excess score beyond GHT score, focal and TGIF2 monomer scores,
GC and CpG content on 19,497 GCM1 GHT-only loci. Each state used the mean of clone-specific log1p exact
200-bp GCM1 ChIP signals as its primary outcome.

| State | Clone Spearman r | Partial R² | Chromosome-bootstrap 95% interval | Screening p | Positive in all folds | Passed |
|---|---:|---:|---:|---:|---:|---:|
| EVT | 0.680 | 0.002471 | 0.001009–0.003741 | 1/101 | Yes | No |
| ST | 0.692 | 0.003471 | 0.001611–0.005342 | 1/101 | Yes | No |

All four clone-specific median coefficients were positive. EVT clone partial R² values were 0.003618 (B31)
and 0.001140 (CT27); ST values were 0.003158 and 0.003038. Both state-level effects nevertheless fell below
the predeclared `0.005` threshold. The replication criterion therefore failed, and the 1,000-permutation final
test was not triggered. The directionality is consistent with a weak transferable sequence association, but
the current evidence does not support a TGIF2–GCM1 cooperative mechanism.

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
Establishing causality requires a time-matched partner or composite-motif perturbation. The official GHT-only
primary analysis and focused external-context test are now complete; neither supports a general CAP grammar
mechanism or GCM1-dependent transcriptional causality.

No alternative profile weighting, promoter window, contact threshold, expression filter or locus aggregator
will be selected after observing these endpoints.
