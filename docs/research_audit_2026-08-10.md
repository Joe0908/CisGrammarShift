# Research audit and publication decision — 2026-08-10

## Research question

Does CAP-SELEX-derived cooperative DNA grammar explain held-out focal-TF genomic occupancy beyond intrinsic
GHT-SELEX binding, focal and available-partner monomer motifs, sequence composition and an independent
accessibility proxy, when the tested loci are selected without using ChIP outcomes?

The question is biologically meaningful and answerable with the assembled data. The current results answer
the broad version negatively. They do not identify an in-cell partner or establish cooperative causality.

## Central data-supported claim

CAP-derived sequence features show small and heterogeneous occupancy increments on outcome-independent
GHT-only loci: FLI1 passed the frozen TF-level criterion, whereas GABPA, GCM1 and RFX5 did not; a focused
TGIF2–GCM1 score was directionally consistent in two trophoblast states but failed its predeclared effect
threshold in both.

This is the strongest defensible claim. “CAP grammar generally determines genomic targeting” and
“TGIF2 cooperates with GCM1 in trophoblast” are not supported.

## Novelty assessment

The potentially useful contribution is methodological and empirical rather than a new confirmed biological
mechanism:

1. It separates locus selection from the ChIP outcome and exposes inflation in an outcome-informed historical
   universe.
2. It controls focal and partner monomer sequence preferences before attributing signal to a composite CAP
   profile.
3. It requires partner availability, replicate direction agreement, chromosome-held-out prediction,
   chromosome-block intervals and spatial nulls.
4. It reports a predeclared negative panel result and a failed focused replication without changing the gate.

A positive mechanistic paper would require direct partner occupancy or perturbation. A methods/resource paper
would require a broader benchmark and comparison with alternative grammar representations.

## Data audit

### Critical problems

| Problem | Consequence | What would resolve it |
|---|---|---|
| Three of four primary TFs failed the frozen effect gate | The general biological hypothesis is rejected | Do not rescue it by retuning; formulate a new, independently frozen question |
| TGIF2–GCM1 was selected using Codebook ChIP and failed the external effect gate | It is not a replicated named-partner mechanism | Partner ChIP/CUT&RUN or TGIF2/composite-site perturbation in the same trophoblast state |
| CAP PWM occurrence is not evidence of simultaneous protein occupancy | Sequence association cannot establish cooperation | Same-cell co-occupancy or perturbational assay |
| Existing transcriptional endpoints are negative | Occupancy associations do not extend to the claimed function | Time-matched motif or partner perturbation with expression readout |

### Important problems

| Problem | Consequence | Required treatment |
|---|---|---|
| Only four expression-evaluable primary TFs | Limited generality and family coverage | Treat as focused feasibility, or expand with matched expression/occupancy contexts |
| HEK293 DNase comes from an independent experiment and hg19 | Imperfect chromatin control | Report mapping and no-accessibility sensitivities; obtain matched hg38 ATAC/DNase if available |
| Partner expression is RNA from a separate wild-type HEK293 dataset | TPM does not guarantee protein availability in each ChIP experiment | Keep as an availability filter, not proof of presence |
| Pair decomposition is outcome-informed and tests many profiles | Pair ranks are exploratory and subject to multiplicity | Use only for prospective candidate freezing; do not assign discovery p-values |
| Only 100 screening permutations were run | Minimum attainable p is 1/101 | Sufficient for the stop/go rule; do not report as high-resolution significance |
| Chromosome bootstrap has few biological blocks | Intervals reflect genomic block uncertainty, not biological replication | Retain replicate-resolved results and avoid treating loci as biological replicates |
| Fixed-genome sensitivity has not been run for the final feature implementation | Remaining locus-selection dependence within GHT-selected regions | Needed only for a methods manuscript; it cannot reverse the failed frozen panel gate |

### Minor problems

- The CAP, GHT, ChIP, DNase and RNA assays are not all measured in one matched biological sample.
- PWM maximum scores compress motif multiplicity and local motif arrangement.
- The expression threshold of 1 TPM is a pragmatic availability rule, not a biological activity threshold.
- The GHT-only universe tests occupancy variation within intrinsic-binding-selected regions, not genome-wide
  recruitment.

## Statistical audit

The primary analysis uses five chromosome folds with nested ridge selection. The effect is out-of-fold
partial R² from adding one CAP aggregate to a frozen baseline. Uncertainty is estimated by resampling
autosomal chromosomes, and spatial specificity by circularly shifting the CAP feature within chromosome.
Both ChIP replicates are modeled separately. These choices address leakage and local sequence correlation
better than random-locus splits.

The p-values are screening quantities. With large locus counts, all four primary TFs had `p=1/101`, including
three effects below threshold. The effect-size gate is therefore the decision criterion. No post hoc threshold
change, feature reweighting or 1,000-permutation test is justified after a failed screen.

## Figure-driven paper decision

A positive biological-mechanism manuscript is not supported by the present evidence. If this work is developed
as a rigorous reanalysis/methods paper, the minimum coherent figure set would be:

1. Outcome-independent study design, data intersections and leakage controls.
2. Official v2 per-TF partial R² with chromosome intervals and replicate-resolved effects.
3. Comparison of historical outcome-informed and GHT-only universes, including fixed-genome sensitivity.
4. Exploratory pair decomposition with an explicit discovery/validation boundary.
5. Frozen TGIF2–GCM1 external-context result and negative functional endpoints.

The paper's value would be the falsification framework and the demonstrated heterogeneity, not a claimed
TGIF2–GCM1 mechanism.

## Strongest rejection argument

The study integrates impressive public resources but does not yet deliver a new biological mechanism: the
broad predeclared claim failed, the named-pair candidate was chosen using an outcome and missed its external
effect threshold, and no partner occupancy or perturbation is available.


**Current status: B — potentially publishable after major additional analyses.**

The code and data provenance are strong enough for a research-experience portfolio now. Submission should
wait for one of two non-overlapping upgrades:

- **Biological route:** direct same-context partner occupancy plus motif/partner perturbation; or
- **Methods route:** a broader cross-TF benchmark, fixed-genome sensitivity, calibrated comparison with
  alternative monomer/composite models, and publication-quality figures.

Until one route is completed, drafting a full positive-result manuscript would create unsupported claims.
