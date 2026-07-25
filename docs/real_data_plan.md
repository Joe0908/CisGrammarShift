# Real-data validation plan

The synthetic benchmark is valuable because the true motif positions and rule are known. It is insufficient
for a biological conclusion. The next stage is a locked, accessioned analysis of pluripotency-factor binding.

## Biological question

Do models that recover periodic syntax in the counterfactual benchmark also detect reproducible,
orientation-aware motif-spacing effects in experimental pluripotency-factor binding data?

## Candidate source

The BPNet study provides a strong validation target because it used base-resolution ChIP-nexus profiles for
Oct4, Sox2, Nanog and Klf4, evaluated on held-out chromosomes, identified Nanog-associated helical
periodicity, and included CRISPR-mutated loci. The study and its supplementary tables should be treated as the
primary source for exact experiment identifiers:

- Avsec Ž, et al. *Nature Genetics* 53, 354–366 (2021).
- DOI: <https://doi.org/10.1038/s41588-021-00782-6>

## Locked split

To permit direct comparison with the published work:

- training chromosomes: all eligible chromosomes except validation and test;
- validation chromosomes: 2, 3, and 4;
- test chromosomes: 1, 8, and 9.

The final implementation must verify this split against the primary paper and record genome assembly,
blacklists, peak unions, and every accession before execution.

## Analysis stages

1. Download accessioned signal and control files; record checksums.
2. Reproduce peak and profile preprocessing with a container or locked environment.
3. Train a sequence-to-profile model and a binary-classification control.
4. Call motif instances on held-out chromosomes without using test labels for thresholds.
5. Quantify motif-pair distance and orientation effects with bootstrap intervals over genomic regions.
6. Test whether the synthetic-benchmark ranking predicts real-data grammar recovery.
7. Where available, evaluate predicted direction and magnitude at CRISPR-mutated loci.

## Required controls

- dinucleotide- or locus-matched negatives;
- mappability and blacklist filtering;
- chromosome-level separation;
- assay-control or bias track;
- replicate-aware evaluation;
- TF-identity and motif-family ambiguity analysis;
- sensitivity to motif caller and score threshold;
- comparison with a PWM scanner and a limited-receptive-field network.

## Go/no-go criteria

Do not claim endogenous syntax unless:

- the effect replicates on held-out chromosomes;
- the direction is consistent across biological replicates;
- the effect remains after GC, accessibility, and motif-strength matching;
- the spacing/orientation analysis was pre-specified;
- at least one perturbational or orthogonal validation supports the inferred interaction.

## Status

This repository does not currently include experimental files or real-data results. That omission is
intentional: synthetic fallback sequences must never be silently relabelled as biological observations.
