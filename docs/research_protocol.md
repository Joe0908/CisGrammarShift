# Research protocol

## Objective

Determine whether neural sequence models recover an explicitly programmed cis-regulatory spacing rule under
counterfactual matching and distribution shift.

The benchmark is designed to separate three capabilities:

1. **motif recognition** — detecting POU5F1 and NANOG instances;
2. **motif integration** — combining evidence from two sites;
3. **grammar recovery** — predicting a label determined only by their relative spacing.

## Pre-specified hypotheses

- **H1:** A PWM-presence baseline will be near chance on matched-pair ranking because both members contain the
  same motif instances.
- **H2:** A local CNN will recognise motifs but underperform models with a sequence-wide receptive field on the
  periodic syntax task.
- **H3:** IID discrimination will overestimate grammar recovery; gap, orientation, GC, or motif-strength shifts
  will reveal architecture-specific failures.
- **H4:** Models that genuinely use the implanted motifs will assign more gradient-times-input mass to the
  ground-truth motif positions than expected from the fraction of sequence occupied by those positions.

These hypotheses are evaluated independently for each seed. They must not be rewritten after observing the
test results.

## Data-generating process

For every pair:

1. Sample one Markov background sequence with pre-specified GC content and persistence.
2. Sample one instance from each TF probability matrix.
3. Sample one orientation for each instance.
4. Construct a positive sequence with a gap satisfying the periodic rule.
5. Construct a negative sequence from the same background and instances with a gap violating the rule.

The rule is:

```text
positive(gap) = (gap mod period) in allowed_phases
```

The default period is 10 bp and the default allowed phases are 0, 1, and 9. This is a controlled mathematical
rule, not an estimate of an endogenous biochemical response.

Pair IDs are unique within a generated dataset. Splits are generated independently with different derived
seeds, so neither a pair nor its background is shared across splits.

## Primary outcomes

1. AUROC on each evaluation condition.
2. Matched-pair accuracy: fraction of pairs for which the positive receives a higher probability than its
   counterfactual negative.

The primary comparison is the mean across five pre-specified seeds. Per-seed values remain available and
should be plotted rather than hidden behind an aggregate.

## Secondary outcomes

- AUPRC
- balanced accuracy at a validation-selected threshold
- Matthews correlation coefficient
- Brier score
- expected calibration error (10 equal-width bins)
- mean positive-minus-negative probability within matched pairs
- gradient-times-input mass fraction inside implanted motif positions
- parameter count and best validation loss

## Model selection

Training uses binary cross-entropy with logits and AdamW. The state with the lowest validation loss is
retained. Early stopping uses the validation set only. A decision threshold is chosen on validation
predictions by maximising balanced accuracy and then frozen for every test condition.

No test condition may be used for early stopping, hyperparameter selection, or threshold selection.

## Robustness conditions

- `iid`: same generator settings as training, independent seed.
- `gap_ood`: non-overlapping gap range while preserving the same periodic label rule.
- `orientation_ood`: reverse-complement instances introduced after forward-only training.
- `gc_ood`: altered background nucleotide composition.
- `strength_ood`: higher matrix sampling temperature, producing less consensus-like instances.

Each condition changes one factor where possible. The condition definitions are versioned in the experiment
configuration.

## Interpretation analysis

Attribution is computed for a fixed, seeded subset of examples; examples are not selected for being correctly
classified. The reported localisation score is:

```text
sum(|gradient × input| within implanted motif positions)
---------------------------------------------------------
sum(|gradient × input| across the full sequence)
```

The sequence-coverage baseline is the fraction of positions occupied by the two motifs. Attribution
localisation is supportive diagnostic evidence, not proof of mechanism.

## Exclusion and failure rules

- A run with non-finite loss is failed, not silently omitted.
- A generated dataset must be exactly balanced and contain exactly two rows per pair.
- Positive gaps must satisfy the rule and negative gaps must violate it.
- Any overlap between motif intervals is rejected.
- If only one class appears in an evaluation set, the run is invalid.
- Missing provenance or synthetic fallback labels are not accepted as real-data validation.

## Reporting

Report all configured models and conditions, including negative results. Do not describe synthetic benchmark
performance as TF-binding accuracy. Any real-data extension must publish accession IDs, genome build,
processing commands, chromosome partitions, and negative-sampling rules.
