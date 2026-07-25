# Synthetic data specification

## Unit of analysis

The independent unit is a counterfactual pair, not an individual sequence. Each pair contains:

- one positive sequence whose motif gap satisfies the rule;
- one negative sequence whose motif gap violates the rule.

Both share a background sequence, POU5F1 instance, NANOG instance, and motif orientations. Only the second
motif's position changes.

## Output tensors

| Field | Shape | Meaning |
|---|---:|---|
| `x` | `(2 × pairs, length, 4)` | A/C/G/T one-hot sequence |
| `y` | `(2 × pairs,)` | Binary grammar label |
| `pair_ids` | `(2 × pairs,)` | Matched-pair identifier |
| `motif_mask` | `(2 × pairs, length)` | Ground-truth implanted positions |
| `records` | one record per row | Gap, phase, positions, orientations, instances, condition |

## Background model

The first base is sampled from the requested GC composition. Later bases either repeat the preceding base
with probability `background_persistence` or are sampled from the GC-controlled base probabilities. This is a
minimal first-order background model, not a substitute for genomic sequence.

## Motif sampling temperature

Each probability row is transformed as:

```text
p_temperature(base) proportional to p(base) ** (1 / temperature)
```

`temperature=1` preserves the supplied matrix. Larger values flatten it and generate weaker,
less-consensus-like instances. The same transformed matrix is used for both members of a pair.

## Provenance

The bundled POU5F1 and NANOG matrices were supplied with the extended project and labelled as CIS-BP 3.00:

- POU5F1: `M05705_3.00`
- NANOG: `M05219_3.00`

The exact numeric matrices and source metadata live in `src/cisgrammar/motifs.py`. A future release should
pin a downloadable motif archive checksum so the matrices can be independently verified.

## Not represented

The synthetic generator does not model chromatin accessibility, nucleosomes, DNA shape, cofactor abundance,
cell state, three-dimensional genome organisation, transcriptional output, or experimental assay bias.
Consequently, the benchmark measures recovery of programmed sequence rules only.
