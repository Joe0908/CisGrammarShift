# Intrinsic specificity versus cooperative DNA grammar in genomic targeting

## Frozen question

Does an independently measured CAP cooperative grammar explain held-out variation in focal-TF ChIP signal
after accounting for intrinsic GHT binding, focal and partner monomer motifs, sequence composition,
accessibility, and genomic context?

The model comparison is nested:

\[
M_0 = f(\mathrm{GHT},\ \mathrm{monomers},\ \mathrm{GC},\ \mathrm{accessibility},\ \mathrm{context})
\]

\[
M_1 = M_0 + \mathrm{CAP\ grammar}.
\]

## Anti-circularity contract

- The universe is the union of fixed 200-bp hg38 tiles containing a ChIP summit or GHT peak midpoint.
- Exact sequence features use a fixed 400-bp context around each tile.
- ChIP intensity never chooses the exact sequence-window centre.
- Source membership may be audited but is not a predictive feature.
- Train/validation/test separation is by chromosome.
- CAP feature multiplicity correction and monomer residualization use training chromosomes only.
- `TOPs`, `CTOPs`, MOODS triple-optimized intersections, and triple-optimized peak intersections are
  permanently forbidden in the primary analysis because they already incorporate outcome information.
- Focal-TF selection and all success criteria were frozen before reading the primary nested result.

The versioned contract is represented by `AnalysisContract` and `configs/capselex_reanalysis.yaml`.

## Panel and coverage

Primary focal TFs are FLI1, GABPA, GCM1, PAX7 and RFX5. They span multiple structural families and have
adequate CAP grammar coverage. MAX lacks a primary-eligible CAP grammar and is retained as an explicit
CAP-null control.

The constructed universe contains 1,204,394 loci, of which 1,035,968 belong to the five-TF primary panel.

## Primary result

| Analysis | Held-out partial R² | Interval / null |
|---|---:|---|
| Full CAP branch | 0.02030 | chromosome-block 95%: 0.01924–0.02131; p=1/101 |
| No DNase proxy | 0.01696 | positive sensitivity |
| Composite only | 0.01999 | retains nearly all signal |
| Spacing only | 0.00060 | minor component |
| Representative profiles only | 0.00565 | passes the frozen 0.005 gate |
| 20×20 monomer-bin residualization | 0.01914 | 0.01812–0.02015; p=1/101 |

All five primary TFs show a positive increment. GABPA replicate-separated fits and omission of the
contaminated GCM1 replicate `THC_0621` preserve the conclusion.

## Interpretation

The result is not “CAP predicts ChIP” in the abstract. It is evidence that CAP-measured composite sequence
grammar contains incremental information about held-out genomic occupancy after strong measured controls.
It does not prove that a named partner is present or that cooperativity causes transcription.

The independent GCM1 validation and functional claim boundary are documented in
[`trophoblast_validation.md`](trophoblast_validation.md).
