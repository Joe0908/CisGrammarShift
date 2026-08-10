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

- The primary universe contains fixed 200-bp hg38 tiles selected by GHT peak midpoints only.
- Every primary window is centred on the fixed tile, never on a ChIP or GHT summit.
- ChIP can annotate a primary tile and supplies the outcome, but cannot determine primary inclusion or centre.
- Exact sequence features use a fixed 400-bp context around each tile.
- Source membership may be audited but is not a predictive feature.
- Train/validation/test separation is by chromosome.
- CAP feature multiplicity correction and monomer residualization use training chromosomes only.
- `TOPs`, `CTOPs`, MOODS triple-optimized intersections, and triple-optimized peak intersections are
  permanently forbidden in the primary analysis because they already incorporate outcome information.
- Focal-TF selection and all success criteria were frozen before reading the primary nested result.

The versioned contract is represented by `AnalysisContract` and `configs/capselex_reanalysis.yaml`.

Three predeclared sensitivity universes test selection dependence:

| Universe | ChIP selects loci? | ChIP centres sequence? | Role |
|---|---:|---:|---|
| `ght-only` | No | No | Primary |
| `assay-union` | Yes | No | Outcome-informed selection sensitivity |
| `fixed-genome` | No | No | Fully selection-independent sensitivity |
| `legacy-assay-union` | Yes | Yes | Historical result reproduction only |

## Panel and coverage

Primary focal TFs are FLI1, GABPA, GCM1, PAX7 and RFX5. They span multiple structural families and have
adequate CAP grammar coverage. MAX is retained outside the primary panel as a predeclared low-coverage
sensitivity analysis. It is not a CAP-null control: the official Supplementary Table 3 contains one
representative MAX–TEAD4 PWM, and Supplementary Table 7 contains five MAX spacing pairs.

The official Nature supplements are now checksummed and parsed directly. Across the six-TF acquisition
panel, the audit found 171, 45, 25, 6, 12 and 15 positive directed pairs for FLI1, GABPA, GCM1, MAX, PAX7
and RFX5, respectively. Every focal TF has raw CAP feature assets. This is a coverage statement only;
partner expression, monomer controls and the frozen feature-selection rule determine model eligibility.

The reported historical universe contains 1,204,394 loci, of which 1,035,968 belong to the five-TF panel.
Primary GHT-only counts have not yet been calculated.

## Historical result requiring primary rerun

| Analysis | Held-out partial R² | Interval / null |
|---|---:|---|
| Full CAP branch | 0.02030 | chromosome-block 95%: 0.01924–0.02131; p=1/101 |
| No DNase proxy | 0.01696 | positive sensitivity |
| Composite only | 0.01999 | retains nearly all signal |
| Spacing only | 0.00060 | minor component |
| Representative profiles only | 0.00565 | passes the frozen 0.005 gate |
| 20×20 monomer-bin residualization | 0.01914 | 0.01812–0.02015; p=1/101 |

These estimates came from `legacy-assay-union`; they are retained for provenance and will not be presented
as primary manuscript evidence. GABPA replicate-separated fits and omission of the contaminated GCM1
replicate `THC_0621` were internally stable, but they do not repair outcome-informed locus construction.

## Interpretation

If the GHT-only and fixed-genome reruns retain the increment, the result would support the narrower statement
that CAP-measured composite sequence grammar contains incremental information about held-out occupancy after
measured controls. The historical analysis alone does not establish that claim. No analysis here proves that
a named partner is present or that cooperativity causes transcription.

The independent GCM1 validation and functional claim boundary are documented in
[`trophoblast_validation.md`](trophoblast_validation.md).
