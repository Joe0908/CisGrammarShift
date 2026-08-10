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

The acquisition panel is FLI1, GABPA, GCM1, PAX7, RFX5 and MAX. Before any ChIP model was run, two
wild-type HEK293 RNA-seq replicates (GSM3611199) were used as an outcome-independent partner-availability
audit. Requiring summed gene TPM >= 1 in each replicate left 16/31 FLI1, 4/5 GABPA, 4/7 GCM1 and 3/5 RFX5
representative pairs. These four TFs form the expression-evaluable primary panel. PAX7–TBX4 is the only
representative PAX7 pair, but TBX4 is 0 TPM in both replicates; PAX7 is retained as a predeclared partner-
availability negative control. MAX–TEAD4 is retained outside the primary panel as a low-coverage
sensitivity and is not a CAP-null control.

The official Nature supplements are now checksummed and parsed directly. Across the six-TF acquisition
panel, the audit found 171, 45, 25, 6, 12 and 15 positive directed pairs for FLI1, GABPA, GCM1, MAX, PAX7
and RFX5, respectively. Every focal TF has raw CAP feature assets. This is a coverage statement only;
partner expression, monomer controls and the frozen feature-selection rule determine model eligibility.
All six focal monomer controls are available from the held-out-ranked Codebook MEX top-1 archive. Frozen
partner controls cover every representative pair except ZBTB20–FLI1; that single profile is excluded before
outcome analysis rather than scored without its matched monomer control.

The reported historical universe contains 1,204,394 loci, of which 1,035,968 belong to the former five-TF panel.
The official Codebook v2 GHT-only call produces 50,250 FLI1, 5,057 GABPA, 19,497 GCM1, 18,485 PAX7 and
38,894 RFX5 autosomal selected intervals; MAX contributes 71,564 sensitivity intervals. These counts use
the frozen `fdr <= 0.05` and positive refined-coefficient rule.

The raw-data fallback has been narrowed without weakening the contract. Author-selected focal target
FASTQs total 2.260 GiB for the five-TF acquisition panel (2.896 GiB including MAX), but an exact corrected-v2
MAGIX rebuild still needs production design/batch-aggregate inputs absent from the accessible public source
snapshot. Supplementary Table 3 resolves the nine focal batch identifiers, but
the production file lists do not state whether aggregates contain all batch experiments (79.476 GiB) or only
approved experiments (26.870 GiB). The public production scripts require a roughly 13-million-bin
library-size fit and request 100 GB RAM. A focal-only substitute would be a different model and is not
accepted as primary evidence. This raw reconstruction is no longer required for the primary analysis:
the official Codebook v2 processed archive has been acquired, and its six focal members pass the frozen
integrity and eligibility audit.

Primary sequence features directly scan author-marked representative CAP PWMs on fixed 400-bp hg38
contexts. Multiple representative models within one pair/mechanism are combined by maximum score. Each raw
CAP score is residualized against focal and partner monomer scores using training-chromosome-derived 20x20
quantile bins. Mechanisms are averaged within a pair before equal-weight averaging across expression-
supported pairs, preventing a pair annotated with both composite and spacing PWMs from receiving double
weight. All-pair aggregates are a required sensitivity.

Accessibility comes from independent HEK293 DNase-seq (GSM2902639). Because the submitted BigWig is hg19,
only included hg38 loci are projected through the frozen UCSC hg38-to-hg19 chain and queried against the
source BigWig with exact interval means. Mapped and unmapped loci are recorded. The model is repeated without
accessibility on both the mapped set and the complete GHT-only set; the DNase track is a cell-line proxy, not
matched chromatin from each Codebook ChIP experiment.

## Official Codebook v2 primary screening

The outcome-independent rerun is complete. A focal TF passed only when its mean-outcome partial R² was at
least 0.005, its one-sided chromosome-shift p-value was at most 0.05, and median standardized CAP
coefficients were positive for the mean outcome and both ChIP replicates.

| TF | Supported pairs | Mapped loci | Partial R² | Chromosome-bootstrap 95% interval | Replicate partial R² | Passed |
|---|---:|---:|---:|---:|---:|---:|
| FLI1 | 16 | 50,060 | 0.010679 | 0.009111–0.012231 | 0.005944 / 0.013159 | Yes |
| GABPA | 4 | 5,041 | 0.002881 | 0.000354–0.005568 | 0.004630 / 0.000462 | No |
| GCM1 | 4 | 19,369 | 0.003826 | 0.001804–0.005478 | 0.002603 / 0.004388 | No |
| RFX5 | 3 | 38,701 | 0.001190 | 0.000442–0.002150 | 0.001393 / 0.000468 | No |

All four one-sided screening tests had `p=1/101`, but only one TF exceeded the predeclared effect-size
threshold. The required four-of-four panel gate therefore failed. The 1,000-permutation final panel test was
not run.

The negative-control and sensitivity cases do not rescue the general claim. PAX7 had no expression-supported
partner and is not interpretable as a primary cooperative test. MAX–TEAD4 produced partial R² 0.002550, below
the effect threshold.

The legacy assay-union estimate of 0.02030 remains available in `reports/result_summary.json` only as
provenance. It will not be presented as primary evidence because ChIP affected the historical locus universe.

## Outcome-informed pair decomposition

After the panel gate failed, each expression-supported pair was evaluated separately for candidate
generation. This analysis is exploratory because the same Codebook ChIP outcomes rank the pairs.

- GATA3–FLI1 composite: partial R² 0.010506, chromosome-bootstrap interval 0.008736–0.012463.
- TGIF2–GCM1 composite: partial R² 0.012558, interval 0.010080–0.015243.
- GCM1–HOXA2 spacing: partial R² 0.005830, interval 0.004038–0.007752.

The TGIF2–GCM1 profile was taken forward only after a separate hTSC RNA audit showed TGIF2 TPM >=1 in all
predeclared EVT/ST clone profiles. The external ChIP result is documented in
[`trophoblast_validation.md`](trophoblast_validation.md).

## Interpretation

The data support heterogeneous, TF-specific occupancy associations rather than a general CAP grammar effect.
FLI1 passed the frozen TF-level criterion; GABPA, GCM1 and RFX5 did not. Pair-level rankings cannot establish
that a named partner is present, and no analysis here demonstrates cooperative causality.

The independent GCM1 validation and functional claim boundary are documented in
[`trophoblast_validation.md`](trophoblast_validation.md).
