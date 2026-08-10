# CisGrammarShift

**Intrinsic TF specificity, cooperative DNA grammar, and genomic targeting across experimental contexts.**

[![CI](https://github.com/Joe0908/CisGrammarShift/actions/workflows/ci.yml/badge.svg)](https://github.com/Joe0908/CisGrammarShift/actions/workflows/ci.yml)
[![Python 3.10+](https://img.shields.io/badge/python-3.10%2B-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)

## Biological question

TF genomic occupancy is often predicted from a focal TF motif and chromatin context. CAP-SELEX nevertheless
reports thousands of directed, DNA-guided TF-pair interactions. This project asks a narrower mechanistic
question:

> Does independently measured cooperative DNA grammar explain held-out variation in focal-TF genomic
> occupancy beyond intrinsic GHT-SELEX binding, monomer motifs, sequence context, and accessibility—and
> does that increment extend to GCM1-dependent transcription?

The repository keeps the original counterfactual synthetic benchmark as a positive-control layer, but the
primary study is now an accessioned CAP-SELEX × GHT-SELEX × ChIP-seq analysis with independent human
trophoblast validation.

## Result in one paragraph

CAP composite grammar contributes a small, reproducible increment to held-out ChIP intensity. In the
five-TF Codebook panel the macro chromosome-held-out partial R² is `0.02030` (chromosome-block 95% interval
`0.01924–0.02131`; empirical `p = 1/101`); all five focal TFs are positive. A stricter training-only 20×20
focal/partner monomer-bin control retains `0.01914`. Independent GCM1 ChIP-seq in hTSC differentiation
replicates a smaller increment in EVT (`0.00153`) and ST (`0.00176`), each at `p = 1/101`. The transferable
signal is a GCM1–ETS-like composite grammar class, not a uniquely identified ETV1 partner. Two binary target
benchmarks and the exact R 4.2.2 / DESeq2 1.36.0 continuous-effect analysis are negative, so downstream
transcriptional causality is not claimed.

## Study layers

| Layer | Role | Leakage control | Endpoint |
|---|---|---|---|
| Synthetic counterfactual | Architecture/implementation positive control | Matched motif instances and pair-level splits | Programmed spacing rule |
| CAP pair feasibility | Diagnose coverage and shortcut learning | Pair-, TF-, family- and paralogue-held-out splits | Directed CAP labels |
| Codebook primary panel | Test endogenous occupancy increment | Fixed assay-union loci and chromosome-held-out fitting | Continuous focal-TF ChIP |
| Independent hTSC | Test cross-context reproducibility | Same frozen GCM1 loci/features; external outcomes only | EVT/ST GCM1 ChIP |
| HiChIP + knockout RNA | Test functional extension | Outcome-blind linkage and nested chromosome folds | Binary targets and DESeq2 LFC |

## Primary results

### Codebook occupancy panel

- Primary focal TFs: `FLI1`, `GABPA`, `GCM1`, `PAX7`, and `RFX5`; `MAX` is retained as a CAP-null control.
- Locus universe: 1,204,394 fixed 200-bp hg38 assay-union loci; 1,035,968 loci in the primary panel.
- M0: continuous GHT, focal/partner monomer scores, GC, accessibility proxy, and genomic context.
- M1: M0 plus CAP composite and spacing/orientation grammar features calibrated on training chromosomes.
- Full CAP partial R²: `0.02030`; without the DNase proxy: `0.01696`.
- Composite-only: `0.01999`; spacing-only: `0.00060`.
- Representative CAP profiles only: `0.00565`.
- Training-only binned-monomer residualization: `0.01914` (`0.01812–0.02015`, `p = 1/101`).

### Independent GCM1 trophoblast occupancy

The frozen 204,916-locus GCM1 universe was projected into `GSE244252` without allowing the external ChIP
signal to choose loci or tune features.

| State | Full CAP partial R² | Binned-monomer sensitivity | Empirical p |
|---|---:|---:|---:|
| EVT | 0.001532 | 0.001538 | 1/101 |
| ST | 0.001759 | 0.001493 | 1/101 |

ETV1 is nearly absent in the relevant hTSC states, while GABPA, ELK3 and other ETS factors are expressed.
The defensible biological unit is therefore an ETS-family sequence grammar, not a literal GCM1–ETV1
protein complex.

### Functional follow-up

| Endpoint | EVT | ST | Interpretation |
|---|---:|---:|---|
| 100-kb nearest-locus binary target benchmark | all metrics slightly worse | all metrics slightly worse | Negative |
| H3K4me3 HiChIP-linked binary benchmark | no coherent gain | all metrics worse | Negative |
| Continuous DESeq2 partial R² | -0.000239 | -0.000460 | Negative |
| Positive grammar coefficients across 5 folds | 1/5 | 2/5 | Direction unstable |

These three tests jointly reject the current claim that the frozen GCM1–ETS grammar independently predicts
GCM1-dependent transcription. Differing ChIP, HiChIP, and RNA time points remain a limitation; a stronger
claim now requires a time-matched motif or partner perturbation rather than additional feature reweighting.

## Reproducible workflow

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -e ".[dev]"

# Freeze the non-circular contract
cisgrammar capselex contract --output results/capselex/analysis_contract.json

# Build the directed CAP pair dataset and run the TF-held-out panel
cisgrammar capselex dataset \
  --interaction-table data/capselex/interaction_matrix.xlsx \
  --construct-table data/capselex/constructs.xlsx \
  --output results/capselex/tf_pair_dataset.tsv

cisgrammar capselex tf-panel \
  --dataset results/capselex/tf_pair_dataset.tsv \
  --target composite_motif \
  --output results/capselex/tf_panel_composite.json

# Build outcome-independent Codebook ChIP/GHT assay-union loci
cisgrammar capselex build-loci \
  --chip-peaks data/codebook/chip/GCM1.narrowPeak \
  --ght-peaks data/codebook/ght/GCM1.narrowPeak \
  --focal-tf GCM1 \
  --output results/capselex/GCM1.loci.tsv

# Export exact count matrices, then run the locked R analysis
cisgrammar capselex trophoblast-deseq2-export \
  --rna-dataset data/trophoblast/pnas.2311372120.sd04.xlsx \
  --raw-rsem-filter-audit data/trophoblast/GSM7810026_V2_2_EVT.genes.results.gz \
  --output-directory results/capselex/deseq2_inputs

Rscript scripts/run_gcm1_deseq2.R \
  results/capselex/deseq2_inputs \
  results/capselex/deseq2_results

cisgrammar capselex trophoblast-deseq2-hichip-benchmark \
  --evt-gene-features results/capselex/EVT.linked_gene_features.tsv.gz \
  --st-gene-features results/capselex/ST.linked_gene_features.tsv.gz \
  --deseq2-directory results/capselex/deseq2_results \
  --output results/capselex/gcm1_continuous_benchmark.json \
  --predictions results/capselex/gcm1_continuous_predictions.tsv.gz

ruff check .
pytest
```

Raw and large processed files are intentionally excluded from Git. Every analysis output records input
paths, sizes, SHA-256 hashes, software versions, and the locked design. See
[`docs/data_access.md`](docs/data_access.md) for public accessions and the small-file download strategy.

## Repository layout

```text
configs/                 Frozen synthetic and CAP reanalysis configurations
docs/                    Feasibility, data, modelling, validation and claim-boundary documents
reports/                 Small versioned result summaries; no raw biological data
scripts/                 Checksum downloader and exact DESeq2 runner
src/cisgrammar/          Synthetic benchmark and complete CAP/GHT/ChIP/hTSC pipeline
tests/                   Unit, leakage-control and synthetic-recovery tests
legacy/original/         Preserved student-stage scripts and figures
```

## Claim boundary

Supported:

> Independently measured CAP composite grammar contains reproducible information about held-out focal-TF
> occupancy beyond intrinsic binding, monomer motifs, sequence context, and the available accessibility
> covariate; the independent GCM1 signal localizes to an ETS-like grammar class.

Not supported:

- a uniquely identified in-cell ETS protein partner;
- CAP grammar as a general predictor of GCM1-dependent transcription;
- cooperative causality without motif or partner perturbation.

## Documentation

- [CAP feasibility and shortcut diagnosis](docs/capselex_feasibility.md)
- [Primary genomic-targeting design](docs/genomic_targeting_phase0.md)
- [Independent trophoblast and functional validation](docs/trophoblast_validation.md)
- [Public data access and anti-circularity rules](docs/data_access.md)
- [Machine-readable result summary](reports/result_summary.json)

## Name

`CisGrammarShift` is retained deliberately: “shift” now refers to transfer across assay and biological
contexts—CAP-SELEX → GHT/ChIP → trophoblast → transcription—not merely neural-network distribution shift.
The correct spelling is **Grammar**, not “Grammer”.

## Citation

Use the metadata in [`CITATION.cff`](CITATION.cff). Until perturbational validation exists, cite this as a
leakage-controlled computational reanalysis of cooperative DNA grammar and genomic occupancy.
