# CisGrammarShift

**Two linked studies of cis-regulatory grammar: controlled computational learnability and genomic relevance.**

[![CI](https://github.com/Joe0908/CisGrammarShift/actions/workflows/ci.yml/badge.svg)](https://github.com/Joe0908/CisGrammarShift/actions/workflows/ci.yml)
[![Python 3.10+](https://img.shields.io/badge/python-3.10%2B-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)

## Two linked studies

### 1. Counterfactual grammar-learning benchmark

Can sequence models learn relative motif syntax rather than exploit
motif-presence shortcuts? This controlled PyTorch study uses matched
POU5F1–NANOG sequence pairs, a PWM-presence baseline, LocalCNN, DilatedCNN, and
Transformer models, five random seeds, and IID plus gap, orientation, GC, and
motif-strength distribution shifts.

### 2. Genomic cooperative-grammar study

Does experimentally measured cooperative sequence grammar explain real focal-TF
occupancy beyond intrinsic binding and conventional sequence features? This is
the current main study and integrates CAP-SELEX, GHT-SELEX, monomer motifs,
sequence composition, accessibility, and held-out genomic outcomes.

The first study asks whether regulatory grammar is computationally learnable
under controlled counterfactual interventions; the second asks whether
experimentally observed grammar carries additional information in endogenous
genomic targeting.

## Current genomic study

On an outcome-independent set of genomic loci, do CAP-SELEX-derived
cooperative sequence features explain held-out focal-TF occupancy beyond
intrinsic GHT-SELEX binding, monomer motifs, sequence composition, and an
independent accessibility proxy?

## Why it matters

Apparent grammar effects can arise when the genomic outcome is used to select
loci, centre sequence windows, tune features, or split nearby loci across train
and test sets. This project separates in-vitro sequence evidence from genomic
outcomes and evaluates incremental signal under chromosome-held-out fitting.

## Genomic study approach

- Freeze GHT-only 200-bp hg38 loci without using ChIP signal for inclusion or
  centring.
- Build focal and partner monomer controls plus CAP composite and
  spacing/orientation features.
- Filter candidate partners using independent HEK293 expression data.
- Compare nested ridge models on held-out chromosomes using out-of-fold partial
  R², chromosome-block intervals, replicate-direction checks, and
  within-chromosome spatial nulls.
- Repeat the frozen specification with an independent ChIP processing pipeline
  and an external trophoblast context.

The primary comparison is:

\[
M_0 = f(\mathrm{GHT},\ \mathrm{monomers},\ \mathrm{GC/CpG},\
\mathrm{accessibility},\ \mathrm{context})
\]

\[
M_1 = M_0 + \mathrm{CAP\ grammar}.
\]

## Counterfactual study

The earlier counterfactual grammar-learning study motivated the later real-data
analysis and now also provides a controlled benchmark within this repository.
It remains a distinct scientific project: matched sequence pairs hold motif
content constant while varying relative syntax, allowing motif-presence
baselines and neural models to be compared under both IID and targeted
distribution shifts. Its results concern recovery of a programmed sequence
rule, not endogenous TF binding.

## Data and provenance

The analysis integrates public CAP-SELEX, GHT-SELEX, reference motifs,
HEK293 DNase/RNA, and trophoblast datasets, together with an author-supplied
McGill-GPHN ChIP bigWig panel. Large source files are not redistributed.

Exact accessions, URLs, genome builds, checksums, file manifests, expected
paths, and access limitations are documented in
[`docs/data_access.md`](docs/data_access.md). Frozen manifests are in
`configs/`, and resolved acquisition/QC records are in `reports/`.

## Repository structure

| Path | Contents |
|---|---|
| `src/cisgrammar/` | Shared infrastructure; counterfactual modules (`data.py`, `models.py`, `baselines.py`, `training.py`) and genomic `capselex_*` modules |
| `scripts/` | Public-data acquisition, QC, feature building, and model entry points |
| `configs/` | Frozen manifests and analysis parameters |
| `reports/` | Machine-readable QC and reference result summaries |
| `tests/` | Unit tests for data contracts, features, models, and leakage controls |
| `docs/data_access.md` | Data provenance and acquisition instructions |
| `docs/data_specification.md` | Synthetic-control data contract |
| `provenance/gpzn_sensitivity/` | Reproduction notes for the public GPZN processing sensitivity |

## Reproduction

Python 3.10 or newer is required.

```bash
git clone https://github.com/Joe0908/CisGrammarShift.git
cd CisGrammarShift
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e '.[dev]'
pytest -q
```

Install the optional CPU or CUDA build of PyTorch to run the counterfactual
grammar-learning study:

```bash
python -m pip install -e '.[dev,ml]'
cisgrammar run --config configs/quick.yaml --output results/quick --device cpu
```

For the real-data workflow, first resolve and verify assets with
`scripts/download_manifest.py` and the audit commands listed in
[`docs/data_access.md`](docs/data_access.md). The primary feature and model
entry points are:

```bash
python scripts/build_capselex_primary_features.py --help
python scripts/run_capselex_primary_models.py --help
```

The full reference run requires large public resources and the checksummed
author-supplied GPHN tracks named in
`configs/codebook_chip_gphn_panel_manifest.json`.

## Headline findings

- Under the frozen GHT-only design, one of four expression-evaluable focal TFs
  exceeded the prespecified partial-R² effect threshold; the panel-level rule
  was therefore not met.
- Repeating the same analysis with the public Toronto-GPZN ChIP processing
  pipeline preserved all focal-TF pass/fail decisions and effect directions.
- A frozen external-context TGIF2-GCM1 sequence score was directionally
  consistent across two trophoblast states but did not exceed the prespecified
  effect threshold in either state.

Detailed per-TF estimates, QC records, and sensitivity outputs remain available
as machine-readable files under `reports/`; this README does not reconstruct a
results or discussion narrative.

## Scope and status

The genomic study supports a heterogeneous, TF-specific feasibility assessment
rather than a general cooperative-mechanism claim. CAP motif scores do not
establish simultaneous protein occupancy, and genomic loci are not independent
biological replicates. The counterfactual study answers a separate question:
whether models recover a programmed sequence rule under matched interventions
and distribution shift; it is not treated as endogenous binding evidence.

The public repository intentionally excludes manuscript planning, publication
audits, speculative extensions, obsolete implementations, and future research
roadmaps.

## Citation and license

Use [`CITATION.cff`](CITATION.cff) until an archival DOI is available. Code is
released under the [MIT License](LICENSE); source datasets retain their original
terms.
