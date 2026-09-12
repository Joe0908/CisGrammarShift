# CisGrammarShift

**Leakage-controlled evaluation of cooperative DNA grammar in genomic TF targeting.**

[![CI](https://github.com/Joe0908/CisGrammarShift/actions/workflows/ci.yml/badge.svg)](https://github.com/Joe0908/CisGrammarShift/actions/workflows/ci.yml)
[![Python 3.10+](https://img.shields.io/badge/python-3.10%2B-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)

The earlier counterfactual grammar-learning study is preserved on the
[`old-version` branch](https://github.com/Joe0908/CisGrammarShift/tree/old-version).

## Question

On an outcome-independent set of genomic loci, do CAP-SELEX-derived
cooperative sequence features explain held-out focal-TF occupancy beyond
intrinsic GHT-SELEX binding, monomer motifs, sequence composition, and an
independent accessibility proxy?

## Why it matters

Apparent grammar effects can arise when the genomic outcome is used to select
loci, centre sequence windows, tune features, or split nearby loci across train
and test sets. This project separates in-vitro sequence evidence from genomic
outcomes and evaluates incremental signal under chromosome-held-out fitting.

## Approach

- Freeze GHT-only 200-bp hg38 loci without using ChIP signal for inclusion or
  centring.
- Build focal and partner monomer controls plus CAP composite and
  spacing/orientation features.
- Filter candidate partners using independent HEK293 expression data.
- Compare nested ridge models on held-out chromosomes using out-of-fold partial
  R², chromosome-block intervals, replicate-direction checks, and
  within-chromosome spatial nulls.
- Test the frozen sequence hypothesis in an external trophoblast context.

The primary comparison is:

\[
M_0 = f(\mathrm{GHT},\ \mathrm{monomers},\ \mathrm{GC/CpG},\
\mathrm{accessibility},\ \mathrm{context})
\]

\[
M_1 = M_0 + \mathrm{CAP\ grammar}.
\]

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
| `src/cisgrammar/` | Reusable CAP/GHT feature, model, metric, and provenance modules |
| `scripts/` | Public-data acquisition, QC, feature building, and model entry points |
| `configs/` | Frozen manifests and analysis parameters |
| `reports/` | Machine-readable QC and reference result summaries |
| `tests/` | Unit tests for data contracts, features, models, and leakage controls |
| `docs/data_access.md` | Data provenance and acquisition instructions |

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

For the genomic workflow, first resolve and verify assets with
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
- A frozen external-context TGIF2-GCM1 sequence score was directionally
  consistent across two trophoblast states but did not exceed the prespecified
  effect threshold in either state.

Detailed per-TF estimates, QC records, and external-context outputs remain available
as machine-readable files under `reports/`; this README does not reconstruct a
results or discussion narrative.

## Scope and status

The completed analysis supports a heterogeneous, TF-specific feasibility
assessment rather than a general cooperative-mechanism claim. CAP motif scores
do not establish simultaneous protein occupancy, and genomic loci are not
independent biological replicates.

The public repository intentionally excludes manuscript planning, publication
audits, speculative extensions, obsolete implementations, and future research
roadmaps.

## Citation and license

Use [`CITATION.cff`](CITATION.cff) until an archival DOI is available. Code is
released under the [MIT License](LICENSE); source datasets retain their original
terms.
