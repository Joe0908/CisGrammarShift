# CisGrammarShift

**Counterfactual benchmarking of cis-regulatory grammar learning under distribution shift.**

[![CI](https://github.com/Joe0908/CisGrammarShift/actions/workflows/ci.yml/badge.svg?branch=old-version)](https://github.com/Joe0908/CisGrammarShift/actions/workflows/ci.yml)
[![Python 3.10+](https://img.shields.io/badge/python-3.10%2B-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)

The later CAP-SELEX × GHT-SELEX genomic study is maintained on the
[`main` branch](https://github.com/Joe0908/CisGrammarShift).

## Question

Can sequence models learn relative cis-regulatory syntax rather than exploit
motif-presence shortcuts?

## Design

The benchmark uses matched POU5F1–NANOG sequence pairs. Within each pair, the
background sequence, motif instances, motif content, and orientation are held
constant; only the position of the second motif changes. Labels are determined
by a periodic spacing rule, so motif presence alone cannot solve the task.

The comparison includes:

- a PWM-presence baseline;
- a local CNN;
- a dilated CNN;
- a Transformer;
- five fixed random seeds;
- IID evaluation and targeted shifts in motif gap, orientation, background GC,
  and motif strength.

Decision thresholds are selected on validation data. Evaluation reports
sequence-level discrimination, matched-pair ranking, and attribution
localisation to the implanted motifs.

## Why it matters

High predictive accuracy does not by itself show that a model has learned
regulatory grammar. Matched counterfactual pairs remove the simplest
motif-presence shortcut, while distribution-shift tests reveal whether the
learned rule transfers beyond the training regime.

## Data specification

All sequences are generated locally from a frozen configuration. The bundled
POU5F1 and NANOG matrices are labelled as CIS-BP 3.00 models
`M05705_3.00` and `M05219_3.00`. The complete pair-level data contract is
documented in [`docs/data_specification.md`](docs/data_specification.md).

## Repository structure

| Path | Contents |
|---|---|
| `src/cisgrammar/data.py` | Matched-pair generator and grammar rule |
| `src/cisgrammar/models.py` | LocalCNN, DilatedCNN, and Transformer |
| `src/cisgrammar/baselines.py` | PWM-presence baseline |
| `src/cisgrammar/training.py` | Deterministic training and prediction |
| `src/cisgrammar/metrics.py` | Sequence- and pair-level evaluation |
| `src/cisgrammar/interpretation.py` | Attribution localisation |
| `configs/` | Quick and five-seed reference configurations |
| `tests/` | Unit tests for generation, baselines, metrics, and models |

## Reproduction

Python 3.10 or newer is required.

```bash
git clone --branch old-version https://github.com/Joe0908/CisGrammarShift.git
cd CisGrammarShift
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e '.[dev,ml]'
pytest -q
cisgrammar run --config configs/quick.yaml --output results/quick --device cpu
```

Run the full five-seed experiment with:

```bash
cisgrammar run --config configs/research.yaml --output results/research --device auto
```

Each run writes the frozen configuration, dataset manifests, model checkpoints,
per-sequence predictions, metrics, summaries, figures, and environment
metadata.

## Scope

This study tests recovery and transfer of a programmed sequence rule. It does
not model endogenous chromatin, cofactor abundance, simultaneous TF occupancy,
or causal transcriptional regulation.

## Citation and license

Use [`CITATION.cff`](CITATION.cff) until an archival DOI is available. Code is
released under the [MIT License](LICENSE).
