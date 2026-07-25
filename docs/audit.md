# Audit of the student-stage project

## What was useful

- The original project asked a sensible introductory question: which architectures can detect an implanted
  motif?
- It used a biologically derived CTCF probability matrix instead of one fixed consensus string.
- It compared several model families and attempted an input-gradient interpretation.
- The extended scripts moved toward two-motif logic, spacing stress tests, multiple seeds, reverse complements,
  AUROC/AUPRC, and simulated methylation.

Those steps provide a good technical foundation. They do not yet define a mature biological study.

## Problems that prevent research-level claims

| Problem | Consequence | Upgrade |
|---|---|---|
| Positives contain an implanted motif; negatives are uniform random DNA | Motif presence is a shortcut | Both pair members contain identical motif instances |
| Synthetic backgrounds are uniformly random | Unrealistic sequence composition and easy negatives | Configurable GC and first-order persistence |
| One split and initially one seed | Results may depend on sampling or initialisation | Pre-specified multi-seed runs |
| Accuracy is the main endpoint | Misleading under imbalance; no calibration | AUROC, AUPRC, MCC, Brier, ECE, paired metrics |
| Architectures have unequal inductive biases and receptive fields | Ranking does not isolate the scientific question | Explicit local versus long-range controls |
| One correctly predicted positive is selected for saliency | Selection bias; no quantitative test | Seeded, unselected attribution set and localisation metric |
| Train and test use nearly identical generators | No evidence of robustness | Gap, orientation, GC, and motif-strength shifts |
| Model files lack a frozen config, data manifest, and software versions | Results cannot be reproduced | Run snapshot, metadata, per-seed metrics, state dicts |
| Matrix provenance is implicit | Motif identity cannot be audited | Named accessions and source metadata in code |

## Invalid “formal test” artifact

The reviewed `formal_test_data(1).csv` has 400 rows:

- 200 positive rows, all with `location=synthetic_fallback`;
- 200 negative rows, all with `location=random_noise`;
- zero experimentally or genomically sourced rows.

The associated retrieval script labels motif-overlap API output as “validated binding sites,” but a motif
annotation is not equivalent to experimentally observed occupancy. When retrieval failed, it generated random
DNA and retained a positive TF label. This means the file cannot evaluate TF binding and must not be called a
real or formal test set.

## Extended spacing script issue

An intermediate generator iterated a fixed number of attempts but discarded some negative examples whenever
both motifs were sampled. It could therefore return fewer examples than requested and alter class balance.
The new generator uses an exact number of matched pairs and asserts its invariants.

## Methylation simulation issue

The extended methylation script changed motif sampling probabilities using an artificial methylation track,
but the label was not derived from a separately specified causal methylation rule. A model could exploit
generator-specific correlations without learning a biological methylation effect.

Methylation is therefore excluded from the core benchmark. It should return only after:

1. the target TF has methylation-sensitive binding evidence;
2. methylated and unmethylated motifs have traceable experimental matrices or readouts;
3. the intervention and label mechanism are pre-specified;
4. real-data validation is available.

## Claim correction

The old result “CNN achieved 96.3% accuracy” is accurately described as:

> On one synthetic split of 100-bp uniform-random sequences, a small CNN distinguished sequences containing a
> sampled CTCF matrix instance from random negatives with 96.3% accuracy.

It is not evidence that the model predicts endogenous CTCF occupancy or that CNNs generally outperform
Transformers in regulatory genomics.
