# Legacy experiments

`original/` preserves the student-stage CTCF motif-detection scripts, console output, and figures from the
repository history. They are retained for provenance and are not imported by the new package.

The original task is an intentionally easy positive control: a sampled CTCF matrix instance is inserted into
uniform random DNA for positive examples, while negatives are uniform random DNA. Its reported accuracy and
single-example saliency map must not be interpreted as evidence of endogenous TF binding or general
cis-regulatory grammar.

The additional local files reviewed during the upgrade were not copied here. In particular, the file named
`formal_test_data` contained only `synthetic_fallback` positives and `random_noise` negatives, so preserving it
as a formal test set would be misleading.
