# CisGrammarShift

**Intrinsic TF specificity, cooperative DNA grammar, and genomic targeting across experimental contexts.**

[![CI](https://github.com/Joe0908/CisGrammarShift/actions/workflows/ci.yml/badge.svg)](https://github.com/Joe0908/CisGrammarShift/actions/workflows/ci.yml)
[![Python 3.10+](https://img.shields.io/badge/python-3.10%2B-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)

## Biological question

TF genomic occupancy is often predicted from a focal TF motif and chromatin context. CAP-SELEX nevertheless
reports thousands of directed, DNA-guided TF-pair interactions. This project asks a narrower mechanistic
question:

> Within an outcome-independent locus universe, does cooperative DNA grammar explain held-out variation in
> focal-TF occupancy beyond intrinsic GHT-SELEX binding, monomer motifs, sequence context, and
> accessibility—and does any increment associate with GCM1-dependent transcription?

The repository keeps the original counterfactual synthetic benchmark as a positive-control layer, but the
primary study is now an accessioned CAP-SELEX × GHT-SELEX × ChIP-seq analysis with human trophoblast
association follow-up.

## Current audit status

The official Codebook v2 GHT-only rerun with the author-supplied McGill-GPHN ChIP tracks is complete. The
predeclared panel claim **failed**: only FLI1, one of four expression-evaluable focal TFs, exceeded the
partial-R² effect threshold. Repeating the frozen analysis with public Toronto-GPZN tracks, without
retuning, preserved every TF-level pass/fail decision and the direction of all five evaluable effects. A
frozen GCM1–TGIF2 follow-up then showed small, directionally consistent increments in both trophoblast
states, but neither reached its predeclared effect threshold. The repository is therefore a complete,
reproducible feasibility and falsification study; it is **not a positive-mechanism manuscript**.

The historical assay-union result is retained only for provenance. In that universe ChIP affected locus
inclusion and could set the sequence-window centre. The official primary universe instead contains
GHT-selected fixed tiles (`ght-only`), for which ChIP affects neither selection nor centering.

## Study layers

| Layer | Role | Leakage control | Endpoint |
|---|---|---|---|
| Synthetic counterfactual | Architecture/implementation positive control | Matched motif instances and pair-level splits | Programmed spacing rule |
| CAP pair feasibility | Diagnose coverage and shortcut learning | Pair-, TF-, family- and paralogue-held-out splits | Directed CAP labels |
| Codebook primary panel | Test induced-expression occupancy increment | GHT-only fixed tiles and chromosome-held-out fitting | Continuous focal-TF ChIP |
| hTSC association replication | Test cross-context reproducibility | Same predeclared GCM1 loci/features; external outcomes only | EVT/ST GCM1 ChIP |
| HiChIP + knockout RNA | Test functional extension | Outcome-blind linkage and nested chromosome folds | Binary targets and DESeq2 LFC |

## Official v2 primary result

All estimates below use chromosome-held-out fitting on official Codebook v2 GHT-only loci. The primary CAP
feature averages only partners detected at >=1 TPM in each of two wild-type HEK293 RNA-seq replicates. The
success rule required partial R² >=0.005, one-sided chromosome-shift p<=0.05, and positive median addition
coefficients for the mean outcome and both biological replicates.

| Focal TF | Expression-supported pairs | Loci with DNase mapping | Partial R² | Chromosome-bootstrap 95% interval | Replicate partial R² | Passed |
|---|---:|---:|---:|---:|---:|---:|
| FLI1 | 16 | 50,060 | 0.011213 | 0.009461–0.012869 | 0.006122 / 0.015334 | Yes |
| GABPA | 4 | 5,041 | 0.004666 | 0.001177–0.008199 | 0.004286 / 0.002268 | No |
| GCM1 | 4 | 19,369 | 0.004497 | 0.002253–0.006312 | 0.003289 / 0.005386 | No |
| RFX5 | 3 | 38,701 | 0.001423 | 0.000729–0.002265 | 0.001408 / 0.001043 | No |

Each screening permutation test gave `p=1/101`; with tens of thousands of loci this is not a substitute for
the frozen effect-size gate. The required four-of-four panel result was one-of-four, so the general CAP
grammar claim is rejected and the 1,000-permutation panel analysis was not run.

Toronto-GPZN sensitivity partial R² values were 0.010679, 0.002881, 0.003826 and 0.001190 for FLI1, GABPA,
GCM1 and RFX5, respectively. MAX remained below threshold in both pipelines (GPHN 0.003005; GPZN
0.002550). This processing concordance supports robustness of direction and threshold decisions, not a
causal cooperative mechanism.

Exploratory pair decomposition identified GATA3–FLI1 and TGIF2–GCM1 as leading pair-specific profiles.
Because this ranking used the Codebook ChIP outcome, it is candidate generation rather than discovery
evidence. TGIF2 passed an outcome-independent trophoblast expression gate before external ChIP was read.
On the same 19,497 GCM1 GHT-only loci, the frozen TGIF2–GCM1 composite score had positive coefficients in
all folds and both clones, but failed the `0.005` threshold in both states:

| External state | Clone correlation | Partial R² | Chromosome-bootstrap 95% interval | Screening p | Passed |
|---|---:|---:|---:|---:|---:|
| EVT | 0.680 | 0.002471 | 0.001009–0.003741 | 1/101 | No |
| ST | 0.692 | 0.003471 | 0.001611–0.005342 | 1/101 | No |

This supports, at most, a weak transferable sequence association. It does not identify TGIF2 co-occupancy
or causal cooperation, and the predeclared 1,000-permutation final test was not triggered.

## Legacy results retained for provenance

### Historical Codebook occupancy panel

- Historical focal TFs: `FLI1`, `GABPA`, `GCM1`, `PAX7`, and `RFX5`; `MAX` was a low-CAP-coverage
  sensitivity analysis, not a CAP-null control.
- Historical locus universe: 1,204,394 200-bp hg38 legacy assay-union loci; 1,035,968 loci in the panel.
- M0: continuous GHT, focal/partner monomer scores, GC, accessibility proxy, and genomic context.
- M1: M0 plus CAP composite and spacing/orientation grammar features calibrated on training chromosomes.
- Full CAP partial R²: `0.02030`; without the DNase proxy: `0.01696`.
- Composite-only: `0.01999`; spacing-only: `0.00060`.
- Representative CAP profiles only: `0.00565`.
- Training-only binned-monomer residualization: `0.01914` (`0.01812–0.02015`, `p = 1/101`).

### Historical GCM1 trophoblast occupancy association

The historical 204,916-locus GCM1 universe was projected into `GSE244252` without allowing the external
hTSC ChIP signal to choose loci or tune features. Because the universe was originally constructed using
Codebook ChIP, these values remain a cross-context association sensitivity rather than a fully
outcome-independent primary result.

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
python -m pip install -e ".[dev,ml]"  # full tests and synthetic training

# For CAP/GHT/ChIP data processing only, the lightweight core is sufficient:
# python -m pip install -e .

# Freeze the non-circular contract
cisgrammar capselex contract --output results/capselex/analysis_contract.json

# Download and verify the small GEO metadata set
python scripts/download_manifest.py \
  --manifest configs/codebook_geo_metadata_manifest.json \
  --output-directory data/codebook/geo_metadata

# Download the formal GHT paper's batch protocol and experiment metadata tables
python scripts/download_manifest.py \
  --manifest configs/ght_nature_supplement_manifest.json \
  --output-directory data/codebook/ght_nature_supplements \
  --resolved-manifest reports/ght_nature_supplement_resolved.json

# Smoke test only: extract and audit the frozen pre-v2 GEO MAGIX panel
python scripts/audit_codebook_magix.py \
  --config configs/codebook_geo_v1_magix_panel.json \
  --metadata-directory data/codebook/geo_metadata \
  --archive data/codebook/geo_peaks/GSE278858_BED_files.tar.gz \
  --extraction-directory data/codebook/geo_peaks/focal_magix \
  --output reports/codebook_geo_v1_magix_asset_audit.json

# Freeze the exact author-selected focal GHT runs and audit v2 rebuild readiness
curl -L --fail \
  'https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJEB76622&result=read_run&fields=run_accession,experiment_title,fastq_ftp,fastq_md5,fastq_bytes&format=tsv&limit=0' \
  -o data/codebook/ena_metadata/PRJEB76622_fastq_filereport_2026-08-10.tsv

python scripts/audit_codebook_ght_rebuild.py \
  --config configs/codebook_ght_v2_rebuild.json \
  --metadata data/codebook/geo_metadata/GSE278858_PRJEB76622_GHT-SELEX_GEO_Fix_with_sra_accessions28012026.xlsx \
  --experiment-metadata data/codebook/ght_nature_supplements/41592_2026_3177_MOESM5_ESM.xlsx \
  --ena-report data/codebook/ena_metadata/PRJEB76622_fastq_filereport_2026-08-10.tsv \
  --output reports/codebook_ght_v2_rebuild_audit.json

# Download and cross-audit the official CAP-SELEX interaction/PWM/spacing supplements
python scripts/download_manifest.py \
  --manifest configs/capselex_nature_supplement_manifest.json \
  --output-directory data/capselex/nature_supplements \
  --resolved-manifest reports/capselex_nature_supplement_resolved.json \
  --jobs 4

python scripts/audit_capselex_nature_supplements.py \
  --interaction-table data/capselex/nature_supplements/41586_2025_8844_MOESM4_ESM.xlsx \
  --pwm-table data/capselex/nature_supplements/41586_2025_8844_MOESM5_ESM.xlsx \
  --spacing-table data/capselex/nature_supplements/41586_2025_8844_MOESM9_ESM.xlsx \
  --output reports/capselex_nature_supplement_audit.json

# Freeze focal and partner monomer controls from held-out-ranked MEX, CAP HT-SELEX and JASPAR 2024
python scripts/audit_monomer_pwm_panel.py \
  --mex-top1 data/codebook/zenodo_pwm/MEX_top1.zip \
  --jaspar data/reference/jaspar2024/JASPAR2024_CORE_vertebrates_non-redundant_pfms_jaspar.txt \
  --cap-pwm-table data/capselex/nature_supplements/41586_2025_8844_MOESM5_ESM.xlsx \
  --cap-audit reports/capselex_nature_supplement_audit.json \
  --output reports/monomer_pwm_panel_audit.json

# Freeze hg38 sequence plus independent HEK293 accessibility and partner-expression controls
python scripts/download_manifest.py \
  --manifest configs/ucsc_hg38_reference_manifest.json \
  --output-directory data/reference/ucsc_hg38 \
  --resolved-manifest reports/ucsc_hg38_reference_resolved.json \
  --jobs 3

python scripts/download_manifest.py \
  --manifest configs/hek293_dnase_manifest.json \
  --output-directory data/reference/hek293_dnase \
  --resolved-manifest reports/hek293_dnase_resolved.json \
  --jobs 2

python scripts/download_manifest.py \
  --manifest configs/hek293_rnaseq_manifest.json \
  --output-directory data/reference/hek293_rnaseq \
  --resolved-manifest reports/hek293_rnaseq_resolved.json \
  --jobs 2

python scripts/audit_hek293_partner_expression.py \
  --replicate-1 data/reference/hek293_rnaseq/GSM3611199_HEK293_rep1_quant.sf.txt.gz \
  --replicate-2 data/reference/hek293_rnaseq/GSM3611199_HEK293_rep2_quant.sf.txt.gz \
  --gencode data/reference/hek293_rnaseq/gencode.v29.annotation.gtf.gz \
  --monomer-audit reports/monomer_pwm_panel_audit.json \
  --output-table results/capselex/hek293_partner_expression.tsv \
  --output-report reports/hek293_partner_expression_audit.json

# Primary ChIP outcome: place the 12 author-supplied McGill GPHN bigWigs here, then verify all hashes,
# bigWig structure and sentinel hg38 chromosome lengths. The 175.39-GB merged archive is not required.
python scripts/audit_gphn_bigwig_panel.py \
  --manifest configs/codebook_chip_gphn_panel_manifest.json \
  --bigwig-directory data/codebook/chip_gphn \
  --output reports/gphn_bigwig_audit.json

# Build per-locus CAP/GHT/ChIP features after the official v2 GHT loci have been frozen
for tf in FLI1 GABPA GCM1 RFX5 PAX7 MAX; do
  python scripts/build_capselex_primary_features.py \
    --focal-tf "$tf" \
    --loci "results/capselex/primary_codebook_v2/${tf}.ght_only.loci.tsv.gz" \
    --magix "data/codebook/ght_v2/Peaks_MAGIX_McGill/${tf}.bed" \
    --chip-resolved-manifest configs/codebook_chip_gphn_panel_manifest.json \
    --chip-directory data/codebook/chip_gphn \
    --reference-resolved-manifest reports/ucsc_hg38_reference_resolved.json \
    --reference-directory data/reference/ucsc_hg38 \
    --cap-pwm-table data/capselex/nature_supplements/41586_2025_8844_MOESM5_ESM.xlsx \
    --mex-top1 data/codebook/zenodo_pwm/MEX_top1.zip \
    --jaspar data/reference/jaspar2024/JASPAR2024_CORE_vertebrates_non-redundant_pfms_jaspar.txt \
    --monomer-audit reports/monomer_pwm_panel_audit.json \
    --output "results/capselex/primary_features_gphn/${tf}.features.tsv.gz" \
    --report "reports/gphn_primary_features/${tf}.json"
done

# Screening model: four expression-evaluable primary TFs, MAX sensitivity, PAX7 availability control
python scripts/run_capselex_primary_models.py \
  --feature-directory results/capselex/primary_features_gphn \
  --chip-processing-pipeline McGill_GPHN_only \
  --focal-tfs FLI1 GABPA GCM1 RFX5 \
  --sensitivity-tfs MAX \
  --availability-negative-controls PAX7 \
  --partner-expression results/capselex/hek293_partner_expression.tsv \
  --dnase-directory data/reference/hek293_dnase \
  --dnase-resolved-manifest reports/hek293_dnase_resolved.json \
  --permutations 100 \
  --partial-r2-threshold 0.005 \
  --minimum-positive-focal-tfs 4 \
  --seed 20260809 \
  --output reports/gphn_capselex_primary_genomic_model_screening.json

# Exploratory pair decomposition after the frozen panel gate failed
python scripts/decompose_capselex_pair_contributions.py \
  --feature-directory results/capselex/primary_features_gphn \
  --focal-tfs FLI1 GABPA GCM1 RFX5 \
  --partner-expression results/capselex/hek293_partner_expression.tsv \
  --dnase-directory data/reference/hek293_dnase \
  --dnase-resolved-manifest reports/hek293_dnase_resolved.json \
  --primary-screening reports/gphn_capselex_primary_genomic_model_screening.json \
  --output reports/gphn_capselex_pair_decomposition.json

# Public Toronto-GPZN acquisition and the no-retuning rerun are retained separately; see
# provenance/legacy_gpzn_fallback/README.md. After that rerun:
python scripts/compare_capselex_chip_pipelines.py \
  --gphn-screening reports/gphn_capselex_primary_genomic_model_screening.json \
  --gpzn-screening reports/capselex_primary_genomic_model_screening.json \
  --output reports/gphn_gpzn_pipeline_comparison.json

# Freeze, acquire and test the outcome-independent TGIF2 expression gate
python scripts/download_manifest.py \
  --manifest configs/trophoblast_gcm1_replication_manifest.json \
  --output-directory data/trophoblast/gcm1_replication \
  --resolved-manifest reports/trophoblast_gcm1_replication_resolved.json \
  --jobs 4

python scripts/audit_trophoblast_tgif2_expression.py \
  --archive data/trophoblast/gcm1_replication/GSE244254_RAW.tar \
  --resolved-manifest reports/trophoblast_gcm1_replication_resolved.json \
  --output reports/trophoblast_tgif2_expression_gate.json

# Frozen external-context screening; 1,000 permutations run only if both states pass
python scripts/run_trophoblast_tgif2_replication.py \
  --features results/capselex/primary_features_gphn/GCM1.features.tsv.gz \
  --asset-directory data/trophoblast/gcm1_replication \
  --resolved-manifest reports/trophoblast_gcm1_replication_resolved.json \
  --expression-gate reports/trophoblast_tgif2_expression_gate.json \
  --pair-decomposition reports/gphn_capselex_pair_decomposition.json \
  --permutations 100 \
  --partial-r2-threshold 0.005 \
  --output reports/trophoblast_tgif2_gcm1_replication_screening.json

# Build the directed CAP pair dataset and run the TF-held-out panel
cisgrammar capselex dataset \
  --interaction-table data/capselex/interaction_matrix.xlsx \
  --construct-table data/capselex/constructs.xlsx \
  --output results/capselex/tf_pair_dataset.tsv

cisgrammar capselex tf-panel \
  --dataset results/capselex/tf_pair_dataset.tsv \
  --target composite_motif \
  --output results/capselex/tf_panel_composite.json

# Primary: GHT chooses tiles; ChIP is outcome/annotation only
cisgrammar capselex build-loci \
  --universe ght-only \
  --chip-peaks data/codebook/chip/GCM1.narrowPeak \
  --ght-magix data/codebook/ght_v2/GCM1_MAGIX.bed \
  --focal-tf GCM1 \
  --output results/capselex/GCM1.ght_only.loci.tsv.gz

# Sensitivity: neither assay chooses the full-genome tiles (written chromosome-by-chromosome)
cisgrammar capselex build-loci \
  --universe fixed-genome \
  --chrom-sizes data/reference/hg38.chrom.sizes \
  --chip-peaks data/codebook/chip/GCM1.narrowPeak \
  --ght-magix data/codebook/ght_v2/GCM1_MAGIX.bed \
  --focal-tf GCM1 \
  --output results/capselex/GCM1.fixed_genome.loci.tsv.gz

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

`build-loci` writes a sidecar manifest containing the declared universe, outcome-dependence flags, row
counts, input/output sizes and SHA-256 hashes. Raw and large processed files are intentionally excluded from
Git. See
[`docs/data_access.md`](docs/data_access.md) for public accessions and the small-file download strategy.

For MAGIX input, the frozen primary rule is Benjamini-Hochberg FDR <= 0.05 and a positive refined
`coefficient.ar`. The manifest records both choices. The current official Codebook v2 web archive was
acquired on 2026-08-10. Its six focal BED members are byte-identical to the corresponding GEO members even
though the tar archive itself has a different SHA-256 and simplified member names. The v2 provenance,
member hashes and eligibility audit are frozen in `configs/codebook_ght_v2_magix_focal_manifest.json` and
`reports/codebook_ght_v2_magix_focal_audit.json`.

Selective raw-data acquisition is feasible: the exact target FASTQs linked to the five-TF acquisition panel
total 2.260 GiB; adding MAX totals 2.896 GiB. This is not yet an exact v2 rebuild. The released
MAGIX production workflow uses batch-aggregate covariates. Supplementary Table 3 now resolves the nine
relevant batch identifiers, but neither the public source tree nor the paper freezes whether those aggregate
lists contained all experiments (79.476 GiB for the relevant batches) or approved experiments only
(26.870 GiB). The repository therefore keeps `exact_author_v2_rebuild_ready=false` for a from-FASTQ
reconstruction; it does not choose the cheaper aggregate definition or replace the missing production
design with an easier model. That reconstruction limitation no longer blocks the primary analysis because
the official processed v2 focal BED members are now available and checksummed.

The primary McGill-GPHN ChIP outcome retains both biological replicates. Its value is the mean of the two
replicate `log1p` exact 200-bp mean signals; both replicates are additionally fit as separate sensitivity
outcomes, and a claimed CAP increment must have the same direction in both. Toronto-GPZN tracks are kept in
a separate no-retuning processing sensitivity and are never mixed with GPHN tracks.

Partner availability is frozen before ChIP modelling from two wild-type HEK293 RNA-seq replicates in
GSM3611199. Transcript TPMs are summed to GENCODE v29 genes; a partner is expression-supported only when it
has at least 1 TPM in each replicate. PAX7's only representative pair is PAX7–TBX4, while TBX4 is 0 TPM in
both replicates. PAX7 is therefore an availability negative control rather than an interpretable primary
cooperativity case. The primary panel is FLI1, GABPA, GCM1 and RFX5; MAX–TEAD4 remains the predeclared
low-CAP-coverage sensitivity.

## Repository layout

```text
configs/                 Frozen synthetic and CAP reanalysis configurations
docs/                    Feasibility, data, modelling, validation and claim-boundary documents
reports/                 Small versioned result summaries; no raw biological data
scripts/                 Checksum downloader and exact DESeq2 runner
provenance/              Retained non-default acquisition and processing workflows
src/cisgrammar/          Synthetic benchmark and complete CAP/GHT/ChIP/hTSC pipeline
tests/                   Unit, leakage-control and synthetic-recovery tests
legacy/original/         Preserved student-stage scripts and figures
```

## Claim boundary

Supported result:

> CAP-derived sequence features showed heterogeneous held-out occupancy increments on outcome-independent
> GHT-only loci; only FLI1 passed the predeclared TF-level effect criterion.

The directionally consistent TGIF2–GCM1 trophoblast result is reported as a weak association that failed its
predeclared effect criterion, not as replication of a mechanism.

Not supported:

- a uniquely identified in-cell ETS protein partner;
- a general CAP-grammar occupancy mechanism across TFs;
- TGIF2 co-occupancy or causal cooperation with GCM1;
- CAP grammar as a general predictor of GCM1-dependent transcription;
- cooperative causality without motif or partner perturbation.

## Documentation

- [CAP feasibility and shortcut diagnosis](docs/capselex_feasibility.md)
- [Primary genomic-targeting design](docs/genomic_targeting_phase0.md)
- [Independent trophoblast and functional validation](docs/trophoblast_validation.md)
- [Research audit and publication decision](docs/research_audit_2026-08-10.md)
- [Public data access and anti-circularity rules](docs/data_access.md)
- [Retained Toronto-GPZN fallback and sensitivity](provenance/legacy_gpzn_fallback/README.md)
- [Machine-readable result summary](reports/result_summary.json)


## Citation

Use the metadata in [`CITATION.cff`](CITATION.cff). Until perturbational validation exists, cite this as a
leakage-controlled computational reanalysis of cooperative DNA grammar and genomic occupancy.
