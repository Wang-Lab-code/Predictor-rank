# Predictor-rank

**A unified R workflow for benchmarking variant effect predictors.**

[![R benchmark checks](https://github.com/Wang-Lab-code/Predictor-rank/actions/workflows/r-tests.yml/badge.svg?branch=codex%2Fbenchmark-v5-upgrade)](https://github.com/Wang-Lab-code/Predictor-rank/actions/workflows/r-tests.yml)
[![R version](https://img.shields.io/badge/R-4.5.1-276DC3?logo=r&logoColor=white)](#quick-start)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

Compare pathogenic and control variants using continuous predictor scores,
fixed decision thresholds and bootstrap uncertainty. The v5 workflow produces
Excel summary tables and SVG figures, including ROC and precision–recall curves,
confusion-matrix heatmaps and comparisons with a supplied dbNSFP reference matrix.

**[Quick start](#quick-start)** · **[Input format](#input-format)** ·
**[Outputs](#outputs)** · **[Methods](#methods-and-interpretation)** ·
**[Validation](#validation)** · **[Project structure](#project-structure)**

## At a glance

| Analysis | What you receive |
| --- | --- |
| Discrimination | Per-tool, combined and category-specific ROC and PR curves |
| Classification | TP, TN, FP, FN and fixed-threshold performance metrics |
| Uncertainty | Bootstrap summaries and percentile confidence intervals |
| Score distributions | Group comparisons and predictor-level visualizations |
| Reference comparison | Historical and latest dbNSFP AUC alongside cohort estimates |
| Reporting | A multi-sheet Excel workbook and SVG figures |

> **v5 entry points:** `ROC.R`, `Confusion Matrix.R` and `Heatmap.R` each run
> the complete shared workflow. Original standalone scripts remain in
> [legacy/](legacy/) for historical analyses.

## Quick start

Use R 4.5.1, the CI target version. Run from the repository root.

### Install dependencies and run

Prepare the workbook using the [input format](#input-format) below.

```sh
Rscript scripts/install_dependencies.R
Rscript ROC.R --input=/path/to/variants.xlsx --sheet=VEPs --gene=GENE
# Alternative entry points: the same complete workflow
Rscript "Confusion Matrix.R" --input=/path/to/variants.xlsx --gene=GENE
Rscript Heatmap.R --input=/path/to/variants.xlsx --gene=GENE
```

### Run options

| Option | Default | Purpose |
| --- | --- | --- |
| `--input=PATH` | Automatic discovery | Input Excel workbook |
| `--sheet=NAME` | `VEPs` | Preferred worksheet |
| `--gene=LABEL` | Input filename stem | Output label |
| `--bootstrap=N` | `1000` | Resampling replicates; `0` disables bootstrap |

Quote the entire `--input=...` argument if the path contains spaces. Use an explicit
input path; automatic discovery only works when exactly one input workbook exists
at the repository root. The preferred worksheet is `VEPs`; the engine searches
for a sheet containing the required identifier and group columns if necessary.
Use a simple gene label suitable for a filename. Results are written under the
repository root; repeated runs with the same gene label replace matching outputs.
Dependencies are checked, never installed implicitly during analysis.

## Input format

Supply an Excel workbook with one row per variant and unchanged predictor column names.

| Column | Required content |
| --- | --- |
| `variant_id_hg38` | Variant identifier |
| `Type` | Positive: `Disease-associated` or `Pathogenic`; negative: `Control` or `Benign` |
| Predictor columns | Numeric scores named exactly as `tool_config$ScoreColumn` in the engine, e.g. `SIFT_score`, `REVEL_score` |

Unknown truth labels are excluded. Missing or insufficient predictor data are
reported in `SkippedTools`; usable observations can differ between tools.
Deduplicate and quality-control variants before analysis: the engine does not
establish independence, resolve transcript choices or perform clinical adjudication.
Do not commit patient workbooks; Excel inputs and generated results are ignored.

## Outputs

```text
GENE_benchmark_v5_results/
├── 01_tables/
│   └── GENE_benchmark_v5_summary.xlsx
└── 02_plots/
    ├── per_tool/
    └── ... combined curves, heatmaps and reference-comparison SVG files
```

### Workbook guide

| Sheet | Contents |
| --- | --- |
| `SampleInfo` | Input, worksheet and cohort information |
| `ThresholdConfig` | Predictor categories, score columns, operators and cutoffs |
| `ThresholdMetrics` | Confusion counts and classification metrics |
| `ROCSummary` | ROC AUC, PR summaries and ROC threshold |
| `BootstrapValidation` | Resampling summaries and confidence intervals |
| `BenchmarkSummary` | Joined classification and discrimination results |
| `DbnsfpPerformanceRaw` | Supplied reference matrix |
| `DbnsfpHistoricalVsBenchmark` | Historical AUC comparisons |
| `DbnsfpLatestVsBenchmark` | Latest-year AUC comparisons |
| `SkippedTools` | Excluded predictors and reasons |

`02_plots/` contains SVG threshold heatmaps, score comparisons, combined and
category-specific ROC/PR plots and dbNSFP comparisons. `per_tool/` contains
individual predictor plots.

## Methods and interpretation

The uploaded v5 configuration is retained: thresholds, strict comparison operators
and predictor categories are recorded in the exported `ThresholdConfig` sheet.
Fixed-threshold outputs include TP, TN, FP, FN, accuracy, balanced accuracy,
sensitivity, specificity, PPV, NPV, MCC and F1. Undefined denominators yield NA.
Score distributions are compared using the v5 statistical procedure and visualized
alongside per-tool and combined ROC and precision–recall curves.

### ROC and precision–recall

ROC AUC uses pROC with automatic direction, as in the supplied v5; PR scores are
oriented according to the configured cutoff direction. AUPRC uses trapezoidal
integration and average precision uses recall increments. They are distinct
summaries. Auto-selected ROC direction can produce optimistic discrimination
estimates; inspect direction and external validation before interpreting rankings.
Smoothed ROC curves are for display; raw and smoothed coordinate frames are
normalized before combining to avoid incompatible-column errors.

### Bootstrap uncertainty

Bootstrap defaults to 1,000 row-resampling replicates, seed 20260606 and 95%
percentile intervals. `--bootstrap=10` is available for smoke testing;
`--bootstrap=0` disables it. This measures within-cohort sampling variability,
not external validation, and does not account for related or clustered samples.
One-class replicates cannot contribute valid ROC/PR estimates.

### Score scales and reference comparisons

Some v5 thresholds are applied to rankscore columns (notably gMVP, LIST-S2 and
AlphaMissense). These settings are retained for fidelity to the uploaded baseline;
confirm the score scale and threshold provenance before scientific use.
The bundled CSV compares historical/latest dbNSFP AUC with the current cohort.
Differences are descriptive: cohorts, curation, prevalence and evaluation strata
may differ, so these are not paired significance tests or evidence of superiority.
Reference data provenance, upstream dbNSFP version and redistribution terms
should be supplied by the dataset owner; they are not inferred from the filename.

## Validation

```sh
Rscript tests/check.R
Rscript ROC.R --self-test
Rscript tests/smoke.R
```

CI parses all maintained scripts, tests cutoff boundaries and mixed ROC coordinates,
checks the reference schema and runs all three entries on synthetic Excel data,
asserting perfect classification/discrimination and generated workbook/SVG outputs.
Synthetic smoke tests use ten bootstrap replicates; they do not validate biological
claims or replace a full run with the study workbook.

See [Validation and submission](VALIDATION_AND_SUBMISSION.md) for local test details.

## Project structure

```text
Predictor-rank/
├── ROC.R                         # Complete v5 workflow entry point
├── Confusion Matrix.R            # Equivalent compatibility entry point
├── Heatmap.R                     # Equivalent compatibility entry point
├── R/benchmark_v5.R              # Shared engine and predictor configuration
├── data/reference/               # Supplied dbNSFP performance matrix
├── legacy/                       # Original standalone scripts
├── scripts/                      # Dependency installer
├── tests/                        # Helper and synthetic integration checks
├── .github/workflows/            # Continuous integration
├── VALIDATION_AND_SUBMISSION.md
└── LICENSE
```


- `ROC.R`, `Confusion Matrix.R`, `Heatmap.R`: compatible entry names; each now runs the complete shared v5 workflow.
- `R/benchmark_v5.R`: analysis engine derived from the uploaded **VEPs benchmark v5.R**.
- `data/reference/dbnsfp_performance_matrix.csv`: uploaded reference matrix, preserved without modification.
- `legacy/`: original three scripts, preserved for historical analyses (including categorical prediction-label workflows and the older standalone heatmap). These require their original dependencies and path settings.
- `scripts/`, `tests/`, `.github/workflows/`: dependency installation, regression checks and CI.

## License and reporting

The existing MIT software license and copyright notice are preserved unchanged.
The bundled reference dataset remains subject to applicable upstream data terms.
For publications, cite the predictor and dbNSFP resources used in your study and
record the repository commit, R/package versions, cohort selection and score scales.
