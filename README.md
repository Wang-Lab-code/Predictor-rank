# Predictor-rank

Reproducible comparison of variant effect predictors (VEPs) in a labelled
case–control variant cohort. The v5 workflow combines continuous-score
discrimination, fixed-threshold classification, bootstrap uncertainty,
score distributions and comparison with a supplied dbNSFP performance matrix.

## Repository layout

- `ROC.R`, `Confusion Matrix.R`, `Heatmap.R`: compatible entry names; each now runs the complete shared v5 workflow.
- `R/benchmark_v5.R`: analysis engine derived from the uploaded **VEPs benchmark v5.R**.
- `data/reference/dbnsfp_performance_matrix.csv`: uploaded reference matrix, preserved without modification.
- `legacy/`: original three scripts, preserved for historical analyses (including categorical prediction-label workflows and the older standalone heatmap). These require their original dependencies and path settings.
- `scripts/`, `tests/`, `.github/workflows/`: dependency installation, regression checks and CI.

## Environment and execution

R 4.5.1 is the CI target. From the repository root:

```sh
Rscript scripts/install_dependencies.R
Rscript ROC.R --input=/path/to/variants.xlsx --sheet=VEPs --gene=GENE
# Equivalent complete workflows:
Rscript "Confusion Matrix.R" --input=/path/to/variants.xlsx --gene=GENE
Rscript Heatmap.R --input=/path/to/variants.xlsx --gene=GENE
```

Quote the entire `--input=...` argument if the path contains spaces. Use an explicit
input path; automatic discovery only works when exactly one input workbook exists
at the repository root. The preferred worksheet is `VEPs`; the engine searches
for a sheet containing the required identifier and group columns if necessary.
Use a simple gene label suitable for a filename. Results are written under the
repository root; repeated runs with the same gene label replace matching outputs.
Dependencies are checked, never installed implicitly during analysis.

## Input contract

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

## Methods and interpretation

The uploaded v5 configuration is retained: thresholds, strict comparison operators
and predictor categories are recorded in the exported `ThresholdConfig` sheet.
Fixed-threshold outputs include TP, TN, FP, FN, accuracy, balanced accuracy,
sensitivity, specificity, PPV, NPV, MCC and F1. Undefined denominators yield NA.
Score distributions are compared using the v5 statistical procedure and visualized
alongside per-tool and combined ROC and precision–recall curves.

ROC AUC uses pROC with automatic direction, as in the supplied v5; PR scores are
oriented according to the configured cutoff direction. AUPRC uses trapezoidal
integration and average precision uses recall increments. They are distinct
summaries. Auto-selected ROC direction can produce optimistic discrimination
estimates; inspect direction and external validation before interpreting rankings.
Smoothed ROC curves are for display; raw and smoothed coordinate frames are
normalized before combining to avoid incompatible-column errors.

Bootstrap defaults to 1,000 row-resampling replicates, seed 20260606 and 95%
percentile intervals. `--bootstrap=10` is available for smoke testing;
`--bootstrap=0` disables it. This measures within-cohort sampling variability,
not external validation, and does not account for related or clustered samples.
One-class replicates cannot contribute valid ROC/PR estimates.

Some v5 thresholds are applied to rankscore columns (notably gMVP, LIST-S2 and
AlphaMissense). These settings are retained for fidelity to the uploaded baseline;
confirm the score scale and threshold provenance before scientific use.
The bundled CSV compares historical/latest dbNSFP AUC with the current cohort.
Differences are descriptive: cohorts, curation, prevalence and evaluation strata
may differ, so these are not paired significance tests or evidence of superiority.
Reference data provenance, upstream dbNSFP version and redistribution terms
should be supplied by the dataset owner; they are not inferred from the filename.

## Outputs

`GENE_benchmark_v5_results/01_tables/GENE_benchmark_v5_summary.xlsx` contains
`SampleInfo`, `ThresholdConfig`, `ThresholdMetrics`, `ROCSummary`,
`BootstrapValidation`, `BenchmarkSummary`, `DbnsfpPerformanceRaw`,
`DbnsfpHistoricalVsBenchmark`, `DbnsfpLatestVsBenchmark` and `SkippedTools`.
`02_plots/` contains SVG threshold heatmaps, score comparisons, combined and
category-specific ROC/PR plots and dbNSFP comparisons; `per_tool/` holds
individual predictor plots.

## Verification

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

## License

The existing MIT software license and copyright notice are preserved unchanged.
The bundled reference dataset remains subject to applicable upstream data terms.
For publications, cite the predictor and dbNSFP resources used in your study and
record the repository commit, R/package versions, cohort selection and score scales.
