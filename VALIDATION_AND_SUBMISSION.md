# Validation and submission

Base commit: 223c79cbd773a4aa324b89a7d445d697d716136d in Wang-Lab-code/Predictor-rank.
Inputs: the previously uploaded VEPs benchmark v5.R and dbnsfp_performance_matrix.csv.

## Changes

The three existing entry points run the shared v5 engine. Original scripts remain
in legacy/ for historical workflows. ROC, PR, bootstrap, fixed-threshold confusion
matrices, group comparisons, heatmaps and dbNSFP comparisons are retained.
The reference CSV is byte-for-byte identical to the upload. The English README,
dependency installer, ignore rules, tests and GitHub Actions are included.
The existing MIT license is unchanged.

## Local validation

R 4.5.1 passed syntax, cutoff-boundary, score-direction, mixed ROC coordinate and
reference-schema checks, plus the v5 self-test.

All three entry points passed synthetic Excel smoke tests with 20 variants, two
predictors and ten bootstrap replicates. Assertions covered TP/TN counts of ten,
AUC/AUPRC of one, expected workbook sheets and at least five SVG outputs.

Existing local dependencies were used. Some packages reported that they were built
under R 4.5.3; analysis also produced warnings without failing the assertions.
Perfect synthetic curves could not be smoothed and fell back to raw curves.
No real study workbook was available. Remote CI has not run. Package versions
are not locked.

## Submission status

Initial connector writes returned HTTP 403 because the GitHub App was not installed
on the repository owner account. After completing installation, repository discovery
and branch creation succeeded. Changes are submitted through a pull request; consult
the pull request checks for current remote CI status. The main branch is unchanged.

## Apply and submit

After restoring authentication, use a repository checkout based on the commit above:

1. Create a new branch.
2. Run git apply --check against Predictor-rank-v5.patch, then apply it.
3. Install dependencies with Rscript scripts/install_dependencies.R.
4. Run Rscript tests/check.R and Rscript tests/smoke.R.
5. Commit the changes and push the new branch.
6. Open a pull request against main and check the remote CI results.

Alternatively, extract Predictor-rank-v5.zip into the checkout on a new branch.
Preserve the existing LICENSE. Restore connector access to this repository before
retrying connector-based submission; support for writing workflow files must also
be available.

## Naming and documentation policy

Project filenames, directory names and standalone documentation use English.
The existing analysis code, including its embedded comments and runtime messages,
was not changed during the documentation-language update.
