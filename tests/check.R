options(encoding = "UTF-8")
files <- c("ROC.R", "Confusion Matrix.R", "Heatmap.R", "R/benchmark_v5.R")
invisible(lapply(files, function(f) parse(f, encoding = "UTF-8")))
# Evaluate function definitions only; never run the analysis during unit checks.
env <- new.env(parent = globalenv())
for (expr in parse("R/benchmark_v5.R", encoding = "UTF-8")) {
  if (is.call(expr) && identical(expr[[1]], as.name("<-")) &&
      is.call(expr[[3]]) && identical(expr[[3]][[1]], as.name("function"))) eval(expr, env)
}
stopifnot(identical(env$apply_cutoff(c(.04,.05,.06), "<", .05), c(1L,0L,0L)) ||
          identical(env$apply_cutoff(c(.04,.05,.06), "<", .05), c(1,0,0)))
stopifnot(identical(env$orient_score_for_positive(c(1,2), "<"), c(-1,-2)))
raw <- data.frame(threshold=c(Inf,0), specificity=c(1,0), sensitivity=c(0,1))
smooth <- raw[c("specificity","sensitivity")]
combined <- env$combine_roc_coordinate_frames(list(raw=raw,smooth=smooth))
stopifnot(nrow(combined)==4, identical(names(combined),c("FPR","TPR","Method")))
db <- read.csv("data/reference/dbnsfp_performance_matrix.csv")
stopifnot(nrow(db)>0, all(c("Score","ClinVar_Year","AUC") %in% names(db)),
          all(db$AUC >= 0 & db$AUC <= 1, na.rm=TRUE))
message("PASS: syntax, cutoff boundary, direction, mixed ROC coordinates, reference schema")
