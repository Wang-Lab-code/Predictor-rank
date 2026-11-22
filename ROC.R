# Compatibility entry point: the complete v5 benchmark is shared by all entries.
args <- commandArgs(trailingOnly = FALSE)
file_arg <- args[startsWith(args, "--file=")]
source_file <- if (length(file_arg)) sub("--file=", "", file_arg[1], fixed = TRUE) else sys.frame(1)$ofile
root <- if (is.null(source_file)) getwd() else dirname(normalizePath(source_file, mustWork = TRUE))
options(predictor.rank.root = root)
source(file.path(root, "R", "benchmark_v5.R"), encoding = "UTF-8", chdir = FALSE)
