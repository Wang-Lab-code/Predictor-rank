################################################################################

################################################################################

## ===== 0. Adjustable Parameters ===============================================
data_dir   <- "C:/Users/wpypo/Desktop/SCN3A/..."   # ← Modify to the actual path
file_name  <- "SCN3A-dedup.xlsx"
sheet      <- 1

id_col     <- 1
group_col  <- 2
grp_pos    <- "Pathogenic"
grp_neg    <- "Control"

# Set to "manual" mode and provide the specified variable list
mode        <- "manual"           # "auto" / "manual"
manual_vars <- c(
  "SIFT_score", "SIFT4G_score", "Polyphen2_HDIV_score", "Polyphen2_HVAR_score",
  "MutationTaster_score", "MutationAssessor_score", "PROVEAN_score", "MetaSVM_score",
  "MetaLR_score", "MetaRNN_score", "M-CAP_score", "MutPred2_score",
  "PrimateAI_score", "DEOGEN2_score", "BayesDel_addAF_score", "BayesDel_noAF_score",
  "ClinPred_score", "LIST-S2_score", "ESM1b_score", "AlphaMissense_score",
  "fathmm-XF_coding_score"
)

out_file    <- "SCN3A_stats_specified_scores.xlsx"

# Set combined ROC mode to "custom" and provide the same list as manual_vars above
combined_mode <- "custom"    # "all" / "rankscore" / "custom"
custom_combined_vars <- c(
  "SIFT_score", "SIFT4G_score", "Polyphen2_HDIV_score", "Polyphen2_HVAR_score",
  "MutationTaster_score", "MutationAssessor_score", "PROVEAN_score", "MetaSVM_score",
  "MetaLR_score", "MetaRNN_score", "M-CAP_score", "MutPred2_score",
  "PrimateAI_score", "DEOGEN2_score", "BayesDel_addAF_score", "BayesDel_noAF_score",
  "ClinPred_score", "LIST-S2_score", "ESM1b_score", "AlphaMissense_score",
  "fathmm-XF_coding_score"
)

# Number of columns for patchwork plots
patchwork_plot_cols       <- 5   # Number of columns for combined comparison/ROC plots

# Fine-tuning parameters for combined ROC plot
combined_roc_width        <- 8  # Combined ROC plot width (inches)
combined_roc_height       <- 4   # Combined ROC plot height (inches)
combined_roc_legend_cols  <- 2   # Number of columns in combined ROC legend

# Plot output settings
plot_folder <- "plots_specified_scores"
compare_plot_width   <- 3   # Single comparison plot width
compare_plot_height  <- 4   # Single comparison plot height
roc_plot_width       <- 4   # Single ROC plot width
roc_plot_height      <- 4   # Single ROC plot height

# Adjustable font sizes
compare_title_size      <- 12
compare_yaxis_text_size <- 12
## ==============================================================================

## ===== 1. Required Packages ===================================================
pkgs <- c("readxl","dplyr","purrr","pROC","writexl",
          "janitor","ggplot2","tibble","svglite","rlang",
          "scales","patchwork")
suppressPackageStartupMessages({
  for (p in pkgs) {
    if (!requireNamespace(p, quietly = TRUE)) {
      install.packages(p)
    }
    library(p, character.only = TRUE)
  }
})

## ===== 2. Read and Prepare Data ===============================================
df_raw <- tryCatch({
  readxl::read_excel(file.path(data_dir, file_name), sheet = sheet)
}, error = function(e) {
  stop("Unable to read the specified Excel file. Please check the path, file name, and sheet settings:\n", e$message)
})


col_name <- function(idx, df){
  if(is.numeric(idx)) names(df)[idx] else idx
}

df <- df_raw |>
  rename(
    id    = !!sym(col_name(id_col,    df_raw)),
    group = !!sym(col_name(group_col, df_raw))
  ) |>
  mutate(group = factor(group, levels = c(grp_neg, grp_pos)))

if (nlevels(df$group) != 2 || !all(c(grp_neg, grp_pos) %in% levels(df$group))) {
  stop("'group' column could not be set correctly to include '", grp_pos, "' and '", grp_neg, "' as two levels. Please check group_col and group name settings.")
}

## ===== 3. Utility Functions ===================================================
num_summary <- function(v){
  v2 <- na.omit(v)
  c(n = length(v2), mean = mean(v2), sd = sd(v2),
    median = median(v2), iqr = IQR(v2))
}
safe_fname <- function(x) gsub("[^A-Za-z0-9._-]","_",x)

## ===== 4. Main Analysis Function (Including Single ROC Plots) =================
analyse_var <- function(df, var, plot_dir){
  if (!var %in% names(df)) {
    warning("Variable '", var, "' does not exist in the data and has been skipped.")
    return(NULL)
  }
  
  vec <- df[[var]]; grp <- df$group
  n_pos <- sum(!is.na(vec[grp == grp_pos]))
  n_neg <- sum(!is.na(vec[grp == grp_neg]))
  
  if(n_pos < 1 || n_neg < 1 || length(unique(na.omit(vec))) < 2){
    is_numeric <- is.numeric(vec)
    invalid_row <- tibble(
      variable = var, 
      type = if (is_numeric) "invalid_numeric" else "invalid_categorical",
      n_pos = n_pos, 
      n_neg = n_neg
    ) |>
      add_column(!!!setNames(rep(list(NA), 28), c(
        "counts_pos_neg", "test_used", "p_main_raw", "mean_pos", "sd_pos",
        "median_pos", "iqr_pos", "mean_neg", "sd_neg", "median_neg", "iqr_neg",
        "auc", "auc_se", "auc_ci_low", "auc_ci_high", "auc_p_raw",
        "threshold_best", "TP", "TN", "FP", "FN", "accuracy", "sensitivity",
        "specificity", "ppv", "npv", "mcc", "f_score"
      )))
    
    return(list(summary = invalid_row, compare_plot = NULL, roc_plot = NULL))
  }
  
  if(is.numeric(vec)){
    sp <- num_summary(vec[grp == grp_pos])
    sn <- num_summary(vec[grp == grp_neg])
    
    test_used <- "Wilcoxon rank-sum"; tst <- wilcox.test(vec ~ grp, exact = FALSE)
    
    p_lab <- if(tst$p.value < 1e-4) "P < 0.0001" else paste0("P = ", signif(tst$p.value, 3))
    pal <- c("black", "red")
    plt <- ggplot(df, aes(x = group, y = .data[[var]], fill = group, colour = group)) +
      geom_violin(trim = FALSE, alpha = .3, width = .9, colour = NA) +
      geom_boxplot(width = .15, outlier.shape = NA, alpha = .5, colour = "black") +
      geom_jitter(width = .12, size = 1.4, alpha = .7) +
      stat_summary(fun = median, geom = "point", shape = 23, size = 3, fill = "white", colour = "black") +
      scale_fill_manual(values = pal, guide = "none") +
      scale_colour_manual(values = pal, guide = "none") +
      labs(x = NULL, y = var, title = var, subtitle = p_lab) +
      theme_classic(base_size = 13) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title  = element_text(size = compare_title_size),
        axis.title.y = element_text(size = compare_yaxis_text_size)
      )
    if(grepl("rankscore", var, fixed = TRUE)){
      plt <- plt + scale_y_continuous(limits = c(0,1), expand = expansion(mult = c(.02,.02)), oob = scales::oob_keep)
    }
    ggsave(file.path(plot_dir, paste0(safe_fname(var), "_compare.svg")), plt, width = compare_plot_width, height = compare_plot_height)
    
    roc_obj <- suppressWarnings(tryCatch(
      pROC::roc(grp, vec, levels = c(grp_neg, grp_pos), direction = "auto", quiet = TRUE),
      error = function(e) NULL
    ))
    
    if(!is.null(roc_obj)){
      auc_v  <- as.numeric(pROC::auc(roc_obj))
      ci     <- as.numeric(pROC::ci.auc(roc_obj, conf.level = 0.95)[c(1,3)])
      auc_se <- (ci[2] - ci[1]) / (2 * 1.96)
      auc_p  <- suppressWarnings(tryCatch(pROC::roc.test(roc_obj, 0.5, method = "delong")$p.value, error = function(e) NA_real_))
      
      coords_best <- pROC::coords(roc_obj, "best", ret = "all", best.method = "youden")
      thr <- coords_best$threshold
      tp  <- coords_best$tp; tn  <- coords_best$tn; fp  <- coords_best$fp; fn  <- coords_best$fn
      acc <- coords_best$accuracy; se  <- coords_best$sensitivity; spc <- coords_best$specificity
      ppv <- coords_best$ppv; npv <- coords_best$npv
      mcc <- (tp * tn - fp * fn) / sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
      f1  <- 2 * tp / (2 * tp + fp + fn)
      
      roc_plot <- ggroc(roc_obj, colour = "#0072B5", size = 1.1, legacy.axes = TRUE) +
        geom_abline(linetype = "dashed", colour = "grey50") +
        labs(title = paste0("ROC - ", var),
             subtitle = sprintf("AUC = %.3f (95%%CI %.3f–%.3f)", auc_v, ci[1], ci[2]),
             x = "1 – Specificity (FPR)",
             y = "Sensitivity (TPR)") +
        theme_classic(13)
      ggsave(file.path(plot_dir, paste0(safe_fname(var), "_ROC.svg")), roc_plot, width = roc_plot_width, height = roc_plot_height)
      
    } else {
      auc_v <- auc_se <- auc_p <- thr <- tp <- tn <- fp <- fn <- acc <- se <- spc <- ppv <- npv <- mcc <- f1 <- NA_real_
      ci <- c(NA_real_, NA_real_)
      roc_plot <- NULL
    }
    
    summary_row <- tibble(
      variable = var, type = "numeric", n_pos, n_neg, counts_pos_neg = NA_character_,
      test_used, p_main_raw = tst$p.value,
      mean_pos = sp["mean"], sd_pos = sp["sd"], median_pos = sp["median"], iqr_pos = sp["iqr"],
      mean_neg = sn["mean"], sd_neg = sn["sd"], median_neg = sn["median"], iqr_neg = sn["iqr"],
      auc = auc_v, auc_se = auc_se, auc_ci_low = ci[1], auc_ci_high = ci[2], auc_p_raw = auc_p,
      threshold_best = thr, TP = tp, TN = tn, FP = fp, FN = fn,
      accuracy = acc, sensitivity = se, specificity = spc, ppv = ppv, npv = npv, mcc = mcc, f_score = f1
    )
    
    return(list(summary = summary_row, compare_plot = plt, roc_plot = roc_plot))
  }
  
  if (is.character(vec) || is.factor(vec)) {
    tab <- table(vec, grp)
    cps <- paste(rownames(tab), tab[, grp_pos], tab[, grp_neg], sep = "/", collapse = "; ")
    
    if(nrow(tab) == 2 && ncol(tab) == 2){
      tst <- fisher.test(tab); tu <- "Fisher exact"
    } else {
      expv <- suppressWarnings(chisq.test(tab)$expected)
      if(any(expv < 5)){
        tst <- fisher.test(tab, simulate.p.value = TRUE, B = 10000); tu <- "Fisher exact (simulated)"
      } else {
        tst <- chisq.test(tab); tu <- "Chi-square"
      }
    }
    
    summary_row <- tibble(
      variable = var, type = "categorical", n_pos, n_neg, counts_pos_neg = cps,
      test_used = tu, p_main_raw = tst$p.value
    ) |>
      add_column(!!!setNames(rep(list(NA), 25), c(
        "mean_pos", "sd_pos", "median_pos", "iqr_pos", "mean_neg", "sd_neg",
        "median_neg", "iqr_neg", "auc", "auc_se", "auc_ci_low", "auc_ci_high",
        "auc_p_raw", "threshold_best", "TP", "TN", "FP", "FN", "accuracy",
        "sensitivity", "specificity", "ppv", "npv", "mcc", "f_score"
      )))
    return(list(summary = summary_row, compare_plot = NULL, roc_plot = NULL))
  }
  
  return(NULL)
}

## ===== 5. Select Columns to Analyze ===========================================
predictor_cols <- if(tolower(mode) == "manual"){
  manual_vars
} else {
  map_lgl(df, ~ is.numeric(.x) || (is.character(.x) || is.factor(.x)) && n_distinct(na.omit(.x)) < 10) |>
    which() |>
    names() |>
    setdiff(c("id", "group"))
}
message("Analyzing the following ", length(predictor_cols), " predictor variables in '", mode, "' mode: \n", paste(predictor_cols, collapse = ", "))

## ===== 6. Create Plot Directory ===============================================
plot_dir <- file.path(data_dir, plot_folder)
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

## ===== 7. Batch Processing ====================================================
res_list <- purrr::map(predictor_cols, ~ analyse_var(df, .x, plot_dir)) |> purrr::compact()

## ===== 8. Summary Table =======================================================
if (length(res_list) > 0) {
  summary_tbl <- dplyr::bind_rows(purrr::map(res_list, "summary"))
  annotation <- tibble(variable = "variable", description = "Column name of the predictor variable in the original data.") |> bind_rows(tibble(variable = "type", description = "Variable type: numeric (numerical); categorical (categorical); invalid_* (skipped due to insufficient data)."),tibble(variable = c("n_pos","n_neg"), description = "Number of valid observations in the case/control group."),tibble(variable = c("mean_pos","sd_pos","median_pos","iqr_pos", "mean_neg","sd_neg","median_neg","iqr_neg"), description = "Distribution summaries for both groups: mean, standard deviation, median, interquartile range."),tibble(variable = "test_used", description = "Statistical test used: Wilcoxon rank-sum, Chi-square, Fisher exact."),tibble(variable = "p_main_raw", description = "Raw p-value for group comparison."),tibble(variable = c("auc","auc_se","auc_ci_low","auc_ci_high","auc_p_raw"), description = "ROC AUC, standard error, 95% CI lower/upper bound, p-value for AUC vs 0.5 test."),tibble(variable = "threshold_best", description = "Optimal threshold corresponding to the maximum Youden index."),tibble(variable = c("TP","TN","FP","FN"), description = "Confusion matrix counts at the optimal threshold."),tibble(variable = c("accuracy","sensitivity","specificity","ppv","npv","mcc","f_score"), description = "Performance metrics: accuracy, recall (sensitivity), specificity, PPV, NPV, MCC, F1-score."),tibble(variable = "counts_pos_neg", description = "Frequencies of each level in cases/controls for categorical variables, formatted as 'LevelA/case count/control count; ...'."))
  writexl::write_xlsx(list(Results = summary_tbl, Annotation = annotation), path = file.path(data_dir, out_file))
} else { message("No valid analysis results to write to file.") }

## ===== 9. Generate Combined Plots =============================================
if (length(res_list) > 0) {
  plots_all <- purrr::map(res_list, function(x) {list(compare_plot = x$compare_plot, roc_plot = x$roc_plot)}) |> purrr::transpose() |> purrr::map(purrr::compact)
  
  save_combined <- function(plot_list, out_name, single_w, single_h, base_dir, n_cols = 6){
    if(length(plot_list) == 0) return(invisible(NULL))
    ncol_grid <- min(length(plot_list), n_cols)
    nrow_grid <- ceiling(length(plot_list) / ncol_grid)
    
    big_w <- single_w * ncol_grid; big_h <- single_h * nrow_grid
    combined <- patchwork::wrap_plots(plotlist = plot_list, ncol = ncol_grid)
    
    ggsave(filename = file.path(base_dir, out_name), plot = combined, width = big_w, height = big_h, device = "svg", limitsize = FALSE) 
  }
  
  save_combined(plots_all$compare_plot, "Combined_compare_specified.svg", compare_plot_width, compare_plot_height, plot_dir, n_cols = patchwork_plot_cols)
  save_combined(plots_all$roc_plot, "Combined_ROC_specified.svg", roc_plot_width, roc_plot_height, plot_dir, n_cols = patchwork_plot_cols)
}

## ===== 10. Combined ROC Curve =================================================
make_combined_roc <- function(df, mode = "rankscore", custom_vars = NULL, save_dir = ".", 
                              plot_w = 12, plot_h = 8, legend_n_col = 3){
  var_list <- if(tolower(mode) == "all"){ names(df)[sapply(df, is.numeric)]
  } else if(tolower(mode) == "rankscore"){ grep("_rankscore$", names(df), value = TRUE)
  } else if(tolower(mode) == "custom"){ custom_vars
  } else stop("Invalid `combined_mode`")
  
  var_list <- intersect(var_list, names(df))
  if(length(var_list) < 1){ message("No valid specified variables found in 'custom' mode; skipping combined ROC."); return(invisible(NULL)) }
  
  roc_list_with_nulls <- purrr::map(var_list, function(v){
    vec <- df[[v]]
    if(!is.numeric(vec) || length(unique(na.omit(vec))) < 2) return(NULL)
    tryCatch(pROC::roc(df$group, vec, levels = c(grp_neg, grp_pos), quiet = TRUE, direction = "auto"), error = function(e) NULL)
  })
  names(roc_list_with_nulls) <- var_list
  roc_list <- purrr::compact(roc_list_with_nulls)
  
  if(length(roc_list) == 0){ message("Insufficient data to generate any ROC curves; skipping combined ROC."); return(invisible(NULL)) }
  
  # [Modified] Generate a highly distinguishable color set for the combined ROC plot
  n_colors <- length(roc_list)
  # Use the "Zissou 1" palette, which provides rich, easily distinguishable colors
  color_palette <- hcl.colors(n_colors, palette = "Zissou 1")
  
  plt <- ggroc(roc_list, aes = "colour", legacy.axes = TRUE) +
    geom_abline(linetype = "dashed", colour = "grey50") +
    theme_classic(13) +
    labs(title = paste0("Combined ROC Curves (", mode, ")"), x = "1 – Specificity (FPR)", y = "Sensitivity (TPR)", colour = "Variable") +
    guides(colour = guide_legend(ncol = legend_n_col)) +
    # [Modified] Apply custom color scheme
    scale_colour_manual(values = color_palette)
  
  ggsave(file.path(save_dir, paste0("Combined_ROC_", mode, ".svg")),
         plt, width = plot_w, height = plot_h, device = "svg", limitsize = FALSE)
  
  invisible(NULL)
}

make_combined_roc(df,
                  mode = combined_mode,
                  custom_vars = custom_combined_vars,
                  save_dir = plot_dir,
                  plot_w = combined_roc_width,
                  plot_h = combined_roc_height,
                  legend_n_col = combined_roc_legend_cols)

message("\n>>> Analysis completed!")
message("    Statistical results written to: ", file.path(data_dir, out_file))
message("    All plots saved in: ", plot_dir)
