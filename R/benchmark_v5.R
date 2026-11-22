################################################################################
# VEPs benchmark v5
# 整合 ROC、Precision-Recall Curve、Bootstrap 验证、分组比较图、
# dbNSFP 与当前基因 AUC 对照、阈值二值化混淆矩阵及热图分析流程
################################################################################

## ===== 0. 参数区：请优先在这里修改分析设置 ====================================

# [0-1] 自动识别脚本所在目录。
# 若您直接在 RStudio / Terminal 中运行本脚本，会优先使用脚本文件所在目录；
# 如果识别失败，则退回到当前工作目录。
get_script_dir <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_flag <- "--file="
  script_path <- sub(file_flag, "", cmd_args[grepl(file_flag, cmd_args)])
  if (length(script_path) > 0) {
    return(dirname(normalizePath(script_path[1], mustWork = FALSE)))
  }

  if (!is.null(sys.frames()[[1]]$ofile)) {
    return(dirname(normalizePath(sys.frames()[[1]]$ofile, mustWork = FALSE)))
  }

  normalizePath(getwd(), mustWork = FALSE)
}

# 统一 ROC 坐标结构后再合并。原始 roc 的 coords 含 threshold 列，而
# smooth.roc 的 coords 不含该列；若把两者直接交给 pROC::ggroc(list)，
# pROC 内部的 base::rbind 会因列数不同而报错。
combine_roc_coordinate_frames <- function(coord_frames,
                                          display_names = names(coord_frames),
                                          curve_scales = rep(1, length(coord_frames))) {
  if (length(coord_frames) == 0) {
    return(data.frame(FPR = numeric(), TPR = numeric(), Method = character()))
  }
  if (length(display_names) != length(coord_frames) ||
      length(curve_scales) != length(coord_frames)) {
    stop("ROC 坐标、显示名称和坐标尺度的数量必须一致。")
  }

  normalized_frames <- lapply(seq_along(coord_frames), function(i) {
    frame <- coord_frames[[i]]
    required <- c("specificity", "sensitivity")
    missing <- setdiff(required, names(frame))
    if (length(missing) > 0) {
      stop("ROC 坐标缺少列: ", paste(missing, collapse = ", "))
    }

    normalized <- data.frame(
      FPR = curve_scales[[i]] - as.numeric(frame$specificity),
      TPR = as.numeric(frame$sensitivity),
      Method = rep(as.character(display_names[[i]]), nrow(frame)),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    normalized[rev(seq_len(nrow(normalized))), , drop = FALSE]
  })

  do.call(base::rbind, normalized_frames)
}

# 无需加载第三方包的快速回归检查：它先确认旧拼接方式会复现本次错误，
# 再确认统一列结构后的新方式可以安全合并。
if ("--self-test" %in% commandArgs(trailingOnly = TRUE)) {
  raw_coords <- data.frame(
    threshold = c(Inf, 0), specificity = c(1, 0), sensitivity = c(0, 1)
  )
  smooth_coords <- data.frame(
    specificity = c(1, 0), sensitivity = c(0, 1)
  )
  old_result <- try(
    do.call(base::rbind, list(raw = raw_coords, smooth = smooth_coords)),
    silent = TRUE
  )
  stopifnot(inherits(old_result, "try-error"))

  combined <- combine_roc_coordinate_frames(
    list(raw = raw_coords, smooth = smooth_coords),
    display_names = c("raw", "smooth"),
    curve_scales = c(1, 1)
  )
  stopifnot(
    identical(names(combined), c("FPR", "TPR", "Method")),
    nrow(combined) == 4,
    all(combined$FPR >= 0 & combined$FPR <= 1),
    all(combined$TPR >= 0 & combined$TPR <= 1)
  )
  message("SELF_TEST_PASS: mixed raw/smoothed ROC coordinates bind safely")
  quit(save = "no", status = 0)
}

get_cli_option <- function(args, option_name) {
  prefix <- paste0("--", option_name, "=")
  matched <- args[startsWith(args, prefix)]
  if (length(matched) == 0) NULL else sub(prefix, "", matched[[1]], fixed = TRUE)
}

resolve_input_path <- function(script_dir, args) {
  cli_input <- get_cli_option(args, "input")
  positional <- args[!startsWith(args, "--")]
  if (is.null(cli_input) && length(positional) > 0) {
    cli_input <- positional[[1]]
  }

  if (!is.null(cli_input) && nzchar(cli_input)) {
    is_absolute <- grepl("^[A-Za-z]:[\\\\/]|^/", cli_input)
    candidates <- if (is_absolute) cli_input else c(file.path(script_dir, cli_input), cli_input)
    existing <- candidates[file.exists(candidates)]
    if (length(existing) == 0) {
      stop("未找到命令行指定的输入文件: ", cli_input)
    }
    return(normalizePath(existing[[1]], winslash = "/", mustWork = TRUE))
  }

  xlsx_files <- list.files(
    script_dir,
    pattern = "\\.xlsx$",
    full.names = TRUE,
    ignore.case = TRUE
  )
  xlsx_files <- xlsx_files[
    !startsWith(basename(xlsx_files), "~$") &
      !grepl("_benchmark_v5_summary\\.xlsx$", basename(xlsx_files), ignore.case = TRUE)
  ]

  if (length(xlsx_files) == 1) {
    return(normalizePath(xlsx_files[[1]], winslash = "/", mustWork = TRUE))
  }
  if (length(xlsx_files) == 0) {
    stop(
      "脚本目录中没有可用的 .xlsx 输入文件。",
      "请把数据文件与脚本放在同一目录，或使用 --input=完整路径。"
    )
  }
  stop(
    "脚本目录中发现多个 .xlsx 文件，无法自动判断输入文件: ",
    paste(basename(xlsx_files), collapse = ", "),
    "。请使用 --input=文件名 指定。"
  )
}

# [0-2] 输入文件设置。
# 默认自动识别脚本同目录下唯一的 .xlsx，并优先读取含关键列的 "VEPs" 工作表。
# 命令行可覆盖：Rscript "VEPs benchmark v5.R" --input=data.xlsx --sheet=VEPs --gene=GENE
cli_args               <- commandArgs(trailingOnly = TRUE)
work_dir               <- getOption("predictor.rank.root", get_script_dir())
input_path             <- resolve_input_path(work_dir, cli_args)
input_file             <- basename(input_path)
input_sheet            <- get_cli_option(cli_args, "sheet")
if (is.null(input_sheet) || !nzchar(input_sheet)) input_sheet <- "VEPs"
gene_name              <- get_cli_option(cli_args, "gene")
if (is.null(gene_name) || !nzchar(gene_name)) {
  gene_name <- tools::file_path_sans_ext(input_file)
}
if ("--config-test" %in% cli_args) {
  message("CONFIG_TEST_PASS")
  message("input_path=", input_path)
  message("input_sheet=", input_sheet)
  message("gene_name=", gene_name)
  quit(save = "no", status = 0)
}
local_r_library_dir   <- file.path(work_dir, "Rlib")

# [0-3] 分组标签设置。
# positive_labels 中列出的标签都视为阳性 / 致病组；
# negative_labels 中列出的标签都视为阴性 / 对照组。
# display_label 用于图表展示顺序和图例文字。
group_col              <- "Type"
id_col                 <- "variant_id_hg38"
positive_labels        <- c("Disease-associated", "Pathogenic")
negative_labels        <- c("Control", "Benign")
positive_display_label <- "Disease-associated"
negative_display_label <- "Control"

# [0-4] 输出目录设置。
# 所有结果会统一整理到 gene_name_benchmark_v5_results 文件夹中，
# 其中 01_tables 保存 Excel 汇总表，02_plots 保存所有图形结果。
output_root_dir <- file.path(work_dir, paste0(gene_name, "_benchmark_v5_results"))
tables_dir      <- file.path(output_root_dir, "01_tables")
plots_dir       <- file.path(output_root_dir, "02_plots")
per_tool_dir    <- file.path(plots_dir, "per_tool")

# [0-5] 混淆矩阵热图参数。
# heatmap_bar_gap_mm 控制混淆矩阵热图与右侧条形图之间的间距，单位为毫米；
# 默认设置为 0，表示两块图默认不额外留白。
# heat_panel_ratio / bar_panel_ratio 控制热图与条形图的宽度占比。
heatmap_bar_gap_mm           <- 0
heat_panel_ratio             <- 5
bar_panel_ratio              <- 2
heatmap_base_font            <- 10
threshold_heatmap_sort_metric <- "BalancedAccuracy"
threshold_heatmap_sort_desc   <- TRUE
threshold_heatmap_metric_digits <- 3
heatmap_tile_radius_mm        <- 1.2

# [0-6] 单软件图与组合小图参数。
# patchwork_col_count 控制组合单图每行展示几个子图。
# compare_plot_width / height 控制单软件分布比较图尺寸。
# roc_plot_width / height 控制单软件 ROC 图尺寸。
patchwork_col_count          <- 7
compare_plot_width           <- 3.3
compare_plot_height          <- 4.1
roc_plot_width               <- 4.0
roc_plot_height              <- 4.0

# [0-7] 总体 Combined ROC 图参数。
# overall_combined_roc_width / overall_combined_roc_height 控制总体 Combined ROC 图像尺寸。
# overall_combined_roc_legend_position 控制总体 Combined ROC 图例放在右侧还是下方，可选 "right" 或 "bottom"。
# overall_combined_roc_legend_ncol 控制总体 Combined ROC 图例分成几列。
# overall_combined_roc_show_auc_in_legend   控制总体 Combined ROC 图例中是否显示每条曲线的 AUC 值；
# TRUE 表示显示，FALSE 表示仅显示软件名称。
# overall_combined_roc_smooth_curves 控制总体 Combined ROC 中的 ROC 折线是否进行平滑处理；
# TRUE 表示使用 pROC 的 binormal smoothing，FALSE 表示保留原始经验 ROC 折线。
overall_combined_roc_width              <- 10
overall_combined_roc_height             <- 5
overall_combined_roc_legend_position    <- "right"
overall_combined_roc_legend_ncol        <- 2
overall_combined_roc_show_auc_in_legend <- FALSE
overall_combined_roc_smooth_curves      <- TRUE

# [0-8] Precision-Recall Curve 与 Bootstrap 参数。
# PR 曲线会依据阈值方向自动调整分数方向，使“越大越偏阳性”。
# bootstrap_n 控制每个软件的自助法重采样次数；bootstrap_seed 保证结果可复现。
# bootstrap_conf_level 控制 bootstrap 置信区间，默认 95%。
pr_plot_width                 <- 4.0
pr_plot_height                <- 4.0
overall_combined_pr_width     <- 10
overall_combined_pr_height    <- 5
overall_combined_pr_legend_position <- "right"
overall_combined_pr_legend_ncol     <- 2
overall_combined_pr_show_auprc_in_legend <- FALSE
category_combined_pr_width    <- 10
category_combined_pr_height   <- 5
category_combined_pr_legend_position <- "right"
category_combined_pr_legend_ncol     <- 2
category_combined_pr_show_auprc_in_legend <- FALSE
bootstrap_n                   <- 1000
bootstrap_arg <- get_cli_option(cli_args, "bootstrap")
if (!is.null(bootstrap_arg)) {
  bootstrap_n <- suppressWarnings(as.integer(bootstrap_arg))
  if (is.na(bootstrap_n) || bootstrap_n < 0) stop("--bootstrap must be a nonnegative integer")
}
bootstrap_seed                <- 20260606
bootstrap_conf_level          <- 0.95

# [0-9] 按类别 Combined ROC 图参数。
# category_combined_roc_types 控制需要单独输出的类别，类别名称来自下方 tool_config 中的 TypeCategory 列。
# 当前默认分别输出 Meta-predictor、Individual predictor、Conservation 三类。
# 其余参数含义与总体 Combined ROC 相同，但仅作用于按类别拆分的 Combined ROC 图。
category_combined_roc_types              <- c("Meta-predictor", "Individual predictor", "Conservation")
category_combined_roc_width              <- 10
category_combined_roc_height             <- 5
category_combined_roc_legend_position    <- "right"
category_combined_roc_legend_ncol        <- 2
category_combined_roc_show_auc_in_legend <- FALSE
category_combined_roc_smooth_curves      <- TRUE

# [0-10] NPG 配色设置。
# 这里直接写入常用 NPG 调色板颜色，保证脚本不额外依赖 ggsci 包，
# 同时让混淆矩阵热图、分组比较图和 ROC 曲线保持统一风格。
npg_palette <- c(
  "#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF", "#F39B7FFF",
  "#8491B4FF", "#91D1C2FF", "#DC0000FF", "#7E6148FF", "#B09C85FF"
)
group_palette <- stats::setNames(
  c(npg_palette[2], npg_palette[1]),
  c(negative_display_label, positive_display_label)
)
heat_warm_palette <- grDevices::colorRampPalette(c("#FFFFFF", npg_palette[5], npg_palette[1]))(100)
heat_cool_palette <- grDevices::colorRampPalette(c("#FFFFFF", npg_palette[7], npg_palette[4]))(100)
heat_acc_palette  <- grDevices::colorRampPalette(c("#FFFFFF", npg_palette[6], npg_palette[3]))(100)
heat_assoc_palette <- grDevices::colorRampPalette(c("#FFFFFF", npg_palette[2], npg_palette[4]))(100)
stack_palette     <- c("Predicted positive" = npg_palette[8], "Predicted negative" = npg_palette[4])

# [0-11] dbNSFP performance matrix 设置。
# dbNSFP 矩阵只用于与当前输入数据集的 ROC AUC 做对照：
# 1) 历年 dbNSFP AUC 与当前 benchmark AUC 对照；
# 2) 最新年份 dbNSFP AUC 与当前 benchmark AUC 的差异图。
dbnsfp_performance_file <- "data/reference/dbnsfp_performance_matrix.csv"
dbnsfp_plot_width       <- 8
dbnsfp_plot_height      <- 8

# 从当前分析目录开始逐级向上查找项目级文件。
# 这样即使输入文件或脚本位于项目子目录，也可以自动识别项目根目录中的
# dbnsfp_performance_matrix.csv；若未找到，则返回原始候选路径并沿用后续跳过逻辑。
find_file_upwards <- function(filename, start_dir = work_dir) {
  current_dir <- normalizePath(start_dir, winslash = "/", mustWork = FALSE)

  repeat {
    candidate_path <- file.path(current_dir, filename)
    if (file.exists(candidate_path)) {
      return(normalizePath(candidate_path, winslash = "/", mustWork = TRUE))
    }

    parent_dir <- dirname(current_dir)
    if (identical(parent_dir, current_dir)) {
      break
    }
    current_dir <- parent_dir
  }

  file.path(normalizePath(start_dir, winslash = "/", mustWork = FALSE), filename)
}

## ===== 1. 载入依赖包 ==========================================================
# 必须先注册脚本同目录的 Rlib，再创建依赖 tibble 的 tool_config；否则移动目录后
# 即使依赖包已经位于本地 Rlib，也可能在这里提前报“没有叫 tibble 的程序包”。
required_pkgs <- c(
  "readxl", "dplyr", "tidyr", "purrr", "pROC", "writexl",
  "ggplot2", "tibble", "rlang", "scales", "patchwork", "ggnewscale", "svglite"
)

if (dir.exists(local_r_library_dir)) {
  .libPaths(c(normalizePath(local_r_library_dir, mustWork = TRUE), .libPaths()))
}

suppressPackageStartupMessages({
  for (pkg in required_pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("缺少 R 包: ", pkg, "。请先安装后再运行本脚本。")
    }
    library(pkg, character.only = TRUE)
  }
})

# [0-12] 软件阈值设置区。
# 这里不再读取 *_pred 的二分类结果，而是直接使用连续型 score / phred 列，
# 依据“软件预测阈值.xlsx”中的 tool 与 cutoff 手动前置到代码中，
# 并将软件类别信息直接写入 TypeCategory 列，便于后续手动调整
# Meta-predictor、Individual predictor、Conservation 的分类绘图范围。
#
# 字段说明：
# Method          : 结果表和图形中显示的软件名称
# TypeCategory    : 软件所属类别，供按类别绘制 Combined ROC 使用
# ScoreColumn     : 数据源中实际用于分析的连续型列名
# CutoffOperator  : 判定方向；">" 表示高于阈值判为致病，"<" 表示低于阈值判为致病
# CutoffValue     : 具体阈值数值
# ThresholdText   : 便于在结果表和图标题中直接展示的阈值字符串
#
# 注意：
# 1. 仅纳入 输入的.xlsx文件中存在且满足“score / phred 连续值 + cutoff”条件的软件。
# 2. 仅有 rankscore 或 raw、但没有对应 score / phred cutoff 的工具不纳入本次混淆矩阵。
tool_config <- tibble::tribble(
  ~Method,              ~TypeCategory,           ~ScoreColumn,               ~CutoffOperator, ~CutoffValue,  ~ThresholdText,
  "SIFT",               "Individual predictor",  "SIFT_score",               "<",              0.05,          "<0.05",
  "SIFT4G",             "Individual predictor",  "SIFT4G_score",             "<",              0.05,          "<0.05",
  "Polyphen2_HDIV",     "Individual predictor",  "Polyphen2_HDIV_score",     ">",              0.5,           ">0.5",
  "Polyphen2_HVAR",     "Individual predictor",  "Polyphen2_HVAR_score",     ">",              0.5,           ">0.5",
  "MutationTaster",     "Individual predictor",  "MutationTaster_score",     ">",              0.5,           ">0.5",
  "MutationAssessor",   "Individual predictor",  "MutationAssessor_score",   ">",              1.935,         ">1.935",
  "PROVEAN",            "Individual predictor",  "PROVEAN_score",            "<",             -2.5,           "<-2.5",
  "VEST4",              "Individual predictor",  "VEST4_score",              ">",              0.5,           ">0.5",
  "MetaSVM",            "Meta-predictor",        "MetaSVM_score",            ">",              0.0,           ">0",
  "MetaLR",             "Meta-predictor",        "MetaLR_score",             ">",              0.5,           ">0.5",
  "MetaRNN",            "Meta-predictor",        "MetaRNN_score",            ">",              0.5,           ">0.5",
  "M-CAP",              "Meta-predictor",        "M-CAP_score",              ">",              0.025,         ">0.025",
  "REVEL",              "Meta-predictor",        "REVEL_score",              ">",              0.5,           ">0.5",
  "MutPred2",           "Individual predictor",  "MutPred2_score",           ">",              0.737,         ">0.737",
  "MVP",                "Meta-predictor",        "MVP_score",                ">",              0.7,           ">0.7",
  "gMVP",               "Meta-predictor",        "gMVP_rankscore",           ">",              0.75,          ">0.75", #不提供原始分数，考虑到尺度一致，使用rankscore分数代替
  #"MPC",                "Meta-predictor",        "MPC_rankscore",            ">",              2.0,           ">2",    #MPC的官方阈值是>2，但是dbNSFP并未提供原始分数，因此采用rankscore,由于尺度不同停止计算
  "PrimateAI",          "Individual predictor",  "PrimateAI_score",          ">",              0.803,         ">0.803",
  "DEOGEN2",            "Individual predictor",  "DEOGEN2_score",            ">",              0.5,           ">0.5",
  "BayesDel_addAF",     "Meta-predictor",        "BayesDel_addAF_score",     ">",              0.0692655,     ">0.0692655",
  "BayesDel_noAF",      "Meta-predictor",        "BayesDel_noAF_score",      ">",             -0.0570105,     ">-0.0570105",
  "ClinPred",           "Meta-predictor",        "ClinPred_score",           ">",              0.5,           ">0.5",
  "LIST-S2",            "Meta-predictor",        "LIST-S2_rankscore",        ">",              0.85,          ">0.85", #不提供原始分数，使用rankscore分数代替
  "VARITY_R",           "Meta-predictor",        "VARITY_R_score",           ">",              0.5,           ">0.5",
  "VARITY_ER",          "Meta-predictor",        "VARITY_ER_score",          ">",              0.5,           ">0.5",
  "VARITY_R_LOO",       "Meta-predictor",        "VARITY_R_LOO_score",       ">",              0.5,           ">0.5",
  "VARITY_ER_LOO",      "Meta-predictor",        "VARITY_ER_LOO_score",      ">",              0.5,           ">0.5",
  "ESM1b",              "Individual predictor",  "ESM1b_score",              "<",             -7.5,           "<-7.5",
  "AlphaMissense",      "Individual predictor",  "AlphaMissense_rankscore",  ">",              0.564,         ">0.564", #不提供原始分数，使用rankscore分数代替
  "PHACTboost",         "Meta-predictor",        "PHACTboost_score",         ">",              0.62,          ">0.62",
  "MutFormer",          "Meta-predictor",        "MutFormer_score",          ">",              0.8838,        ">0.8838",
  "MutScore",           "Meta-predictor",        "MutScore_score",           ">",              0.5,           ">0.5",
  #"popEVE",             "Individual predictor",  "popEVE_converted_rankscore", "<",            -5.056,         "<-5.056", #popEVE的官方阈值是<-5.056，但是dbNSFP并未提供原始分数，因此采用rankscore,由于尺度不同停止计算
  "CADD_phred",         "Meta-predictor",        "CADD_phred",               ">",             20.0,           ">20",
  "DANN",               "Meta-predictor",        "DANN_score",               ">",              0.96,          ">0.96",
  "fathmm-XF",          "Individual predictor",  "fathmm-XF_coding_score",   ">",              0.5,           ">0.5"
)

## ===== 2. 工具函数 ============================================================
safe_fname <- function(x) {
  gsub("[^A-Za-z0-9._-]", "_", x)
}

ensure_dir <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
}

copy_to_ascii_temp <- function(src_path) {
  safe_name <- safe_fname(basename(src_path))
  temp_path <- file.path(tempdir(), paste0(format(Sys.time(), "%Y%m%d%H%M%S"), "_", safe_name))
  ok <- file.copy(src_path, temp_path, overwrite = TRUE)
  if (!ok) {
    stop("无法将输入文件复制到临时目录: ", src_path)
  }
  temp_path
}

read_excel_safe <- function(path, sheet = 1, ...) {
  temp_path <- copy_to_ascii_temp(path)
  readxl::read_excel(temp_path, sheet = sheet, ...)
}

resolve_input_sheet <- function(path, preferred_sheet, required_columns) {
  temp_path <- copy_to_ascii_temp(path)
  sheet_names <- readxl::excel_sheets(temp_path)
  if (length(sheet_names) == 0) {
    stop("输入 Excel 不包含任何工作表: ", path)
  }

  header_names <- lapply(sheet_names, function(sheet_name) {
    names(readxl::read_excel(temp_path, sheet = sheet_name, n_max = 0))
  })
  has_required <- vapply(
    header_names,
    function(cols) all(required_columns %in% cols),
    logical(1)
  )

  if (preferred_sheet %in% sheet_names) {
    preferred_index <- match(preferred_sheet, sheet_names)
    if (has_required[[preferred_index]]) {
      return(preferred_sheet)
    }
    message(
      ">>> 指定工作表 ", preferred_sheet,
      " 缺少关键列，将自动查找可用工作表。"
    )
  }

  eligible <- which(has_required)
  if (length(eligible) == 0) {
    stop(
      "没有工作表同时包含关键列: ", paste(required_columns, collapse = ", "),
      "。现有工作表: ", paste(sheet_names, collapse = ", ")
    )
  }

  eligible_names <- sheet_names[eligible]
  scores <- integer(length(eligible_names))
  scores <- scores + ifelse(eligible_names == "VEPs", 100L, 0L)
  scores <- scores + ifelse(grepl("^VEPs", eligible_names, ignore.case = TRUE), 10L, 0L)
  scores <- scores - ifelse(grepl("未排除|raw|原始", eligible_names, ignore.case = TRUE), 5L, 0L)
  best <- which(scores == max(scores))
  if (length(best) > 1) {
    stop(
      "多个工作表都包含所需关键列，无法自动选择: ",
      paste(eligible_names[best], collapse = ", "),
      "。请使用 --sheet=工作表名 指定。"
    )
  }

  selected <- eligible_names[[best]]
  message(">>> 自动选择输入工作表: ", selected)
  selected
}

map_truth_label <- function(x) {
  dplyr::case_when(
    x %in% positive_labels ~ 1,
    x %in% negative_labels ~ 0,
    TRUE ~ NA_real_
  )
}

map_group_display <- function(y_true) {
  dplyr::case_when(
    y_true == 1 ~ positive_display_label,
    y_true == 0 ~ negative_display_label,
    TRUE ~ NA_character_
  )
}

apply_cutoff <- function(score, operator, cutoff) {
  if (operator == ">") {
    return(ifelse(is.na(score), NA_real_, ifelse(score > cutoff, 1, 0)))
  }
  if (operator == ">=") {
    return(ifelse(is.na(score), NA_real_, ifelse(score >= cutoff, 1, 0)))
  }
  if (operator == "<") {
    return(ifelse(is.na(score), NA_real_, ifelse(score < cutoff, 1, 0)))
  }
  if (operator == "<=") {
    return(ifelse(is.na(score), NA_real_, ifelse(score <= cutoff, 1, 0)))
  }
  stop("不支持的阈值方向: ", operator)
}

num_summary <- function(v) {
  v2 <- stats::na.omit(v)
  c(
    n = length(v2),
    mean = if (length(v2) > 0) mean(v2) else NA_real_,
    sd = if (length(v2) > 1) stats::sd(v2) else NA_real_,
    median = if (length(v2) > 0) stats::median(v2) else NA_real_,
    iqr = if (length(v2) > 0) stats::IQR(v2) else NA_real_
  )
}

scale_01 <- function(x) {
  if (all(is.na(x))) {
    return(rep(NA_real_, length(x)))
  }
  rng <- range(x, na.rm = TRUE)
  if (diff(rng) == 0) {
    return(rep(0.5, length(x)))
  }
  (x - rng[1]) / diff(rng)
}

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

stroke_pt <- 72.27 / 25.4

GeomRoundedTile <- ggplot2::ggproto(
  "GeomRoundedTile",
  ggplot2::GeomTile,
  draw_panel = function(data, panel_params, coord, radius = grid::unit(1.2, "mm")) {
    coords <- coord$transform(data, panel_params)
    if (nrow(coords) == 0) {
      return(grid::nullGrob())
    }

    grobs <- lapply(seq_len(nrow(coords)), function(i) {
      grid::roundrectGrob(
        x = (coords$xmin[i] + coords$xmax[i]) / 2,
        y = (coords$ymin[i] + coords$ymax[i]) / 2,
        width = abs(coords$xmax[i] - coords$xmin[i]),
        height = abs(coords$ymax[i] - coords$ymin[i]),
        r = radius,
        default.units = "native",
        gp = grid::gpar(
          col = if (is.na(coords$colour[i])) NA else scales::alpha(coords$colour[i], coords$alpha[i]),
          fill = scales::alpha(coords$fill[i], coords$alpha[i]),
          lwd = (coords$linewidth[i] %||% 0.5) * stroke_pt,
          lty = coords$linetype[i]
        )
      )
    })

    do.call(grid::grobTree, grobs)
  }
)

geom_rounded_tile <- function(mapping = NULL, data = NULL, stat = "identity",
                              position = "identity", ..., radius = grid::unit(1.2, "mm"),
                              na.rm = FALSE, show.legend = NA, inherit.aes = TRUE) {
  ggplot2::layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomRoundedTile,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(radius = radius, na.rm = na.rm, ...)
  )
}

calc_binary_metrics_raw <- function(y_true, y_pred) {
  tp <- sum(y_true == 1 & y_pred == 1)
  tn <- sum(y_true == 0 & y_pred == 0)
  fp <- sum(y_true == 0 & y_pred == 1)
  fn <- sum(y_true == 1 & y_pred == 0)

  acc  <- (tp + tn) / (tp + tn + fp + fn)
  sens <- ifelse(tp + fn == 0, NA_real_, tp / (tp + fn))
  spec <- ifelse(tn + fp == 0, NA_real_, tn / (tn + fp))
  bal  <- ifelse(is.na(sens) | is.na(spec), NA_real_, (sens + spec) / 2)
  ppv  <- ifelse(tp + fp == 0, NA_real_, tp / (tp + fp))
  npv  <- ifelse(tn + fn == 0, NA_real_, tn / (tn + fn))

  mcc_denom <- sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
  mcc <- ifelse(mcc_denom == 0, NA_real_, (tp * tn - fp * fn) / mcc_denom)
  f1  <- ifelse(2 * tp + fp + fn == 0, NA_real_, 2 * tp / (2 * tp + fp + fn))

  tibble::tibble(
    TP = tp,
    TN = tn,
    FP = fp,
    FN = fn,
    Accuracy = acc,
    BalancedAccuracy = bal,
    Sensitivity = sens,
    Specificity = spec,
    PPV = ppv,
    NPV = npv,
    MCC = mcc,
    F1_score = f1
  )
}

calc_binary_metrics <- function(y_true, y_pred) {
  calc_binary_metrics_raw(y_true, y_pred) %>%
    dplyr::mutate(
      dplyr::across(
        c(Accuracy, BalancedAccuracy, Sensitivity, Specificity, PPV, NPV, MCC, F1_score),
        ~ round(.x, 4)
      )
    )
}

orient_score_for_positive <- function(score, operator) {
  if (operator %in% c(">", ">=")) {
    return(score)
  }
  if (operator %in% c("<", "<=")) {
    return(-score)
  }
  stop("不支持的阈值方向: ", operator)
}

calc_pr_curve <- function(y_true, predictor) {
  pr_df <- tibble::tibble(
    truth = as.integer(y_true),
    predictor = as.numeric(predictor)
  ) %>%
    dplyr::filter(!is.na(truth), !is.na(predictor))

  n_positive <- sum(pr_df$truth == 1)
  n_negative <- sum(pr_df$truth == 0)

  if (nrow(pr_df) == 0 || n_positive == 0 || n_negative == 0 || dplyr::n_distinct(pr_df$predictor) < 2) {
    return(list(points = NULL, summary = NULL))
  }

  threshold_tbl <- pr_df %>%
    dplyr::arrange(dplyr::desc(predictor)) %>%
    dplyr::group_by(predictor) %>%
    dplyr::summarise(
      PositivesAtThreshold = sum(truth == 1),
      NegativesAtThreshold = sum(truth == 0),
      .groups = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(predictor)) %>%
    dplyr::mutate(
      TP = cumsum(PositivesAtThreshold),
      FP = cumsum(NegativesAtThreshold),
      Recall = TP / n_positive,
      Precision = TP / (TP + FP)
    )

  pr_points <- dplyr::bind_rows(
    tibble::tibble(
      predictor = Inf,
      PositivesAtThreshold = 0,
      NegativesAtThreshold = 0,
      TP = 0,
      FP = 0,
      Recall = 0,
      Precision = 1
    ),
    threshold_tbl
  )

  recall_diff <- diff(pr_points$Recall)
  precision_mid <- (utils::head(pr_points$Precision, -1) + utils::tail(pr_points$Precision, -1)) / 2
  auprc_trapezoid <- sum(recall_diff * precision_mid, na.rm = TRUE)
  average_precision <- sum(recall_diff * utils::tail(pr_points$Precision, -1), na.rm = TRUE)
  prevalence <- n_positive / (n_positive + n_negative)

  list(
    points = pr_points,
    summary = tibble::tibble(
      AUPRC = round(auprc_trapezoid, 4),
      AveragePrecision = round(average_precision, 4),
      PR_baseline = round(prevalence, 4)
    )
  )
}

make_single_pr <- function(df_tool, method, operator) {
  oriented_score <- orient_score_for_positive(df_tool$score, operator)
  pr_res <- calc_pr_curve(df_tool$truth, oriented_score)

  if (is.null(pr_res$points)) {
    return(list(points = NULL, plot = NULL, summary = NULL))
  }

  auprc <- pr_res$summary$AUPRC
  ap <- pr_res$summary$AveragePrecision
  baseline <- pr_res$summary$PR_baseline

  pr_plot <- ggplot2::ggplot(pr_res$points, ggplot2::aes(x = Recall, y = Precision)) +
    ggplot2::geom_step(colour = npg_palette[3], linewidth = 1.1, direction = "hv") +
    ggplot2::geom_hline(yintercept = baseline, linetype = "dashed", colour = "grey55") +
    ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    ggplot2::labs(
      title = paste0("Precision-Recall - ", method),
      subtitle = sprintf("AUPRC = %.3f | AP = %.3f | baseline = %.3f", auprc, ap, baseline),
      x = "Recall",
      y = "Precision"
    ) +
    ggplot2::theme_classic(base_size = 13) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 12),
      plot.subtitle = ggplot2::element_text(size = 10)
    )

  list(points = pr_res$points, plot = pr_plot, summary = pr_res$summary)
}

summarise_bootstrap_vector <- function(x, prefix, conf_level) {
  x <- x[!is.na(x)]
  alpha <- (1 - conf_level) / 2
  low_prob <- alpha
  high_prob <- 1 - alpha

  tibble::tibble(
    Metric = prefix,
    BootstrapN_valid = length(x),
    Mean = if (length(x) > 0) mean(x) else NA_real_,
    SD = if (length(x) > 1) stats::sd(x) else NA_real_,
    CI_low = if (length(x) > 0) as.numeric(stats::quantile(x, low_prob, names = FALSE, na.rm = TRUE)) else NA_real_,
    CI_high = if (length(x) > 0) as.numeric(stats::quantile(x, high_prob, names = FALSE, na.rm = TRUE)) else NA_real_
  )
}

calc_bootstrap_validation <- function(df_tool, operator, cutoff, n_boot, conf_level) {
  if (n_boot <= 0) {
    return(NULL)
  }

  n <- nrow(df_tool)
  if (n == 0) {
    return(NULL)
  }

  metric_names <- c(
    "AUC", "AUPRC", "AveragePrecision", "Accuracy", "BalancedAccuracy",
    "Sensitivity", "Specificity", "PPV", "NPV", "MCC", "F1_score"
  )

  boot_mat <- matrix(NA_real_, nrow = n_boot, ncol = length(metric_names))
  colnames(boot_mat) <- metric_names
  oriented_score <- orient_score_for_positive(df_tool$score, operator)

  for (b in seq_len(n_boot)) {
    idx <- sample.int(n, size = n, replace = TRUE)
    y_b <- df_tool$truth[idx]
    score_b <- df_tool$score[idx]
    oriented_b <- oriented_score[idx]
    pred_b <- apply_cutoff(score_b, operator, cutoff)

    binary_b <- calc_binary_metrics_raw(y_b, pred_b)
    boot_mat[b, c("Accuracy", "BalancedAccuracy", "Sensitivity", "Specificity", "PPV", "NPV", "MCC", "F1_score")] <-
      as.numeric(binary_b[1, c("Accuracy", "BalancedAccuracy", "Sensitivity", "Specificity", "PPV", "NPV", "MCC", "F1_score")])

    if (length(unique(stats::na.omit(y_b))) == 2 && length(unique(stats::na.omit(score_b))) >= 2) {
      roc_b <- suppressWarnings(
        tryCatch(
          pROC::roc(
            response = factor(
              dplyr::if_else(y_b == 1, positive_display_label, negative_display_label),
              levels = c(negative_display_label, positive_display_label)
            ),
            predictor = score_b,
            levels = c(negative_display_label, positive_display_label),
            direction = "auto",
            quiet = TRUE
          ),
          error = function(e) NULL
        )
      )
      if (!is.null(roc_b)) {
        boot_mat[b, "AUC"] <- as.numeric(pROC::auc(roc_b))
      }

      pr_b <- calc_pr_curve(y_b, oriented_b)
      if (!is.null(pr_b$summary)) {
        boot_mat[b, "AUPRC"] <- pr_b$summary$AUPRC
        boot_mat[b, "AveragePrecision"] <- pr_b$summary$AveragePrecision
      }
    }
  }

  purrr::map_dfr(metric_names, ~ summarise_bootstrap_vector(boot_mat[, .x], .x, conf_level)) %>%
    dplyr::mutate(
      Mean = round(Mean, 4),
      SD = round(SD, 4),
      CI_low = round(CI_low, 4),
      CI_high = round(CI_high, 4)
    )
}

format_p_value <- function(p_value) {
  if (is.na(p_value)) {
    return("P = NA")
  }
  if (p_value < 1e-4) {
    return("P < 0.0001")
  }
  paste0("P = ", signif(p_value, 3))
}

save_plot_if_not_null <- function(plot_obj, path, width, height) {
  if (is.null(plot_obj)) {
    return(invisible(NULL))
  }
  ggplot2::ggsave(
    filename = path,
    plot = plot_obj,
    width = width,
    height = height,
    device = "svg",
    limitsize = FALSE
  )
}

make_compare_plot <- function(df_tool, method, score_col, threshold_text) {
  p_value <- tryCatch(
    stats::wilcox.test(score ~ group, data = df_tool, exact = FALSE)$p.value,
    error = function(e) NA_real_
  )

  plt <- ggplot2::ggplot(
    df_tool,
    ggplot2::aes(x = group, y = score, fill = group, colour = group)
  ) +
    ggplot2::geom_violin(trim = FALSE, alpha = 0.28, width = 0.92, colour = NA) +
    ggplot2::geom_boxplot(width = 0.16, outlier.shape = NA, alpha = 0.55, colour = "black") +
    ggplot2::geom_jitter(width = 0.12, size = 1.5, alpha = 0.75) +
    ggplot2::stat_summary(
      fun = stats::median,
      geom = "point",
      shape = 23,
      size = 3,
      fill = "white",
      colour = "black"
    ) +
    ggplot2::scale_fill_manual(values = group_palette, guide = "none") +
    ggplot2::scale_colour_manual(values = group_palette, guide = "none") +
    ggplot2::labs(
      x = NULL,
      y = score_col,
      title = method,
      subtitle = paste0("Threshold ", threshold_text, " | ", format_p_value(p_value))
    ) +
    ggplot2::theme_classic(base_size = 13) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 35, hjust = 1),
      plot.title = ggplot2::element_text(face = "bold", size = 12),
      plot.subtitle = ggplot2::element_text(size = 10)
    )

  list(plot = plt, p_value = p_value)
}

make_single_roc <- function(df_tool, method) {
  if (length(unique(stats::na.omit(df_tool$score))) < 2) {
    return(list(roc_obj = NULL, plot = NULL, summary = NULL))
  }

  roc_obj <- suppressWarnings(
    tryCatch(
      pROC::roc(
        response = df_tool$group,
        predictor = df_tool$score,
        levels = c(negative_display_label, positive_display_label),
        direction = "auto",
        quiet = TRUE
      ),
      error = function(e) NULL
    )
  )

  if (is.null(roc_obj)) {
    return(list(roc_obj = NULL, plot = NULL, summary = NULL))
  }

  auc_value <- as.numeric(pROC::auc(roc_obj))
  auc_ci    <- as.numeric(pROC::ci.auc(roc_obj, conf.level = 0.95)[c(1, 3)])
  best_info <- pROC::coords(
    roc_obj,
    x = "best",
    ret = c("threshold", "accuracy", "sensitivity", "specificity"),
    best.method = "youden",
    transpose = FALSE
  )

  roc_plot <- pROC::ggroc(roc_obj, colour = npg_palette[1], size = 1.1, legacy.axes = TRUE) +
    ggplot2::geom_abline(linetype = "dashed", colour = "grey55") +
    ggplot2::labs(
      title = paste0("ROC - ", method),
      subtitle = sprintf("AUC = %.3f (95%% CI %.3f-%.3f)", auc_value, auc_ci[1], auc_ci[2]),
      x = "1 - Specificity (FPR)",
      y = "Sensitivity (TPR)"
    ) +
    ggplot2::theme_classic(base_size = 13) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = 12),
      plot.subtitle = ggplot2::element_text(size = 10)
    )

  roc_summary <- tibble::tibble(
    AUC = round(auc_value, 4),
    AUC_low = round(auc_ci[1], 4),
    AUC_high = round(auc_ci[2], 4),
    ROC_best_threshold = round(as.numeric(best_info["threshold"]), 6),
    ROC_best_accuracy = round(as.numeric(best_info["accuracy"]), 4),
    ROC_best_sensitivity = round(as.numeric(best_info["sensitivity"]), 4),
    ROC_best_specificity = round(as.numeric(best_info["specificity"]), 4)
  )

  list(roc_obj = roc_obj, plot = roc_plot, summary = roc_summary)
}

analyse_tool <- function(df, config_row, save_dir) {
  method <- config_row$Method
  type_category <- config_row$TypeCategory
  score_col <- config_row$ScoreColumn
  operator <- config_row$CutoffOperator
  cutoff <- config_row$CutoffValue
  threshold_text <- config_row$ThresholdText

  if (!score_col %in% names(df)) {
    return(list(
      classification = NULL,
      roc_summary = NULL,
      bootstrap_summary = NULL,
      compare_plot = NULL,
      roc_plot = NULL,
      pr_plot = NULL,
      roc_obj = NULL,
      pr_points = NULL,
      skipped = tibble::tibble(Method = method, ScoreColumn = score_col, Reason = "Score/phred column not found in input file")
    ))
  }

  score_vec <- suppressWarnings(as.numeric(df[[score_col]]))
  keep_idx <- !is.na(df$truth) & !is.na(score_vec)

  if (sum(keep_idx) == 0) {
    return(list(
      classification = NULL,
      roc_summary = NULL,
      bootstrap_summary = NULL,
      compare_plot = NULL,
      roc_plot = NULL,
      pr_plot = NULL,
      roc_obj = NULL,
      pr_points = NULL,
      skipped = tibble::tibble(Method = method, ScoreColumn = score_col, Reason = "No usable score/phred values after NA filtering")
    ))
  }

  df_tool <- df[keep_idx, , drop = FALSE] %>%
    dplyr::mutate(score = score_vec[keep_idx])

  y_true <- df_tool$truth
  y_pred <- apply_cutoff(df_tool$score, operator, cutoff)

  class_metrics <- calc_binary_metrics(y_true, y_pred) %>%
    dplyr::mutate(
      Method = method,
      TypeCategory = type_category,
      ScoreColumn = score_col,
      ThresholdRule = threshold_text,
      N_used = length(y_true),
      N_positive = sum(y_true == 1),
      N_negative = sum(y_true == 0),
      .before = 1
    )

  compare_res <- make_compare_plot(df_tool, method, score_col, threshold_text)
  roc_res <- make_single_roc(df_tool, method)
  pr_res <- make_single_pr(df_tool, method, operator)

  stats_pos <- num_summary(df_tool$score[df_tool$group == positive_display_label])
  stats_neg <- num_summary(df_tool$score[df_tool$group == negative_display_label])

  roc_summary <- tibble::tibble(
    Method = method,
    TypeCategory = type_category,
    ScoreColumn = score_col,
    ThresholdRule = threshold_text,
    N_used = nrow(df_tool),
    N_positive = sum(df_tool$group == positive_display_label),
    N_negative = sum(df_tool$group == negative_display_label),
    P_value = round(compare_res$p_value, 6),
    Mean_positive = round(stats_pos["mean"], 6),
    Mean_negative = round(stats_neg["mean"], 6),
    Median_positive = round(stats_pos["median"], 6),
    Median_negative = round(stats_neg["median"], 6)
  )

  if (!is.null(roc_res$summary)) {
    roc_summary <- dplyr::bind_cols(roc_summary, roc_res$summary)
  } else {
    roc_summary <- dplyr::bind_cols(
      roc_summary,
      tibble::tibble(
        AUC = NA_real_,
        AUC_low = NA_real_,
        AUC_high = NA_real_,
        ROC_best_threshold = NA_real_,
        ROC_best_accuracy = NA_real_,
        ROC_best_sensitivity = NA_real_,
        ROC_best_specificity = NA_real_
      )
    )
  }

  if (!is.null(pr_res$summary)) {
    roc_summary <- dplyr::bind_cols(roc_summary, pr_res$summary)
  } else {
    roc_summary <- dplyr::bind_cols(
      roc_summary,
      tibble::tibble(
        AUPRC = NA_real_,
        AveragePrecision = NA_real_,
        PR_baseline = round(mean(df_tool$truth == 1), 4)
      )
    )
  }

  bootstrap_tbl <- calc_bootstrap_validation(
    df_tool = df_tool,
    operator = operator,
    cutoff = cutoff,
    n_boot = bootstrap_n,
    conf_level = bootstrap_conf_level
  )

  if (!is.null(bootstrap_tbl)) {
    bootstrap_tbl <- bootstrap_tbl %>%
      dplyr::mutate(
        Method = method,
        TypeCategory = type_category,
        ScoreColumn = score_col,
        ThresholdRule = threshold_text,
        BootstrapN_requested = bootstrap_n,
        ConfidenceLevel = bootstrap_conf_level,
        .before = 1
      )
  }

  compare_path <- file.path(save_dir, paste0(safe_fname(method), "_compare.svg"))
  roc_path     <- file.path(save_dir, paste0(safe_fname(method), "_ROC.svg"))
  pr_path      <- file.path(save_dir, paste0(safe_fname(method), "_PR.svg"))
  save_plot_if_not_null(compare_res$plot, compare_path, compare_plot_width, compare_plot_height)
  save_plot_if_not_null(roc_res$plot, roc_path, roc_plot_width, roc_plot_height)
  save_plot_if_not_null(pr_res$plot, pr_path, pr_plot_width, pr_plot_height)

  list(
    classification = class_metrics,
    roc_summary = roc_summary,
    bootstrap_summary = bootstrap_tbl,
    compare_plot = compare_res$plot,
    roc_plot = roc_res$plot,
    pr_plot = pr_res$plot,
    roc_obj = roc_res$roc_obj,
    pr_points = pr_res$points,
    skipped = NULL
  )
}

make_threshold_heatmap <- function(classification_tbl, save_path) {
  sort_metric <- threshold_heatmap_sort_metric
  if (!sort_metric %in% names(classification_tbl)) {
    stop("threshold_heatmap_sort_metric not found in classification table: ", sort_metric)
  }

  arrange_fun <- if (isTRUE(threshold_heatmap_sort_desc)) dplyr::desc else identity

  plot_tbl <- classification_tbl %>%
    dplyr::mutate(
      PredictedPositive = TP + FP,
      PredictedNegative = TN + FN
    ) %>%
    dplyr::filter(!is.na(.data[[sort_metric]])) %>%
    dplyr::arrange(arrange_fun(.data[[sort_metric]]), dplyr::desc(Accuracy))

  if (nrow(plot_tbl) == 0) {
    return(invisible(NULL))
  }

  plot_tbl$Method <- factor(plot_tbl$Method, levels = rev(plot_tbl$Method))
  metric_label_fmt <- paste0("%.", threshold_heatmap_metric_digits, "f")

  tile_width <- 2
  within_pair_gap <- 0
  between_pair_gap <- 1

  centers <- numeric(7)
  centers[1] <- 0
  centers[2] <- centers[1] + tile_width + within_pair_gap
  centers[3] <- centers[2] + tile_width + between_pair_gap
  centers[4] <- centers[3] + tile_width + within_pair_gap
  centers[5] <- centers[4] + tile_width + between_pair_gap
  centers[6] <- centers[5] + tile_width + within_pair_gap
  centers[7] <- centers[6] + tile_width + within_pair_gap
  x_pos <- stats::setNames(centers, c("FN", "TP", "FP", "TN", "BalancedAccuracy", "MCC", "F1_score"))

  heat_tbl <- plot_tbl %>%
    dplyr::select(Method, FN, TP, FP, TN, BalancedAccuracy, MCC, F1_score) %>%
    tidyr::pivot_longer(-Method, names_to = "Metric", values_to = "Value") %>%
    dplyr::mutate(
      X = x_pos[Metric],
      Group = dplyr::case_when(
        Metric %in% c("FN", "TP") ~ "warm",
        Metric %in% c("FP", "TN") ~ "cool",
        Metric == "BalancedAccuracy" ~ "acc",
        TRUE ~ "assoc"
      )
    ) %>%
    dplyr::group_by(Group) %>%
    dplyr::mutate(Norm = scale_01(Value)) %>%
    dplyr::ungroup() %>%
  dplyr::mutate(
      Label = dplyr::case_when(
        Metric %in% c("BalancedAccuracy", "MCC", "F1_score") ~ sprintf(metric_label_fmt, Value),
        TRUE ~ as.character(Value)
      )
    )

  bar_tbl <- plot_tbl %>%
    dplyr::select(Method, PredictedPositive, PredictedNegative) %>%
    tidyr::pivot_longer(-Method, names_to = "Type", values_to = "Count") %>%
    dplyr::mutate(
      Type = dplyr::recode(
        Type,
        PredictedPositive = "Predicted positive",
        PredictedNegative = "Predicted negative"
      )
    )

  guide_bottom <- function(order) {
    ggplot2::guide_colourbar(
      order = order,
      direction = "horizontal",
      title.position = "bottom",
      title.hjust = 0.5,
      barwidth = grid::unit(3.4, "cm"),
      barheight = grid::unit(0.35, "cm")
    )
  }

  p_heat <- ggplot2::ggplot() +
    geom_rounded_tile(
      data = dplyr::filter(heat_tbl, Group == "warm"),
      ggplot2::aes(X, Method, fill = Norm),
      width = tile_width,
      height = 0.9,
      radius = grid::unit(heatmap_tile_radius_mm, "mm")
    ) +
    ggplot2::scale_fill_gradientn(
      colours = heat_warm_palette,
      limits = c(0, 1),
      na.value = "grey95",
      name = "TP / FN scale",
      guide = guide_bottom(1)
    ) +
    ggnewscale::new_scale_fill() +
    geom_rounded_tile(
      data = dplyr::filter(heat_tbl, Group == "cool"),
      ggplot2::aes(X, Method, fill = Norm),
      width = tile_width,
      height = 0.9,
      radius = grid::unit(heatmap_tile_radius_mm, "mm")
    ) +
    ggplot2::scale_fill_gradientn(
      colours = heat_cool_palette,
      limits = c(0, 1),
      na.value = "grey95",
      name = "TN / FP scale",
      guide = guide_bottom(2)
    ) +
    ggnewscale::new_scale_fill() +
    geom_rounded_tile(
      data = dplyr::filter(heat_tbl, Group == "acc"),
      ggplot2::aes(X, Method, fill = Norm),
      width = tile_width,
      height = 0.9,
      radius = grid::unit(heatmap_tile_radius_mm, "mm")
    ) +
    ggplot2::scale_fill_gradientn(
      colours = heat_acc_palette,
      limits = c(0, 1),
      na.value = "grey95",
      name = "Balanced accuracy",
      guide = guide_bottom(3)
    ) +
    ggnewscale::new_scale_fill() +
    geom_rounded_tile(
      data = dplyr::filter(heat_tbl, Group == "assoc"),
      ggplot2::aes(X, Method, fill = Norm),
      width = tile_width,
      height = 0.9,
      radius = grid::unit(heatmap_tile_radius_mm, "mm")
    ) +
    ggplot2::scale_fill_gradientn(
      colours = heat_assoc_palette,
      limits = c(0, 1),
      na.value = "grey95",
      name = "MCC / F1 scale",
      guide = guide_bottom(4)
    ) +
    ggplot2::geom_text(
      data = heat_tbl,
      ggplot2::aes(X, Method, label = Label, colour = Norm > 0.6),
      size = 3,
      show.legend = FALSE
    ) +
    ggplot2::scale_colour_manual(values = c("TRUE" = "white", "FALSE" = "black")) +
    ggplot2::scale_x_continuous(
      breaks = centers,
      labels = c("FN", "TP", "FP", "TN", "Balanced Acc.", "MCC", "F1"),
      expand = ggplot2::expansion(mult = 0, add = 0)
    ) +
    ggplot2::scale_y_discrete(expand = c(0, 0)) +
    ggplot2::theme_minimal(base_size = heatmap_base_font) +
    ggplot2::theme(
      axis.title = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 35, hjust = 1, vjust = 1),
      panel.grid = ggplot2::element_blank(),
      legend.position = "top",
      legend.box = "horizontal",
      legend.box.just = "left",
      legend.spacing.x = grid::unit(5, "pt"),
      legend.title = ggplot2::element_text(size = heatmap_base_font - 1),
      plot.margin = ggplot2::margin(t = 4, r = 0, b = 4, l = 4)
    )

  p_bar <- ggplot2::ggplot(bar_tbl, ggplot2::aes(Count, Method, fill = Type)) +
    ggplot2::geom_col(width = 0.85) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = 0, add = 0)) +
    ggplot2::scale_fill_manual(values = stack_palette, name = NULL) +
    ggplot2::labs(x = "Count") +
    ggplot2::theme_minimal(base_size = heatmap_base_font) +
    ggplot2::theme(
      axis.title.y = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      legend.position = "top",
      plot.margin = ggplot2::margin(t = 4, r = 4, b = 4, l = 0)
    )

  if (heatmap_bar_gap_mm > 0) {
    full_plot <- p_heat + patchwork::plot_spacer() + p_bar +
      patchwork::plot_layout(
        widths = grid::unit(
          c(heat_panel_ratio, heatmap_bar_gap_mm, bar_panel_ratio),
          c("null", "mm", "null")
        )
      )
  } else {
    full_plot <- p_heat + p_bar +
      patchwork::plot_layout(widths = c(heat_panel_ratio, bar_panel_ratio))
  }

  ggplot2::ggsave(
    filename = save_path,
    plot = full_plot,
    width = 15,
    height = max(5, nlevels(plot_tbl$Method) * 0.33),
    device = "svg",
    limitsize = FALSE
  )
}

save_combined_patchwork <- function(plot_list, save_path, item_width, item_height) {
  plot_list <- purrr::compact(plot_list)
  if (length(plot_list) == 0) {
    return(invisible(NULL))
  }

  ncol_grid <- min(length(plot_list), patchwork_col_count)
  combined_plot <- patchwork::wrap_plots(plotlist = plot_list, ncol = ncol_grid)

  ggplot2::ggsave(
    filename = save_path,
    plot = combined_plot,
    width = item_width * ncol_grid,
    height = item_height * ceiling(length(plot_list) / ncol_grid),
    device = "svg",
    limitsize = FALSE
  )
}

make_combined_roc_overlay <- function(roc_objects,
                                      save_path,
                                      plot_width,
                                      plot_height,
                                      legend_position,
                                      legend_ncol,
                                      show_auc_in_legend,
                                      smooth_curves) {
  roc_objects <- roc_objects[!vapply(roc_objects, is.null, logical(1))]
  if (length(roc_objects) == 0) {
    return(invisible(NULL))
  }

  legend_position <- match.arg(legend_position, c("right", "bottom"))

  auc_tbl <- tibble::tibble(
    Method = names(roc_objects),
    AUC = purrr::map_dbl(roc_objects, ~ as.numeric(pROC::auc(.x)))
  ) %>%
    dplyr::arrange(dplyr::desc(AUC))

  roc_objects <- roc_objects[auc_tbl$Method]
  display_names <- if (isTRUE(show_auc_in_legend)) {
    paste0(auc_tbl$Method, " (AUC=", sprintf("%.3f", auc_tbl$AUC), ")")
  } else {
    auc_tbl$Method
  }

  plot_roc_objects <- roc_objects
  if (isTRUE(smooth_curves)) {
    plot_roc_objects <- purrr::imap(
      roc_objects,
      function(roc_obj, method) suppressWarnings(
        tryCatch(
          pROC::smooth(roc_obj, method = "binormal"),
          error = function(e) {
            message(
              ">>> ROC 平滑跳过: ", method,
              "（", conditionMessage(e), "）"
            )
            roc_obj
          }
        )
      )
    )
  }
  names(plot_roc_objects) <- display_names

  overlay_colors <- grDevices::colorRampPalette(npg_palette)(length(roc_objects))
  curve_scales <- purrr::map_dbl(
    plot_roc_objects,
    ~ if (isTRUE(.x$percent)) 100 else 1
  )
  if (length(unique(curve_scales)) != 1) {
    stop("ROC 曲线的 percent 设置不一致，无法绘制在同一坐标系中。")
  }

  coord_frames <- purrr::map(
    plot_roc_objects,
    ~ pROC::coords(
      .x,
      x = "all",
      ret = c("specificity", "sensitivity"),
      transpose = FALSE,
      ignore.partial.auc = TRUE
    )
  )
  overlay_tbl <- combine_roc_coordinate_frames(
    coord_frames = coord_frames,
    display_names = names(plot_roc_objects),
    curve_scales = curve_scales
  )
  overlay_tbl$Method <- factor(overlay_tbl$Method, levels = display_names)
  axis_max <- unique(curve_scales)[[1]]

  overlay_plot <- ggplot2::ggplot(
    overlay_tbl,
    ggplot2::aes(x = FPR, y = TPR, colour = Method, group = Method)
  ) +
    ggplot2::geom_line(linewidth = 0.95) +
    ggplot2::geom_abline(linetype = "dashed", colour = "grey55") +
    ggplot2::coord_cartesian(xlim = c(0, axis_max), ylim = c(0, axis_max)) +
    ggplot2::scale_colour_manual(values = stats::setNames(overlay_colors, display_names)) +
    ggplot2::labs(
      title = NULL,
      x = "1 - Specificity (FPR)",
      y = "Sensitivity (TPR)",
      colour = NULL
    ) +
    ggplot2::guides(
      colour = ggplot2::guide_legend(
        ncol = legend_ncol,
        byrow = TRUE
      )
    ) +
    ggplot2::theme_classic(base_size = 13) +
    ggplot2::theme(
      legend.position = legend_position,
      legend.text = ggplot2::element_text(size = 8),
      legend.title = ggplot2::element_blank(),
      plot.title = ggplot2::element_blank()
    )

  ggplot2::ggsave(
    filename = save_path,
    plot = overlay_plot,
    width = plot_width,
    height = plot_height,
    device = "svg",
    limitsize = FALSE
  )
}

make_combined_pr_overlay <- function(pr_point_objects,
                                     save_path,
                                     plot_width,
                                     plot_height,
                                     legend_position,
                                     legend_ncol,
                                     show_auprc_in_legend) {
  pr_point_objects <- pr_point_objects[!vapply(pr_point_objects, is.null, logical(1))]
  if (length(pr_point_objects) == 0) {
    return(invisible(NULL))
  }

  legend_position <- match.arg(legend_position, c("right", "bottom"))

  pr_auc_tbl <- tibble::tibble(
    Method = names(pr_point_objects),
    AUPRC = purrr::map_dbl(pr_point_objects, function(pr_tbl) {
      recall_diff <- diff(pr_tbl$Recall)
      precision_mid <- (utils::head(pr_tbl$Precision, -1) + utils::tail(pr_tbl$Precision, -1)) / 2
      sum(recall_diff * precision_mid, na.rm = TRUE)
    }),
    Baseline = purrr::map_dbl(pr_point_objects, function(pr_tbl) {
      max(pr_tbl$TP, na.rm = TRUE) / (max(pr_tbl$TP, na.rm = TRUE) + max(pr_tbl$FP, na.rm = TRUE))
    })
  ) %>%
    dplyr::arrange(dplyr::desc(AUPRC))

  display_lookup <- if (isTRUE(show_auprc_in_legend)) {
    stats::setNames(
      paste0(pr_auc_tbl$Method, " (AUPRC=", sprintf("%.3f", pr_auc_tbl$AUPRC), ")"),
      pr_auc_tbl$Method
    )
  } else {
    stats::setNames(pr_auc_tbl$Method, pr_auc_tbl$Method)
  }

  plot_tbl <- purrr::imap_dfr(pr_point_objects[pr_auc_tbl$Method], function(pr_tbl, method) {
    pr_tbl %>%
      dplyr::mutate(Method = display_lookup[[method]])
  })

  display_names <- unname(display_lookup[pr_auc_tbl$Method])
  overlay_colors <- grDevices::colorRampPalette(npg_palette)(length(display_names))
  baseline <- stats::weighted.mean(pr_auc_tbl$Baseline, w = rep(1, nrow(pr_auc_tbl)), na.rm = TRUE)

  overlay_plot <- ggplot2::ggplot(plot_tbl, ggplot2::aes(x = Recall, y = Precision, colour = Method)) +
    ggplot2::geom_step(linewidth = 0.95, direction = "hv") +
    ggplot2::geom_hline(yintercept = baseline, linetype = "dashed", colour = "grey55") +
    ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    ggplot2::scale_colour_manual(values = stats::setNames(overlay_colors, display_names)) +
    ggplot2::labs(
      title = NULL,
      x = "Recall",
      y = "Precision",
      colour = NULL
    ) +
    ggplot2::guides(
      colour = ggplot2::guide_legend(
        ncol = legend_ncol,
        byrow = TRUE
      )
    ) +
    ggplot2::theme_classic(base_size = 13) +
    ggplot2::theme(
      legend.position = legend_position,
      legend.text = ggplot2::element_text(size = 8),
      legend.title = ggplot2::element_blank(),
      plot.title = ggplot2::element_blank()
    )

  ggplot2::ggsave(
    filename = save_path,
    plot = overlay_plot,
    width = plot_width,
    height = plot_height,
    device = "svg",
    limitsize = FALSE
  )
}

normalize_score_name <- function(x) {
  x %>%
    tolower() %>%
    gsub("_ensembl.*$", "", .) %>%
    gsub("_v[0-9].*$", "", .) %>%
    gsub("_coding$", "", .) %>%
    gsub("_rankscore$", "", .) %>%
    gsub("_phred$", "", .) %>%
    gsub("[^a-z0-9]", "", .)
}

read_dbnsfp_performance <- function(path) {
  if (!file.exists(path)) {
    message(">>> dbNSFP performance matrix 未找到，跳过: ", path)
    return(NULL)
  }

  db_tbl <- utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  required_cols <- c(
    "Score", "Score_Category", "ClinVar_Year", "Mode", "Stratum",
    "AUC", "AUC_Low", "AUC_High", "pAUC_FPR_05", "pAUC_FPR_05_Low",
    "pAUC_FPR_05_High", "pAUC_TPR_95", "pAUC_TPR_95_Low",
    "pAUC_TPR_95_High", "Coverage", "N_Path", "N_Benign"
  )

  missing_cols <- setdiff(required_cols, names(db_tbl))
  if (length(missing_cols) > 0) {
    stop("dbNSFP performance matrix 缺少列: ", paste(missing_cols, collapse = ", "))
  }

  db_tbl %>%
    dplyr::mutate(
      ClinVar_Year = as.integer(ClinVar_Year),
      dplyr::across(
        c(AUC, AUC_Low, AUC_High, pAUC_FPR_05, pAUC_FPR_05_Low, pAUC_FPR_05_High,
          pAUC_TPR_95, pAUC_TPR_95_Low, pAUC_TPR_95_High, Coverage, N_Path, N_Benign),
        ~ suppressWarnings(as.numeric(.x))
      ),
      ScoreKey = normalize_score_name(Score),
      ScoreCategory = tools::toTitleCase(Score_Category),
      TotalN = N_Path + N_Benign
    )
}

build_dbnsfp_benchmark_comparison <- function(db_tbl, roc_summary_tbl) {
  benchmark_auc <- roc_summary_tbl %>%
    dplyr::mutate(ScoreKey = normalize_score_name(Method)) %>%
    dplyr::select(
      ScoreKey,
      Method,
      TypeCategory,
      ScoreColumn,
      Benchmark_AUC = AUC,
      Benchmark_AUC_low = AUC_low,
      Benchmark_AUC_high = AUC_high,
      Benchmark_AUPRC = AUPRC
    )

  latest_year <- max(db_tbl$ClinVar_Year, na.rm = TRUE)

  db_tbl %>%
    dplyr::inner_join(benchmark_auc, by = "ScoreKey") %>%
    dplyr::mutate(
      dbNSFP_Score = Score,
      dbNSFP_Year = ClinVar_Year,
      dbNSFP_AUC = AUC,
      dbNSFP_AUC_Low = AUC_Low,
      dbNSFP_AUC_High = AUC_High,
      dbNSFP_Coverage = Coverage,
      IsLatestDbnsfpYear = ClinVar_Year == latest_year,
      AUC_Delta_Benchmark_minus_dbNSFP = Benchmark_AUC - dbNSFP_AUC,
      DisplayLabel = paste0(Method, " / ", dbNSFP_Score)
    ) %>%
    dplyr::select(
      Method, dbNSFP_Score, ScoreKey, TypeCategory, Score_Category, ScoreCategory,
      ScoreColumn, dbNSFP_Year, IsLatestDbnsfpYear,
      Benchmark_AUC, Benchmark_AUC_low, Benchmark_AUC_high, Benchmark_AUPRC,
      dbNSFP_AUC, dbNSFP_AUC_Low, dbNSFP_AUC_High, dbNSFP_Coverage,
      pAUC_FPR_05, pAUC_FPR_05_Low, pAUC_FPR_05_High,
      pAUC_TPR_95, pAUC_TPR_95_Low, pAUC_TPR_95_High,
      N_Path, N_Benign, TotalN, AUC_Delta_Benchmark_minus_dbNSFP, DisplayLabel
    ) %>%
    dplyr::arrange(Method, dbNSFP_Year)
}

save_dbnsfp_historical_vs_benchmark_plot <- function(history_tbl, save_path) {
  if (is.null(history_tbl) || nrow(history_tbl) == 0) {
    return(invisible(NULL))
  }

  latest_delta_order <- history_tbl %>%
    dplyr::group_by(Method) %>%
    dplyr::slice_max(dbNSFP_Year, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(TypeCategory, AUC_Delta_Benchmark_minus_dbNSFP) %>%
    dplyr::pull(Method)

  plot_tbl <- history_tbl %>%
    dplyr::mutate(
      Method = factor(Method, levels = latest_delta_order),
      Year = factor(dbNSFP_Year),
      DeltaLabel = sprintf("%+.2f", AUC_Delta_Benchmark_minus_dbNSFP)
    )

  delta_limit <- max(abs(plot_tbl$AUC_Delta_Benchmark_minus_dbNSFP), na.rm = TRUE)
  delta_limit <- max(0.05, ceiling(delta_limit * 20) / 20)

  plt <- ggplot2::ggplot(plot_tbl, ggplot2::aes(x = Year, y = Method)) +
    ggplot2::geom_tile(ggplot2::aes(fill = AUC_Delta_Benchmark_minus_dbNSFP), colour = "white", linewidth = 0.35) +
    ggplot2::geom_text(ggplot2::aes(label = DeltaLabel), size = 2.7, colour = "black") +
    ggplot2::facet_grid(TypeCategory ~ ., scales = "free_y", space = "free_y", switch = "y") +
    ggplot2::scale_fill_gradient2(
      low = npg_palette[1],
      mid = "white",
      high = npg_palette[3],
      midpoint = 0,
      limits = c(-delta_limit, delta_limit),
      oob = scales::squish,
      name = paste0(gene_name, " AUC - dbNSFP AUC")
    ) +
    ggplot2::labs(
      x = "ClinVar year in dbNSFP matrix",
      y = NULL,
      caption = paste0("Cell value is AUC difference: positive values mean the ", gene_name, " benchmark AUC is higher; negative values mean dbNSFP AUC is higher.")
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      strip.placement = "outside",
      strip.text.y.left = ggplot2::element_text(face = "bold", angle = 0, margin = ggplot2::margin(r = 6)),
      axis.text.x = ggplot2::element_text(angle = 35, hjust = 1),
      axis.text.y = ggplot2::element_text(size = 8),
      legend.position = "bottom",
      plot.caption = ggplot2::element_text(size = 8, hjust = 0),
      plot.margin = ggplot2::margin(t = 8, r = 8, b = 8, l = 8)
    ) +
    ggplot2::guides(
      fill = ggplot2::guide_colourbar(
        direction = "horizontal",
        title.position = "top",
        barwidth = grid::unit(8, "cm"),
        barheight = grid::unit(0.35, "cm")
      )
    )

  ggplot2::ggsave(
    save_path,
    plt,
    width = 14,
    height = max(dbnsfp_plot_height, dplyr::n_distinct(plot_tbl$Method) * 0.28),
    device = "svg",
    limitsize = FALSE
  )
}

save_dbnsfp_latest_dumbbell_plot <- function(latest_tbl, save_path) {
  if (is.null(latest_tbl) || nrow(latest_tbl) == 0) {
    return(invisible(NULL))
  }

  plot_tbl <- latest_tbl %>%
    dplyr::arrange(AUC_Delta_Benchmark_minus_dbNSFP) %>%
    dplyr::mutate(
      Method = factor(Method, levels = Method),
      DeltaDirection = dplyr::if_else(AUC_Delta_Benchmark_minus_dbNSFP >= 0, "Benchmark higher", "dbNSFP higher"),
      DeltaLabel = sprintf("%+.3f", AUC_Delta_Benchmark_minus_dbNSFP),
      Delta_low = Benchmark_AUC_low - dbNSFP_AUC_High,
      Delta_high = Benchmark_AUC_high - dbNSFP_AUC_Low,
      LabelX = dplyr::if_else(AUC_Delta_Benchmark_minus_dbNSFP >= 0, AUC_Delta_Benchmark_minus_dbNSFP + 0.012, AUC_Delta_Benchmark_minus_dbNSFP - 0.012),
      LabelHjust = dplyr::if_else(AUC_Delta_Benchmark_minus_dbNSFP >= 0, 0, 1)
    )

  latest_year <- max(plot_tbl$dbNSFP_Year, na.rm = TRUE)
  x_min <- min(plot_tbl$Delta_low, plot_tbl$LabelX, na.rm = TRUE) - 0.02
  x_max <- max(plot_tbl$Delta_high, plot_tbl$LabelX, na.rm = TRUE) + 0.02

  plt <- ggplot2::ggplot(plot_tbl, ggplot2::aes(y = Method)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "grey45", linewidth = 0.6) +
    ggplot2::geom_col(ggplot2::aes(x = AUC_Delta_Benchmark_minus_dbNSFP, fill = DeltaDirection), width = 0.7, alpha = 0.88) +
    ggplot2::geom_errorbar(ggplot2::aes(xmin = Delta_low, xmax = Delta_high), orientation = "y", width = 0.18, colour = "grey20", alpha = 0.65) +
    ggplot2::geom_text(ggplot2::aes(x = LabelX, label = DeltaLabel, hjust = LabelHjust), size = 3.0) +
    ggplot2::facet_grid(TypeCategory ~ ., scales = "free_y", space = "free_y", switch = "y") +
    ggplot2::scale_fill_manual(values = c("Benchmark higher" = npg_palette[3], "dbNSFP higher" = npg_palette[1]), name = NULL) +
    ggplot2::coord_cartesian(xlim = c(x_min, x_max), clip = "off") +
    ggplot2::labs(
      x = paste0(gene_name, " benchmark AUC - dbNSFP ", latest_year, " AUC"),
      y = NULL,
      caption = "Bars show AUC difference. Error bars use the outer bounds of the two 95% CIs."
    ) +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(
      legend.position = "bottom",
      legend.box = "vertical",
      legend.title = ggplot2::element_text(size = 9),
      legend.text = ggplot2::element_text(size = 9),
      strip.background = ggplot2::element_rect(fill = "grey92", colour = NA),
      strip.placement = "outside",
      strip.text.y.left = ggplot2::element_text(face = "bold", angle = 0, margin = ggplot2::margin(r = 6)),
      plot.caption = ggplot2::element_text(size = 8, hjust = 0),
      plot.margin = ggplot2::margin(t = 8, r = 14, b = 8, l = 8)
    ) +
    ggplot2::guides(
      fill = ggplot2::guide_legend(nrow = 1)
    )

  ggplot2::ggsave(save_path, plt, width = 10.5, height = max(7.5, nrow(plot_tbl) * 0.30), device = "svg", limitsize = FALSE)
}

## ===== 3. 建立输出目录 ========================================================
ensure_dir(output_root_dir)
ensure_dir(tables_dir)
ensure_dir(plots_dir)
ensure_dir(per_tool_dir)

## ===== 4. 读取并预处理数据 ====================================================
if (!file.exists(input_path)) {
  stop("未找到输入文件: ", input_path)
}

message(">>> 正在读取输入文件: ", input_path)
input_sheet <- resolve_input_sheet(
  path = input_path,
  preferred_sheet = input_sheet,
  required_columns = c(id_col, group_col)
)
message(">>> 使用工作表: ", input_sheet)
df_raw <- read_excel_safe(input_path, sheet = input_sheet)

if (!all(c(id_col, group_col) %in% names(df_raw))) {
  stop("输入文件缺少关键列，请检查是否包含: ", paste(c(id_col, group_col), collapse = ", "))
}

df <- df_raw %>%
  dplyr::mutate(
    variant_id = .data[[id_col]],
    truth = map_truth_label(.data[[group_col]]),
    group = factor(
      map_group_display(truth),
      levels = c(negative_display_label, positive_display_label)
    )
  ) %>%
  dplyr::filter(!is.na(truth))

message(">>> 数据集分组分布:")
print(table(df$group, useNA = "ifany"))

sample_info_tbl <- tibble::tibble(
  Item = c("Gene", "Input file", "Input sheet", "Total variants used", "Negative group", "Positive group"),
  Value = c(
    gene_name,
    input_file,
    as.character(input_sheet),
    as.character(nrow(df)),
    paste(negative_labels, collapse = ", "),
    paste(positive_labels, collapse = ", ")
  )
)

## ===== 5. 按软件循环分析 ======================================================
message(">>> 开始分析 ", nrow(tool_config), " 个软件...")
set.seed(bootstrap_seed)

analysis_results <- purrr::map(seq_len(nrow(tool_config)), function(i) {
  analyse_tool(df = df, config_row = tool_config[i, ], save_dir = per_tool_dir)
})

classification_tbl <- purrr::map(analysis_results, "classification") %>% dplyr::bind_rows()
roc_summary_tbl    <- purrr::map(analysis_results, "roc_summary") %>% dplyr::bind_rows()
bootstrap_summary_tbl <- purrr::map(analysis_results, "bootstrap_summary") %>% dplyr::bind_rows()
skipped_tbl        <- purrr::map(analysis_results, "skipped") %>% dplyr::bind_rows()

if (nrow(classification_tbl) == 0) {
  stop("没有生成任何有效的分类指标结果，请检查输入表格和参数设置。")
}

benchmark_summary_tbl <- classification_tbl %>%
  dplyr::left_join(
    roc_summary_tbl %>% dplyr::select(Method, ScoreColumn, AUC, AUC_low, AUC_high, AUPRC, AveragePrecision, PR_baseline, ROC_best_threshold),
    by = c("Method", "ScoreColumn")
  ) %>%
  dplyr::arrange(dplyr::desc(BalancedAccuracy), dplyr::desc(AUC), dplyr::desc(AUPRC))

## ===== 6. dbNSFP performance matrix 比较 =====================================
dbnsfp_path <- find_file_upwards(dbnsfp_performance_file, work_dir)
dbnsfp_tbl <- read_dbnsfp_performance(dbnsfp_path)
dbnsfp_historical_vs_benchmark_tbl <- tibble::tibble()
dbnsfp_latest_vs_benchmark_tbl <- tibble::tibble()

if (!is.null(dbnsfp_tbl)) {
  message(">>> 正在分析 dbNSFP performance matrix: ", dbnsfp_path)
  dbnsfp_historical_vs_benchmark_tbl <- build_dbnsfp_benchmark_comparison(dbnsfp_tbl, roc_summary_tbl)
  dbnsfp_latest_vs_benchmark_tbl <- dbnsfp_historical_vs_benchmark_tbl %>%
    dplyr::filter(IsLatestDbnsfpYear) %>%
    dplyr::arrange(dplyr::desc(abs(AUC_Delta_Benchmark_minus_dbNSFP)))
}

## ===== 7. 导出 Excel 汇总表 ===================================================
excel_output_path <- file.path(tables_dir, paste0(gene_name, "_benchmark_v5_summary.xlsx"))
writexl::write_xlsx(
  list(
    SampleInfo = sample_info_tbl,
    ThresholdConfig = tool_config,
    ThresholdMetrics = classification_tbl,
    ROCSummary = roc_summary_tbl,
    BootstrapValidation = if (nrow(bootstrap_summary_tbl) > 0) bootstrap_summary_tbl else tibble::tibble(),
    BenchmarkSummary = benchmark_summary_tbl,
    DbnsfpPerformanceRaw = if (!is.null(dbnsfp_tbl)) dbnsfp_tbl else tibble::tibble(),
    DbnsfpHistoricalVsBenchmark = dbnsfp_historical_vs_benchmark_tbl,
    DbnsfpLatestVsBenchmark = dbnsfp_latest_vs_benchmark_tbl,
    SkippedTools = if (nrow(skipped_tbl) > 0) skipped_tbl else tibble::tibble(Method = character(), ScoreColumn = character(), Reason = character())
  ),
  path = excel_output_path
)
message(">>> Excel 汇总表已保存: ", excel_output_path)

## ===== 8. 绘制混淆矩阵热图 ====================================================
heatmap_output_path <- file.path(plots_dir, paste0(gene_name, "_threshold_metrics_heatmap.svg"))
make_threshold_heatmap(classification_tbl, heatmap_output_path)
message(">>> 阈值混淆矩阵热图已保存: ", heatmap_output_path)

## ===== 9. 生成 ROC / PR 组合图 ================================================
compare_plot_list <- purrr::map(analysis_results, "compare_plot")
roc_plot_list     <- purrr::map(analysis_results, "roc_plot")
pr_plot_list      <- purrr::map(analysis_results, "pr_plot")

save_combined_patchwork(
  compare_plot_list,
  file.path(plots_dir, paste0(gene_name, "_Combined_Compare.svg")),
  compare_plot_width,
  compare_plot_height
)

save_combined_patchwork(
  roc_plot_list,
  file.path(plots_dir, paste0(gene_name, "_Combined_Individual_ROC.svg")),
  roc_plot_width,
  roc_plot_height
)

save_combined_patchwork(
  pr_plot_list,
  file.path(plots_dir, paste0(gene_name, "_Combined_Individual_PR.svg")),
  pr_plot_width,
  pr_plot_height
)

roc_objects <- purrr::map(analysis_results, "roc_obj")
roc_objects <- stats::setNames(roc_objects, tool_config$Method)
make_combined_roc_overlay(
  roc_objects = roc_objects,
  save_path = file.path(plots_dir, paste0(gene_name, "_Combined_ROC_Overlay.svg")),
  plot_width = overall_combined_roc_width,
  plot_height = overall_combined_roc_height,
  legend_position = overall_combined_roc_legend_position,
  legend_ncol = overall_combined_roc_legend_ncol,
  show_auc_in_legend = overall_combined_roc_show_auc_in_legend,
  smooth_curves = overall_combined_roc_smooth_curves
)
message(">>> 总体 Combined ROC 已保存: ", file.path(plots_dir, paste0(gene_name, "_Combined_ROC_Overlay.svg")))

pr_point_objects <- purrr::map(analysis_results, "pr_points")
pr_point_objects <- stats::setNames(pr_point_objects, tool_config$Method)
make_combined_pr_overlay(
  pr_point_objects = pr_point_objects,
  save_path = file.path(plots_dir, paste0(gene_name, "_Combined_PR_Overlay.svg")),
  plot_width = overall_combined_pr_width,
  plot_height = overall_combined_pr_height,
  legend_position = overall_combined_pr_legend_position,
  legend_ncol = overall_combined_pr_legend_ncol,
  show_auprc_in_legend = overall_combined_pr_show_auprc_in_legend
)
message(">>> 总体 Combined PR 已保存: ", file.path(plots_dir, paste0(gene_name, "_Combined_PR_Overlay.svg")))

available_type_config <- tool_config %>%
  dplyr::select(Method, TypeCategory) %>%
  dplyr::distinct()

for (type_name in category_combined_roc_types) {
  method_subset <- available_type_config %>%
    dplyr::filter(TypeCategory == type_name) %>%
    dplyr::pull(Method)

  category_roc_objects <- roc_objects[method_subset]
  category_roc_objects <- category_roc_objects[!vapply(category_roc_objects, is.null, logical(1))]

  if (length(category_roc_objects) == 0) {
    message(">>> 分类 Combined ROC 跳过: ", type_name, "（当前无可用软件）")
    next
  }

  category_save_path <- file.path(
    plots_dir,
    paste0(gene_name, "_Combined_ROC_", safe_fname(type_name), ".svg")
  )

  make_combined_roc_overlay(
    roc_objects = category_roc_objects,
    save_path = category_save_path,
    plot_width = category_combined_roc_width,
    plot_height = category_combined_roc_height,
    legend_position = category_combined_roc_legend_position,
    legend_ncol = category_combined_roc_legend_ncol,
    show_auc_in_legend = category_combined_roc_show_auc_in_legend,
    smooth_curves = category_combined_roc_smooth_curves
  )
  message(">>> 分类 Combined ROC 已保存: ", category_save_path)

  category_pr_objects <- pr_point_objects[method_subset]
  category_pr_objects <- category_pr_objects[!vapply(category_pr_objects, is.null, logical(1))]

  if (length(category_pr_objects) == 0) {
    message(">>> 分类 Combined PR 跳过: ", type_name, "（当前无可用软件）")
    next
  }

  category_pr_save_path <- file.path(
    plots_dir,
    paste0(gene_name, "_Combined_PR_", safe_fname(type_name), ".svg")
  )

  make_combined_pr_overlay(
    pr_point_objects = category_pr_objects,
    save_path = category_pr_save_path,
    plot_width = category_combined_pr_width,
    plot_height = category_combined_pr_height,
    legend_position = category_combined_pr_legend_position,
    legend_ncol = category_combined_pr_legend_ncol,
    show_auprc_in_legend = category_combined_pr_show_auprc_in_legend
  )
  message(">>> 分类 Combined PR 已保存: ", category_pr_save_path)
}

## ===== 10. 生成 dbNSFP 与当前 benchmark AUC 对照图 ===========================
if (nrow(dbnsfp_historical_vs_benchmark_tbl) > 0) {
  obsolete_dbnsfp_plot_files <- file.path(
    plots_dir,
    paste0(
      gene_name,
      c(
        "_dbNSFP_AUC_year_heatmap.svg",
        "_dbNSFP_top_AUC_trends.svg",
        "_dbNSFP_latest_metric_heatmap.svg",
        "_dbNSFP_AUC_year_delta.svg",
        "_dbNSFP_vs_GRIN1_AUC.svg"
      )
    )
  )
  unlink(obsolete_dbnsfp_plot_files[file.exists(obsolete_dbnsfp_plot_files)])

  dbnsfp_history_path <- file.path(plots_dir, paste0(gene_name, "_dbNSFP_historical_AUC_vs_Benchmark.svg"))
  save_dbnsfp_historical_vs_benchmark_plot(dbnsfp_historical_vs_benchmark_tbl, dbnsfp_history_path)
  message(">>> dbNSFP 历年 AUC 与当前 benchmark AUC 对照图已保存: ", dbnsfp_history_path)

  dbnsfp_latest_path <- file.path(plots_dir, paste0(gene_name, "_dbNSFP_latest_AUC_delta_vs_Benchmark.svg"))
  save_dbnsfp_latest_dumbbell_plot(dbnsfp_latest_vs_benchmark_tbl, dbnsfp_latest_path)
  message(">>> dbNSFP 最新年份 AUC 与当前 benchmark AUC 差异图已保存: ", dbnsfp_latest_path)
}

## ===== 11. 输出完成信息 ========================================================
message("\n>>> 所有分析已完成！")
message(">>> 结果文件夹: ", output_root_dir)
