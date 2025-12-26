############################################################
############################################################

# install.packages(c("dplyr", "caret", "openxlsx", "readxl")) # if needed
library(readxl)
library(dplyr)
library(caret)
library(openxlsx)

# ---------- 1. Load ----------
df <- read_xlsx("C:/Users/wpypo/Desktop/SCN3A/...")   # <- adjust path

# ---------- 2. Clean abnormal header ----------
names(df) <- sub(".*fathmm-XF_coding_pred", "fathmm-XF_coding_pred", names(df))

# ---------- 3. Create fathmm-XF_coding_pred if only numeric score present ----------
if ("fathmm-XF_coding" %in% names(df) && !"fathmm-XF_coding_pred" %in% names(df)) {
  df <- df %>% mutate(`fathmm-XF_coding_pred` =
                        ifelse(`fathmm-XF_coding` >= 0.5, "D", "N"))
}

# ---------- 4. Mapping helpers ----------
label_true <- function(x) ifelse(x == "Pathogenic", 1, 0)

label_pred <- function(pred, col) {
  pred <- as.character(pred)
  switch(col,
         "Polyphen2_HDIV_pred"   = ifelse(pred == "D", 1, ifelse(pred %in% c("B","P"), 0, NA)),
         "Polyphen2_HVAR_pred"   = ifelse(pred == "D", 1, ifelse(pred %in% c("B","P"), 0, NA)),
         "MutationTaster_pred"   = ifelse(pred %in% c("A","D"), 1, ifelse(pred %in% c("N","P"), 0, NA)),
         "MutationAssessor_pred" = ifelse(pred %in% c("H","M"), 1, ifelse(pred %in% c("N","L"), 0, NA)),
         "M-CAP_pred"            = ifelse(pred == "D", 1, ifelse(pred == "T", 0, NA)),
         "AlphaMissense_pred"    = ifelse(pred == "P", 1, ifelse(pred %in% c("B","A"), 0, NA)),
         "fathmm-XF_coding_pred" = ifelse(pred == "D", 1, ifelse(pred == "N", 0, NA)),
         # default binary
         { ifelse(pred %in% c("D"), 1,
                  ifelse(pred %in% c("T"), 0, NA)) }
  )
}

# ---------- 5. Metrics ----------
pred_cols <- grep("_pred$", names(df), value = TRUE)
results <- lapply(pred_cols, function(col) {
  y_true <- label_true(df$Type)
  y_pred <- label_pred(df[[col]], col)
  idx <- which(!is.na(y_true) & !is.na(y_pred))
  y_true <- y_true[idx]; y_pred <- y_pred[idx]
  
  TP <- sum(y_true==1 & y_pred==1)
  TN <- sum(y_true==0 & y_pred==0)
  FP <- sum(y_true==0 & y_pred==1)
  FN <- sum(y_true==1 & y_pred==0)
  
  acc   <- (TP+TN)/(TP+TN+FP+FN)
  sens  <- ifelse(TP+FN==0, NA, TP/(TP+FN))
  spec  <- ifelse(TN+FP==0, NA, TN/(TN+FP))
  bal   <- ifelse(is.na(sens) | is.na(spec), NA, (sens+spec)/2)
  ppv   <- ifelse(TP+FP==0, NA, TP/(TP+FP))
  npv   <- ifelse(TN+FN==0, NA, TN/(TN+FN))
  mcc_d <- sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
  mcc   <- ifelse(mcc_d==0, NA, (TP*TN - FP*FN)/mcc_d)
  f1    <- ifelse(2*TP+FP+FN==0, NA, 2*TP/(2*TP+FP+FN))
  
  data.frame(
    Method          = col,
    TP, TN, FP, FN,
    Accuracy        = round(acc,4),
    BalancedAccuracy= round(bal,4),
    Sensitivity     = round(sens,4),
    Specificity     = round(spec,4),
    PPV             = round(ppv,4),
    NPV             = round(npv,4),
    MCC             = round(mcc,4),
    F1_score        = round(f1,4)
  )
})

result <- bind_rows(results)

# ---------- 6. Export ----------
write.xlsx(result, "SCN3A_pred_metrics.xlsx", rowNames = FALSE)
print(result)
Script 3: Plotting Confusion Matrix (English Comments Version)
language-R# ==========================================================
#  SCN3A Heatmap + Stacked Bar Chart
#  – FN-TP Adjacent, FP-TN Adjacent, BalancedAccuracy Separate
#  – legend_left_shift / heat_left_shift / panel_ratio Adjustable
# ==========================================================

## ====== Adjustable Parameters ==========================================
tile_width         <- 2       # Width of each tile
within_pair_gap    <- 0.0     # Increment in center distance for FN-TP and FP-TN pairs
between_pair_gap   <- 1.0     # Gap between groups (≥ tile_width)
heat_panel_ratio   <- 5       # Relative width of heatmap panel in patchwork
bar_panel_ratio    <- 2       # Relative width of bar chart panel
legend_left_shift  <- 10      # Left shift of legend box (pt; negative = leftward)
heat_left_shift    <- 50      # Right shift of entire heatmap (pt; positive = rightward)
bar_left_shift     <- -55     # Left shift of bar chart (pt; negative = toward heatmap)
bar_width          <- 0.85    # Width of stacked bar columns
axis_title_gap     <- -20     # Distance between "Count" and x-axis
base_font          <- 10      # Uniform font size (pt)
## ==========================================================

library(readxl);    library(dplyr);  library(tidyr);   library(ggplot2)
library(ggnewscale);                 library(patchwork);  library(scales)

# ---------- 1. Read and Prepare Data --------------------------------------
data_dir  <- "C:/Users/wpypo/Desktop/SCN3A/赵嘉兴专属文件夹/重新分析"
xlsx_file <- file.path(data_dir, "SCN3A_pred_metrics.xlsx")

stats <- read_xlsx(xlsx_file) |>
  rename(variable = Method) |>
  mutate(
    DEL = TP + FP,
    TOL = TN + FN
  ) |>
  filter(!is.na(BalancedAccuracy)) |>
  arrange(desc(BalancedAccuracy))

stats$variable <- factor(stats$variable, levels = rev(stats$variable))

# ---------- 2. Calculate Column Center Coordinates --------------------------------------
centers <- numeric(5)
centers[1] <- 0
centers[2] <- centers[1] + tile_width + within_pair_gap
centers[3] <- centers[2] + tile_width + between_pair_gap
centers[4] <- centers[3] + tile_width + within_pair_gap
centers[5] <- centers[4] + tile_width + between_pair_gap
x_pos <- setNames(
  centers,
  c("FN_SCN3A","TP_SCN3A","FP_SCN3A","TN_SCN3A","BalancedAccuracy")
)

# ---------- 3. Long Format for Heatmap --------------------------------------------
heat <- stats |>
  select(variable,
         FN_SCN3A = FN, TP_SCN3A = TP,
         FP_SCN3A = FP, TN_SCN3A = TN,
         BalancedAccuracy) |>
  pivot_longer(-variable, names_to = "metric", values_to = "val") |>
  mutate(
    x = x_pos[metric],
    grp = case_when(
      metric %in% c("FN_SCN3A","TP_SCN3A") ~ "warm",
      metric %in% c("FP_SCN3A","TN_SCN3A") ~ "cool",
      TRUE                                 ~ "acc"
    )
  ) |>
  group_by(grp) |>
  mutate(
    norm = if (all(is.na(val))) NA else
      (val - min(val, na.rm = TRUE)) /
      (max(val, na.rm = TRUE) - min(val, na.rm = TRUE))
  ) |>
  ungroup() |>
  mutate(
    label = ifelse(
      metric == "BalancedAccuracy",
      percent(val, 0.1),
      ifelse(!is.na(val), as.character(val), "")
    )
  )

# ---------- 4. Data for Bar Chart ------------------------------------------
bar <- stats |>
  select(variable, DEL, TOL) |>
  pivot_longer(-variable, names_to = "type", values_to = "cnt")

# ---------- 5. Color Palettes & Legend -------------------------------------
pal_warm <- colorRampPalette(c("#ffffff","#ffd27f","#d73027"))(100)
pal_cool <- colorRampPalette(c("#ffffff","#9ecae1","#08519c"))(100)
pal_acc  <- colorRampPalette(c("#ffffff","#8bc98b","#006400"))(100)

guide_bottom <- function(order){
  guide_colourbar(order = order, title.position = "bottom", title.hjust = 0.5)
}

# ---------- 6. Heatmap -------------------------------------------------
p_heat <- ggplot() +
  geom_tile(data = filter(heat, grp == "warm"),
            aes(x, variable, fill = norm),
            width = tile_width, height = 0.9) +
  scale_fill_gradientn(colours = pal_warm, limits = c(0,1),
                       na.value = "grey95", name = "TP & FN density scale",
                       guide = guide_bottom(1)) +
  ggnewscale::new_scale_fill() +
  geom_tile(data = filter(heat, grp == "cool"),
            aes(x, variable, fill = norm),
            width = tile_width, height = 0.9) +
  scale_fill_gradientn(colours = pal_cool, limits = c(0,1),
                       na.value = "grey95", name = "TN & FP density scale",
                       guide = guide_bottom(2)) +
  ggnewscale::new_scale_fill() +
  geom_tile(data = filter(heat, grp == "acc"),
            aes(x, variable, fill = norm),
            width = tile_width, height = 0.9) +
  scale_fill_gradientn(colours = pal_acc, limits = c(0,1),
                       na.value = "grey95", name = "Balanced accuracy",
                       guide = guide_bottom(3)) +
  geom_text(data = filter(heat, !is.na(val)),
            aes(x, variable, label = label, colour = norm > 0.6),
            size = 3, show.legend = FALSE) +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "black")) +
  scale_x_continuous(
    breaks = centers,
    labels = c("FN-SCN3A","TP-SCN3A",
               "FP-SCN3A","TN-SCN3A",
               "Balanced Acc."),
    expand = expansion(add = between_pair_gap / 2)
  ) +
  scale_y_discrete(expand = c(0, 0)) +
  coord_fixed() +
  theme_minimal(base_size = base_font) +
  theme(
    plot.margin    = margin(l = heat_left_shift, r = 0, t = 5, b = 5),
    axis.title     = element_blank(),
    axis.text.x    = element_text(angle = 40, hjust = 1, vjust = 1),
    panel.grid     = element_blank(),
    legend.position= "top",
    legend.box.just= "left",
    legend.margin  = margin(l = legend_left_shift),
    legend.title   = element_text(size = base_font - 1)
  )

# ---------- 7. Stacked Bar Chart --------------------------------------
p_bar <- ggplot(bar, aes(cnt, variable, fill = type)) +
  geom_col(width = bar_width) +
  scale_fill_manual(values = c(DEL = "#e41a1c", TOL = "#377eb8"),
                    name = "Variant effect") +
  theme_minimal(base_size = base_font) +
  theme(
    plot.margin       = margin(l = bar_left_shift, r = 5, t = 5, b = 5),
    axis.title.y      = element_blank(),
    axis.text.y       = element_blank(),
    axis.ticks.y      = element_blank(),
    panel.grid.major.y= element_blank(),
    legend.position   = "top",
    axis.title.x      = element_text(margin = margin(t = axis_title_gap))
  ) +
  labs(x = "Count")

# ---------- 8. Combine & Save -----------------------------------------
full <- p_heat + p_bar +
  plot_layout(widths = c(heat_panel_ratio, bar_panel_ratio))

out_file <- file.path(data_dir, "SCN3A_pred_metrics_heatmap.svg")
ggsave(out_file, full,
       width  = 12,
       height = max(4, nlevels(stats$variable) * 0.3),
       device = "svg")
message(">>> Generated: ", out_file)