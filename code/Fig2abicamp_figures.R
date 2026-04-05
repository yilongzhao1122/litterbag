############################################################
# Fig. 2a and Fig. 2b: iCAMP process importance comparison
# Project: litterbag
#
# Expected project structure:
# litterbag/
# ├── code/
# │   └── Fig2ab_iCAMP_process_importance.R
# ├── data/
# │   └── icamp/
# │       ├── litterbag.16S.iCAMP.BootSummary.csv
# │       ├── litterbag.16S.iCAMP.Compare.csv
# │       ├── litterbag.ITS..iCAMP.BootSummary.csv
# │       └── litterbag.ITS..iCAMP.Compare.csv
# └── output/
#
# This script:
# 1) reads iCAMP boot summary and pairwise comparison tables
#    for bacterial and fungal communities,
# 2) prepares process-importance summaries for half, control,
#    and double precipitation treatments,
# 3) draws the paired Fig. 2a and Fig. 2b panels,
# 4) saves the combined figure.
#
# Input files used:
# - data/icamp/litterbag.16S.iCAMP.BootSummary.csv
# - data/icamp/litterbag.16S.iCAMP.Compare.csv
# - data/icamp/litterbag.ITS..iCAMP.BootSummary.csv
# - data/icamp/litterbag.ITS..iCAMP.Compare.csv
#
# Output files:
# - output/icamp/Fig2ab_iCAMP_process_importance.pdf
# - output/icamp/Fig2ab_iCAMP_process_importance.png
############################################################

# Locate project root automatically
root_candidates <- c(
  getwd(),
  dirname(getwd()),
  dirname(dirname(getwd()))
)

root_candidates <- unique(normalizePath(root_candidates, winslash = "/", mustWork = FALSE))

project_root <- NULL
for (x in root_candidates) {
  if (dir.exists(file.path(x, "code")) && dir.exists(file.path(x, "data"))) {
    project_root <- x
    break
  }
}

if (is.null(project_root)) {
  stop("Project root not found. Please run this script from the project root or from the code/ directory.")
}

data_dir <- file.path(project_root, "data", "icamp")
output_dir <- file.path(project_root, "output", "icamp")

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

boot_bac_file <- file.path(data_dir, "litterbag.16S.iCAMP.BootSummary.csv")
cmp_bac_file  <- file.path(data_dir, "litterbag.16S.iCAMP.Compare.csv")
boot_fun_file <- file.path(data_dir, "litterbag.ITS..iCAMP.BootSummary.csv")
cmp_fun_file  <- file.path(data_dir, "litterbag.ITS..iCAMP.Compare.csv")

required_files <- c(
  boot_bac_file, cmp_bac_file,
  boot_fun_file, cmp_fun_file
)

missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop(
    "Missing required input file(s):\n",
    paste(missing_files, collapse = "\n")
  )
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(ggpubr)
  library(stringr)
})

## =========================
## Global settings
## =========================
proc_levels <- c(
  "Heterogeneous.Selection",
  "Homogeneous.Selection",
  "Dispersal.Limitation",
  "Homogenizing.Dispersal",
  "Drift.and.Others"
)

proc_labels <- c(
  "Heterogeneous\nSelection",
  "Homogeneous\nSelection",
  "Dispersal\nLimitation",
  "Homogenizing\nDispersal",
  "Drift"
)

groups_3 <- c("half", "normal", "double")
group_colors <- c(
  "half"   = "#f39b7f",
  "normal" = "#3c5488",
  "double" = "#4dbbd5"
)

pair_levels <- c("half-normal", "normal-double", "half-double")
offsets <- c(half = -0.27, normal = 0.00, double = 0.27)
pair_height <- c(
  `half-normal`   = 0.90,
  `normal-double` = 0.95,
  `half-double`   = 1.00
)
label_lift <- 0.02

## ==========================================
## Data preparation
## - Remove column 13 from BootSummary to avoid
##   duplicated-column issues from some exports
## - Use pairwise effect-size labels and p-value stars
## ==========================================
prepare_icamp_data <- function(boot_summary_file, compare_file, community) {
  # Boot summary
  boot_raw <- read.csv(boot_summary_file, check.names = FALSE)
  if (ncol(boot_raw) >= 13) {
    boot_raw <- boot_raw[, -13]
  }
  
  boot_summary <- boot_raw %>%
    dplyr::filter(
      Group %in% groups_3,
      Process %in% proc_levels
    ) %>%
    dplyr::transmute(
      Group = factor(Group, levels = groups_3),
      Process = factor(Process, levels = proc_levels, labels = proc_labels),
      Observed = as.numeric(Observed),
      Stdev = as.numeric(Stdev),
      Community = community
    )
  
  # Pairwise comparison results
  cmp_raw <- read.csv(compare_file, check.names = FALSE)
  
  cmp_long <- cmp_raw %>%
    dplyr::mutate(
      Pair = dplyr::case_when(
        (Group1 == "half"   & Group2 == "normal") | (Group1 == "normal" & Group2 == "half")   ~ "half-normal",
        (Group1 == "normal" & Group2 == "double") | (Group1 == "double" & Group2 == "normal") ~ "normal-double",
        (Group1 == "half"   & Group2 == "double") | (Group1 == "double" & Group2 == "half")   ~ "half-double",
        TRUE ~ NA_character_
      )
    ) %>%
    dplyr::filter(!is.na(Pair)) %>%
    tidyr::pivot_longer(
      cols = dplyr::matches("(_Effect.Size|_P.value)$"),
      names_to = c("ProcessRaw", ".value"),
      names_sep = "_"
    ) %>%
    dplyr::mutate(
      Process = factor(ProcessRaw, levels = proc_levels, labels = proc_labels),
      Effect.Size = as.character(Effect.Size),
      P.value = suppressWarnings(as.numeric(P.value)),
      Significance = dplyr::case_when(
        !is.na(P.value) & P.value < 0.001 ~ "***",
        !is.na(P.value) & P.value < 0.01  ~ "**",
        !is.na(P.value) & P.value < 0.05  ~ "*",
        TRUE ~ ""
      ),
      Pair = factor(Pair, levels = pair_levels)
    ) %>%
    dplyr::group_by(Pair, Process) %>%
    dplyr::summarise(
      EffectLabel = dplyr::first(Effect.Size[!is.na(Effect.Size)]),
      Significance = dplyr::first(Significance),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      Label = stringr::str_trim(paste0(EffectLabel, " ", Significance))
    )
  
  list(boot = boot_summary, cmp = cmp_long)
}

## =======================================================
## Plotting function
## - Three treatment bars per process
## - Three pairwise brackets per process
## - Labels use "effect size class + significance stars"
## =======================================================
create_process_plot <- function(boot_df, cmp_df, community) {
  plot_data <- boot_df %>% dplyr::filter(Community == community)
  
  y_max <- max(plot_data$Observed + plot_data$Stdev, na.rm = TRUE) * 1.25
  
  base_df <- expand.grid(
    Process = levels(plot_data$Process),
    Pair = pair_levels,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  ) %>%
    dplyr::as_tibble() %>%
    dplyr::left_join(cmp_df, by = c("Process", "Pair")) %>%
    dplyr::mutate(
      x_center = as.numeric(factor(Process, levels = levels(plot_data$Process))),
      x1 = dplyr::case_when(
        Pair == "half-normal"   ~ x_center + offsets["half"],
        Pair == "normal-double" ~ x_center + offsets["normal"],
        Pair == "half-double"   ~ x_center + offsets["half"]
      ),
      x2 = dplyr::case_when(
        Pair == "half-normal"   ~ x_center + offsets["normal"],
        Pair == "normal-double" ~ x_center + offsets["double"],
        Pair == "half-double"   ~ x_center + offsets["double"]
      ),
      y = y_max * unname(pair_height[Pair]),
      ytick = y - 0.01 * y_max,
      label_x = (x1 + x2) / 2,
      label_y = y + label_lift * y_max
    )
  
  ggplot(plot_data, aes(x = Process, y = Observed, fill = Group)) +
    geom_bar(
      stat = "identity",
      position = position_dodge(width = 0.8),
      width = 0.7,
      color = "black",
      size = 0.5
    ) +
    geom_errorbar(
      aes(
        ymin = pmax(0, Observed - Stdev),
        ymax = Observed + Stdev
      ),
      position = position_dodge(width = 0.8),
      width = 0.25,
      size = 0.8,
      color = "black"
    ) +
    geom_segment(
      data = base_df,
      aes(x = x1, xend = x2, y = y, yend = y),
      inherit.aes = FALSE,
      size = 0.8,
      color = "black"
    ) +
    geom_segment(
      data = base_df,
      aes(x = x1, xend = x1, y = ytick, yend = y),
      inherit.aes = FALSE,
      size = 0.8,
      color = "black"
    ) +
    geom_segment(
      data = base_df,
      aes(x = x2, xend = x2, y = ytick, yend = y),
      inherit.aes = FALSE,
      size = 0.8,
      color = "black"
    ) +
    geom_text(
      data = base_df,
      aes(x = label_x, y = label_y, label = Label),
      inherit.aes = FALSE,
      size = 4.8,
      fontface = "bold"
    ) +
    scale_fill_manual(values = group_colors, breaks = groups_3) +
    labs(
      title = paste("Process Importance in", community),
      subtitle = "Pairwise comparisons: half-normal, normal-double, half-double",
      x = "Ecological process",
      y = "Relative importance",
      fill = "Precipitation"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 0.8),
      axis.text.x = element_text(angle = 0, hjust = 0.5, face = "bold", size = 12),
      axis.title = element_text(face = "bold", size = 12),
      plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
      plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = 12),
      legend.position = "top",
      legend.title = element_text(face = "bold")
    ) +
    scale_y_continuous(
      limits = c(0, y_max),
      expand = expansion(mult = c(0, 0.1))
    )
}

#########################################################
# 3. Load data: Bacteria / Fungi
#########################################################
bac <- prepare_icamp_data(
  boot_summary_file = boot_bac_file,
  compare_file = cmp_bac_file,
  community = "Bacteria"
)

fun <- prepare_icamp_data(
  boot_summary_file = boot_fun_file,
  compare_file = cmp_fun_file,
  community = "Fungi"
)

#########################################################
# 4. Generate plots
#########################################################
plot_bac <- create_process_plot(bac$boot, bac$cmp, "Bacteria")
plot_fun <- create_process_plot(fun$boot, fun$cmp, "Fungi")

combined_plot <- ggarrange(
  plot_bac, plot_fun,
  ncol = 2,
  labels = c("A", "B"),
  font.label = list(size = 18, face = "bold"),
  common.legend = TRUE,
  legend = "top"
)

final_plot <- annotate_figure(
  combined_plot,
  top = text_grob(
    "Comparative Analysis of Community Assembly Processes",
    face = "bold",
    size = 18,
    color = "black"
  ),
  bottom = text_grob(
    paste(
      "Error bars: standard deviation |",
      "effect size labels from iCAMP Compare |",
      "Significance: *p<0.05, **p<0.01, ***p<0.001"
    ),
    color = "grey40",
    size = 11,
    face = "italic"
  )
)

final_plot

#########################################################
# 5. Save figure
#########################################################
ggsave(
  filename = file.path(output_dir, "Fig2ab_iCAMP_process_importance.pdf"),
  plot = final_plot,
  width = 16,
  height = 9,
  dpi = 300
)

ggsave(
  filename = file.path(output_dir, "Fig2ab_iCAMP_process_importance.png"),
  plot = final_plot,
  width = 16,
  height = 9,
  dpi = 300
)