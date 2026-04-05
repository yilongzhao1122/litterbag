############################################################
# Fig. 2c: Partial Mantel correlation plot
# Project: litterbag
#
# Expected project structure:
# litterbag/
# ├── code/
# │   └── Fig2c_partial_mantel_correlation.R
# ├── data/
# │   ├── treat2.csv
# │   ├── env_clear.csv
# │   └── partial_mantel.csv
# └── output/
#
# This script:
# 1) reads litterbag treatment metadata and environmental data,
# 2) builds the lower-triangle environmental Pearson correlation matrix,
# 3) reads partial Mantel results from partial_mantel.csv,
# 4) reconstructs the curved links used in Fig. 2c,
# 5) saves the final figure.
#
# Input files used:
# - data/treat2.csv
# - data/env_clear.csv
# - data/partial_mantel.csv
#
# Output files:
# - output/correlation/Fig2c_partial_mantel_one_matrix_FINAL.pdf
# - output/correlation/Fig2c_partial_mantel_one_matrix_FINAL.png
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

data_dir <- file.path(project_root, "data")
output_dir <- file.path(project_root, "output", "correlation")

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

treat_file <- file.path(data_dir, "treat2.csv")
env_file <- file.path(data_dir, "env_clear.csv")
partial_mantel_file <- file.path(data_dir, "partial_mantel.csv")

required_files <- c(treat_file, env_file, partial_mantel_file)
missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
  stop(
    "Missing required input file(s):\n",
    paste(missing_files, collapse = "\n")
  )
}

suppressPackageStartupMessages({
  library(linkET)
  library(ggplot2)
  library(dplyr)
  library(tibble)
  library(RColorBrewer)
  library(cowplot)
  library(grid)
})

############################################################
# 0. Read and process environmental data for the lower triangle
############################################################
treat <- read.csv(
  treat_file,
  row.names = 1,
  check.names = FALSE
)

if (!("sample" %in% colnames(treat))) {
  treat$sample <- rownames(treat)
}

treat <- subset(treat, Clip == "uncliping")
treat <- treat[, -c(6)]
treat <- subset(treat, type == "litter")

env <- read.csv(
  env_file,
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)

if (!("sample" %in% colnames(env))) {
  env$sample <- rownames(env)
}

env <- subset(env, Clip == "uncliping")
env <- env[, -c(22)]  # keep consistent with the original analysis
rownames(env) <- env$sample

# Keep only litterbag samples present in treat
env <- env[rownames(env) %in% treat$sample, , drop = FALSE]

# Variables used in Table 3
env_vars_raw <- c(
  "NO3.N",
  "NH4.N",
  "TN",
  "TC",
  "pH",
  "plant.richness",
  "moisture_samplingmonth",
  "annual_moisture",
  "temperature_annual",
  "FlTotl",
  "FlC4",
  "FlC3"
)

env_labels <- c(
  "NO3-N",
  "NH4-N",
  "Soil TN",
  "Soil TC",
  "pH",
  "Plant richness",
  "Sampling-month\nmoisture",
  "Annual\nmoisture",
  "Annual\ntemperature",
  "Total plant\nbiomass",
  "C4 biomass",
  "C3 biomass"
)

env_heat <- env[, env_vars_raw, drop = FALSE]
env_heat <- env_heat[complete.cases(env_heat), , drop = FALSE]
colnames(env_heat) <- env_labels

############################################################
# 1. Read partial Mantel results from CSV
############################################################
table3_long <- read.csv(
  partial_mantel_file,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  quote = "\""
)

# Remove the first unnamed index column if present
first_col_name <- colnames(table3_long)[1]
if (
  is.na(first_col_name) ||
  first_col_name == "" ||
  grepl("^X$", first_col_name) ||
  grepl("^\\.\\.\\.", first_col_name)
) {
  table3_long <- table3_long[, -1, drop = FALSE]
}

required_pm_cols <- c("env", "spec", "r", "p")
if (!all(required_pm_cols %in% colnames(table3_long))) {
  stop(
    "partial_mantel.csv must contain the columns: ",
    paste(required_pm_cols, collapse = ", ")
  )
}

table3_long <- table3_long %>%
  dplyr::select(env, spec, r, p) %>%
  mutate(
    env = as.character(env),
    spec = as.character(spec),
    r = as.numeric(r),
    p = as.numeric(p)
  )

# Normalize environmental labels so they exactly match env_labels used in plotting
table3_long$env <- gsub("\r\n|\r|\n", "\n", table3_long$env)

# Optional safety check
missing_env_labels <- setdiff(unique(table3_long$env), env_labels)
if (length(missing_env_labels) > 0) {
  warning(
    "The following environmental labels in partial_mantel.csv do not match env_labels and may not be plotted correctly:\n",
    paste(missing_env_labels, collapse = "\n")
  )
}

############################################################
# 2. Build a mantel_test() skeleton for geom_couple
############################################################
set.seed(123)

dummy_spec <- data.frame(
  Bacteria_half    = rnorm(nrow(env_heat)) + 10,
  Bacteria_control = rnorm(nrow(env_heat)) + 10,
  Bacteria_double  = rnorm(nrow(env_heat)) + 10,
  Fungi_half       = rnorm(nrow(env_heat)) + 10,
  Fungi_control    = rnorm(nrow(env_heat)) + 10,
  Fungi_double     = rnorm(nrow(env_heat)) + 10
)
rownames(dummy_spec) <- rownames(env_heat)

mantel_skeleton <- mantel_test(
  spec = dummy_spec,
  env = env_heat,
  spec_select = list(
    "Bacteria · Half"    = "Bacteria_half",
    "Bacteria · Control" = "Bacteria_control",
    "Bacteria · Double"  = "Bacteria_double",
    "Fungi · Half"       = "Fungi_half",
    "Fungi · Control"    = "Fungi_control",
    "Fungi · Double"     = "Fungi_double"
  )
)

print(names(mantel_skeleton))
print(head(mantel_skeleton))

############################################################
# 3. Overwrite skeleton r and p with partial_mantel.csv values
############################################################
mantel_plot_df <- mantel_skeleton %>%
  left_join(table3_long, by = c("spec", "env"), suffix = c(".old", "")) %>%
  mutate(
    r = r,
    p = p,
    rd = cut(
      abs(r),
      breaks = c(-Inf, 0.05, 0.10, 0.20, Inf),
      labels = c("< 0.05", "0.05–0.10", "0.10–0.20", "≥ 0.20")
    ),
    pd = ifelse(p < 0.05, "P < 0.05", "P ≥ 0.05")
  )

if (any(is.na(mantel_plot_df$r)) || any(is.na(mantel_plot_df$p))) {
  warning("Some spec-env combinations in the plotting skeleton were not matched to partial_mantel.csv.")
}

############################################################
# 3.5 Diagonal-center coordinates for environmental labels
############################################################
n_env <- length(env_labels)

env_endpoint_df <- tibble(
  env = env_labels,
  x = seq(1.5, by = 1, length.out = n_env),
  y = rev(seq(1.5, by = 1, length.out = n_env))
) %>%
  mutate(
    x_lab = x + case_when(
      env == "Sampling-month\nmoisture" ~ 0.40,
      env == "Annual\nmoisture"         ~ 0.20,
      env == "Annual\ntemperature"      ~ 0.25,
      env == "Total plant\nbiomass"     ~ 0.30,
      TRUE                              ~ 0.10
    ),
    y_lab = y,
    label_txt = env
  )

print(env_endpoint_df)

############################################################
# 4. Plot
# First use geom_couple to obtain spec-side anchor positions,
# then redraw the curves with corrected environmental endpoints
############################################################

# 4.1 Temporary plot for extracting spec-side anchor positions
p_tmp_anchor <- qcorrplot(
  correlate(env_heat, method = "pearson"),
  type = "lower",
  diag = FALSE
) +
  geom_square() +
  geom_couple(
    aes(colour = spec, size = rd, linetype = pd),
    data = mantel_plot_df,
    curvature = nice_curvature(),
    node.colour = c("grey25", NA),
    node.fill = c("white", NA),
    node.size = c(3.0, 0.01)
  )

gb <- ggplot_build(p_tmp_anchor)

candidate_layers <- which(vapply(
  gb$data,
  function(z) all(c("x", "y", "xend", "yend") %in% names(z)),
  logical(1)
))

if (length(candidate_layers) == 0) {
  stop("No layer with x/y/xend/yend found in ggplot_build output.")
}

layer_nrows <- vapply(gb$data[candidate_layers], nrow, integer(1))
curve_layer_id <- candidate_layers[which.min(abs(layer_nrows - nrow(mantel_plot_df)))]

couple_curve_df <- as_tibble(gb$data[[curve_layer_id]])

if (nrow(couple_curve_df) != nrow(mantel_plot_df)) {
  stop(
    paste0(
      "Extracted curve layer row count (", nrow(couple_curve_df),
      ") does not match mantel_plot_df row count (", nrow(mantel_plot_df),
      "). Please inspect gb$data manually."
    )
  )
}

# 4.2 Build the manual curve table
curve_df <- bind_cols(
  mantel_plot_df %>%
    dplyr::select(spec, env, rd, pd),
  couple_curve_df %>%
    dplyr::select(x, y)
) %>%
  left_join(
    env_endpoint_df %>% dplyr::select(env, x_lab, y_lab),
    by = "env"
  )

spec_node_df <- curve_df %>%
  group_by(spec) %>%
  summarise(
    x = median(x, na.rm = TRUE),
    y = median(y, na.rm = TRUE),
    .groups = "drop"
  )

print(spec_node_df)
print(head(curve_df))

############################################################
# 3.6 Spec-side text labels
############################################################
spec_text_df <- spec_node_df %>%
  mutate(
    x_lab = x - 0.22,
    y_lab = y,
    label_txt = spec
  )

print(spec_text_df)

# 4.3 Final plot
p_pm_one <- qcorrplot(
  correlate(env_heat, method = "pearson"),
  type = "lower",
  diag = FALSE
) +
  geom_square() +
  geom_curve(
    data = curve_df,
    aes(
      x = x, y = y,
      xend = x_lab, yend = y_lab,
      colour = spec,
      linewidth = rd,
      linetype = pd
    ),
    curvature = 0.12,
    lineend = "round",
    inherit.aes = FALSE
  ) +
  geom_point(
    data = spec_node_df,
    aes(x = x, y = y, colour = spec),
    inherit.aes = FALSE,
    size = 3.0,
    shape = 16,
    show.legend = FALSE
  ) +
  geom_text(
    data = spec_text_df,
    aes(x = x_lab, y = y_lab, label = label_txt, colour = spec),
    inherit.aes = FALSE,
    hjust = 1,
    vjust = 0.5,
    size = 3.0,
    lineheight = 0.95,
    show.legend = FALSE
  ) +
  geom_text(
    data = env_endpoint_df,
    aes(x = x_lab, y = y_lab, label = label_txt),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 0.5,
    size = 3.0,
    lineheight = 0.9,
    colour = "black"
  ) +
  scale_fill_gradientn(
    colours = RColorBrewer::brewer.pal(11, "RdBu"),
    limits = c(-1, 1),
    breaks = seq(-1, 1, by = 0.2),
    guide = guide_colorbar(
      title = "Pearson's r\n(env–env)",
      title.position = "top",
      barheight = unit(92, "mm"),
      barwidth = unit(5.5, "mm"),
      nbin = 300
    )
  ) +
  scale_linewidth_manual(
    values = c(
      "< 0.05"    = 0.35,
      "0.05–0.10" = 0.70,
      "0.10–0.20" = 1.20,
      "≥ 0.20"    = 1.80
    )
  ) +
  scale_linetype_manual(
    values = c("P < 0.05" = "solid", "P ≥ 0.05" = "dashed")
  ) +
  scale_colour_manual(
    values = c(
      "Bacteria · Half"    = "#f39b7f",
      "Bacteria · Control" = "#3c5488",
      "Bacteria · Double"  = "#4dbbd5",
      "Fungi · Half"       = "#f39b7f",
      "Fungi · Control"    = "#3c5488",
      "Fungi · Double"     = "#4dbbd5"
    )
  ) +
  guides(
    colour = guide_legend(order = 1, title = NULL),
    linetype = guide_legend(order = 2, title = "Partial Mantel significance"),
    linewidth = guide_legend(order = 3, title = "|Partial Mantel r|"),
    fill = guide_colorbar(order = 4)
  ) +
  expand_limits(
    x = max(env_endpoint_df$x_lab, na.rm = TRUE) + 0.8
  ) +
  coord_cartesian(clip = "off") +
  theme(
    legend.position = "right",
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 8),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    plot.margin = margin(20, 100, 20, 140)
  )

p_pm_one

############################################################
# 5. Re-layout legends and axis labels
############################################################

# 5.1 Main matrix panel: no legends, no left-side y labels
p_matrix <- p_pm_one +
  guides(
    fill = "none",
    colour = "none",
    linetype = "none",
    linewidth = "none"
  ) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 8),
    plot.margin = margin(10, 0, 10, 0)
  )

# 5.2 Fill legend only
p_fill_only <- p_pm_one +
  guides(
    colour = "none",
    linetype = "none",
    linewidth = "none",
    fill = guide_colorbar(
      title = "Pearson's r\n(env–env)",
      title.position = "top",
      barheight = unit(92, "mm"),
      barwidth = unit(5.5, "mm"),
      nbin = 300
    )
  ) +
  theme(
    legend.position = "left",
    legend.justification = c(1, 0.5),
    legend.box.margin = margin(0, 0, 0, 0),
    legend.margin = margin(0, 0, 0, 0)
  )

leg_fill <- cowplot::get_legend(p_fill_only)

# 5.3 Right-side legends only
p_right_only <- p_pm_one +
  guides(
    fill = "none",
    colour = guide_legend(order = 1, title = NULL),
    linetype = guide_legend(order = 2, title = "Partial Mantel significance"),
    linewidth = guide_legend(order = 3, title = "|Partial Mantel r|")
  ) +
  theme(
    legend.position = "right",
    legend.box.margin = margin(0, 0, 0, 0),
    legend.margin = margin(0, 0, 0, 0)
  )

leg_right <- cowplot::get_legend(p_right_only)

# 5.4 Assemble final figure
p_pm_final <- cowplot::plot_grid(
  cowplot::ggdraw(leg_fill),
  p_matrix,
  cowplot::ggdraw(leg_right),
  nrow = 1,
  rel_widths = c(0.10, 1, 0.34),
  align = "h",
  axis = "tb"
)

p_pm_final

############################################################
# 6. Save figure
############################################################
ggsave(
  filename = file.path(output_dir, "Fig2c_partial_mantel_one_matrix_FINAL.pdf"),
  plot = p_pm_final,
  width = 10.4,
  height = 7.8,
  units = "in",
  dpi = 300,
  device = cairo_pdf
)

ggsave(
  filename = file.path(output_dir, "Fig2c_partial_mantel_one_matrix_FINAL.png"),
  plot = p_pm_final,
  width = 10.4,
  height = 7.8,
  units = "in",
  dpi = 300
)