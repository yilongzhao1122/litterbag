############################################################
# Fig. 1a standardized coefficient plot
# Project: litterbag
#
# Expected project structure:
# litterbag/
# ├── code/
# ├── data/
# │   ├── treat2.csv
# │   ├── env_clear.csv
# │   └── LMM_picture_WM.RData
# └── output/
#
# Variable order follows the manuscript text.
# This version does NOT calculate, save, or plot significance labels.
############################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(ggplot2)
  library(lmerTest)
  library(patchwork)
})

############################################################
# 0. Helper functions
############################################################

locate_project_root <- function() {
  candidates <- c(
    getwd(),
    dirname(getwd()),
    dirname(dirname(getwd()))
  )
  
  candidates <- unique(normalizePath(candidates, winslash = "/", mustWork = FALSE))
  
  for (x in candidates) {
    if (dir.exists(file.path(x, "code")) && dir.exists(file.path(x, "data"))) {
      return(x)
    }
  }
  
  stop(
    "Project root not found. Please run this script from the project root ",
    "or from the code/ directory."
  )
}

check_required_files <- function(files) {
  missing_files <- files[!file.exists(files)]
  if (length(missing_files) > 0) {
    stop(
      "Missing required input file(s):\n",
      paste(missing_files, collapse = "\n")
    )
  }
}

theme_coef_panel <- function(base_size = 9) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
      strip.background = element_blank(),
      legend.title = element_blank(),
      legend.key = element_blank(),
      axis.title = element_text(face = "bold"),
      plot.title = element_text(face = "bold", hjust = 0.5)
    )
}

extract_rate_df <- function(df, obj_name = "rate_obj") {
  if (!("sample" %in% colnames(df))) {
    stop(paste0(obj_name, " lacks the 'sample' column."))
  }
  if (!("loss" %in% colnames(df))) {
    stop(paste0(obj_name, " lacks the 'loss' column."))
  }
  
  df %>%
    dplyr::select(sample, loss) %>%
    mutate(
      sample = as.character(sample),
      loss = as.numeric(loss)
    )
}

fit_std_lmm <- function(var_name, dat) {
  sub_df <- dat %>% filter(variable == var_name)
  
  mod <- lmer(
    value_z ~ Precipitation + (1 | block) + (1 | year),
    data = sub_df
  )
  
  sm <- coef(summary(mod))
  
  tibble(
    variable        = var_name,
    estimate_half   = sm["Precipitationhalf",   "Estimate"],
    se_half         = sm["Precipitationhalf",   "Std. Error"],
    estimate_double = sm["Precipitationdouble", "Estimate"],
    se_double       = sm["Precipitationdouble", "Std. Error"]
  )
}

make_coef_panel <- function(df_sub, var_order, panel_title, show_legend = FALSE) {
  df_sub <- df_sub %>%
    filter(label %in% var_order) %>%
    mutate(label = factor(label, levels = rev(var_order)))
  
  xmin_raw <- min(df_sub$conf.low, na.rm = TRUE)
  xmax_raw <- max(df_sub$conf.high, na.rm = TRUE)
  
  x_span <- xmax_raw - xmin_raw
  if (is.na(x_span) || x_span == 0) {
    x_span <- 0.1
  }
  
  x_abs_max <- max(
    abs(df_sub$conf.low),
    abs(df_sub$conf.high),
    abs(df_sub$estimate),
    na.rm = TRUE
  )
  
  x_limit <- x_abs_max * 1.08
  x_breaks <- pretty(c(-x_limit, x_limit), n = 5)
  
  ggplot(df_sub, aes(x = estimate, y = label, colour = contrast)) +
    geom_vline(
      xintercept = 0,
      linetype = 2,
      linewidth = 0.4,
      colour = "grey40"
    ) +
    geom_errorbarh(
      aes(xmin = conf.low, xmax = conf.high),
      height = 0.18,
      linewidth = 0.45,
      position = position_dodge(width = 0.55)
    ) +
    geom_point(
      size = 2.2,
      position = position_dodge(width = 0.55)
    ) +
    scale_colour_manual(values = c(
      "Half vs Control"   = "#F39B7FFF",
      "Double vs Control" = "#4DBBD5FF"
    )) +
    scale_x_continuous(breaks = x_breaks) +
    labs(
      title = panel_title,
      x = "Standardized coefficient",
      y = NULL
    ) +
    coord_cartesian(
      xlim = c(-x_limit, x_limit),
      clip = "off",
      expand = FALSE
    ) +
    theme_coef_panel(9) +
    theme(
      legend.position = if (show_legend) "top" else "none",
      plot.margin = margin(5.5, 20, 5.5, 5.5)
    )
}

############################################################
# 1. Locate directories and input files
############################################################

root_dir <- locate_project_root()
data_dir <- file.path(root_dir, "data")
output_dir <- file.path(root_dir, "output/coefficient")

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

treat_file <- file.path(data_dir, "treat2.csv")
env_file <- file.path(data_dir, "env_clear.csv")
rdata_file <- file.path(data_dir, "LMM_picture_WM.RData")

check_required_files(c(treat_file, env_file, rdata_file))

############################################################
# 2. Read metadata and environmental data
############################################################

treat <- read.csv(
  treat_file,
  row.names = 1,
  check.names = FALSE
)

if (!("sample" %in% colnames(treat))) {
  treat$sample <- rownames(treat)
}

treat <- treat %>%
  filter(Clip == "uncliping")

env <- read.csv(
  env_file,
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)

env <- env %>%
  filter(Clip == "uncliping") %>%
  dplyr::select(-any_of(c(
    "precipitation_Oct_cm",
    "precipitation_annual_cm",
    "root_biomass"
  )))

env$year <- as.factor(env$year)

# Unwarming, precipitation-only analysis
env2_P <- env %>%
  filter(Warm == "unwarming")

############################################################
# 3. Load processed litterbag decomposition objects from .RData
############################################################

rdata_env <- new.env(parent = emptyenv())
load(rdata_file, envir = rdata_env)

required_objects <- c("rate1", "rate2", "rate3")
missing_objects <- required_objects[!vapply(required_objects, exists, logical(1), envir = rdata_env)]

if (length(missing_objects) > 0) {
  stop(
    "The following required objects are missing from LMM_picture_WM.RData:\n",
    paste(missing_objects, collapse = "\n")
  )
}

rate1 <- get("rate1", envir = rdata_env)
rate2 <- get("rate2", envir = rdata_env)
rate3 <- get("rate3", envir = rdata_env)

############################################################
# 4. Variable groups and display labels
############################################################

edaphic_vars <- c(
  "moisture_samplingmonth",
  "annual_moisture",
  "temperature_annual",
  "NO3.N",
  "NH4.N",
  "TN",
  "TC",
  "pH"
)

veg_prod_vars <- c(
  "FlTotl",
  "FlC4",
  "FlC3",
  "plant.richness"
)

flux_vars <- c(
  "ER_annualmean",
  "Autotrophic",
  "Heterotrophic",
  "total_soil_respiration",
  "GPP_annualmean",
  "NEE_annualmean"
)

plot_vars <- c(edaphic_vars, veg_prod_vars, flux_vars)

pretty_labels <- c(
  "moisture_samplingmonth" = "Sampling-month moisture",
  "annual_moisture"        = "Annual moisture",
  "temperature_annual"     = "Annual temperature",
  "NO3.N"                  = "Soil nitrate nitrogen",
  "NH4.N"                  = "Soil ammonium nitrogen",
  "TN"                     = "Soil total nitrogen",
  "TC"                     = "Soil total carbon",
  "pH"                     = "Soil pH",
  "FlTotl"                 = "Total aboveground plant biomass",
  "FlC4"                   = "C4 plant biomass",
  "FlC3"                   = "C3 plant biomass",
  "plant.richness"         = "Plant richness",
  "GPP_annualmean"         = "Gross primary productivity",
  "NEE_annualmean"         = "Net ecosystem exchange",
  "ER_annualmean"          = "Ecosystem respiration",
  "Autotrophic"            = "Autotrophic respiration",
  "Heterotrophic"          = "Heterotrophic respiration",
  "total_soil_respiration" = "Total Rs",
  "loss"                   = "Cellulose mass loss"
)

group_labels <- c(
  "moisture_samplingmonth" = "Edaphic variables",
  "annual_moisture"        = "Edaphic variables",
  "temperature_annual"     = "Edaphic variables",
  "NO3.N"                  = "Edaphic variables",
  "NH4.N"                  = "Edaphic variables",
  "TN"                     = "Edaphic variables",
  "TC"                     = "Edaphic variables",
  "pH"                     = "Edaphic variables",
  "FlTotl"                 = "Vegetation and productivity",
  "FlC4"                   = "Vegetation and productivity",
  "FlC3"                   = "Vegetation and productivity",
  "plant.richness"         = "Vegetation and productivity",
  "GPP_annualmean"         = "Carbon fluxes",
  "NEE_annualmean"         = "Carbon fluxes",
  "ER_annualmean"          = "Carbon fluxes",
  "Autotrophic"            = "Carbon fluxes",
  "Heterotrophic"          = "Carbon fluxes",
  "total_soil_respiration" = "Carbon fluxes",
  "loss"                   = "Decomposition"
)

############################################################
# 5. Reshape environmental data
############################################################

env_long <- env2_P %>%
  pivot_longer(
    cols = all_of(plot_vars),
    names_to = "variable",
    values_to = "value"
  ) %>%
  mutate(
    Precipitation = factor(Precipitation, levels = c("normal", "half", "double")),
    block = factor(block),
    year = factor(year)
  ) %>%
  group_by(variable) %>%
  mutate(value_z = as.numeric(scale(value))) %>%
  ungroup()

############################################################
# 6. Fit LMMs for environmental variables
############################################################

coef_env <- purrr::map_dfr(plot_vars, fit_std_lmm, dat = env_long)

############################################################
# 7. Prepare decomposition data and fit LMM
############################################################

rate1_use <- extract_rate_df(rate1, "rate1")
rate2_use <- extract_rate_df(rate2, "rate2")
rate3_use <- extract_rate_df(rate3, "rate3")

rate_raw <- bind_rows(rate1_use, rate2_use, rate3_use)

treat_litter <- treat %>%
  filter(type == "litter") %>%
  dplyr::select(sample, Precipitation, Warm, block, year, type)

rate <- rate_raw %>%
  left_join(treat_litter, by = "sample") %>%
  mutate(
    Precipitation = factor(Precipitation, levels = c("normal", "half", "double")),
    Warm = factor(Warm),
    block = factor(block),
    year = factor(year)
  )

rate2016 <- rate %>%
  filter(Warm == "unwarming", year == "2016") %>%
  mutate(loss_z = as.numeric(scale(loss)))

mod_rate2016 <- lmer(loss_z ~ Precipitation + (1 | block), data = rate2016)
sm_rate2016 <- coef(summary(mod_rate2016))

coef_rate <- tibble(
  variable        = "loss",
  estimate_half   = sm_rate2016["Precipitationhalf",   "Estimate"],
  se_half         = sm_rate2016["Precipitationhalf",   "Std. Error"],
  estimate_double = sm_rate2016["Precipitationdouble", "Estimate"],
  se_double       = sm_rate2016["Precipitationdouble", "Std. Error"]
)

############################################################
# 8. Combine coefficients
############################################################

coef_all <- bind_rows(coef_env, coef_rate) %>%
  pivot_longer(
    cols = c(
      estimate_half, se_half,
      estimate_double, se_double
    ),
    names_to = c(".value", "contrast_key"),
    names_pattern = "(estimate|se)_(half|double)"
  ) %>%
  mutate(
    contrast = recode(
      contrast_key,
      "half" = "Half vs Control",
      "double" = "Double vs Control"
    ),
    conf.low = estimate - 1.96 * se,
    conf.high = estimate + 1.96 * se,
    label = pretty_labels[variable],
    group = group_labels[variable]
  ) %>%
  mutate(
    group = factor(
      group,
      levels = c(
        "Edaphic variables",
        "Vegetation and productivity",
        "Carbon fluxes",
        "Decomposition"
      )
    )
  )

############################################################
# 9. Variable order in the figure
############################################################

coef_order <- c(
  "Annual moisture",
  "Annual temperature",
  "Soil pH",
  "Soil nitrate nitrogen",
  "Soil ammonium nitrogen",
  "Soil total nitrogen",
  "Soil total carbon",
  "Sampling-month moisture",
  "Total aboveground plant biomass",
  "C4 plant biomass",
  "C3 plant biomass",
  "Plant richness",
  "Gross primary productivity",
  "Net ecosystem exchange",
  "Total Rs",
  "Autotrophic respiration",
  "Heterotrophic respiration",
  "Ecosystem respiration",
  "Cellulose mass loss"
)

coef_all <- coef_all %>%
  mutate(label = factor(label, levels = rev(coef_order)))

group_edaphic <- c(
  "Annual moisture",
  "Annual temperature",
  "Soil pH",
  "Soil nitrate nitrogen",
  "Soil ammonium nitrogen",
  "Soil total nitrogen",
  "Soil total carbon",
  "Sampling-month moisture"
)

group_plant <- c(
  "Total aboveground plant biomass",
  "C4 plant biomass",
  "C3 plant biomass",
  "Plant richness"
)

group_flux <- c(
  "Gross primary productivity",
  "Ecosystem respiration",
  "Total Rs",
  "Autotrophic respiration",
  "Heterotrophic respiration",
  "Net ecosystem exchange"
)

group_loss <- c(
  "Cellulose mass loss"
)

############################################################
# 10. Draw panels
############################################################

p_coef_edaphic <- make_coef_panel(
  df_sub = coef_all,
  var_order = group_edaphic,
  panel_title = "Edaphic variables",
  show_legend = TRUE
)

p_coef_plant <- make_coef_panel(
  df_sub = coef_all,
  var_order = group_plant,
  panel_title = "Vegetation and productivity"
)

p_coef_flux <- make_coef_panel(
  df_sub = coef_all,
  var_order = group_flux,
  panel_title = "Carbon fluxes"
)

p_coef_loss <- make_coef_panel(
  df_sub = coef_all,
  var_order = group_loss,
  panel_title = "Decomposition"
)

############################################################
# 11. Combine and save
############################################################

fig_coef_split <- (
  p_coef_edaphic /
    p_coef_plant /
    p_coef_flux /
    p_coef_loss
) + plot_layout(heights = c(1.6, 1.2, 0.9, 0.55))

ggsave(
  filename = file.path(output_dir, "Fig1a_standardized_coefficient_plot.pdf"),
  plot = fig_coef_split,
  width = 7.6,
  height = 10.8,
  units = "in",
  dpi = 300
)

ggsave(
  filename = file.path(output_dir, "Fig1a_standardized_coefficient_plot.png"),
  plot = fig_coef_split,
  width = 7.6,
  height = 10.8,
  units = "in",
  dpi = 300
)

message("Done. Figure files were saved to: ", output_dir)