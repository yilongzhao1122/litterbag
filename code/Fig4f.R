############################################################
# Fig. 4f: gMEND-predicted C-pool and C-flux variables
# Project: litterbag
#
# Expected project structure:
# litterbag/
# ├── code/
# │   └── Fig4f_CPool.R
# ├── model/
# │   ├── Model_half/
# │   ├── Model_Control/
# │   └── Model_double/
# └── output/
#
# This script reads daily gMEND output files for half, control, and
# double precipitation treatments, summarizes predicted C-related
# variables, and generates the bar plots used for Fig. 4f.
#
# Input files used:
# - model/Model_half/userio/MEND3/gMEND/H_VAR_day.out
# - model/Model_half/userio/MEND3/gMEND/H_rate_day.out
# - model/Model_half/userio/MEND3/gMEND/H_FLX_day.out
# - model/Model_Control/userio/MEND3/gMEND/H_VAR_day.out
# - model/Model_Control/userio/MEND3/gMEND/H_rate_day.out
# - model/Model_Control/userio/MEND3/gMEND/H_FLX_day.out
# - model/Model_double/userio/MEND3/gMEND/D_VAR_day.out
# - model/Model_double/userio/MEND3/gMEND/D_rate_day.out
# - model/Model_double/userio/MEND3/gMEND/D_FLX_day.out
#
# Output files:
# - output/mend/Fig4f_gMEND_Rh.pdf/png
# - output/mend/Fig4f_gMEND_rate.pdf/png
# - output/mend/Fig4f_gMEND_active_biomass_fraction.pdf/png
# - output/mend/Fig4f_gMEND_CUE.pdf/png
# - output/mend/Fig4f_gMEND_MB_MBA.pdf/png
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
  if (dir.exists(file.path(x, "code")) && dir.exists(file.path(x, "model"))) {
    project_root <- x
    break
  }
}

if (is.null(project_root)) {
  stop("Project root not found. Please run this script from the project root or from the code/ directory.")
}

model_dir <- file.path(project_root, "model")
output_dir <- file.path(project_root, "output", "mend")

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

suppressPackageStartupMessages({
  library(reshape2)
  library(ggplot2)
  library(Rmisc)
})

######################
# 1. input
######################

# double
file_var_D <- file.path(model_dir, "Model_double", "userio", "MEND3", "gMEND", "D_VAR_day.out")
file_rate_D <- file.path(model_dir, "Model_double", "userio", "MEND3", "gMEND", "D_rate_day.out")
file_flx_D <- file.path(model_dir, "Model_double", "userio", "MEND3", "gMEND", "D_FLX_day.out")

if (!file.exists(file_var_D)) stop("Missing file: ", file_var_D)
if (!file.exists(file_rate_D)) stop("Missing file: ", file_rate_D)
if (!file.exists(file_flx_D)) stop("Missing file: ", file_flx_D)

var_D <- read.delim(file_var_D, sep = "", skip = 1)
colnames(var_D)
var_D <- var_D[, c(1, 3:7, 9:16)]

rate_D <- read.delim(file_rate_D, sep = "", skip = 1)
colnames(rate_D)
rate_D <- rate_D[, c(1:6, 10, 13, 14)]

Rh_D <- read.delim(file_flx_D, sep = "", skip = 1)
Rh_D <- Rh_D[, c(1, 31)]

# control
file_var_C <- file.path(model_dir, "Model_Control", "userio", "MEND3", "gMEND", "H_VAR_day.out")
file_rate_C <- file.path(model_dir, "Model_Control", "userio", "MEND3", "gMEND", "H_rate_day.out")
file_flx_C <- file.path(model_dir, "Model_Control", "userio", "MEND3", "gMEND", "H_FLX_day.out")

if (!file.exists(file_var_C)) stop("Missing file: ", file_var_C)
if (!file.exists(file_rate_C)) stop("Missing file: ", file_rate_C)
if (!file.exists(file_flx_C)) stop("Missing file: ", file_flx_C)

var_C <- read.delim(file_var_C, sep = "", skip = 1)
colnames(var_C)
var_C <- var_C[, c(1, 3:7, 9:16)]

rate_C <- read.delim(file_rate_C, sep = "", skip = 1)
rate_C <- rate_C[, c(1:6, 10, 13, 14)]

Rh_C <- read.delim(file_flx_C, sep = "", skip = 1)
colnames(Rh_C)
Rh_C <- Rh_C[, c(1, 31)]

# half
file_var_H <- file.path(model_dir, "Model_half", "userio", "MEND3", "gMEND", "H_VAR_day.out")
file_rate_H <- file.path(model_dir, "Model_half", "userio", "MEND3", "gMEND", "H_rate_day.out")
file_flx_H <- file.path(model_dir, "Model_half", "userio", "MEND3", "gMEND", "H_FLX_day.out")

if (!file.exists(file_var_H)) stop("Missing file: ", file_var_H)
if (!file.exists(file_rate_H)) stop("Missing file: ", file_rate_H)
if (!file.exists(file_flx_H)) stop("Missing file: ", file_flx_H)

var_H <- read.delim(file_var_H, sep = "", skip = 1)
colnames(var_H)
var_H <- var_H[, c(1, 3:7, 9:16)]

rate_H <- read.delim(file_rate_H, sep = "", skip = 1)
rate_H <- rate_H[, c(1:6, 10, 13, 14)]

Rh_H <- read.delim(file_flx_H, sep = "", skip = 1)
Rh_H <- Rh_H[, c(1, 31)]

mean(Rh_D$CO2_gmo)
mean(Rh_C$CO2_gmo)
mean(Rh_H$CO2_gmo)

data_C <- merge(var_C, rate_C, by = "Day")
data_C <- merge(data_C, Rh_C, by = "Day")

data_D <- merge(var_D, rate_D, by = "Day")
data_D <- merge(data_D, Rh_D, by = "Day")

data_H <- merge(var_H, rate_H, by = "Day")
data_H <- merge(data_H, Rh_H, by = "Day")

data_C$Precipitation <- "control"
data_D$Precipitation <- "double"
data_H$Precipitation <- "half"

data <- rbind(data_C, data_D, data_H)
colnames(data)

data <- data[, c(1, 23, 15, 16, 17, 18, 21, 22, 8:10, 24)]
colnames(data) <- c("Day", "Rh", "POM1", "POM2", "MOM", "DOM", "Active biomass", "CUE", "MB", "MBA", "MBD", "Precipitation")

data$Rh <- data$Rh * 24 * 10000 * 100 / 1000
data$POM1 <- data$POM1 * 1000000 / 1000
data$POM2 <- data$POM2 * 1000000 / 1000
data$MOM <- data$MOM * 1000000 / 1000
data$DOM <- data$DOM * 1000000 / 1000
data$`Active biomass` <- data$`Active biomass` * 100

# double
out_D_gMEND <- read.delim(file.path(model_dir, "Model_double", "userio", "MEND3", "gMEND", "D_SIM_obs.out"), sep = "", skip = 1)
out_D_gMEND_R2 <- out_D_gMEND[1:22, 3:6]
out_D_gMEND_R2$precipitation <- "double"
out_D_gMEND_R2$OBS_avg <- as.numeric(as.character(out_D_gMEND_R2$OBS_avg))
out_D_gMEND_R2$SIM_avg <- as.numeric(as.character(out_D_gMEND_R2$SIM_avg))

# control
out_C_gMEND <- read.delim(file.path(model_dir, "Model_Control", "userio", "MEND3", "gMEND", "H_SIM_obs.out"), sep = "", skip = 1)
out_C_gMEND_R2 <- out_C_gMEND[1:37, 3:6]
out_C_gMEND_R2$precipitation <- "control"
out_C_gMEND_R2$OBS_avg <- as.numeric(as.character(out_C_gMEND_R2$OBS_avg))
out_C_gMEND_R2$SIM_avg <- as.numeric(as.character(out_C_gMEND_R2$SIM_avg))

# half
out_H_gMEND <- read.delim(file.path(model_dir, "Model_half", "userio", "MEND3", "gMEND", "H_SIM_obs.out"), sep = "", skip = 1)
out_H_gMEND_R2 <- out_H_gMEND[1:22, 3:6]
out_H_gMEND_R2$precipitation <- "half"
out_H_gMEND_R2$OBS_avg <- as.numeric(as.character(out_H_gMEND_R2$OBS_avg))
out_H_gMEND_R2$SIM_avg <- as.numeric(as.character(out_H_gMEND_R2$SIM_avg))

mean(out_D_gMEND_R2$SIM_avg)
mean(out_C_gMEND_R2$SIM_avg)
mean(out_H_gMEND_R2$SIM_avg)

out_gMEND_R2 <- rbind(out_D_gMEND_R2, out_C_gMEND_R2, out_H_gMEND_R2)
out_gMEND_R2$OBS_avg <- out_gMEND_R2$OBS_avg * 24 * 10000 * 100 / 1000
out_gMEND_R2$SIM_avg <- out_gMEND_R2$SIM_avg * 24 * 10000 * 100 / 1000

model <- lm(OBS_avg ~ precipitation, out_gMEND_R2)
test1 <- aov(model)
TukeyHSD(test1)

(mean(out_D_gMEND_R2$SIM_avg) - mean(out_C_gMEND_R2$SIM_avg)) / mean(out_C_gMEND_R2$SIM_avg)

######################
# 2. picture
######################
data2 <- melt(data, id = c("Day", "Precipitation"))

data_filt <- subset(data, CUE > -0.2 & Rh >= 1 & Rh <= 9)
data2_filt <- subset(
  data2,
  (variable == "CUE" & value > -0.2) |
    (variable == "Rh" & value >= 1 & value <= 9) |
    !(variable %in% c("Rh", "CUE"))
)

rh_control <- subset(data2_filt, Precipitation == "control" & variable == "Rh")

data3 <- summarySE(data2_filt, measurevar = "value", groupvars = c("variable", "Precipitation"))
data3$Precipitation <- factor(data3$Precipitation, levels = c("half", "control", "double"))

ggplot(data3, aes(variable, value, fill = Precipitation)) +
  facet_wrap(~variable, scales = "free", nrow = 3) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7) +
  geom_errorbar(aes(ymin = value - se, ymax = value + se),
                width = 0.25, linewidth = 0.3, position = position_dodge(0.9)) +
  scale_fill_manual(values = c("#F39B7FFF", "#3C5488FF", "#4DBBD5FF")) +
  labs(title = NULL, x = NULL, y = NULL, fill = NULL, size = 0.3) +
  scale_y_continuous(expand = expansion(mult = c(0, .15))) +
  theme(
    panel.grid = element_blank(),
    legend.position = "right",
    panel.background = element_rect(fill = "transparent", color = "black"),
    strip.background = element_rect(color = "transparent", fill = "transparent"),
    axis.text.y = element_text(size = 7)
  ) +
  theme(strip.text = element_text(size = 7), axis.ticks.x = element_blank())

p_Rh <- ggplot(subset(data3, variable == "Rh"), aes(variable, value, fill = Precipitation)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7) +
  geom_errorbar(aes(ymin = value - se, ymax = value + se),
                width = 0.25, linewidth = 0.3, position = position_dodge(0.75)) +
  scale_fill_manual(values = c("#F39B7FFF", "#3C5488FF", "#4DBBD5FF")) +
  labs(title = NULL, x = NULL, y = expression(R[h] * " " * (g * " " * C * " " * m^-2 * d^-1)), fill = NULL, size = 0.3) +
  scale_y_continuous(expand = expansion(mult = c(0, .15))) +
  theme(
    panel.grid = element_blank(),
    legend.position = "right",
    panel.background = element_rect(fill = "transparent", color = "black"),
    strip.background = element_rect(color = "transparent", fill = "transparent"),
    axis.text.y = element_text(size = 7)
  ) +
  theme(strip.text = element_text(size = 7), axis.ticks.x = element_blank()) +
  scale_x_discrete(labels = c(expression(italic(R)[h])))
p_Rh

p_POM2 <- ggplot(subset(data3, variable == "POM2"), aes(variable, value, fill = Precipitation)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7) +
  geom_errorbar(aes(ymin = value - se, ymax = value + se),
                width = 0.25, linewidth = 0.3, position = position_dodge(0.75)) +
  scale_fill_manual(values = c("#F39B7FFF", "#3C5488FF", "#4DBBD5FF")) +
  labs(title = NULL, x = NULL, y = expression(R[h] * " " * (g * " " * C * " " * m^-2 * d^-1)), fill = NULL, size = 0.3) +
  scale_y_continuous(expand = expansion(mult = c(0, .15))) +
  theme(
    panel.grid = element_blank(),
    legend.position = "right",
    panel.background = element_rect(fill = "transparent", color = "black"),
    strip.background = element_rect(color = "transparent", fill = "transparent"),
    axis.text.y = element_text(size = 7)
  ) +
  theme(strip.text = element_text(size = 7), axis.ticks.x = element_blank()) +
  scale_x_discrete(labels = c(expression(POM[2])))
p_POM2

p_rate <- ggplot(subset(data3, variable %in% c("POM1", "MOM")), aes(variable, value, fill = Precipitation)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7) +
  geom_errorbar(aes(ymin = value - se, ymax = value + se),
                width = 0.25, linewidth = 0.3, position = position_dodge(0.75)) +
  scale_fill_manual(values = c("#F39B7FFF", "#3C5488FF", "#4DBBD5FF")) +
  labs(title = NULL, x = NULL, y = expression("Decomposition rate" * " " * (g * " " * C * " " * m^-3 * d^-1)), fill = NULL, size = 0.3) +
  scale_y_continuous(expand = expansion(mult = c(0, .15))) +
  theme(
    panel.grid = element_blank(),
    legend.position = "none",
    panel.background = element_rect(fill = "transparent", color = "black"),
    strip.background = element_rect(color = "transparent", fill = "transparent"),
    axis.text.y = element_text(size = 7)
  ) +
  theme(strip.text = element_text(size = 7), axis.ticks.x = element_blank()) +
  scale_x_discrete(labels = c(expression(POM[1]), expression(MAOM)))
p_rate

p_MB <- ggplot(subset(data3, variable %in% c("Active biomass")), aes(variable, value, fill = Precipitation)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7) +
  geom_errorbar(aes(ymin = value - se, ymax = value + se),
                width = 0.25, linewidth = 0.3, position = position_dodge(0.75)) +
  scale_fill_manual(values = c("#F39B7FFF", "#3C5488FF", "#4DBBD5FF")) +
  labs(title = NULL, x = NULL, y = expression("Percentage" * " " * "(%)"), fill = NULL, size = 0.3) +
  scale_y_continuous(expand = expansion(mult = c(0, .15))) +
  theme(
    panel.grid = element_blank(),
    legend.position = "none",
    panel.background = element_rect(fill = "transparent", color = "black"),
    strip.background = element_rect(color = "transparent", fill = "transparent"),
    axis.text.y = element_text(size = 7)
  ) +
  theme(strip.text = element_text(size = 7), axis.ticks.x = element_blank())
p_MB

p_CUE <- ggplot(subset(data3, variable %in% c("CUE")), aes(variable, value, fill = Precipitation)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7) +
  geom_errorbar(aes(ymin = value - se, ymax = value + se),
                width = 0.25, linewidth = 0.3, position = position_dodge(0.75)) +
  scale_fill_manual(values = c("#F39B7FFF", "#3C5488FF", "#4DBBD5FF")) +
  labs(title = NULL, x = NULL, y = expression("Carbon use efficiency"), fill = NULL, size = 0.3) +
  scale_y_continuous(expand = expansion(mult = c(0, .15))) +
  theme(
    panel.grid = element_blank(),
    legend.position = "none",
    panel.background = element_rect(fill = "transparent", color = "black"),
    strip.background = element_rect(color = "transparent", fill = "transparent"),
    axis.text.y = element_text(size = 7)
  ) +
  theme(strip.text = element_text(size = 7), axis.ticks.x = element_blank())
p_CUE

p_MB2 <- ggplot(subset(data3, variable %in% c("MB", "MBA")), aes(variable, value, fill = Precipitation)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7) +
  geom_errorbar(aes(ymin = value - se, ymax = value + se),
                width = 0.25, linewidth = 0.3, position = position_dodge(0.75)) +
  scale_fill_manual(values = c("#F39B7FFF", "#3C5488FF", "#4DBBD5FF")) +
  labs(title = NULL, x = NULL, y = expression("mg C" * " " * cm^-3), fill = NULL, size = 0.3) +
  scale_y_continuous(expand = expansion(mult = c(0, .15))) +
  theme(
    panel.grid = element_blank(),
    legend.position = "none",
    panel.background = element_rect(fill = "transparent", color = "black"),
    strip.background = element_rect(color = "transparent", fill = "transparent"),
    axis.text.y = element_text(size = 7)
  ) +
  theme(strip.text = element_text(size = 7), axis.ticks.x = element_blank())
p_MB2

######################
# 3. Save figures
######################
ggsave(
  plot = p_Rh,
  filename = file.path(output_dir, "Fig4f_gMEND_Rh.pdf"),
  width = 2.1, height = 2.6, units = "in", device = "pdf"
)

ggsave(
  plot = p_Rh,
  filename = file.path(output_dir, "Fig4f_gMEND_Rh.png"),
  width = 2.1, height = 2.6, units = "in", dpi = 300
)

ggsave(
  plot = p_rate,
  filename = file.path(output_dir, "Fig4f_gMEND_rate.pdf"),
  width = 4.6, height = 2.6, units = "in", device = "pdf"
)

ggsave(
  plot = p_rate,
  filename = file.path(output_dir, "Fig4f_gMEND_rate.png"),
  width = 4.6, height = 2.6, units = "in", dpi = 300
)

ggsave(
  plot = p_MB,
  filename = file.path(output_dir, "Fig4f_gMEND_active_biomass_fraction.pdf"),
  width = 2.15, height = 2.6, units = "in", device = "pdf"
)

ggsave(
  plot = p_MB,
  filename = file.path(output_dir, "Fig4f_gMEND_active_biomass_fraction.png"),
  width = 2.15, height = 2.6, units = "in", dpi = 300
)

ggsave(
  plot = p_CUE,
  filename = file.path(output_dir, "Fig4f_gMEND_CUE.pdf"),
  width = 2.15, height = 2.6, units = "in", device = "pdf"
)

ggsave(
  plot = p_CUE,
  filename = file.path(output_dir, "Fig4f_gMEND_CUE.png"),
  width = 2.15, height = 2.6, units = "in", dpi = 300
)

ggsave(
  plot = p_MB2,
  filename = file.path(output_dir, "Fig4f_gMEND_MB_MBA.pdf"),
  width = 3.5, height = 2.6, units = "in", device = "pdf"
)

ggsave(
  plot = p_MB2,
  filename = file.path(output_dir, "Fig4f_gMEND_MB_MBA.png"),
  width = 3.5, height = 2.6, units = "in", dpi = 300
)