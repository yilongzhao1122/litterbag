############################################################
# Fig. 4d: gMEND model performance for heterotrophic respiration
# Project: litterbag
#
# Expected project structure:
# litterbag/
# ├── code/
# │   └── Fig4d_gMEND_R2.R
# ├── model/
# │   ├── Model_half/
# │   ├── Model_Control/
# │   └── Model_double/
# └── output/
#
# This script reads simulated-vs-observed MEND output files from the
# gMEND model, compares model performance across half, control, and
# double precipitation treatments, and generates Fig. 4d.
#
# Input files used:
# - model/Model_half/userio/MEND3/gMEND/H_SIM_obs.out
# - model/Model_Control/userio/MEND3/gMEND/H_SIM_obs.out
# - model/Model_double/userio/MEND3/gMEND/D_SIM_obs.out
#
# Output files:
# - output/mend/Fig4d_gMEND_R2.pdf
# - output/mend/Fig4d_gMEND_R2.png
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
  library(ggplot2)
})

###############################
# 1. Read gMEND simulation results
###############################

# double
file_D <- file.path(model_dir, "Model_double", "userio", "MEND3", "gMEND", "D_SIM_obs.out")
if (!file.exists(file_D)) stop("Missing file: ", file_D)

out_D_gMEND <- read.delim(file_D, sep = "", skip = 1)
out_D_gMEND_R2 <- out_D_gMEND[1:22, 3:6]
out_D_gMEND_R2$precipitation <- "double"
out_D_gMEND_R2$OBS_avg <- as.numeric(as.character(out_D_gMEND_R2$OBS_avg))
out_D_gMEND_R2$SIM_avg <- as.numeric(as.character(out_D_gMEND_R2$SIM_avg))

# control
file_C <- file.path(model_dir, "Model_Control", "userio", "MEND3", "gMEND", "H_SIM_obs.out")
if (!file.exists(file_C)) stop("Missing file: ", file_C)

out_C_gMEND <- read.delim(file_C, sep = "", skip = 1)
out_C_gMEND_R2 <- out_C_gMEND[1:37, 3:6]
out_C_gMEND_R2$precipitation <- "control"
out_C_gMEND_R2$OBS_avg <- as.numeric(as.character(out_C_gMEND_R2$OBS_avg))
out_C_gMEND_R2$SIM_avg <- as.numeric(as.character(out_C_gMEND_R2$SIM_avg))

# half
file_H <- file.path(model_dir, "Model_half", "userio", "MEND3", "gMEND", "H_SIM_obs.out")
if (!file.exists(file_H)) stop("Missing file: ", file_H)

out_H_gMEND <- read.delim(file_H, sep = "", skip = 1)
out_H_gMEND_R2 <- out_H_gMEND[1:22, 3:6]
out_H_gMEND_R2$precipitation <- "half"
out_H_gMEND_R2$OBS_avg <- as.numeric(as.character(out_H_gMEND_R2$OBS_avg))
out_H_gMEND_R2$SIM_avg <- as.numeric(as.character(out_H_gMEND_R2$SIM_avg))

mean(out_D_gMEND_R2$SIM_avg)
mean(out_C_gMEND_R2$SIM_avg)
mean(out_H_gMEND_R2$SIM_avg)

t.test(out_C_gMEND_R2$SIM_avg, out_H_gMEND_R2$SIM_avg)
t.test(out_C_gMEND_R2$SIM_avg[1:22], out_H_gMEND_R2$SIM_avg, paired = TRUE)
wilcox.test(out_C_gMEND_R2$SIM_avg, out_H_gMEND_R2$SIM_avg)

t.test(out_C_gMEND_R2$SIM_avg, out_D_gMEND_R2$SIM_avg)
t.test(out_C_gMEND_R2$SIM_avg[1:22], out_D_gMEND_R2$SIM_avg, paired = TRUE)
wilcox.test(out_C_gMEND_R2$SIM_avg, out_D_gMEND_R2$SIM_avg)
t.test(out_C_gMEND_R2$OBS_avg, out_D_gMEND_R2$OBS_avg)

t.test(out_H_gMEND_R2$SIM_avg, out_D_gMEND_R2$SIM_avg)
t.test(out_H_gMEND_R2$SIM_avg, out_D_gMEND_R2$SIM_avg, paired = TRUE)
wilcox.test(out_H_gMEND_R2$SIM_avg, out_D_gMEND_R2$SIM_avg)

###############################
# 2. Prepare combined data
###############################

out_gMEND_R2 <- rbind(out_D_gMEND_R2, out_C_gMEND_R2, out_H_gMEND_R2)  # mg C cm^-3 h^-1
out_gMEND_R2$OBS_avg <- out_gMEND_R2$OBS_avg * 24 * 10000 * 100 / 1000  # g C m^-2 d^-1
out_gMEND_R2$SIM_avg <- out_gMEND_R2$SIM_avg * 24 * 10000 * 100 / 1000  # g C m^-2 d^-1

fit_D <- lm(data = out_D_gMEND_R2, SIM_avg ~ OBS_avg)
fit_C <- lm(data = out_C_gMEND_R2, SIM_avg ~ OBS_avg)
fit_H <- lm(data = out_H_gMEND_R2, SIM_avg ~ OBS_avg)

out_gMEND_R2$precipitation <- factor(
  out_gMEND_R2$precipitation,
  levels = c("half", "control", "double")
)

###############################
# 3. Plot Fig. 4d
###############################

p1 <- ggplot(out_gMEND_R2, aes(x = OBS_avg, y = SIM_avg, color = precipitation)) +
  geom_point(size = 1.5) +
  scale_color_manual(values = c("#F39B7FFF", "#3C5488FF", "#4DBBD5FF")) +
  geom_smooth(method = "lm", se = FALSE) +
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(color = "black", fill = "transparent"),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5),
    axis.title = element_text(size = 9),
    title = element_text(size = 10),
    strip.text = element_text(size = 7),
    strip.background = element_rect(color = "transparent", fill = "transparent")
  ) +
  scale_y_continuous(limit = c(0, 7.8)) +
  scale_x_continuous(limit = c(0, 7.8)) +
  labs(
    x = expression("Observed " * italic(R)[h] * " " * (g * " " * C * " " * m^-2 * " " * d^-1)),
    y = expression("Simulated " * italic(R)[h] * " " * (g * " " * C * " " * m^-2 * " " * d^-1)),
    title = ""
  )

label_fit <- data.frame(
  formula_H = sprintf(paste0("half (", "italic(R^2) == %.2f", ")"), summary(fit_H)$adj.r.squared),
  formula_C = sprintf(paste0("control (", "italic(R^2) == %.2f", ")"), summary(fit_C)$adj.r.squared),
  formula_D = sprintf(paste0("double (", "italic(R^2) == %.2f", ")"), summary(fit_D)$adj.r.squared)
)

p1 <- p1 +
  annotate(geom = "point", x = 0, y = 7.5, colour = "#F39B7FFF", size = 1.5) +
  annotate(geom = "point", x = 0, y = 6.8, colour = "#3C5488FF", size = 1.5) +
  annotate(geom = "point", x = 0, y = 6.1, colour = "#4DBBD5FF", size = 1.5) +
  geom_text(
    x = 0.5, y = 7.5, aes(label = formula_H), data = label_fit, parse = TRUE,
    hjust = 0, color = "black", show.legend = FALSE, size = 2.8
  ) +
  geom_text(
    x = 0.5, y = 6.8, aes(label = formula_C), data = label_fit, parse = TRUE,
    hjust = 0, color = "black", show.legend = FALSE, size = 2.8
  ) +
  geom_text(
    x = 0.5, y = 6.1, aes(label = formula_D), data = label_fit, parse = TRUE,
    hjust = 0, color = "black", show.legend = FALSE, size = 2.8
  )

p1

###############################
# 4. Save figure
###############################

ggsave(
  plot = p1,
  filename = file.path(output_dir, "Fig4d_gMEND_R2.pdf"),
  width = 2.8,
  height = 2.8,
  units = "in",
  device = "pdf"
)

ggsave(
  plot = p1,
  filename = file.path(output_dir, "Fig4d_gMEND_R2.png"),
  width = 2.8,
  height = 2.8,
  units = "in",
  dpi = 300
)

