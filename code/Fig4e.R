############################################################
# Fig. 4e: uncertainty reduction from tMEND to gMEND
# Project: litterbag
#
# Expected project structure:
# litterbag/
# ├── code/
# │   └── Fig4e_uncertainty.R
# ├── model/
# │   ├── Model_half/
# │   ├── Model_Control/
# │   └── Model_double/
# └── output/
#
# This script compares parameter uncertainty between tMEND and gMEND
# under half, control, and double precipitation treatments, and
# generates the boxplot figure used for Fig. 4e.
#
# Input files used:
# - model/Model_half/userio/MEND3/uncertainty/tMEND/H_UQpar.out
# - model/Model_half/userio/MEND3/uncertainty/gMEND/H_UQpar.out
# - model/Model_Control/userio/MEND3/uncertainty/tMEND/H_UQpar.out
# - model/Model_Control/userio/MEND3/uncertainty/gMEND/H_UQpar.out
# - model/Model_double/userio/MEND3/uncertainty/tMEND/D_UQpar.out
# - model/Model_double/userio/MEND3/uncertainty/gMEND/D_UQpar.out
#
# Output files:
# - output/mend/Fig4e_double_uncertainty.pdf
# - output/mend/Fig4e_control_uncertainty.pdf
# - output/mend/Fig4e_half_uncertainty.pdf
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
})

########################
# 1. calculation
########################
cal_cv <- function(x) {
  y <- na.omit(x)
  return(sd(y) / mean(y))
}

###################
# 2. double
###################
file_D_tMEND <- file.path(model_dir, "Model_double", "userio", "MEND3", "uncertainty", "tMEND", "D_UQpar.out")
file_D_gMEND <- file.path(model_dir, "Model_double", "userio", "MEND3", "uncertainty", "gMEND", "D_UQpar.out")

if (!file.exists(file_D_tMEND)) stop("Missing file: ", file_D_tMEND)
if (!file.exists(file_D_gMEND)) stop("Missing file: ", file_D_gMEND)

UQ_tMEND <- read.delim(file_D_tMEND, sep = "", skip = 1)
UQ_gMEND <- read.delim(file_D_gMEND, sep = "", skip = 1)

UQ_tMEND <- UQ_tMEND[, c(18:26, 28)]
UQ_gMEND <- UQ_gMEND[, c(18:26, 28)]

stat_tMEND <- as.data.frame(apply(UQ_tMEND, 2, cal_cv))
stat_gMEND <- as.data.frame(apply(UQ_gMEND, 2, cal_cv))

colnames(stat_tMEND) <- "tMEND"
colnames(stat_gMEND) <- "gMEND"

stat_tMEND$parameter <- rownames(stat_tMEND)
stat_gMEND$parameter <- rownames(stat_gMEND)

mean(apply(UQ_tMEND, 2, cal_cv))
mean(apply(UQ_gMEND, 2, cal_cv))

mean(apply(UQ_tMEND, 2, cal_cv)) - mean(apply(UQ_gMEND, 2, cal_cv))

t.test(apply(UQ_tMEND, 2, cal_cv), apply(UQ_gMEND, 2, cal_cv), paired = TRUE)
t.test(apply(UQ_tMEND, 2, cal_cv), apply(UQ_gMEND, 2, cal_cv))
wilcox.test(apply(UQ_tMEND, 2, cal_cv), apply(UQ_gMEND, 2, cal_cv))

stat <- merge(stat_tMEND, stat_gMEND, by = "parameter")
stat <- melt(stat, id = "parameter")

stat$parameter <- factor(stat$parameter, levels = stat_tMEND$parameter)

p1 <- ggplot(stat, aes(x = parameter, y = value * 100, fill = variable)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.7) +
  scale_fill_manual(values = c("#00A087FF", "#91D1C1FF")) +
  labs(title = NULL, x = NULL, y = "Coefficient of variation (%)", fill = NULL) +
  theme(
    panel.grid = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.position = c(0.9, 0.85)
  ) +
  scale_x_discrete(labels = c(
    expression(italic(r)[E]), expression(italic(p)[EP]), expression(italic(fp)[EP]),
    expression(italic(f)[D]), expression(italic(g)[D]), expression(italic(V)[g]),
    expression(italic(alpha)), expression(italic(K)[D]), expression(italic(Y)[g]),
    expression(italic(Q)[10])
  ))

p2 <- ggplot(stat, aes(x = variable, y = value * 100)) +
  geom_boxplot(outlier.size = 0.45, aes(color = variable)) +
  geom_jitter(aes(x = variable, y = value * 100, color = variable),
              position = position_jitterdodge(), size = 3, alpha = 1) +
  scale_color_manual(values = c("#00A087FF", "#91D1C1FF")) +
  labs(title = NULL, x = NULL, y = "Coefficient of variation (%)", fill = NULL) +
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(color = "black", fill = "transparent"),
    legend.position = "none",
    axis.line = element_line(colour = "black")
  )

###################
# 3. control
###################
file_C_tMEND <- file.path(model_dir, "Model_Control", "userio", "MEND3", "uncertainty", "tMEND", "H_UQpar.out")
file_C_gMEND <- file.path(model_dir, "Model_Control", "userio", "MEND3", "uncertainty", "gMEND", "H_UQpar.out")

if (!file.exists(file_C_tMEND)) stop("Missing file: ", file_C_tMEND)
if (!file.exists(file_C_gMEND)) stop("Missing file: ", file_C_gMEND)

UQ_tMEND_C <- read.delim(file_C_tMEND, sep = "", skip = 1)
UQ_gMEND_C <- read.delim(file_C_gMEND, sep = "", skip = 1)

UQ_tMEND_C <- UQ_tMEND_C[, c(18:26, 28)]
UQ_gMEND_C <- UQ_gMEND_C[, c(18:26, 28)]

stat_tMEND_C <- as.data.frame(apply(UQ_tMEND_C, 2, cal_cv))
stat_gMEND_C <- as.data.frame(apply(UQ_gMEND_C, 2, cal_cv))

colnames(stat_tMEND_C) <- "tMEND"
colnames(stat_gMEND_C) <- "gMEND"

stat_tMEND_C$parameter <- rownames(stat_tMEND_C)
stat_gMEND_C$parameter <- rownames(stat_gMEND_C)

mean(apply(UQ_tMEND_C, 2, cal_cv))
mean(apply(UQ_gMEND_C, 2, cal_cv))

mean(apply(UQ_tMEND_C, 2, cal_cv)) - mean(apply(UQ_gMEND_C, 2, cal_cv))

t.test(apply(UQ_tMEND_C, 2, cal_cv), apply(UQ_gMEND_C, 2, cal_cv), paired = TRUE)
t.test(apply(UQ_tMEND_C, 2, cal_cv), apply(UQ_gMEND_C, 2, cal_cv))
wilcox.test(apply(UQ_tMEND_C, 2, cal_cv), apply(UQ_gMEND_C, 2, cal_cv))

stat_C <- merge(stat_tMEND_C, stat_gMEND_C, by = "parameter")
stat_C <- melt(stat_C, id = "parameter")

stat_C$parameter <- factor(stat_C$parameter, levels = stat_tMEND_C$parameter)

p3 <- ggplot(stat_C, aes(x = parameter, y = value * 100, fill = variable)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.7) +
  scale_fill_manual(values = c("#00A087FF", "#91D1C1FF")) +
  labs(title = NULL, x = NULL, y = "Coefficient of variation (%)", fill = NULL) +
  theme(
    panel.grid = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.position = c(0.9, 0.85)
  ) +
  scale_x_discrete(labels = c(
    expression(italic(r)[E]), expression(italic(p)[EP]), expression(italic(fp)[EP]),
    expression(italic(f)[D]), expression(italic(g)[D]), expression(italic(V)[g]),
    expression(italic(alpha)), expression(italic(K)[D]), expression(italic(Y)[g]),
    expression(italic(Q)[10])
  ))

p4 <- ggplot(stat_C, aes(x = variable, y = value * 100)) +
  geom_boxplot(outlier.size = 0.45, aes(color = variable)) +
  geom_jitter(aes(x = variable, y = value * 100, color = variable),
              position = position_jitterdodge(), size = 3, alpha = 1) +
  scale_color_manual(values = c("#00A087FF", "#91D1C1FF")) +
  labs(title = NULL, x = NULL, y = "Coefficient of variation (%)", fill = NULL) +
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(color = "black", fill = "transparent"),
    legend.position = "none",
    axis.line = element_line(colour = "black")
  )

###################
# 4. half
###################
file_H_tMEND <- file.path(model_dir, "Model_half", "userio", "MEND3", "uncertainty", "tMEND", "H_UQpar.out")
file_H_gMEND <- file.path(model_dir, "Model_half", "userio", "MEND3", "uncertainty", "gMEND", "H_UQpar.out")

if (!file.exists(file_H_tMEND)) stop("Missing file: ", file_H_tMEND)
if (!file.exists(file_H_gMEND)) stop("Missing file: ", file_H_gMEND)

UQ_tMEND <- read.delim(file_H_tMEND, sep = "", skip = 1)
UQ_gMEND <- read.delim(file_H_gMEND, sep = "", skip = 1)

UQ_tMEND <- UQ_tMEND[, c(18:26, 28)]
UQ_gMEND <- UQ_gMEND[, c(18:26, 28)]

stat_tMEND <- as.data.frame(apply(UQ_tMEND, 2, cal_cv))
stat_gMEND <- as.data.frame(apply(UQ_gMEND, 2, cal_cv))

colnames(stat_tMEND) <- "tMEND"
colnames(stat_gMEND) <- "gMEND"

stat_tMEND$parameter <- rownames(stat_tMEND)
stat_gMEND$parameter <- rownames(stat_gMEND)

mean(apply(UQ_tMEND, 2, cal_cv))
mean(apply(UQ_gMEND, 2, cal_cv))

mean(apply(UQ_tMEND, 2, cal_cv)) - mean(apply(UQ_gMEND, 2, cal_cv))

t.test(apply(UQ_tMEND, 2, cal_cv), apply(UQ_gMEND, 2, cal_cv), paired = TRUE)
t.test(apply(UQ_tMEND, 2, cal_cv), apply(UQ_gMEND, 2, cal_cv))
wilcox.test(apply(UQ_tMEND, 2, cal_cv), apply(UQ_gMEND, 2, cal_cv))

stat <- merge(stat_tMEND, stat_gMEND, by = "parameter")
stat <- melt(stat, id = "parameter")

stat$parameter <- factor(stat$parameter, levels = stat_tMEND$parameter)

p5 <- ggplot(stat, aes(x = parameter, y = value * 100, fill = variable)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.7) +
  scale_fill_manual(values = c("#00A087FF", "#91D1C1FF")) +
  labs(title = NULL, x = NULL, y = "Coefficient of variation (%)", fill = NULL) +
  theme(
    panel.grid = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.position = c(0.9, 0.85)
  ) +
  scale_x_discrete(labels = c(
    expression(italic(r)[E]), expression(italic(p)[EP]), expression(italic(fp)[EP]),
    expression(italic(f)[D]), expression(italic(g)[D]), expression(italic(V)[g]),
    expression(italic(alpha)), expression(italic(K)[D]), expression(italic(Y)[g]),
    expression(italic(Q)[10])
  ))

p6 <- ggplot(stat, aes(x = variable, y = value * 100)) +
  geom_boxplot(outlier.size = 0.45, aes(color = variable)) +
  geom_jitter(aes(x = variable, y = value * 100, color = variable),
              position = position_jitterdodge(), size = 3, alpha = 1) +
  scale_color_manual(values = c("#00A087FF", "#91D1C1FF")) +
  labs(title = NULL, x = NULL, y = "Coefficient of variation (%)", fill = NULL) +
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(color = "black", fill = "transparent"),
    legend.position = "none",
    axis.line = element_line(colour = "black")
  )

###################
# 5. Save figures
###################
ggsave(
  plot = p2,
  filename = file.path(output_dir, "Fig4e_double_uncertainty.pdf"),
  width = 3.1, height = 2.6, units = "in", device = "pdf"
)

ggsave(
  plot = p2,
  filename = file.path(output_dir, "Fig4e_double_uncertainty.png"),
  width = 3.1, height = 2.6, units = "in", dpi = 300
)

ggsave(
  plot = p4,
  filename = file.path(output_dir, "Fig4e_control_uncertainty.pdf"),
  width = 3.1, height = 2.6, units = "in", device = "pdf"
)

ggsave(
  plot = p4,
  filename = file.path(output_dir, "Fig4e_control_uncertainty.png"),
  width = 3.1, height = 2.6, units = "in", dpi = 300
)

ggsave(
  plot = p6,
  filename = file.path(output_dir, "Fig4e_half_uncertainty.pdf"),
  width = 3.1, height = 2.6, units = "in", device = "pdf"
)

ggsave(
  plot = p6,
  filename = file.path(output_dir, "Fig4e_half_uncertainty.png"),
  width = 3.1, height = 2.6, units = "in", dpi = 300
)