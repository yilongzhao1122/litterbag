############################################################
# Fig. 1b and Fig. 1c: Alpha diversity boxplots
# Project: litterbag
#
# Expected project structure:
# litterbag/
# ├── code/
# │   └── Fig1bc_alpha_diversity_boxplot.R
# ├── data/
# │   ├── 16S_30236.txt
# │   ├── ITS_10391.txt
# │   ├── 16S_tree.nwk
# │   ├── treat2.csv
# └── output/
#     └── alpha_diversity/
#
# This script:
# 1) reads bacterial and fungal amplicon tables, phylogenetic tree,
#    and metadata,
# 2) calculates alpha diversity metrics,
# 3) keeps only unwarming samples for the precipitation comparison,
# 4) draws the bacterial and fungal alpha diversity boxplots,
# 5) saves Fig. 1b and Fig. 1c.
############################################################

############################################################
# 0. Locate project root and define paths
############################################################
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
output_dir <- file.path(project_root, "output", "alpha_diversity")

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

otu16_file <- file.path(data_dir, "16S_30236.txt")
its_file <- file.path(data_dir, "ITS_10391.txt")
tree_file <- file.path(data_dir, "16S_tree.nwk")

# Use treat2.csv if present
treat_file <- file.path(data_dir, "treat2.csv")


required_files <- c(otu16_file, its_file, tree_file, treat_file)
missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
  stop(
    "Missing required input file(s):\n",
    paste(missing_files, collapse = "\n")
  )
}

############################################################
# 1. Load packages
############################################################
suppressPackageStartupMessages({
  library(ieggr)
  library(vegan)
  library(picante)
  library(ape)
  library(ggplot2)
  library(reshape2)
  library(dplyr)
  library(showtext)
  library(grid)
})

showtext_auto()

############################################################
# 2. Read input data
############################################################
otu.file <- read.table(otu16_file, header = TRUE, row.names = 1, check.names = FALSE)
ITS.file <- read.table(its_file, header = TRUE, row.names = 1, check.names = FALSE)
tree <- ape::read.tree(tree_file, quote = "\"'")

treat <- read.csv(treat_file, row.names = 1, check.names = FALSE)
if (!("sample" %in% colnames(treat))) {
  treat$sample <- rownames(treat)
}

treat <- subset(treat, Clip == "uncliping")

# Keep the original metadata handling logic
if (ncol(treat) >= 6) {
  treat <- treat[, -6]
}

otu <- t(otu.file)
ITS <- t(ITS.file)

otu <- subset(otu, rownames(otu) %in% treat$sample)
ITS <- subset(ITS, rownames(ITS) %in% treat$sample)

otu <- otu[, colSums(otu) != 0, drop = FALSE]
ITS <- ITS[, colSums(ITS) != 0, drop = FALSE]

############################################################
# 3. Calculate alpha diversity
############################################################

#################################
# 3.1 Bacterial alpha diversity
#################################
# Taxonomic alpha diversity
alpha.tax <- alpha.g(otu, trace = TRUE)
alpha.tax <- as.data.frame(alpha.tax)
alpha.tax$sample <- rownames(alpha.tax)

# Phylogenetic alpha diversity
Faith.PD <- picante::pd(otu, tree, include.root = FALSE)
Faith.PD$sample <- rownames(Faith.PD)

alpha_16S <- merge(Faith.PD[, c(1, 3)], alpha.tax, by = "sample")
alpha_16S$shannon_diversity <- exp(1)^alpha_16S$shannon

#################################
# 3.2 Fungal alpha diversity
#################################
alpha_ITS <- alpha.g(ITS, trace = TRUE)
alpha_ITS <- as.data.frame(alpha_ITS)
alpha_ITS$shannon_diversity <- exp(1)^alpha_ITS$shannon
alpha_ITS$sample <- rownames(alpha_ITS)

############################################################
# 4. Keep only unwarming samples for precipitation plots
############################################################
alpha_16S <- merge(treat, alpha_16S, by = "sample")
alpha_16S <- subset(alpha_16S, Warm == "unwarming")

alpha_ITS <- merge(treat, alpha_ITS, by = "sample")
alpha_ITS <- subset(alpha_ITS, Warm == "unwarming")

############################################################
# 5. Reshape data for plotting
############################################################
alpha_16S <- melt(alpha_16S, id = c("sample", "year", "block", "Warm", "Precipitation", "type"))
alpha_ITS <- melt(alpha_ITS, id = c("sample", "year", "block", "Warm", "Precipitation", "type"))

############################################################
# 6. Plot styling
############################################################

# Match the original display labels
alpha_16S[alpha_16S$type == "litter", "type"] <- "litterbag"
alpha_16S$type <- as.factor(alpha_16S$type)

alpha_16S$Precipitation <- as.character(alpha_16S$Precipitation)
alpha_16S[alpha_16S$Precipitation == "half", "Precipitation"] <- "Half"
alpha_16S[alpha_16S$Precipitation == "normal", "Precipitation"] <- "Control"
alpha_16S[alpha_16S$Precipitation == "double", "Precipitation"] <- "Double"
alpha_16S$Precipitation <- factor(alpha_16S$Precipitation, levels = c("Half", "Control", "Double"))

alpha_ITS[alpha_ITS$type == "litter", "type"] <- "litterbag"
alpha_ITS$type <- as.factor(alpha_ITS$type)

alpha_ITS$Precipitation <- as.character(alpha_ITS$Precipitation)
alpha_ITS[alpha_ITS$Precipitation == "half", "Precipitation"] <- "Half"
alpha_ITS[alpha_ITS$Precipitation == "normal", "Precipitation"] <- "Control"
alpha_ITS[alpha_ITS$Precipitation == "double", "Precipitation"] <- "Double"
alpha_ITS$Precipitation <- factor(alpha_ITS$Precipitation, levels = c("Half", "Control", "Double"))

############################################################
# 7. Fig. 1b: Bacterial alpha diversity boxplot
############################################################
alpha_16S1 <- subset(alpha_16S, variable %in% c("richness", "shannon"))
alpha_16S1$variable <- factor(alpha_16S1$variable, levels = c("richness", "shannon"))

p1 <- ggplot(alpha_16S1, aes(x = type, y = value)) +
  facet_wrap(~variable, scales = "free", nrow = 1) +
  geom_boxplot(outlier.size = 0.45, aes(color = Precipitation)) +
  geom_jitter(
    aes(x = type, y = value, color = Precipitation),
    position = position_jitterdodge(),
    size = 0.45,
    alpha = 1
  ) +
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "transparent", color = "black"),
    legend.title = element_blank(),
    legend.key = element_blank(),
    legend.position = "right",
    strip.background = element_rect(color = "transparent", fill = "transparent")
  ) +
  theme(
    strip.switch.pad.grid = unit(2, "inch"),
    axis.text.x = element_text(color = "black"),
    axis.text.y = element_text(color = "black")
  ) +
  labs(x = "", y = "") +
  scale_color_manual(values = c("#F39B7FFF", "#3C5488FF", "#4DBBD5FF"))

############################################################
# 8. Fig. 1c: Fungal alpha diversity boxplot
############################################################
alpha_ITS1 <- subset(alpha_ITS, variable %in% c("richness", "shannon"))
alpha_ITS1$variable <- factor(alpha_ITS1$variable, levels = c("richness", "shannon"))

p2 <- ggplot(alpha_ITS1, aes(x = type, y = value)) +
  facet_wrap(~variable, scales = "free", nrow = 1) +
  geom_boxplot(outlier.size = 0.45, aes(color = Precipitation)) +
  geom_jitter(
    aes(x = type, y = value, color = Precipitation),
    position = position_jitterdodge(),
    size = 0.45,
    alpha = 1
  ) +
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "transparent", color = "black"),
    legend.title = element_blank(),
    legend.key = element_blank(),
    legend.position = "right",
    strip.background = element_rect(color = "transparent", fill = "transparent")
  ) +
  theme(
    strip.switch.pad.grid = unit(2, "inch"),
    axis.text.x = element_text(color = "black"),
    axis.text.y = element_text(color = "black")
  ) +
  labs(x = "", y = "") +
  scale_color_manual(values = c("#F39B7FFF", "#3C5488FF", "#4DBBD5FF"))

############################################################
# 9. Save figures
############################################################
ggsave(
  plot = p1,
  filename = file.path(output_dir, "Fig1b_alpha_diversity_16S.pdf"),
  width = 5.4,
  height = 2.4,
  units = "in",
  dpi = 300
)

ggsave(
  plot = p1,
  filename = file.path(output_dir, "Fig1b_alpha_diversity_16S.png"),
  width = 5.4,
  height = 2.4,
  units = "in",
  dpi = 300
)

ggsave(
  plot = p2,
  filename = file.path(output_dir, "Fig1c_alpha_diversity_ITS.pdf"),
  width = 5.4,
  height = 2.4,
  units = "in",
  dpi = 300
)

ggsave(
  plot = p2,
  filename = file.path(output_dir, "Fig1c_alpha_diversity_ITS.png"),
  width = 5.4,
  height = 2.4,
  units = "in",
  dpi = 300
)