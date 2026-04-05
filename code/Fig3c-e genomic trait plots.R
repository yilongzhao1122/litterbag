############################################################
# Fig. 3c-e: Genomic trait plots under precipitation treatments
# Project: litterbag
#
# Expected project structure:
# litterbag/
# ├── code/
# │   └── Fig3cde_genomic_traits.R
# ├── data/
# │   ├── treat2.csv
# │   ├── 16S_30236.txt
# │   └── genometrait.txt
# └── output/
#
# This script:
# 1) reads treatment metadata, bacterial OTU table, and genomic trait table,
# 2) calculates community-level copy number, abundance-weighted genome size,
#    and abundance-weighted GC content,
# 3) fits precipitation LMMs under unwarming conditions,
# 4) generates the precipitation genomic trait plots used for Fig. 3c-e.
#
# Input files used:
# - data/treat2.csv
# - data/16S_30236.txt
# - data/genometrait.txt
#
# Output files:
# - output/genomic_traits/Fig3cde_genomic_traits_precipitation.pdf
# - output/genomic_traits/Fig3cde_genomic_traits_precipitation.png
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
output_dir <- file.path(project_root, "output", "genomic_traits")

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

treat_file <- file.path(data_dir, "treat2.csv")
otu_file <- file.path(data_dir, "16S_30236.txt")
trait_file <- file.path(data_dir, "genometrait.txt")

required_files <- c(treat_file, otu_file, trait_file)
missing_files <- required_files[!file.exists(required_files)]

if (length(missing_files) > 0) {
  stop(
    "Missing required input file(s):\n",
    paste(missing_files, collapse = "\n")
  )
}

suppressPackageStartupMessages({
  library(readr)
  library(ieggr)
  library(dplyr)
  library(reshape2)
  library(ggplot2)
  library(lmerTest)
  library(emmeans)
})

################################
# 1. read in data
################################
treat <- read.csv(treat_file, row.names = 1, check.names = FALSE)
treat <- subset(treat, treat$Clip == "uncliping")
treat <- treat[, -c(6)]
rownames(treat) <- treat$sample

otu.file <- read.table(otu_file, header = TRUE, row.names = 1, check.names = FALSE)
otu <- as.data.frame(t(otu.file))
otu <- subset(otu, rownames(otu) %in% treat$sample)
otu <- otu[, colSums(otu) != 0]

rrn <- read.csv(trait_file, check.names = FALSE, stringsAsFactors = FALSE)

# Remove the first index column if present
if (colnames(rrn)[1] == "" || grepl("^X$", colnames(rrn)[1]) || grepl("^\\.\\.\\.", colnames(rrn)[1])) {
  rrn <- rrn[, -1, drop = FALSE]
}

required_cols <- c("OTU.ID", "hit", "size", "GC")
if (!all(required_cols %in% colnames(rrn))) {
  stop("genometrait.txt must contain the columns: OTU.ID, hit, size, GC")
}

if (any(duplicated(colnames(rrn)))) {
  stop("genometrait.txt contains duplicated column names: ",
       paste(unique(colnames(rrn)[duplicated(colnames(rrn))]), collapse = ", "))
}

samp <- match.name(rn.list = list(otu = otu, treat = treat))
otu <- samp$otu
treat <- samp$treat

################################################
# 2. calculate community-level rrn copy number, genome size and GC
################################################
otu <- as.data.frame(t(otu))
otu$OTU.ID <- rownames(otu)
otu <- merge(otu, rrn, by = "OTU.ID")

sample_cols <- colnames(otu)[2:(ncol(otu) - 3)]
n_samples <- length(sample_cols)

otu_rrn <- data.frame(sample = sample_cols, stringsAsFactors = FALSE)

for (i in seq_len(n_samples)) {
  otu_rrn[i, "CopyNumber"] <- sum(otu[, i + 1]) / sum(otu[, i + 1] / otu$hit)
}

for (i in seq_len(n_samples)) {
  otu_rrn[i, "GenomeSize"] <- sum(otu[, i + 1] * otu$size) / sum(otu[, i + 1])
}

for (i in seq_len(n_samples)) {
  otu_rrn[i, "GC"] <- sum(otu[, i + 1] * otu$GC) / sum(otu[, i + 1])
}

################################################
# 3. mean copy number, genome size and GC
################################################
for (i in seq_len(n_samples)) {
  otui <- otu[, c("OTU.ID", sample_cols[i], "hit", "size", "GC")]
  otui <- subset(otui, otui[, 2] > 0)
  otu_rrn[i, "CopyNumber_mean"] <- mean(otui$hit)
  otu_rrn[i, "GenomeSize_mean"] <- mean(otui$size)
  otu_rrn[i, "GC_mean"] <- mean(otui$GC)
}

otu_rrn <- merge(treat, otu_rrn, by = "sample")
otu_rrn$Precipitation <- factor(otu_rrn$Precipitation, levels = c("normal", "half", "double"))

otu_rrn %>% group_by(type) %>% summarise(mean(CopyNumber), .groups = "drop")

################################################
# 4. precipitation models
################################################
rrn_P <- subset(otu_rrn, Warm == "unwarming")

anova(lmer(CopyNumber ~ Precipitation * type + (1 | block) + (1 | year), data = rrn_P))
summary(lmer(CopyNumber ~ Precipitation + type + (1 | block) + (1 | year), data = rrn_P))
emmeans(
  lmer(CopyNumber ~ Precipitation * type + (1 | block) + (1 | year), data = rrn_P),
  pairwise ~ Precipitation | type
) %>% .$contr

anova(lmer(GenomeSize ~ Precipitation * type + (1 | block) + (1 | year), data = rrn_P))
summary(lmer(GenomeSize ~ Precipitation + type + (1 | block) + (1 | year), data = rrn_P))
emmeans(
  lmer(GenomeSize ~ Precipitation * type + (1 | block) + (1 | year), data = rrn_P),
  pairwise ~ Precipitation | type
) %>% .$contr

anova(lmer(GC ~ Precipitation * type + (1 | block) + (1 | year), data = rrn_P))
summary(lmer(GC ~ Precipitation + type + (1 | block) + (1 | year), data = rrn_P))

anova(lmer(CopyNumber_mean ~ Precipitation * type + (1 | block) + (1 | year), data = rrn_P))
summary(lmer(CopyNumber_mean ~ Precipitation + type + (1 | block) + (1 | year), data = rrn_P))
emmeans(
  lmer(CopyNumber_mean ~ Precipitation * type + (1 | block) + (1 | year), data = rrn_P),
  pairwise ~ Precipitation | type
) %>% .$contr

anova(lmer(GenomeSize_mean ~ Precipitation * type + (1 | block) + (1 | year), data = rrn_P))
summary(lmer(GenomeSize_mean ~ Precipitation + type + (1 | block) + (1 | year), data = rrn_P))

anova(lmer(GC_mean ~ Precipitation * type + (1 | block) + (1 | year), data = rrn_P))
summary(lmer(GC_mean ~ Precipitation + type + (1 | block) + (1 | year), data = rrn_P))

################################################
# 5. picture for Fig. 3c-e
################################################
rrn_P2 <- rrn_P
rrn_P2 <- melt(rrn_P2, id = c("sample", "year", "block", "Precipitation", "Warm", "type"))
rrn_P2$variable <- as.character(rrn_P2$variable)

rrn_P2[rrn_P2$variable == "GenomeSize", "variable"] <- "GenomeSize (Mbp)"
rrn_P2[rrn_P2$variable == "GC", "variable"] <- "GC (%)"

rrn_P2 <- subset(rrn_P2, variable %in% c("CopyNumber", "GenomeSize (Mbp)", "GC (%)"))
rrn_P2$variable <- factor(rrn_P2$variable, levels = c("CopyNumber", "GenomeSize (Mbp)", "GC (%)"))
rrn_P2$Precipitation <- factor(rrn_P2$Precipitation, levels = c("half", "normal", "double"))

p2 <- ggplot(rrn_P2, aes(x = type, y = value)) +
  facet_wrap(~variable, scales = "free", nrow = 1) +
  geom_boxplot(outlier.size = 0.45, aes(color = Precipitation), width = 0.8) +
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
    legend.position = "none",
    strip.background = element_rect(color = "transparent", fill = "transparent")
  ) +
  theme(
    strip.text = element_text(color = "black"),
    axis.text.x = element_text(color = "black"),
    axis.text.y = element_text(color = "black")
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.25))) +
  labs(x = "", y = "", title = "") +
  scale_color_manual(values = c("#F39B7FFF", "#3C5488FF", "#4DBBD5FF"))
p2

################################################
# 6. save figure
################################################
ggsave(
  plot = p2,
  filename = file.path(output_dir, "Fig3cde_genomic_traits_precipitation.pdf"),
  width = 6,
  height = 2.5,
  units = "in",
  dpi = 300
)

ggsave(
  plot = p2,
  filename = file.path(output_dir, "Fig3cde_genomic_traits_precipitation.png"),
  width = 6,
  height = 2.5,
  units = "in",
  dpi = 300
)