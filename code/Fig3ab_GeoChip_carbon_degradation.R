############################################################
# Fig. 3a and Fig. 3b: GeoChip carbon degradation response ratio
# Project: litterbag
#
# Expected project structure:
# litterbag/
# ├── code/
# │   └── Fig3ab_GeoChip_response_ratio.R
# ├── data/
# │   ├── treat2.csv
# │   └── group_.cut.treat.0.25.log.RA.txt
# └── output/
#
# This script:
# 1) reads the transformed GeoChip relative abundance table,
# 2) keeps carbon degradation probes,
# 3) screens significant probes under precipitation treatments
#    in litterbag samples under unwarming conditions,
# 4) aggregates significant probes to gene level,
# 5) fits pairwise LMMs for half vs control and double vs control,
# 6) generates Fig. 3a and Fig. 3b.
#
# Input files used:
# - data/treat2.csv
# - data/group_.cut.treat.0.25.log.RA.txt
#
# Output files:
# - output/geochip/LMM_L_PHN_EN.pdf
# - output/geochip/LMM_L_PDN_EN.pdf
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
output_dir <- file.path(project_root, "output", "geochip")

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

treat_file <- file.path(data_dir, "treat2.csv")
geochip_file <- file.path(data_dir, "group_.cut.treat.0.25.log.RA.txt")

if (!file.exists(treat_file)) stop("Missing file: ", treat_file)
if (!file.exists(geochip_file)) stop("Missing file: ", geochip_file)

save.wd <- output_dir

suppressPackageStartupMessages({
  library(readr)
  library(ieggr)
  library(dplyr)
  library(ggplot2)
  library(ggsci)
  library(reshape2)
})

############################
# 1. read in data
############################
treat <- read.csv(treat_file, row.names = 1, check.names = FALSE)

treat <- subset(treat, treat$Clip == "uncliping")
treat <- treat[, -c(6)]
rownames(treat) <- treat$sample
treat_L <- subset(treat, type == "litter")
treat_S <- subset(treat, type == "soil")

GeoChip <- read_delim(geochip_file, show_col_types = FALSE)
GeoChip <- as.data.frame(GeoChip[, -1])
rownames(GeoChip) <- GeoChip$uniqueID

GeoChip_C <- subset(GeoChip, subcategory1 == "Carbon Degradation" | subcategory1 == "Carbon degradation")
unique(GeoChip_C$gene)
colnames(GeoChip_C)[11:ncol(GeoChip_C)] <- substring(colnames(GeoChip_C)[11:ncol(GeoChip_C)], 2, nchar(colnames(GeoChip_C)[11:ncol(GeoChip_C)]))

GeoChip_C_anno <- unique(GeoChip_C[, c(4, 8, 10)])

############################
# 2. P-litterbag
############################
# 2.1 all probes
treat_LP <- subset(treat_L, Warm == "unwarming")
GeoChip_C_L <- as.data.frame(t(GeoChip_C[, colnames(GeoChip_C) %in% treat_LP$sample]))
GeoChip_C_L <- GeoChip_C_L[, colSums(GeoChip_C_L, na.rm = TRUE) > 0]

sampc = match.name(rn.list = list(GeoChip_C_L = GeoChip_C_L, treat_LP = treat_LP))
GeoChip_C_L = sampc$GeoChip_C_L
treat_LP = sampc$treat_LP

GeoChip_C_L <- cbind(treat_LP, GeoChip_C_L)

GeoChip_C_L <- melt(GeoChip_C_L, id = c("sample", "year", "block", "Precipitation", "Warm", "type"))

unique(GeoChip_C_L$variable)
GeoChip_C_L$Precipitation <- factor(GeoChip_C_L$Precipitation, levels = c("normal", "half", "double"))

GeoChip_C_L[is.na(GeoChip_C_L) == TRUE] <- 0

LMM_P <- as.data.frame(matrix(nrow = length(unique(GeoChip_C_L$variable)), ncol = 2))
rownames(LMM_P) <- unique(GeoChip_C_L$variable)
colnames(LMM_P) <- c("F value", "Pr(>F)")

for (i in 1:length(unique(GeoChip_C_L$variable))) {
  subLMM <- subset(GeoChip_C_L, variable == unique(GeoChip_C_L$variable)[i])
  lm <- lmerTest::lmer(value ~ Precipitation + (1 | block) + (1 | year), data = subLMM)
  lm_summ <- anova(lm)
  LMM_P[i, "F value"] <- lm_summ$`F value`
  LMM_P[i, "Pr(>F)"] <- lm_summ$`Pr(>F)`
}

LMM_P_sig <- subset(LMM_P, `Pr(>F)` <= 0.1)

# 2.2 significant probes again
LMM_P_sig2 <- LMM_P_sig
colnames(LMM_P_sig2)
LMM_P_sig2$uniqueID <- rownames(LMM_P_sig2)

LMM_P_sig2 <- merge(LMM_P_sig2, GeoChip_C, by = "uniqueID")

LMM_P_sig2 <- aggregate(LMM_P_sig2[, 14:ncol(LMM_P_sig2)], by = list(LMM_P_sig2$gene), FUN = sum, na.rm = TRUE)
rownames(LMM_P_sig2) <- LMM_P_sig2$Group.1
LMM_P_sig2 <- LMM_P_sig2[, -1]

LMM_P_sig2 <- as.data.frame(t(LMM_P_sig2))

sampc = match.name(rn.list = list(LMM_P_sig2 = LMM_P_sig2, treat_LP = treat_LP))
LMM_P_sig2 = sampc$LMM_P_sig2
treat_LP = sampc$treat_LP

LMM_P_sig2 <- cbind(treat_LP, LMM_P_sig2)

LMM_P_sig2 <- melt(LMM_P_sig2, id = c("sample", "year", "block", "Precipitation", "Warm", "type"))

# 2.3 half/normal
LMM_P_sig2_HN <- subset(LMM_P_sig2, Precipitation != "double")
unique(LMM_P_sig2_HN$variable)
LMM_P_sig2_HN$Precipitation <- factor(LMM_P_sig2_HN$Precipitation, levels = c("normal", "half"))

LMM_P_2 <- as.data.frame(matrix(nrow = length(unique(LMM_P_sig2_HN$variable)), ncol = 3))
rownames(LMM_P_2) <- unique(LMM_P_sig2_HN$variable)
colnames(LMM_P_2) <- c("estimate", "t value", "Pr(>|t|)")

for (i in 1:length(unique(LMM_P_sig2_HN$variable))) {
  subLMM <- subset(LMM_P_sig2_HN, variable == unique(LMM_P_sig2_HN$variable)[i])
  lm <- lmerTest::lmer(value ~ Precipitation + (1 | block) + (1 | year), data = subLMM)
  lm_summ <- summary(lm)
  LMM_P_2[i, "estimate"] <- lm_summ$coefficients[2, 1]
  LMM_P_2[i, "t value"] <- lm_summ$coefficients[2, 4]
  LMM_P_2[i, "Pr(>|t|)"] <- lm_summ$coefficients[2, 5]
}

# 2.4 double/normal
LMM_P_sig2_DN <- subset(LMM_P_sig2, Precipitation != "half")
unique(LMM_P_sig2_DN$variable)
LMM_P_sig2_DN$Precipitation <- factor(LMM_P_sig2_DN$Precipitation, levels = c("normal", "double"))

LMM_P_3 <- as.data.frame(matrix(nrow = length(unique(LMM_P_sig2_DN$variable)), ncol = 3))
rownames(LMM_P_3) <- unique(LMM_P_sig2_DN$variable)
colnames(LMM_P_3) <- c("estimate", "t value", "Pr(>|t|)")

for (i in 1:length(unique(LMM_P_sig2_DN$variable))) {
  subLMM <- subset(LMM_P_sig2_DN, variable == unique(LMM_P_sig2_DN$variable)[i])
  lm <- lmerTest::lmer(value ~ Precipitation + (1 | block) + (1 | year), data = subLMM)
  lm_summ <- summary(lm)
  LMM_P_3[i, "estimate"] <- lm_summ$coefficients[2, 1]
  LMM_P_3[i, "t value"] <- lm_summ$coefficients[2, 4]
  LMM_P_3[i, "Pr(>|t|)"] <- lm_summ$coefficients[2, 5]
}

############################
# 3. picture
############################

LMM_P_4 <- LMM_P_2
LMM_P_4[LMM_P_4$estimate > 0, "direction"] <- "+"
LMM_P_4[LMM_P_4$estimate < 0, "direction"] <- "-"

LMM_P_4$gene <- rownames(LMM_P_4)
LMM_P_4 <- merge(LMM_P_4, GeoChip_C_anno, by = "gene")

unique(LMM_P_4$subcategory2)
LMM_P_4 <- subset(LMM_P_4, LMM_P_4$subcategory2 != "Phospholipids")
LMM_P_4$subcategory2 <- factor(
  LMM_P_4$subcategory2,
  levels = c("Glyoxylate cycle", "Lactose", "Inulin",
             "Starch", "Hemicellulose", "Pectin",
             "Cellulose", "Cutin", "Chitin",
             "Vanillin/Lignin", "Lignin", "Tannins",
             "Terpenes", "Other")
)
LMM_P_4 <- LMM_P_4[order(LMM_P_4$subcategory2), ]
LMM_P_4$gene <- factor(LMM_P_4$gene, levels = c(LMM_P_4$gene))

LMM_P_4$sig = ""
LMM_P_4[LMM_P_4$`Pr(>|t|)` < 0.1, "sig"] = "+"
LMM_P_4[LMM_P_4$`Pr(>|t|)` < 0.05, "sig"] = "*"
LMM_P_4[LMM_P_4$`Pr(>|t|)` < 0.01, "sig"] = "**"
LMM_P_4[LMM_P_4$`Pr(>|t|)` < 0.001, "sig"] = "***"

LMM_P_4$stary = ""
for (i in 1:nrow(LMM_P_4)) {
  if (LMM_P_4$estimate[i] > 0) {
    LMM_P_4$stary[i] = LMM_P_4$estimate[i] + 1.5
  } else {
    LMM_P_4$stary[i] = 1.5
  }
}
LMM_P_4$stary <- as.numeric(LMM_P_4$stary)

FDR = 0
p1 <- ggplot(data = LMM_P_4, aes(x = gene, y = estimate, fill = subcategory2)) +
  geom_bar(aes(x = gene, y = estimate, fill = subcategory2), alpha = 0.8, stat = "identity") +
  geom_hline(yintercept = FDR, linewidth = 0.5, color = "black") +
  labs(x = "", y = "effect") +
  scale_fill_d3(palette = "category20") +
  geom_text(aes(x = gene, y = stary, label = sig), size = 4, vjust = 0) +
  theme_bw() +
  ggtitle("") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "transparent", color = "black"),
    text = element_text(size = 14),
    axis.text = element_text(colour = "black"),
    axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1),
    legend.position = "none"
  ) +
  geom_vline(xintercept = 4.5, lwd = 0.5, colour = "black", linetype = "dashed") +
  geom_vline(xintercept = 5.5, lwd = 0.5, colour = "black", linetype = "dashed") +
  geom_vline(xintercept = 6.5, lwd = 0.5, colour = "black", linetype = "dashed") +
  geom_vline(xintercept = 12.5, lwd = 0.5, colour = "black", linetype = "dashed") +
  geom_vline(xintercept = 15.5, lwd = 0.5, colour = "black", linetype = "dashed") +
  geom_vline(xintercept = 25.5, lwd = 0.5, colour = "black", linetype = "dashed") +
  geom_vline(xintercept = 29.5, lwd = 0.5, colour = "black", linetype = "dashed") +
  geom_vline(xintercept = 30.5, lwd = 0.5, colour = "black", linetype = "dashed") +
  geom_vline(xintercept = 33.5, lwd = 0.5, colour = "black", linetype = "dashed") +
  geom_vline(xintercept = 35.5, lwd = 0.5, colour = "black", linetype = "dashed") +
  geom_vline(xintercept = 39.5, lwd = 0.5, colour = "black", linetype = "dashed") +
  geom_vline(xintercept = 40.5, lwd = 0.5, colour = "black", linetype = "dashed") +
  geom_vline(xintercept = 43.5, lwd = 0.5, colour = "black", linetype = "dashed") +
  annotate("text", x = 2.5, y = -55, label = "Glyoxylate cycle", colour = "black", size = 3.5) +
  annotate("text", x = 5, y = -55, label = "Lactose", colour = "black", size = 3.5, angle = 90) +
  annotate("text", x = 6, y = -55, label = "Inulin", colour = "black", size = 3.5, angle = 90) +
  annotate("text", x = 9.5, y = -55, label = "Starch", colour = "black", size = 3.5) +
  annotate("text", x = 14, y = -55, label = "Hemicellulose", colour = "black", size = 3.5) +
  annotate("text", x = 20, y = -55, label = "Pectin", colour = "black", size = 3.5) +
  annotate("text", x = 27.5, y = -55, label = "Cellulose", colour = "black", size = 3.5) +
  annotate("text", x = 30, y = -55, label = "Cutin", colour = "black", size = 3.5, angle = 90) +
  annotate("text", x = 32, y = -55, label = "Chitin", colour = "black", size = 3.5) +
  annotate("text", x = 34.5, y = -55, label = "Vanillin/Lignin", colour = "black", size = 3.5, angle = 90) +
  annotate("text", x = 37.5, y = -55, label = "Lignin", colour = "black", size = 3.5) +
  annotate("text", x = 40, y = -55, label = "Tannins", colour = "black", size = 3.5, angle = 90) +
  annotate("text", x = 42, y = -55, label = "Terpenes", colour = "black", size = 3.5) +
  annotate("text", x = 44, y = -55, label = "Other", colour = "black", size = 3.5, angle = 90)
p1

LMM_P_5 <- LMM_P_3
LMM_P_5[LMM_P_5$estimate > 0, "direction"] <- "+"
LMM_P_5[LMM_P_5$estimate < 0, "direction"] <- "-"

LMM_P_5$gene <- rownames(LMM_P_5)
LMM_P_5 <- merge(LMM_P_5, GeoChip_C_anno, by = "gene")

unique(LMM_P_5$subcategory2)
LMM_P_5 <- subset(LMM_P_5, LMM_P_5$subcategory2 != "Phospholipids")
LMM_P_5$subcategory2 <- factor(
  LMM_P_5$subcategory2,
  levels = c("Glyoxylate cycle", "Lactose", "Inulin",
             "Starch", "Hemicellulose", "Pectin",
             "Cellulose", "Cutin", "Chitin",
             "Vanillin/Lignin", "Lignin", "Tannins",
             "Terpenes", "Other")
)
LMM_P_5 <- LMM_P_5[order(LMM_P_5$subcategory2), ]
LMM_P_5$gene <- factor(LMM_P_5$gene, levels = c(LMM_P_5$gene))

LMM_P_5$sig = ""
LMM_P_5[LMM_P_5$`Pr(>|t|)` < 0.1, "sig"] = "+"
LMM_P_5[LMM_P_5$`Pr(>|t|)` < 0.05, "sig"] = "*"
LMM_P_5[LMM_P_5$`Pr(>|t|)` < 0.01, "sig"] = "**"
LMM_P_5[LMM_P_5$`Pr(>|t|)` < 0.001, "sig"] = "***"

LMM_P_5$stary = ""
for (i in 1:nrow(LMM_P_5)) {
  if (LMM_P_5$estimate[i] > 0) {
    LMM_P_5$stary[i] = LMM_P_5$estimate[i] + 1.5
  } else {
    LMM_P_5$stary[i] = 1.5
  }
}
LMM_P_5$stary <- as.numeric(LMM_P_5$stary)

FDR = 0
p2 <- ggplot(data = LMM_P_5, aes(x = gene, y = estimate, fill = subcategory2)) +
  geom_bar(aes(x = gene, y = estimate, fill = subcategory2), alpha = 0.8, stat = "identity") +
  geom_hline(yintercept = FDR, linewidth = 0.5, color = "black") +
  labs(x = "", y = "effect") +
  scale_fill_d3(palette = "category20") +
  geom_text(aes(x = gene, y = stary, label = sig), size = 4, vjust = 0) +
  theme_bw() +
  ggtitle("") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "transparent", color = "black"),
    text = element_text(size = 14),
    axis.text = element_text(colour = "black"),
    axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1),
    legend.position = "none"
  ) +
  geom_vline(xintercept = 4.5, lwd = 0.5, colour = "black", linetype = "dashed") +
  geom_vline(xintercept = 5.5, lwd = 0.5, colour = "black", linetype = "dashed") +
  geom_vline(xintercept = 6.5, lwd = 0.5, colour = "black", linetype = "dashed") +
  geom_vline(xintercept = 12.5, lwd = 0.5, colour = "black", linetype = "dashed") +
  geom_vline(xintercept = 15.5, lwd = 0.5, colour = "black", linetype = "dashed") +
  geom_vline(xintercept = 25.5, lwd = 0.5, colour = "black", linetype = "dashed") +
  geom_vline(xintercept = 29.5, lwd = 0.5, colour = "black", linetype = "dashed") +
  geom_vline(xintercept = 30.5, lwd = 0.5, colour = "black", linetype = "dashed") +
  geom_vline(xintercept = 33.5, lwd = 0.5, colour = "black", linetype = "dashed") +
  geom_vline(xintercept = 35.5, lwd = 0.5, colour = "black", linetype = "dashed") +
  geom_vline(xintercept = 39.5, lwd = 0.5, colour = "black", linetype = "dashed") +
  geom_vline(xintercept = 40.5, lwd = 0.5, colour = "black", linetype = "dashed") +
  geom_vline(xintercept = 43.5, lwd = 0.5, colour = "black", linetype = "dashed") +
  annotate("text", x = 2.5, y = 60, label = "Glyoxylate cycle", colour = "black", size = 3.5) +
  annotate("text", x = 5, y = 60, label = "Lactose", colour = "black", size = 3.5, angle = 90) +
  annotate("text", x = 6, y = 60, label = "Inulin", colour = "black", size = 3.5, angle = 90) +
  annotate("text", x = 9.5, y = 60, label = "Starch", colour = "black", size = 3.5) +
  annotate("text", x = 14, y = 60, label = "Hemicellulose", colour = "black", size = 3.5) +
  annotate("text", x = 20, y = 60, label = "Pectin", colour = "black", size = 3.5) +
  annotate("text", x = 27.5, y = 60, label = "Cellulose", colour = "black", size = 3.5) +
  annotate("text", x = 30, y = 60, label = "Cutin", colour = "black", size = 3.5, angle = 90) +
  annotate("text", x = 32, y = 60, label = "Chitin", colour = "black", size = 3.5) +
  annotate("text", x = 34.5, y = 60, label = "Vanillin/Lignin", colour = "black", size = 3.5, angle = 90) +
  annotate("text", x = 37.5, y = 60, label = "Lignin", colour = "black", size = 3.5) +
  annotate("text", x = 40, y = 60, label = "Tannins", colour = "black", size = 3.5, angle = 90) +
  annotate("text", x = 42, y = 60, label = "Terpenes", colour = "black", size = 3.5) +
  annotate("text", x = 44, y = 60, label = "Other", colour = "black", size = 3.5, angle = 90)
p2

library(showtext)
showtext_auto()

ggsave(
  plot = p1,
  dpi = 300,
  filename = file.path(save.wd, "LMM_L_PHN_EN.pdf"),
  width = 11,
  height = 5.5,
  units = "in"
)

ggsave(
  plot = p2,
  dpi = 300,
  filename = file.path(save.wd, "LMM_L_PDN_EN.pdf"),
  width = 10.9,
  height = 5.5,
  units = "in"
)