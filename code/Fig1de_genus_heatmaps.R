############################################################
# Fig. 1d and Fig. 1e: Basic genus-level heatmaps for bacteria and fungi
# Project: litterbag
#
# Expected project structure:
# litterbag/
# ├── code/
# │   └── Fig1de_heatmap_basic.R
# ├── data/
# │   ├── 16S_30236.txt
# │   ├── ITS_10391.txt
# │   ├── taxonomy_16S2.txt
# │   ├── taxonomy_ITS2.txt
# │   └── treat2.csv
# └── output/
#     └── heatmap/
#
# This script:
# 1) reads bacterial/fungal amplicon tables, taxonomy tables, and metadata,
# 2) aggregates abundance to the genus level,
# 3) calculates mean relative abundance under half, control, and double precipitation,
# 4) log-transforms and row-standardizes the abundance matrix,
# 5) draws basic heatmaps without phylum annotation or significance labels,
# 6) saves Fig. 1d and Fig. 1e.
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
output_dir <- file.path(project_root, "output", "heatmap")

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

otu16_file <- file.path(data_dir, "16S_30236.txt")
its_file <- file.path(data_dir, "ITS_10391.txt")
tax16_file <- file.path(data_dir, "taxonomy_16S2.txt")
taxits_file <- file.path(data_dir, "taxonomy_ITS2.txt")
treat_file <- file.path(data_dir, "treat2.csv")

required_files <- c(otu16_file, its_file, tax16_file, taxits_file, treat_file)
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
  library(RColorBrewer)
  library(ComplexHeatmap)
  library(showtext)
  library(grid)
})

showtext_auto()

############################################################
# 2. Helper function
############################################################
make_basic_heatmap_matrix <- function(mat) {
  mat2 <- mat
  for (i in 1:nrow(mat2)) {
    mat2[i, ] <- scale(log(unlist(mat2[i, ]) + 1, 2))
  }
  as.matrix(mat2)
}

############################################################
# 3. Read data and preprocess
############################################################
otu.file <- read.table(otu16_file, header = TRUE, row.names = 1, check.names = FALSE)
ITS.file <- read.table(its_file, header = TRUE, row.names = 1, check.names = FALSE)

# Keep the taxonomy import style identical to the original workflow
taxonomy_16S <- read.table(tax16_file)
taxonomy_ITS <- read.table(taxits_file)

treat <- read.csv(treat_file, row.names = 1, check.names = FALSE)

treat <- subset(treat, treat$Clip == "uncliping" & treat$Warm == "unwarming")
treat <- treat[, -c(5, 6)]

# Soil samples only, following the original Fig. 1d/e workflow
treat <- subset(treat, type == "soil")
rownames(treat) <- treat$sample

otu <- as.data.frame(t(otu.file))
ITS <- as.data.frame(t(ITS.file))

otu <- subset(otu, rownames(otu) %in% treat$sample)
ITS <- subset(ITS, rownames(ITS) %in% treat$sample)

otu <- otu[, colSums(otu) != 0, drop = FALSE]
ITS <- ITS[, colSums(ITS) != 0, drop = FALSE]

############################################################
# 4. Bacterial heatmap
############################################################

##########################################
# 4.1 Genus-level bacterial abundance table
##########################################
rownames(taxonomy_16S) <- taxonomy_16S$Feature.ID
taxonomy_16S <- taxonomy_16S[, -1]

otu2 <- as.data.frame(t(otu))
samp <- match.name(rn.list = list(otu2 = otu2, taxonomy_16S = taxonomy_16S))
otu2 <- samp$otu2
taxonomy_16S <- samp$taxonomy_16S
otu2 <- cbind(taxonomy_16S, otu2)

otu2 <- aggregate(otu2[, 9:ncol(otu2)], by = list(otu2$Genus), FUN = sum)
rownames(otu2) <- otu2$Group.1
otu2 <- otu2[, -1]

##########################################
# 4.2 Mean abundance by precipitation treatment
##########################################
otu3 <- as.data.frame(t(otu2))
spc <- match.name(rn.list = list(treat = treat, otu3 = otu3))
otu3 <- cbind(spc$treat, spc$otu3)

otu3 <- aggregate(otu3[, 6:ncol(otu3)], by = list(otu3$Precipitation), FUN = mean)
rownames(otu3) <- otu3$Group.1
otu3 <- as.data.frame(t(otu3[, -1]))

otu3.ra <- decostand(otu3, "total", MARGIN = 2)
otu3.ra <- otu3.ra[rowMeans(otu3.ra) > 0.001, ]
otu3.ra <- otu3.ra[, c(2, 3, 1)]

# Remove ambiguous taxa
otu3.ra <- otu3.ra[
  !grepl("uncultured", rownames(otu3.ra), ignore.case = TRUE) &
    rownames(otu3.ra) != "Unassigned" &
    rownames(otu3.ra) != "Ambiguous_taxa" &
    rownames(otu3.ra) != "metagenome",
]

##########################################
# 4.3 Transform for heatmap
##########################################
otu3.ra2 <- make_basic_heatmap_matrix(otu3.ra)

##########################################
# 4.4 Draw and save bacterial heatmap
##########################################
heat1 <- Heatmap(
  otu3.ra2,
  col = colorRampPalette(brewer.pal(9, "RdBu"))(100),
  heatmap_legend_param = list(
    title = NULL,
    grid_height = unit(10, "mm"),
    labels_gp = gpar(fontsize = 8)
  ),
  show_row_names = TRUE,
  column_names_rot = 0,
  top_annotation = NULL,
  column_order = c("half", "normal", "double"),
  width = ncol(otu3.ra2) * unit(20, "mm"),
  height = nrow(otu3.ra2) * unit(1.5, "mm"),
  column_names_centered = TRUE,
  column_names_gp = gpar(fontsize = 8),
  row_names_gp = gpar(fontsize = 8)
)

pdf(file.path(output_dir, "Fig1d_bacterial_heatmap.pdf"), width = 5.5, height = 5.5)
draw(heat1)
dev.off()

png(file.path(output_dir, "Fig1d_bacterial_heatmap.png"), width = 5.5, height = 5.5, units = "in", res = 300)
draw(heat1)
dev.off()

############################################################
# 5. Fungal heatmap
############################################################

##########################################
# 5.1 Genus-level fungal abundance table
##########################################
rownames(taxonomy_ITS) <- taxonomy_ITS$Feature.ID
taxonomy_ITS <- taxonomy_ITS[, -1]

ITS2 <- as.data.frame(t(ITS))
samp <- match.name(rn.list = list(ITS2 = ITS2, taxonomy_ITS = taxonomy_ITS))
ITS2 <- samp$ITS2
taxonomy_ITS <- samp$taxonomy_ITS
ITS2 <- cbind(taxonomy_ITS, ITS2)

ITS2 <- aggregate(ITS2[, 9:ncol(ITS2)], by = list(ITS2$Genus), FUN = sum)
rownames(ITS2) <- ITS2$Group.1
ITS2 <- ITS2[, -1]

##########################################
# 5.2 Mean abundance by precipitation treatment
##########################################
ITS3 <- as.data.frame(t(ITS2))
spc <- match.name(rn.list = list(treat = treat, ITS3 = ITS3))
ITS3 <- cbind(spc$treat, spc$ITS3)

ITS3 <- aggregate(ITS3[, 6:ncol(ITS3)], by = list(ITS3$Precipitation), FUN = mean)
rownames(ITS3) <- ITS3$Group.1
ITS3 <- as.data.frame(t(ITS3[, -1]))

ITS3.ra <- decostand(ITS3, "total", MARGIN = 2)
ITS3.ra <- ITS3.ra[rowMeans(ITS3.ra) > 0.001, ]
ITS3.ra <- ITS3.ra[, c(2, 3, 1)]

ITS3.ra <- ITS3.ra[
  !grepl("uncultured", rownames(ITS3.ra), ignore.case = TRUE) &
    rownames(ITS3.ra) != "Unassigned" &
    rownames(ITS3.ra) != "Ambiguous_taxa" &
    rownames(ITS3.ra) != "metagenome",
]

##########################################
# 5.3 Transform for heatmap
##########################################
ITS3.ra2 <- make_basic_heatmap_matrix(ITS3.ra)

##########################################
# 5.4 Draw and save fungal heatmap
##########################################
heat2 <- Heatmap(
  ITS3.ra2,
  col = colorRampPalette(brewer.pal(9, "RdBu"))(100),
  heatmap_legend_param = list(
    title = NULL,
    grid_height = unit(10, "mm"),
    labels_gp = gpar(fontsize = 8)
  ),
  show_row_names = TRUE,
  column_names_rot = 0,
  top_annotation = NULL,
  column_order = c("half", "normal", "double"),
  width = ncol(ITS3.ra2) * unit(20, "mm"),
  height = nrow(ITS3.ra2) * unit(1.5, "mm"),
  column_names_centered = TRUE,
  column_names_gp = gpar(fontsize = 8),
  row_names_gp = gpar(fontsize = 8)
)

pdf(file.path(output_dir, "Fig1e_fungal_heatmap.pdf"), width = 5.5, height = 5.5)
draw(heat2)
dev.off()

png(file.path(output_dir, "Fig1e_fungal_heatmap.png"), width = 5.5, height = 5.5, units = "in", res = 300)
draw(heat2)
dev.off()