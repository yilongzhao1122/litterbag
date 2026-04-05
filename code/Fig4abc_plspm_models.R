############################################################
# Fig. 4a-c: PLS-PM models based on community ordination,
# GeoChip carbon degradation genes, and ecosystem variables
# Project: litterbag
#
# Expected project structure:
# litterbag/
# ├── code/
# │   └── Fig4abc_plspm_gene.R
# ├── data/
# │   ├── 16S_30236.txt
# │   ├── ITS_10391.txt
# │   ├── treat2.csv
# │   ├── env_clear.csv
# │   ├── group_.cut.treat.0.25.log.RA.txt
# │   └── GeoChip_C_sig.txt
# └── output/
#
# This script:
# 1) reads bacterial, fungal, environmental, and GeoChip data,
# 2) prepares treatment-specific datasets for half, control, and
#    double precipitation conditions,
# 3) builds three PLS-PM models corresponding to Fig. 4a-c,
# 4) saves the model plots and summary tables.
#
# Input files used:
# - data/16S_30236.txt
# - data/ITS_10391.txt
# - data/treat2.csv
# - data/env_clear.csv
# - data/group_.cut.treat.0.25.log.RA.txt
# - data/GeoChip_C_sig.txt
#
# Output files:
# - output/plspm/Fig4a_plspm_half.pdf/png
# - output/plspm/Fig4b_plspm_control.pdf/png
# - output/plspm/Fig4c_plspm_double.pdf/png
# - output/plspm/Fig4a_path_coefs.csv
# - output/plspm/Fig4a_inner_model.csv
# - output/plspm/Fig4b_path_coefs.csv
# - output/plspm/Fig4b_inner_model.csv
# - output/plspm/Fig4c_path_coefs.csv
# - output/plspm/Fig4c_inner_model.csv
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
output_dir <- file.path(project_root, "output", "plspm")

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

save.wd <- output_dir

otu_file_path <- file.path(data_dir, "16S_30236.txt")
its_file_path <- file.path(data_dir, "ITS_10391.txt")
treat_file_path <- file.path(data_dir, "treat2.csv")
env_file_path <- file.path(data_dir, "env_clear.csv")
geochip_file_path <- file.path(data_dir, "group_.cut.treat.0.25.log.RA.txt")
geochip_sig_file_path <- file.path(data_dir, "GeoChip_C_sig.txt")

required_files <- c(
  otu_file_path,
  its_file_path,
  treat_file_path,
  env_file_path,
  geochip_file_path,
  geochip_sig_file_path
)

missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop(
    "Missing required input file(s):\n",
    paste(missing_files, collapse = "\n")
  )
}

suppressPackageStartupMessages({
  library(readr)
  library(stringr)
  library(ieggr)
  library(vegan)
  library(lavaan)
  library(plspm)
})

######################
# 1. Read in data
######################
otu.file <- read.table(otu_file_path, header = TRUE, row.names = 1, check.names = FALSE)
ITS.file <- read.table(its_file_path, header = TRUE, row.names = 1, check.names = FALSE)

treat <- read.csv(treat_file_path, row.names = 1, check.names = FALSE)
treat <- subset(treat, treat$Clip == "uncliping")
treat <- treat[, -c(6)]
treat <- subset(treat, type == "litter")

otu = as.data.frame(t(otu.file))
ITS = as.data.frame(t(ITS.file))

otu <- subset(otu, rownames(otu) %in% treat$sample)
ITS <- subset(ITS, rownames(ITS) %in% treat$sample)

otu <- otu[, colSums(otu) != 0]
ITS <- ITS[, colSums(ITS) != 0]

# Environmental data
env <- read.csv(env_file_path, header = TRUE, row.names = 1, check.names = FALSE)
env <- subset(env, Clip == "uncliping")
colnames(env)
env <- env[, -c(22)]

colnames(env)
rownames(env) <- env$sample
eco <- env[, c(1:13, 28:33)]
env <- env[, c(1:24, 27)]

# GeoChip transformed table
GeoChip <- read_delim(geochip_file_path, show_col_types = FALSE)
GeoChip <- as.data.frame(GeoChip)

if (!("uniqueID" %in% colnames(GeoChip))) {
  if (ncol(GeoChip) >= 2 && "uniqueID" %in% colnames(GeoChip)[-1]) {
    GeoChip <- GeoChip[, -1]
  } else {
    stop("GeoChip transformed file must contain a 'uniqueID' column.")
  }
}

rownames(GeoChip) <- GeoChip$uniqueID
GeoChip <- subset(GeoChip, subcategory1 == "Carbon Degradation" | subcategory1 == "Carbon degradation")
GeoChip_C_anno <- unique(GeoChip[, c(4, 8, 10)])

# Significant GeoChip gene table used in the original workflow
GeoChip_C <- read.table(geochip_sig_file_path, header = TRUE, row.names = 1, check.names = FALSE)
GeoChip_C$sum <- rowSums(GeoChip_C)
GeoChip_C$sum_cellulose <- rowSums(GeoChip_C[, c(10, 13, 17, 19)])

######################
# 2. half
######################
env_half <- subset(env, Precipitation == "half")

sampc = match.name(rn.list = list(env_half = env_half, otu = otu, ITS = ITS, GeoChip_C = GeoChip_C, eco = eco))
otu_half = sampc$otu
env_half = sampc$env_half
ITS_half = sampc$ITS
GeoChip_C_half = sampc$GeoChip_C
eco_half = sampc$eco

otu_half <- otu_half[, colSums(otu_half) > 0]
ITS_half <- ITS_half[, colSums(ITS_half) > 0]
GeoChip_C_half <- GeoChip_C_half[, colSums(GeoChip_C_half) > 0]

# PCoA
bray_otu_half <- vegdist(otu_half, method = "bray", binary = TRUE)
pcoa_otu_half <- cmdscale(bray_otu_half, k = (nrow(otu_half) - 1), eig = TRUE)
pcoa_otu_half$eig / sum(pcoa_otu_half$eig)
site_otu_half <- scores(pcoa_otu_half)
colnames(site_otu_half)[1] <- "otuPC1"

bray_ITS_half <- vegdist(ITS_half, method = "bray", binary = TRUE)
pcoa_ITS_half <- cmdscale(bray_ITS_half, k = (nrow(ITS_half) - 1), eig = TRUE)
pcoa_ITS_half$eig / sum(pcoa_ITS_half$eig)
site_ITS_half <- scores(pcoa_ITS_half)
colnames(site_ITS_half)[1] <- "ITSPC1"

data_half <- cbind(env_half, eco_half[, 14:19], site_otu_half, site_ITS_half, GeoChip_C_half)
data_half <- scale(data_half[, 14:ncol(data_half)])

dat_blocks_half <- list(
  tem = c("temperature_annual"),
  moi = c("annual_moisture", "moisture_samplingmonth"),
  soil = c("TN", "TC"),
  plant = c("FlC4", "FlC3", "plant.richness"),
  gene = c("Cdh", "Vdh"),
  rs = c("Autotrophic", "Heterotrophic")
)

tem <- c(0, 0, 0, 0, 0, 0)
moi <- c(0, 0, 0, 0, 0, 0)
soil <- c(0, 1, 0, 0, 0, 0)
plant <- c(1, 1, 0, 0, 0, 0)
gene <- c(1, 0, 0, 1, 0, 0)
rs <- c(0, 0, 0, 1, 1, 0)

dat_path_half <- rbind(tem, moi, soil, plant, gene, rs)
colnames(dat_path_half) <- rownames(dat_path_half)
dat_path_half

dat_modes_half <- rep("A", 6)
dat_modes_half

dat_pls_half <- plspm(data_half, dat_path_half, dat_blocks_half, modes = dat_modes_half)
dat_pls_half
summary(dat_pls_half)

dat_pls_half$path_coefs
dat_pls_half$inner_model

innerplot(dat_pls_half, colpos = "red", colneg = "blue", show.values = TRUE, lcol = "gray", box.lwd = 0)
dat_pls_half$gof

pdf(file.path(save.wd, "Fig4a_plspm_half.pdf"), width = 7, height = 5.5)
innerplot(dat_pls_half, colpos = "red", colneg = "blue", show.values = TRUE, lcol = "gray", box.lwd = 0)
dev.off()

png(file.path(save.wd, "Fig4a_plspm_half.png"), width = 7, height = 5.5, units = "in", res = 300)
innerplot(dat_pls_half, colpos = "red", colneg = "blue", show.values = TRUE, lcol = "gray", box.lwd = 0)
dev.off()

######################
# 3. normal
######################
env_normal <- subset(env, Precipitation == "normal")

sampc = match.name(rn.list = list(env_normal = env_normal, otu = otu, ITS = ITS, GeoChip_C = GeoChip_C, eco = eco))
otu_normal = sampc$otu
env_normal = sampc$env_normal
ITS_normal = sampc$ITS
GeoChip_C_normal = sampc$GeoChip_C
eco_normal = sampc$eco

otu_normal <- otu_normal[, colSums(otu_normal) > 0]
ITS_normal <- ITS_normal[, colSums(ITS_normal) > 0]
GeoChip_C_normal <- GeoChip_C_normal[, colSums(GeoChip_C_normal) > 0]

# PCoA
bray_otu_normal <- vegdist(otu_normal, method = "bray", binary = TRUE)
pcoa_otu_normal <- cmdscale(bray_otu_normal, k = (nrow(otu_normal) - 1), eig = TRUE)
pcoa_otu_normal$eig / sum(pcoa_otu_normal$eig)
site_otu_normal <- scores(pcoa_otu_normal)
colnames(site_otu_normal)[1] <- "otuPC1"

bray_ITS_normal <- vegdist(ITS_normal, method = "bray", binary = TRUE)
pcoa_ITS_normal <- cmdscale(bray_ITS_normal, k = (nrow(ITS_normal) - 1), eig = TRUE)
pcoa_ITS_normal$eig / sum(pcoa_ITS_normal$eig)
site_ITS_normal <- scores(pcoa_ITS_normal)
colnames(site_ITS_normal)[1] <- "ITSPC1"

data_normal <- cbind(env_normal, eco_normal[, 14:19], site_otu_normal, site_ITS_normal, GeoChip_C_normal)
data_normal <- scale(data_normal[, 14:ncol(data_normal)])

dat_blocks_normal <- list(
  tem_normal = c("temperature_annual"),
  moi_normal = c("annual_moisture", "moisture_samplingmonth"),
  soil_normal = c("TC"),
  plant_normal = c("FlC4"),
  gene_normal = c("Alpha_galactosidase_fungi", "Lmo", "Limeh", "Glx", "Chitin_deacetylase_fungi",
                  "Endoglucanase", "Cellobiase", "AmyA", "Lactase_fungi", "AceA"),
  rs_normal = c("total_soil_respiration")
)

tem_normal <- c(0, 0, 0, 0, 0, 0)
moi_normal <- c(0, 0, 0, 0, 0, 0)
soil_normal <- c(0, 1, 0, 0, 0, 0)
plant_normal <- c(0, 1, 0, 0, 0, 0)
gene_normal <- c(1, 0, 1, 0, 0, 0)
rs_normal <- c(0, 0, 0, 1, 1, 0)

dat_path_normal <- rbind(tem_normal, moi_normal, soil_normal, plant_normal, gene_normal, rs_normal)
colnames(dat_path_normal) <- rownames(dat_path_normal)
dat_path_normal
dat_modes_normal <- rep("A", 6)
dat_modes_normal

dat_pls_normal <- plspm(data_normal, dat_path_normal, dat_blocks_normal, modes = dat_modes_normal)
dat_pls_normal
summary(dat_pls_normal)

dat_pls_normal$path_coefs
dat_pls_normal$inner_model

innerplot(dat_pls_normal, colpos = "red", colneg = "blue", show.values = TRUE, lcol = "gray", box.lwd = 0)
dat_pls_normal$gof

pdf(file.path(save.wd, "Fig4b_plspm_control.pdf"), width = 7, height = 5.5)
innerplot(dat_pls_normal, colpos = "red", colneg = "blue", show.values = TRUE, lcol = "gray", box.lwd = 0)
dev.off()

png(file.path(save.wd, "Fig4b_plspm_control.png"), width = 7, height = 5.5, units = "in", res = 300)
innerplot(dat_pls_normal, colpos = "red", colneg = "blue", show.values = TRUE, lcol = "gray", box.lwd = 0)
dev.off()

######################
# 4. double
######################
env_double <- subset(env, Precipitation == "double")

sampc = match.name(rn.list = list(env_double = env_double, otu = otu, ITS = ITS, GeoChip_C = GeoChip_C, eco = eco))
otu_double = sampc$otu
env_double = sampc$env_double
ITS_double = sampc$ITS
GeoChip_C_double = sampc$GeoChip_C
eco_double = sampc$eco

otu_double <- otu_double[, colSums(otu_double) > 0]
ITS_double <- ITS_double[, colSums(ITS_double) > 0]
GeoChip_C_double <- GeoChip_C_double[, colSums(GeoChip_C_double) > 0]

# PCoA
bray_otu_double <- vegdist(otu_double, method = "bray", binary = TRUE)
pcoa_otu_double <- cmdscale(bray_otu_double, k = (nrow(otu_double) - 1), eig = TRUE)
pcoa_otu_double$eig / sum(pcoa_otu_double$eig)
site_otu_double <- scores(pcoa_otu_double)
colnames(site_otu_double)[1] <- "otuPC1"

bray_ITS_double <- vegdist(ITS_double, method = "bray", binary = TRUE)
pcoa_ITS_double <- cmdscale(bray_ITS_double, k = (nrow(ITS_double) - 1), eig = TRUE)
pcoa_ITS_double$eig / sum(pcoa_ITS_double$eig)
site_ITS_double <- scores(pcoa_ITS_double)
colnames(site_ITS_double)[1] <- "ITSPC1"

data_double <- cbind(env_double, eco_double[, 14:19], site_otu_double, site_ITS_double, GeoChip_C_double)
data_double <- scale(data_double[, 14:ncol(data_double)])

dat_blocks_double <- list(
  tem_double = c("temperature_annual"),
  moi_double = c("moisture_samplingmonth"),
  soil_double = c("TC"),
  plant_double = c("FlC4", "FlC3"),
  gene_double = "sum_cellulose",
  rs_double = c("Autotrophic", "Heterotrophic")
)

tem_double <- c(0, 0, 0, 0, 0, 0)
moi_double <- c(0, 0, 0, 0, 0, 0)
soil_double <- c(1, 0, 0, 0, 0, 0)
plant_double <- c(0, 1, 1, 0, 0, 0)
gene_double <- c(0, 1, 1, 0, 0, 0)
rs_double <- c(0, 0, 0, 1, 1, 0)

dat_path_double <- rbind(tem_double, moi_double, soil_double, plant_double, gene_double, rs_double)
colnames(dat_path_double) <- rownames(dat_path_double)
dat_path_double
dat_modes_double <- rep("A", 6)
dat_modes_double

dat_pls_double <- plspm(data_double, dat_path_double, dat_blocks_double, modes = dat_modes_double)
dat_pls_double
summary(dat_pls_double)

dat_pls_double$path_coefs
dat_pls_double$inner_model

innerplot(dat_pls_double, colpos = "red", colneg = "blue", show.values = TRUE, lcol = "gray", box.lwd = 0)
dat_pls_double$gof

pdf(file.path(save.wd, "Fig4c_plspm_double.pdf"), width = 7, height = 5.5)
innerplot(dat_pls_double, colpos = "red", colneg = "blue", show.values = TRUE, lcol = "gray", box.lwd = 0)
dev.off()

png(file.path(save.wd, "Fig4c_plspm_double.png"), width = 7, height = 5.5, units = "in", res = 300)
innerplot(dat_pls_double, colpos = "red", colneg = "blue", show.values = TRUE, lcol = "gray", box.lwd = 0)
dev.off()