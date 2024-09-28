# ---------------------------------------------------------------------------- #
# Script Name: <name_of_the_script>.R
# Author: Eleftherios Zormpas
# Date: <date_of_creation_or_modification>
#
# Description:
#   - Analysis of human lung pulmonary fibrosis data
#   - This script is for pre-processing the spatial data.
#   - Brief summary of the script’s purpose.
#   - Highlight any specific biological or computational goals, e.g. data preprocessing,
#     spatial clustering, differential gene expression analysis, etc.
#
# Input:
#   - <List of input files, formats, and paths (e.g. .csv, .rds)>
#   - Specify any external datasets being used (single-cell, spatial omics, etc.)
#
# Output:
#   - <Describe the output of the script: plots, processed data, or analysis results>
#
# Dependencies:
#   - <List of required R packages (e.g. Seurat, scran, spatialDE)>
#   - <Any external scripts or tools used>
#
# Instructions:
#   - How to run the script, including any required parameters or settings.
#   - Example command or execution method if applicable.
#
# Notes:
#   - <Any assumptions, considerations, or specific details to keep in mind>
#   - <References or links to relevant research articles or documentation>
#
# Version History:
#   - <Version updates, e.g. added new functionality, fixed bugs>
# ---------------------------------------------------------------------------- #

####################### The preprocessing starts here ##########################

# ---------------------------------------------------------------------------- #
# ----Load packages ------------------------------------------------------------
## Data manipulation/ plotting
library(ggplot2)
library(tidyverse)
library(patchwork)
library(ggpubr)
library(cols4all)
library(tidyquant)
library(ggdist)
library(ggthemes)

## Analysis
library(GSVA)
library(harmony)
library(Seurat)
library(STExplorer)


# ---------------------------------------------------------------------------- #
# ----Set parameters -----------------------------------------------------------
set.seed(1)
DIR_PROJECT <- getwd()
DIR_DATA <- paste0(DIR_PROJECT, "/data/spatial/hs_visium_spaceranger_output/")
DIR_METADATA <- paste0(DIR_PROJECT, "/data/spatial/hs_misc/")
DIR_SAMPLES <- list.dirs(DIR_PROJECT, recursive = FALSE)
DIR_RES <- file.path(DIR_ROOT, "graphics_out")
fig_res <- 300


# ---------------------------------------------------------------------------- #
# ----Create log file ----------------------------------------------------------
#' Create a log file to store parameter settings and intermediate information
analysis_txt_fname <- file.path(DIR_RES, paste0("logs/log.hs_visium_preprocessing_A.txt"))
writeLines("Analysis of human lung visium data generated 2021 by script hs_visium_preprocessing_A.R.", analysis_txt_fname)


# ---------------------------------------------------------------------------- #
# ----Source functions ---------------------------------------------------------
source(file.path(DIR_ROOT, "scripts", "custom_functions.R"))


# ---------------------------------------------------------------------------- #
# ----Load Metadata ------------------------------------------------------------
## The metadata in this section come from the original publication's available
## data to download
##
## Load the metadata
metaDataFile <- paste0(projectFolder, DIR_METADATA, "hs_visium_metadata.tsv")
metadata <- read.table(metaDataFile, header = TRUE) %>%
  mutate(fibrotic_extent_score_by_pathologist_0.3 = as.character(fibrotic_extent_score_by_pathologist_0.3))
rownames(metadata) <- metadata$sample_id

## Load spot annotations
annotationFile <- paste0(projectFolder, DIR_METADATA, "hs_visium_merged_histo_annotations.csv")
annotation <- read.table(annotationFile, header = TRUE, sep = ",") %>%
  rename(Barcode = barcode, sample_id = sample_slide_id)

## Load cell density-related information as generated using Cell2location
## --> see script analysis_cellDensities.py
cellAbundanceFolder <- paste0(projectFolder, "/data/cell2location_habermann2020/")
cellGroupsFile <- paste0(projectFolder, "/data/hs_misc/habermann_cell_type_groups.csv")
cellGroups <- read.table(cellGroupsFile, header = TRUE, sep = ";")

# ---------------------------------------------------------------------------- #
# ----Set samples --------------------------------------------------------------
sampleNames <- list.files(DIR_DATA)
sampleNames <- sampleNames[sampleNames != "README.md"]
names(sampleDir) <- sampleNames

## Update metadata to keep only the selected samples
## This is because the dataset includes 25 samples but we are working only with
## the 16 that have pathologist annotations provided by the authors.
metadata <- metadata[metadata$sample_id %in% sampleNames,]



# ---------------------------------------------------------------------------- #
# ----Load data ----------------------------------------------------------------
## Create the MSFE object
msfe <- MetaSpatialFeatureExperiment()

## Import samples
for (i in seq_along(sampleNames)) {
  id <- sampleNames[i]
  message("Adding sample: ", sampleNames[i])
  msfe <- addSFE(msfe,
                 read10xVisiumSFE(samples = sampleDir[id],
                                  sample_id = sampleNames[i],
                                  type = "HDF5",
                                  data = "filtered",
                                  images = "lowres",
                                  style = "W",
                                  zero.policy = TRUE))
}

## Import annotations
gTruth_list <- list()
for (id in sampleNames) {
  gTruth_list[[id]] <- annotation[annotation$sample_id == id, c("Barcode", "sample_id", "annotation")]
}

## Import cell type abundances
cellAbundance <- vector(mode = "list", length = length(sampleNames))
names(cellAbundance) <- sampleNames
for (s in sampleNames) {
  file <- paste0(cellAbundanceFolder, s, "_spot_cell_abundances_5pc.csv")
  df <- read.csv(file) %>%
    rename(Barcode = spot_id) %>%
    mutate(Barcode = gsub(paste0(s, "_"), "", Barcode),
           sample_id = s, .after = Barcode)
  cellAbundance[[s]] <- df
}


# ---------------------------------------------------------------------------- #
# ----Select genes -------------------------------------------------------------
## Get annotations
biomart <- createBiomart()
annotation_gene <- data.frame(gene_name = rowData(msfe@sfe_data[[1]])[["symbol"]])
rownames(annotation_gene) <- rownames(msfe@sfe_data[[1]])
annotation_gene <- annotateDataFrame(annotation_gene,
                                     biomart = biomart,
                                     add = c("biotype", "chromosome", "description"))

## Senescence Markers
senescenceMarkers <- c("CDKN2A", "CDKN1A", "GLB1", "LMNB1")
senescenceMarkers_ENSG <- setNames(
  annotation_gene$id[annotation_gene$gene_name %in% senescenceMarkers],
  annotation_gene$gene_name[annotation_gene$gene_name %in% senescenceMarkers]
)

## Genes Related to senescence and fibrosis
otherGenes <- c("IL6", "IL6R", "IL6ST", "GDF15", "CCL2", "TGFB1", "FAP", "VIM",
                "ACTA2", "COL1A1", "COL1A2", "FN1", "CDH1", "CDH5")
otherGenes_ENSG <- setNames(
  annotation_gene$id[annotation_gene$gene_name %in% otherGenes],
  annotation_gene$gene_name[annotation_gene$gene_name %in% otherGenes]
)

## SenMayo senescence markers group
senMayoGenes <- read_csv(paste0(DIR_METADATA), "geneset_senmayo.txt")
senMayoGenes_ENSG <- setNames(
  annotation_gene$id[annotation_gene$gene_name %in% senMayoGenes],
  annotation_gene$gene_name[annotation_gene$gene_name %in% senMayoGenes]
)

## Epithelial to mesenchymal transition gene set from MSigDB HALLMARK GENESET
emtGenes <- read_csv(paste0(DIR_METADATA), "geneset_emt")
emtGenes_ENSG <- setNames(
  annotation_gene$id[annotation_gene$gene_name %in% emtGenes],
  annotation_gene$gene_name[annotation_gene$gene_name %in% emtGenes]
)

## GO Extracellular Matrix gene set from MSigDB GO::CC GENESET
ecmGenes <- read_csv(paste0(DIR_METADATA), "geneset_ecm.txt")
ecmGenes_ENSG <- setNames(
  annotation_gene$id[annotation_gene$gene_name %in% ecmGenes],
  annotation_gene$gene_name[annotation_gene$gene_name %in% ecmGenes]
) # Removed "ANXA2P2" "SSPOP" because are marked as Pseudogenes

## TGF-beta & IL4 gene set as found in : DOI:
tgfbGenes <- read_csv(paste0(DIR_METADATA), "geneset_tgfb.txt")
tgfbGenes_ENSG <- setNames(
  annotation_gene$id[annotation_gene$gene_name %in% tgfbGenes],
  annotation_gene$gene_name[annotation_gene$gene_name %in% tgfbGenes]
) # Removed "AC083837.2" "GBP1P1" because are marked as Pseudogenes


# ---------------------------------------------------------------------------- #
# ----Analysis -----------------------------------------------------------------


# ---------------------------------------------------------------------------- #
## ----QC ----------------------------------------------------------------------
### ----Calculate QC Metrics ---------------------------------------------------
## Mark a subset of mitochondrial genes
is_mito <- getSubset(msfe,
                     sample_id = TRUE,
                     subset = "(^MT-)|(^mt-)|(^MTRNR)",
                     set = "rowData",
                     col_name = "symbol")
## Mark a subset of haemoglobin genes
is_hemo <- getSubset(msfe,
                     sample_id = TRUE,
                     subset = "^HB",
                     set = "rowData",
                     col_name = "symbol")
## Mark a subset of Ribosomal genes
is_ribo <- getSubset(msfe,
                     sample_id = TRUE,
                     subset = "^RPS|^RPL",
                     set = "rowData",
                     col_name = "symbol")

for (id in sampleNames) {
  message("Working on sample: ", id)
  ## Add location-related statistics
  msfe <- addPerLocQC(msfe,
                      sample_id = id,
                      gTruth = gTruth_list[[id]],
                      assay = "counts",
                      MARGIN = 2,
                      subsets = list(mito = is_mito[[id]],
                                     hemo = is_hemo[[id]],
                                     ribo = is_ribo[[id]]))
  message("\tAdded location-related statistics")

  ## Add geometries
  msfe <- addGeometries(msfe,
                        samples = sampleDir[id],
                        sample_id = id,
                        res = "fullres",
                        flipped = FALSE,
                        geoms = "both")
  message("\tAdded geometries")

  ## Add gene/feature-related statistics
  msfe <- addPerGeneQC(msfe,
                       sample_id = id,
                       assay = "counts",
                       version = NULL,
                       mirror = NULL,
                       add = "none")
  message("\tAdded gene/feature-related statistics")
}


# ---------------------------------------------------------------------------- #
### ----Quality Control: Spots--------------------------------------------------
#### ----Library size ----------------------------------------------------------
## Create folder
graphics_out <- paste0(projectFolder, "/graphics_out/QC_Spots/QC_library_size/")
dir.create(graphics_out, recursive = TRUE)

for (id in sampleNames) {
  message("Working on sample: ", id)

  ## Extract sample
  sfe <- getSFE(msfe, sample_id = id)

  ## Select filter
  min <- 350 # at least 350 unique molecules per spot

  ## Density and histogram of library sizes
  p1 <- plotQC_hist(sfe, metric = "libsize",
                    limits = c(min, NA),
                    hist_args = list(bins = 100)) +
    theme(axis.text.x = element_text(size = 12, angle = 45))

  ## Map the library sizes
  p2 <- plotQC_map(sfe, metric = "libsize", type = "hex", size = 0.5)

  ## Select threshold
  sfe <- setQCthresh_LibSize(sfe, sample_id = TRUE,
                             min_t = min,
                             max_t = NA)

  ## Check putative spatial patterns of removed spots
  p3 <- plotQC_filtered(sfe, metric = "libsize", sample_id = TRUE, type = "hex", size = 0.5)

  ## Plot tissue image
  p4 <- plotQC_tissueImg(sfe, res = "lowres", sample_id = id, type = "none")

  ## Patch, print and save
  print((p1 + p2) / (p3 + p4))

  prfx <- id
  main <- "_QC_LibSize"
  sfx <- ""
  other <- ""

  ggsave(paste0(graphics_out, prfx, main, sfx, other, ".png"),
         device = "png",
         width = grDevices::dev.size(units = "in")[1],
         height = grDevices::dev.size(units = "in")[2],
         units = "in",
         dpi = 300)

  ## Update msfe object
  msfe@sfe_data[[id]] <- sfe

  ## Housekeeping
  rm(sfe, p1, p2, p3, p4)
}


# ---------------------------------------------------------------------------- #
#### ----Number of expressed genes ---------------------------------------------
## Create folder
graphics_out <- paste0(projectFolder, "/graphics_out/QC_Spots/QC_gene_number/")
dir.create(graphics_out, recursive = TRUE)

for (id in sampleNames) {
  message("Working on sample: ", id)

  ## Extract sample
  sfe <- getSFE(msfe, sample_id = id)

  ## Select filter
  min <- 100 # at least 100 genes per spot

  ## Density and histogram of expressed genes
  p1 <- plotQC_hist(sfe, metric = "detected", c(min, NA))

  ## Map the library sizes
  p2 <- plotQC_map(sfe, metric = "detected", type = "hex", size = 0.5)

  ## Select threshold
  sfe <- setQCthresh_GenesExpr(sfe, sample_id = TRUE,
                               min_t = min,
                               max_t = NA)

  ## Check putative spatial patterns of removed spots
  p3 <- plotQC_filtered(sfe, metric = "detected", sample_id = TRUE, type = "hex", size = 0.5)

  ## Patch, print and save
  print(p1 + p2 + p3)

  prfx <- id
  main <- "_QC_GeneNumber"
  sfx <- ""
  other <- ""

  ggsave(paste0(graphics_out, prfx, main, sfx, other, ".png"),
         device = "png",
         width = grDevices::dev.size(units = "in")[1],
         height = grDevices::dev.size(units = "in")[2],
         units = "in",
         dpi = 300)

  ## Update msfe object
  msfe@sfe_data[[id]] <- sfe

  ## Housekeeping
  rm(sfe, p1, p2, p3)
}


# ---------------------------------------------------------------------------- #
#### ----Percent of mito expression --------------------------------------------
## Create folder
graphics_out <- paste0(projectFolder, "/graphics_out/QC_Spots/QC_mito_percent/")
dir.create(graphics_out, recursive = TRUE)

for (id in sampleNames) {
  message("Working on sample: ", id)

  ## Extract sample
  sfe <- getSFE(msfe, sample_id = id)

  ## Density and histogram of percentage of mitochondrial expression
  p1 <- plotQC_hist(sfe, metric = "mito", limits = c(NA, 15))

  ## Map the library sizes
  p2 <- plotQC_map(sfe, metric = "mito", type = "hex", size = 0.5)

  ## Select threshold
  sfe <- setQCthresh_Mito(sfe, sample_id = TRUE, min_t = NA, max_t = 15)

  ## Check putative spatial patterns of removed spots
  p3 <- plotQC_filtered(sfe, metric = "mito", sample_id = TRUE, type = "hex", size = 0.5)

  ## Patch, print and save
  print(p1 + p2 + p3)

  prfx <- id
  main <- "_QC_MitoPercent"
  sfx <- ""
  other <- ""

  ggsave(paste0(graphics_out, prfx, main, sfx, other, ".png"),
         device = "png",
         width = grDevices::dev.size(units = "in")[1],
         height = grDevices::dev.size(units = "in")[2],
         units = "in",
         dpi = 300)

  ## Update msfe object
  msfe@sfe_data[[id]] <- sfe

  ## Housekeeping
  rm(sfe, p1, p2, p3)
}


# ---------------------------------------------------------------------------- #
#### ----Percent of haemoglobin expression -------------------------------------
## Create folder
graphics_out <- paste0(projectFolder, "/graphics_out/QC_Spots/QC_hemo_percent/")
dir.create(graphics_out, recursive = TRUE)

for (id in sampleNames) {
  message("Working on sample: ", id)

  ## Extract sample
  sfe <- getSFE(msfe, sample_id = id)

  ## Get the vector
  hemo_out <- list(id = sfe@colData$subsets_hemo_percent > 30)
  names(hemo_out) <- id

  ## Map the library sizes
  p1 <- plotQC_map(sfe, metric = "custom", type = "hex", size = 0.5, metric_name = "subsets_hemo_percent")

  ## Select threshold
  sfe <- setQCthresh_custom(sfe, sample_id = TRUE, MARGIN = 2, qcMetric = hemo_out)

  ## Check putative spatial patterns of removed spots
  p2 <- plotQC_filtered(sfe,
                        metric = "custom",
                        sample_id = TRUE,
                        type = "hex",
                        size = 0.5,
                        metric_name = "qc_hemo_out",
                        metric_lab = "Filtered for % Hemoglobin Expression")

  ## Patch, print and save
  print(p1 + p2)

  prfx <- id
  main <- "_QC_HemoPercent"
  sfx <- ""
  other <- ""

  ggsave(paste0(graphics_out, prfx, main, sfx, other, ".png"),
         device = "png",
         width = grDevices::dev.size(units = "in")[1],
         height = grDevices::dev.size(units = "in")[2],
         units = "in",
         dpi = 300)

  ## Update msfe object
  msfe@sfe_data[[id]] <- sfe

  ## Housekeeping
  rm(sfe, p1, p2, hemo_out)
}


# ---------------------------------------------------------------------------- #
#### ----Remove low-quality spots ----------------------------------------------
## Create folder
graphics_out <- paste0(projectFolder, "/graphics_out/QC_Spots/QC_filter_spots/")
dir.create(graphics_out, recursive = TRUE)

for (id in sampleNames) {
  message("Working on sample: ", id)

  ## Extract sample
  sfe <- getSFE(msfe, sample_id = id)

  ## Set the combined filtering threshold using the QC metrics
  sfe <- setQCtoDiscard_loc(sfe, sample_id = TRUE, filters = TRUE)

  ## Check putative spatial patterns of removed spots
  plotQC_filtered(sfe, metric = "discard", sample_id = TRUE, type = "hex")

  prfx <- id
  main <- "_QC_SpotsFiltered"
  sfx <- ""
  other <- ""

  ggsave(paste0(graphics_out, prfx, main, sfx, other, ".png"),
         device = "png",
         width = grDevices::dev.size(units = "in")[1],
         height = grDevices::dev.size(units = "in")[2],
         units = "in",
         dpi = 300)

  ## Remove combined set of low-quality spots
  sfe <- applyQCthresh_loc(sfe, sample_id = TRUE)

  ## Update msfe object
  msfe@sfe_data[[id]] <- sfe

  ## Housekeeping
  rm(sfe)
}


# ---------------------------------------------------------------------------- #
### ----Quality Control: Features ----------------------------------------------

# ---------------------------------------------------------------------------- #
#### ----Mark specific gene sets -----------------------------------------------
## Genes in the XY chromosomes
xychr_genes <- unique(subset(annotation_gene, annotation_gene$chromosome %in% c("X", "Y"))$gene_name)
is_xy <- lapply(sampleNames, function(x){rowData(msfe@sfe_data[[x]])[["gene_name"]] %in% xychr_genes})
names(is_xy) <- sampleNames

## Keep only protein coding, IG and TR genes
keep_biotype <- c("protein_coding", grep("_gene", unique(annotation_gene$biotype), value = TRUE))
coding_genes <- annotation_gene[annotation_gene$biotype %in% keep_biotype, "gene_name"]
is_ncGene <- lapply(sampleNames, function(x){!rowData(msfe@sfe_data[[x]])[["gene_name"]] %in% coding_genes$gene_name})
names(is_ncGene) <- sampleNames

## Remove mitochondrial, hemoglobin, ribosomal, XY, and anything NOT protein coding or IG or TR.
for (id in sampleNames) {
  msfe <- setQCthresh_custom(msfe, MARGIN = 1, qcMetric = is_mito[[id]])
  msfe <- setQCthresh_custom(msfe, MARGIN = 1, qcMetric = is_hemo[[id]])
  msfe <- setQCthresh_custom(msfe, MARGIN = 1, qcMetric = is_ribo[[id]])
  msfe <- setQCthresh_custom(msfe, MARGIN = 1, qcMetric = is_xy[[id]])
  msfe <- setQCthresh_custom(msfe, MARGIN = 1, qcMetric = is_ncGene[[id]])
}


# ---------------------------------------------------------------------------- #
#### ----Remove pre-specified gene sets ----------------------------------------
## QC discard Features
## Set the combined filtering threshold using the QC metrics
msfe <- setQCtoDiscard_feat(msfe, filters = TRUE)

## Apply gene-level QC threshold
msfe <- applyQCthresh_feat(msfe)


# ---------------------------------------------------------------------------- #
## ----Normalisation of Counts -------------------------------------------------
### ----Get MultiSFE object ----------------------------------------------------
sfe_multi <- getMultiSFE(msfe, sampleNames)

## Add the subject alias in the colData to regress it out downstream
sfe_multi@colData$subject_alias <- left_join(as.data.frame(sfe_multi@colData["sample_id"]),
                                             metadata,
                                             by = "sample_id")[["subject_alias"]]

## Remove genes with zero expression in all samples
sfe_multi <- setQCthresh_ZeroExpr(sfe_multi)
sfe_multi <- setQCtoDiscard_feat(sfe_multi)
sfe_multi <- applyQCthresh_feat(sfe_multi)


### ----Convert to Seurat object -----------------------------------------------
seu <- as.Seurat(sfe_multi, data = NULL)


### ----SCTransform counts -----------------------------------------------------
## Normalise and regress out sampleID and donor, to remove major effects from
## technical and interindividual differences.
vars_to_reg <- c("sample_id", "subject_alias")

options(future.globals.maxSize = 10000 * 1024^2)

seu <- SCTransform(seu, vars.to.regress = vars_to_reg, assay = "originalexp", conserve.memory = TRUE)


## ----Generate Module Scores --------------------------------------------------
### ----Senescence -------------------------------------------------------------
seu <- AddModuleScore(seu, features = list(senMayo = senMayoGenes_ENSG), assay = "SCT", name = "senMayo")

### ----EMT --------------------------------------------------------------------
seu <- AddModuleScore(seu, features = list(emt = emtGenes_ENSG), assay = "SCT", name = "emt")

### ----ECM --------------------------------------------------------------------
seu <- AddModuleScore(seu, features = list(ecm = ecmGenes_ENSG), assay = "SCT", name = "ecm")

### ----TGFB --------------------------------------------------------------------
seu <- AddModuleScore(seu, features = list(tgfb = tgfbGenes_ENSG), assay = "SCT", name = "tgfb")

## Some genes from the above feature sets are not present in the data set.
## These genes have been filtered out. However, it is good to know which are
## these genes:
cat("Genes from SenMayo gene set not present in the dataset:", sort(names(senMayoGenes_ENSG[!senMayoGenes_ENSG %in% rownames(seu)])), "\n")
cat("Genes from EMT gene set not present in the dataset:", sort(names(emtGenes_ENSG[!emtGenes_ENSG %in% rownames(seu)])), "\n")
cat("Genes from ECM gene set not present in the dataset:", sort(names(ecmGenes_ENSG[!ecmGenes_ENSG %in% rownames(seu)])), "\n")
cat("Genes from TGFB gene set not present in the dataset:", sort(names(tgfbGenes_ENSG[!tgfbGenes_ENSG %in% rownames(seu)])), "\n")


### ----Convert back to SFE ----------------------------------------------------
## Move the SCT counts and the senMayo/emt/ecm/tgfb scores into the sfe_multi object.
## Then split the sfe_multi object into multiple SFE objects inside a MetaSFE

## 1. First remove the genes that have been filtered out after SCTransform
to_keep <- rownames(seu)
sfe_multi <- sfe_multi[to_keep,]

## 2. Add the SCTransformed counts
assay(sfe_multi, "sct") <- seu[["SCT"]]$data

## 3. Add the senMayo/emt/ecm/tgfb
colData(sfe_multi)$senMayo <- as.vector(seu$senMayo1)
colData(sfe_multi)$emt <- as.vector(seu$emt1)
colData(sfe_multi)$ecm <- as.vector(seu$ecm1)
colData(sfe_multi)$tgfb <- as.vector(seu$tgfb1)

## 4. Split into single SFEs inside an MSFE object
msfe <- addMultiSFE(msfe, sfe_multi)


### ----Plots ------------------------------------------------------------------



# ---------------------------------------------------------------------------- #
# ----Wrap up -----------------------------------------------------------------
saveRDS(sfe_multi, file = file.path(DIR_PROJECT, "analysis_objs/spatial/sfe_multi.rds"))
saveRDS(seu, file = file.path(DIR_PROJECT, "nalysis_objs/spatial/seu.rds"))
saveRDS(msfe, file = file.path(DIR_PROJECT, "analysis_objs/spatial/msfe.rds"))

