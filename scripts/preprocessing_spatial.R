## Analysis of human lung pulmonary fibrosis data
#
# This script is for pre-processing the data.
#
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
analysis_txt_fname <- file.path(DIR_RES, paste0("log.hs_visium_preprocessing_A.txt"))
writeLines("Analysis of human lung visium data generated 2021 by script hs_visium_preprocessing_A.R.", analysis_txt_fname)


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

## Load gene annotations from biomart (as created by original authors)
annotation_gene <- read.table(file = file.path(projectFolder, "data/hs_misc", "hs_gene_biomart.2022-02-23.tsv"),
                              sep = "\t",
                              header = T,
                              stringsAsFactors = F) %>%
  rename(id = ensembl_gene_id, gene_name = hgnc_symbol)

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
# ----Custom functions to use --------------------------------------------------
# A custom function to save plots



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

