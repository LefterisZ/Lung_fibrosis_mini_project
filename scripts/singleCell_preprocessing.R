# ---------------------------------------------------------------------------- #
# Script Name: <name_of_the_script>.R
# Author: Eleftherios Zormpas
# Date: <date_of_creation_or_modification>
#
# Description:
#   - Analysis of human lung pulmonary fibrosis data
#   - This script is for pre-processing the single-cell data.
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
library(readr)
library(tidyverse)

## Analysis
library(Seurat)
library(SeuratDisk)


# ---------------------------------------------------------------------------- #
# ----Set parameters -----------------------------------------------------------
set.seed(1)
DIR_PROJECT <- getwd()
DIR_DATA <- paste0(DIR_PROJECT, "/data/single_cell/GSE135893/")
DIR_METADATA <- paste0(DIR_PROJECT, "/data/single_cell/misc/")
DIR_RES <- file.path(DIR_ROOT, "graphics_out")


# ---------------------------------------------------------------------------- #
# ----Create log file ----------------------------------------------------------
#' Create a log file to store parameter settings and intermediate information
analysis_txt_fname <- file.path(DIR_RES, paste0("logs/log.hs_visium_preprocessing_A.txt"))
writeLines("Analysis of human lung visium data generated 2021 by script hs_visium_preprocessing_A.R.", analysis_txt_fname)


# ---------------------------------------------------------------------------- #
# ----Source functions ---------------------------------------------------------
source(file.path(DIR_ROOT, "scripts", "custom_functions.R"))


# ---------------------------------------------------------------------------- #
# ----Load metadata ------------------------------------------------------------
metadata <- read_csv(file.path(DIR_METADATA, "IPF_metadata.csv")) %>%
  column_to_rownames(var = "...1") %>%
  mutate(Barcode = rownames(.), .before = "orig.ident") %>%
  arrange(Barcode)


# ---------------------------------------------------------------------------- #
# ----Load data ----------------------------------------------------------------
## Read counts matrix
seu <- Read10X(DIR_DATA,
               gene.column = 1,
               cell.column = 1,
               unique.features = TRUE,
               strip.suffix = FALSE)
## Create Seurat object
seu <- CreateSeuratObject(counts = seu)

## Add metadata from the authors
seu <- AddMetaData(seu, metadata = metadata)


# ---------------------------------------------------------------------------- #
# ----QC -----------------------------------------------------------------------
## Identify percentage of mitochondrial reads in each cell
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")

## Visualize QC metrics as a violin plot
p1 <- VlnPlot(seu, features = "nFeature_RNA", pt.size=0, group.by="Sample_Name", split.by = "Status")
p2 <- VlnPlot(seu, features = "nCount_RNA", pt.size=0, group.by="Sample_Name", split.by = "Status")
p3 <- VlnPlot(seu, features = "percent.mt", pt.size=0, group.by="Sample_Name", split.by = "Status")

savePlot((p1 / p2 / p3),
         DIR_RES,
         prfx = "",
         main = "",
         sfx = "",
         other = "")



## Visualise correlation between QC metrics
p4 <- FeatureScatter(seu, feature1 = "nCount_RNA", feature2 = "percent.mt")
p5 <- FeatureScatter(seu, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

## Subset data
## Remove cells with less than 1000 UMIs or more than 25% of reads arising from
## mitochondrial genes.
seu <- subset(seu, subset = nFeature_RNA > 1000 & percent.mt < 25)


# ---------------------------------------------------------------------------- #
# ----Normalisation ------------------------------------------------------------
## To avoid the memory error
options(future.globals.maxSize = 10000 * 1024^2)

seu <- SCTransform(seu,
                   variable.features.n = 3000,
                   cells = 5000,
                   do.scale = TRUE,
                   do.center = TRUE,
                   verbose = TRUE)


# ---------------------------------------------------------------------------- #
# ----Dimensionality Reduction -------------------------------------------------
## Run PCA
seu <- RunPCA(seu, assay = "SCT", npcs = 50)

## Run UMAP
seu <- RunUMAP(seu, dims = 1:20, reduction = "pca")

## View the explained variability
# Plots the standard deviations of the principle components for easy
# identification of an elbow in the graph. This elbow often corresponds well
# with the significant dims and is much faster to run than Jackstraw.
ElbowPlot(seu)


# ---------------------------------------------------------------------------- #
# ----Clustering ---------------------------------------------------------------
## Find Neighbours using reduced dimensions
# Use the first 20 PCs since the variability
seu <- FindNeighbors(object = seu, dims = 1:20, verbose = F)

## Find Clusters
seu <- FindClusters(object = seu, resolution = 0.01, verbose = F)


# ---------------------------------------------------------------------------- #
# ----Wrap up ------------------------------------------------------------------
saveRDS(seu, file = file.path(DIR_PROJECT, "analysis_objs/single_cell/seu_sc.rds"))

