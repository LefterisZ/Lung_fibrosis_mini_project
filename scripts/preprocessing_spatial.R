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


# ---------------------------------------------------------------------------- #
# ----Custom functions to use --------------------------------------------------
# A custom function to save plots



# ---------------------------------------------------------------------------- #
# ----Set samples --------------------------------------------------------------
sampleNames <- list.files(DIR_DATA)
sampleNames <- sampleNames[sampleNames != "README.txt"]
names(sampleDir) <- sampleNames

## Update metadata to keep only the selected samples
## This is because the dataset includes 25 samples but we are working only with
## the 16 that have pathologist annotations provided by the authors.
metadata <- metadata[metadata$sample_id %in% sampleNames_selected,]

