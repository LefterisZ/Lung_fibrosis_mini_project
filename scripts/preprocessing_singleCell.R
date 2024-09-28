library(Seurat)
library(SeuratDisk)

library(readr)
library(tidyverse)

# ----Load metadata ------------------------------------------------------------
metadata <- read_csv("data/single_cell/IPF_metadata.csv") %>%
  column_to_rownames(var = "...1") %>%
  mutate(Barcode = rownames(.), .before = "orig.ident") %>%
  arrange(Barcode)

# ----Load data ----------------------------------------------------------------
dataDir <- "./data/single_cell/GSE135893/"
## Read counts matrix
seu <- Read10X(dataDir,
               gene.column = 1,
               cell.column = 1,
               unique.features = TRUE,
               strip.suffix = FALSE)
## Create Seurat object
seu <- CreateSeuratObject(counts = seu)

## Add metadata from the authors
seu <- AddMetaData(seu, metadata = metadata)

# ----QC -----------------------------------------------------------------------
## Identify percentage of mitochondrial reads in each cell
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")

## Visualize QC metrics as a violin plot
VlnPlot(seu, features = "nFeature_RNA", pt.size=0)
VlnPlot(seu, features = "nCount_RNA", pt.size=0)
VlnPlot(seu, features = "percent.mt", pt.size=0)

## Visualise correlation between QC metrics
FeatureScatter(seu, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(seu, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

## Subset data
## Remove cells with less than 1000 UMIs or more than 25% of reads arising from
## mitochondrial genes.
seu <- subset(seu,
              subset = nFeature_RNA > 1000 & percent.mt < 25)

# ----Normalisation ------------------------------------------------------------
## To avoid the memory error
options(future.globals.maxSize = 8000 * 1024^2)

seu <- SCTransform(seu,
                   variable.features.n = 3000,
                   cells = 5000,
                   do.scale = TRUE,
                   do.center = TRUE,
                   verbose = TRUE)

# ----Dimensionality Reduction -------------------------------------------------
# PCA
seu <- RunPCA(seu, assay = "SCT", npcs = 50)

# UMAP
seu <- RunUMAP(seu, dims = 1:20, reduction = "pca")

# View the explained variability
ElbowPlot(seu)

# ----Clustering ---------------------------------------------------------------
seu <- FindNeighbors(object = seu, dims = 1:20, verbose = F)
seu <- FindClusters(object = seu, resolution = 0.01, verbose = F)



