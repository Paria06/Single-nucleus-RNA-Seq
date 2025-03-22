# Load necessary libraries
library(Seurat)
library(future)
library(dplyr)
library(Matrix)
library(monocle)
library(pheatmap)
library(gplots)
library(ggplot2)
library(cowplot)
library(reshape)
library(RColorBrewer)
library(tidyverse)
library(scales)
library(RCurl)
library(patchwork)
library(metap)

# Define color palette
col <- colorRampPalette(rev(brewer.pal(11, "Spectral")))(10)[-c(5, 6)]

# Set options for future processing
options(future.globals.maxSize = 5000 * 1024^2)

# Load the filtered Seurat object
filtered_seurat <- readRDS("filtered_seurat.rds")

# Update future options for memory size
options(future.globals.maxSize = 4000 * 1024^2)

# Split Seurat object by condition for cell cycle scoring and SCT
split_seurat <- SplitObject(filtered_seurat, split.by = "orig.ident")

# Normalize data and find variable features for each split
split_seurat <- lapply(X = split_seurat, FUN = function(x) {
  x <- NormalizeData(x, verbose = FALSE)
  x <- FindVariableFeatures(x, verbose = FALSE)
  return(x)
})

# Data integration
integ_features <- SelectIntegrationFeatures(object.list = split_seurat, nfeatures = 3000)

# Scale data and perform PCA for each split
split_seurat <- lapply(X = split_seurat, FUN = function(x) {
  x <- ScaleData(x, features = integ_features, verbose = FALSE)
  x <- RunPCA(x, features = integ_features, npcs = 80, verbose = FALSE)
  return(x)
})

# Find integration anchors
integ_anchors <- FindIntegrationAnchors(
  object.list = split_seurat, 
  anchor.features = integ_features, 
  reference = c(28, 30, 13, 19, 20, 11, 3, 4), 
  reduction = "rpca", 
  dims = 1:80, 
  k.anchor = 20
)

# Integrate data
seurat_integrated <- IntegrateData(anchorset = integ_anchors, dims = 1:80)

# Scale integrated data and run PCA
seurat_integrated <- ScaleData(seurat_integrated, verbose = FALSE)
seurat_integrated <- RunPCA(seurat_integrated, verbose = FALSE, npcs = 80)

# Run UMAP for dimensionality reduction
seurat_integrated <- RunUMAP(seurat_integrated, dims = 1:80)

# Save the integrated Seurat object
saveRDS(seurat_integrated, "ocu_sc_med_rpca.rds")

sessionInfo()
