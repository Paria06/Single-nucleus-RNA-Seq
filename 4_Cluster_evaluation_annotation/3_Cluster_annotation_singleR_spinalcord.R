# Load required libraries
library(SingleR)
library(Seurat)
library(future)
library(dplyr)
library(Matrix)
library(monocle)
library(pheatmap)
library(gplots)
library(ggplot2)
library(cowplot)
library(SingleCellExperiment)
library(EnsDb.Hsapiens.v75)
library(EnsDb.Hsapiens.v79)
library(EnsDb.Hsapiens.v86)
library(org.Hs.eg.db)
library(tidyverse)

# Read in the test dataset
ALS <- readRDS("seurat_integrated_DB_removed.rds")
DefaultAssay(ALS) <- "RNA"
als <- as.SingleCellExperiment(ALS)

# Read in the reference dataset
counts <- read.csv("counts.csv")
rownames(counts) <- counts[, 1]
counts <- counts[, -1]
metadata <- read.csv("metadata.csv")
sc <- SingleCellExperiment(assays = list(counts = counts), colData = metadata)
sc <- logNormCounts(sc)

# Run SingleR function
new_label <- SingleR(test = als, ref = sc, labels = sc$subtype_annotation, de.method = "wilcox")

# Save the results
write.csv(x = new_label, file = "label_MD.csv")

# Visualization
ALS <- readRDS("seurat_integrated_DB_removed.rds")
DefaultAssay(ALS) <- "integrated"

# Calculate pairwise Rand index
pairwiseRand(ALS$seurat_clusters, new_label$labels, mode = "index")

# Create a contingency table
tab <- table(cluster = ALS$seurat_clusters, label = new_label$labels)

# Generate heatmap
pdf('pheatmap_sc_MD_Mode.pdf', height = 30, width = 30)
pheatmap::pheatmap(log10(tab + 10))
dev.off()

# Add cell type labels to the ALS object
ALS$cell.type.label <- new_label$labels

# Create UMAP plot
pdf("DimPlot_spinalcord_MD_Mode.pdf", height = 20, width = 40)
DimPlot(ALS, reduction = "umap", group.by = "cell.type.label", label = TRUE, pt.size = 6, label.size = 5)
dev.off()

sessionInfo()