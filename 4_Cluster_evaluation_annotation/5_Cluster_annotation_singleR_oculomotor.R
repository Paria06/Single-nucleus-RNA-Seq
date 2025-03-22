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
library(EnsDb.Hsapiens.v79)
library(org.Hs.eg.db)
library(tidyverse)
library(biomaRt)

# Read in the test dataset
ALS <- readRDS("seurat_integrated_DB_removed.rds")
DefaultAssay(ALS) <- "RNA"
als <- as.SingleCellExperiment(ALS)

# Read in the reference dataset
reference <- readRDS("/home/paria/scratch/als-project/analysis/OCU_MED/ocu-med3/singleR/reference/Midbrain_M_Periaqueductal_gray_and_nearby_nuclei_PAG.rds")

# Create a vector mapping Ensembl IDs to SYMBOLs
ens <- mapIds(EnsDb.Hsapiens.v79, 
              keys = reference@assays$RNA@data@Dimnames[[1]], 
              column = 'SYMBOL', 
              keytype = 'GENEID')

# Remove NA values from the ens vector
keep <- !is.na(ens)
ens <- ens[keep]

# Update the gene names in the reference object
reference <- reference[keep, ]
reference@assays$RNA@data@Dimnames[[1]] <- ens
names(reference@assays$RNA@data@Dimnames[[1]]) <- NULL  # Clear names

# Convert the reference object to SingleCellExperiment format
ref <- as.SingleCellExperiment(reference)

# Run SingleR to label the test dataset
new_label <- SingleR(test = als, 
                     ref = ref, 
                     assay.type.test = 1, 
                     labels = reference$supercluster_term)

# Save the results to a CSV file
write.csv(x = new_label, file = "label_supercluster.csv")

# Visualization
ALS <- readRDS("seurat_integrated_DB_removed.rds")
DefaultAssay(ALS) <- "integrated"

# Calculate pairwise Rand index
pairwiseRand(ALS$seurat_clusters, new_label$labels, mode = "index")

# Create a contingency table
tab <- table(cluster = ALS$seurat_clusters, label = new_label$labels)

# Generate heatmap
pdf('pheatmap_PAG_supercluster.pdf', height = 30, width = 30)
pheatmap::pheatmap(log10(tab + 10))
dev.off()

# Add cell type labels to the ALS object
ALS$cell.type.label <- new_label$labels

# Create UMAP plot
pdf("DimPlot_PAG_supercluster.pdf", height = 20, width = 40)
DimPlot(ALS, reduction = "umap", group.by = "cell.type.label", label = TRUE, pt.size = 6, label.size = 5)
dev.off()

sessionInfo()