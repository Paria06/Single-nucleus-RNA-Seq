# Load necessary libraries
library(Seurat)
library(presto)

# Load the Seurat object
seurat_integrated <- readRDS("seurat_integrated_DB_removed.rds")

# Set the default assay to RNA
DefaultAssay(seurat_integrated) <- "RNA"

# Normalize RNA data for visualization purposes
seurat_integrated <- NormalizeData(seurat_integrated, verbose = FALSE)

# Perform differential expression analysis using presto's wilcoxauc function
presto_results <- wilcoxauc(seurat_integrated, 'cluster_id', seurat_assay = 'RNA')

# Filter markers based on specified criteria
filtered_markers <- subset(presto_results, logFC > log(2) & padj < 0.05 & (pct_in - pct_out) > 10)

# Extract top 100 markers for each cluster
top_markers <- top_markers(filtered_markers, n = 100)

# Save the top markers to a CSV file
write.csv(top_markers, file = "top100_markers_presto.csv", row.names = FALSE)

sessionInfo()