# Load necessary libraries
library(Seurat)
library(future)
library(scDblFinder)
library(scater)
library(scran)
library(BiocSingular)

# Read the integrated Seurat object from an RDS file
seurat_integrated <- readRDS("full_sample_k_100_resolution_0.85_PC_30.rds")

# Set the default assay to "RNA"
DefaultAssay(seurat_integrated) <- "RNA"

# Convert Seurat object to SingleCellExperiment (SCE) object
als <- as.SingleCellExperiment(seurat_integrated)

# Transfer PCA, UMAP, and t-SNE embeddings from Seurat to SCE object
# Note: Reduced dimensions do not automatically transfer during conversion
reducedDim(als, "PCA", withDimnames=TRUE) <- seurat_integrated[['pca']]@cell.embeddings
reducedDim(als, "UMAP", withDimnames=TRUE) <- seurat_integrated[['umap']]@cell.embeddings
reducedDim(als, "TSNE", withDimnames=TRUE) <- seurat_integrated[['TSNE']]@cell.embeddings

# Set the default assay to "integrated" for further analysis
DefaultAssay(seurat_integrated) <- "integrated"

# Run doublet detection using scDblFinder
# Specify the sample identifier and use multicore processing for efficiency
sce <- scDblFinder(als, samples="orig.ident", BPPARAM=MulticoreParam(37))

# Display the counts of detected doublets and singlets
table(sce$scDblFinder.class)

# Add the doublet classification results back to the Seurat object
seurat_integrated$scDblFinder.class <- sce$scDblFinder.class 

# Subset the Seurat object to retain only singlet cells
seurat_integrated <- subset(seurat_integrated, subset = scDblFinder.class == "singlet")

# Save the filtered Seurat object to an RDS file
saveRDS(seurat_integrated, file = "seurat_integrated_DB_removed.rds")

# Display session information for reproducibility
sessionInfo()
