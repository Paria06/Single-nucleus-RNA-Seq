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
library(SingleCellExperiment)
library(tidyverse)
library(scales)
library(RCurl)
library(patchwork)
library(metap)

# Set color palette
col <- colorRampPalette(rev(brewer.pal(11, "Spectral")))(10)[-c(5, 6)]
options(future.globals.maxSize = 5000 * 1024^2)

# Define sample names
samples <- c("T01496", "T01493", "T01494", "T01495", "T01471", "T01472", 
             "T01473", "T01474", "T01497", "T01498", "T01489", "T01480", 
             "T01477", "T01478", "T01513", "T01514", "T01499", "T01500", 
             "T01486", "T01487", "T01488", "T01490", "T01491", "T01492", 
             "T01479", "T01515", "T01516", "T01551", "T01552", "T01553", 
             "T01554", "T01555", "T01556", "T01557", "T01558", "T01559", 
             "T01560")

# Create Seurat objects and add mitochondrial percentage
seurat_objects <- lapply(samples, function(sample) {
    als.data <- Read10X(data.dir = paste0("/home/paria/scratch/als-project/filtered_feature_bc_matrix/", sample))
    seurat_obj <- CreateSeuratObject(als.data, min.features = 100, project = sample)
    seurat_obj <- AddMetaData(seurat_obj, metadata = PercentageFeatureSet(seurat_obj, pattern = "^MT-"), col.name = "percent.chrM")
    return(seurat_obj)
})

# Merge all Seurat objects
merged_seurat <- Reduce(function(x, y) merge(x, y), seurat_objects)

# Remove mitochondrial genes
MTgenes <- "^MT-"
merged_seurat <- subset(merged_seurat, features = setdiff(rownames(merged_seurat), grep(MTgenes, rownames(merged_seurat), value = TRUE)))

# Recalculate percentage of mitochondrial genes
merged_seurat[["percent.chrM"]] <- PercentageFeatureSet(merged_seurat, pattern = "^MT-")

# Add metadata
metadata <- merged_seurat@meta.data
metadata$percent.chrM <- merged_seurat$percent.chrM

# Define status
status_mapping <- list(
    als = c("T01499", "T01500", "T01487", "T01488", "T01491", "T01492", "T01479", "T01515", "T01516", "T01497", "T01489", "T01480", "T01498", "T01486", "T01490", "T01477", "T01478", "T01513", "T01514"),
    con = c("T01471", "T01473", "T01551", "T01552", "T01553", "T01554", "T01555", "T01556", "T01557", "T01558", "T01559", "T01560", "T01493", "T01494", "T01495", "T01472", "T01474", "T01496")
)

for (status in names(status_mapping)) {
    metadata[metadata$orig.ident %in% status_mapping[[status]], "status"] <- status
}

# Define individual mappings
individual_mapping <- list(
    ind1454 = c("T01493", "T01494"),
    ind1888 = c("T01495", "T01496"),
    ind1895 = c("T01471", "T01472"),
    ind1990 = c("T01473", "T01474"),
    ind2158 = c("T01477", "T01480", "T01478", "T01479"),
    ind2177 = c("T01499", "T01500", "T01497", "T01498"),
    ind2166 = c("T01488", "T01487", "T01486"),
    ind2172 = c("T01491", "T01492", "T01489", "T01490"),
    ind2187 = c("T01515", "T01516", "T01513", "T01514"),
    A158_14 = c("T01551", "T01552"),
    A309_99 = c("T01553", "T01554"),
    A127_11 = c("T01555", "T01556"),
    A200_13 = c("T01557", "T01558"),
    A346_10 = c("T01559", "T01560")
)

for (individual in names(individual_mapping)) {
    metadata[metadata$orig.ident %in% individual_mapping[[individual]], "individual"] <- individual
}

# Define sex mappings
sex_mapping <- list(
    male = c("T01479", "T01555", "T01556", "T01473", "T01477", "T01480", "T01474", "T01478"),
    female = c("T01493", "T01495", "T01497", "T01489", "T01513", "T01499", "T01500", "T01487", "T01488", "T01491", "T01492", "T01515", "T01516", "T01551", "T01552", "T01553", "T01554", "T01557", "T01558", "T01559", "T01560", "T01471", "T01472", "T01486", "T01490", "T01494", "T01498", "T01514", "T01496")
)

for (sex in names(sex_mapping)) {
    metadata[metadata$orig.ident %in% sex_mapping[[sex]], "sex"] <- sex
}

# Define batchlib mappings
batchlib_mapping <- list(
    batch1 = c("T01500", "T01479", "T01473", "T01489", "T01493", "T01472", "T01486", "T01496"),
    batch2 = c("T01491", "T01497"),
    batch3 = c("T01487", "T01477", "T01494"),
    batch4 = c("T01492", "T01495", "T01498"),
    batch5 = c("T01499", "T01474", "T01478", "T01490"),
    batch6 = c("T01471", "T01515", "T01516", "T01488", "T01513", "T01514"),
    batch7 = c("T01551", "T01552", "T01553", "T01554", "T01555", "T01556", "T01557", "T01558", "T01559", "T01560", "T01480")
)

for (batch in names(batchlib_mapping)) {
    metadata[metadata$orig.ident %in% batchlib_mapping[[batch]], "batchlib"] <- gsub("batch", "", batch)
}

# Define region mappings
region_mapping <- list(
    sc = c("T01499", "T01487", "T01491", "T01479", "T01515", "T01552", "T01554", "T01556", 
           "T01558", "T01560", "T01500", "T01488", "T01492", "T01516", "T01551", "T01553", 
           "T01555", "T01557", "T01559", "T01480"),
    ocu = c("T01493", "T01495", "T01471", "T01473", "T01497", "T01489", "T01477", "T01513"),
    med = c("T01494", "T01496", "T01472", "T01474", "T01498", "T01486", "T01490", "T01478", "T01514")
)

for (region in names(region_mapping)) {
    metadata[metadata$orig.ident %in% region_mapping[[region]], "region"] <- region
}

# Define group_id mappings
group_id_mapping <- list(
    alsocu = c("T01497", "T01489", "T01477", "T01513"),
    conocu = c("T01493", "T01495", "T01471", "T01473"),
    alssc = c("T01499", "T01500", "T01487", "T01488", "T01491", "T01492", "T01479", "T01480", "T01515", "T01516"),
    consc = c("T01551", "T01552", "T01553", "T01554", "T01555", "T01556", "T01557", "T01558", "T01559", "T01560"),
    conmed = c("T01494", "T01496", "T01472", "T01474"),
    alsmed = c("T01498", "T01486", "T01490", "T01478", "T01514")
)

for (group in names(group_id_mapping)) {
    metadata[metadata$orig.ident %in% group_id_mapping[[group]], "group_id"] <- group
}

# Define PMI mappings
pmi_mapping <- list(
    "23.6" = c("T01493", "T01494"),
    "37.42" = c("T01495", "T01496"),
    "21.07" = c("T01471", "T01472"),
    "7.48" = c("T01473", "T01474"),
    "24.28" = c("T01497", "T01498", "T01499", "T01500"),
    "25.27" = c("T01486", "T01487", "T01488"),
    "25.22" = c("T01489", "T01490", "T01491", "T01492"),
    "24.17" = c("T01477", "T01478", "T01479", "T01480"),
    "35.4" = c("T01513", "T01514", "T01515", "T01516"),
    "27" = c("T01551", "T01552"),
    "21" = c("T01553", "T01554"),
    "23" = c("T01555", "T01556"),
    "30" = c("T01557", "T01558"),
    "34" = c("T01559", "T01560")
)

for (pmi in names(pmi_mapping)) {
    metadata[metadata$orig.ident %in% pmi_mapping[[pmi]], "PMI"] <- pmi
}

# Add metadata back to Seurat object
merged_seurat@meta.data <- metadata

# Filter the Seurat object
filtered_seurat <- subset(merged_seurat, subset = (nFeature_RNA >= 350) & (percent.chrM < 10) & (nCount_RNA < 126000))

# Gene-level filtering
counts <- GetAssayData(filtered_seurat, slot = "counts")
keep_genes <- Matrix::rowSums(counts > 1) >= 10
filtered_counts <- counts[keep_genes, ]

# Create a new filtered Seurat object
filtered_seurat <- CreateSeuratObject(filtered_counts, metadata = filtered_seurat@meta.data)

# Save the filtered Seurat object
saveRDS(filtered_seurat, file = "filtered_seurat.rds")

sessionInfo()