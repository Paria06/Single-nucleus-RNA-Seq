# The analysis of cell type proportion using Annova and ttest were done. I used the below tutorials
# https://github.com/phipsonlab/speckle/blob/main/vignettes/speckle.Rmd
# https://phipsonlab.github.io/propeller-paper-analysis/RealDataAnalysis.html#Heart_development_analysis
# https://rdrr.io/github/Oshlack/speckle/f/vignettes/speckle.Rmd

# Load necessary libraries for analysis and visualization
library(speckle)
library(limma)
library(edgeR)
library(pheatmap)
library(gt)
library(SingleCellExperiment)
library(CellBench)
library(limma)
library(ggplot2)
library(scater)
library(patchwork)
library(edgeR)
library(statmod)
library(Seurat)

# Read in the integrated Seurat object from the specified path
seurat_integrated <- readRDS("/home/paria/scratch/als-project/analysis/ocu_med_sc/rpca/wicoxon_test/speckle/cluster/seurat_integrated_labelled.rds")

# Update cluster information in the Seurat object
seurat_integrated$subcluster <- seurat_integrated$cluster_id
seurat_integrated$cluster_id <- seurat_integrated$main
seurat_integrated$main <- NULL

# Data exploration: Calculate transformed proportions of cell types
props <- getTransformedProps(clusters = seurat_integrated$cluster_id, sample = seurat_integrated$sample_id)

# Create a PDF to visualize broad cell type proportions
pdf("1_Broad_cell_type_proportions.pdf", height = 30, width = 40)
p2 <- plotCellTypeProps(clusters = seurat_integrated$cluster_id, sample = seurat_integrated$sample_id)
p2 + theme(axis.text.x = element_text(angle = 45)) + ggtitle("Broad cell type proportions")
dev.off()

# Calculate baseline proportions of cell types across samples
counts <- table(seurat_integrated$cluster_id, seurat_integrated$sample_id)
baselineN <- rowSums(counts)
N <- sum(baselineN)
baselineprops <- baselineN / N

# Factorize the final cell type variable based on baseline proportions
seurat_integrated$final_ct <- factor(seurat_integrated$cluster_id, levels = names(sort(baselineprops, decreasing = TRUE)))

# Get transformed proportions for the final cell types
props <- getTransformedProps(clusters = seurat_integrated$final_ct, sample = seurat_integrated$sample_id)
cols <- ggplotColors(nrow(props$Proportions))
m <- match(rownames(props$Proportions), levels(factor(seurat_integrated$cluster_id)))

# Create a PDF to visualize cell type proportions estimates
pdf("2_broad_Cell_type_proportions_estimates.pdf", height = 6, width = 10)
par(mfrow = c(1, 1))
par(mar = c(7, 5, 2, 2))
plot(jitter(props$Proportions[, 1]), col = cols[m], pch = 16, ylim = c(0, max(props$Proportions)),
     xaxt = "n", xlab = "", ylab = "Cell type proportion", cex.lab = 1.5, cex.axis = 1.5)
for (i in 2:ncol(props$Proportions)) {
  points(jitter(1:nrow(props$Proportions)), props$Proportions[, i], col = cols[m], pch = 18)
}
axis(side = 1, at = 1:nrow(props$Proportions), las = 2, labels = rownames(props$Proportions))
title("Cell type proportions estimates for 37 tissues")
dev.off()

# Create PDFs to visualize mean-variance relationships
pdf("3_broad_Mean_variance_relationship.pdf", height = 6, width = 6)
plotCellTypeMeanVar(counts)  # Plot mean-variance relationship for cell types
dev.off()

pdf("4_broad_Mean_variance_relationship.pdf", height = 6, width = 6)
plotCellTypePropsMeanVar(counts)  # Plot mean-variance relationship for cell type proportions
dev.off()

# Convert Seurat object to SingleCellExperiment object for further analysis
counts <- seurat_integrated@assays$RNA@counts
metadata <- seurat_integrated@meta.data
metadata$cluster_id <- factor(seurat_integrated$cluster_id)
sce <- SingleCellExperiment(assays = list(counts = counts), colData = metadata)

# Define variables as factors for analysis
factor_vars <- c("sample_id", "cluster_id", "condition", "chemistry", "region", "group_id", "batchlib", "individual", "sex")
sce[factor_vars] <- lapply(sce[factor_vars], as.factor)  # Convert specified columns to factors
sce$PMI <- scale(as.numeric(sce$PMI))  # Scale PMI variable

# Load experimental information from a CSV file
ei <- read.csv("ei.csv")

# Define design matrix for ANOVA
design <- model.matrix(~ 0 + condition + batch + individual + sex, data = ei)

# Perform ANOVA comparisons for specified groups
anova_comparisons <- list(
  all_groups = c(1, 2, 3, 4, 5, 6, 7, 8),
  control_regions = c(2, 4, 6, 8),
  als_regions = c(1, 3, 5, 7)
)

# Loop through each comparison and perform ANOVA
for (comp in names(anova_comparisons)) {
  cat("Performing ANOVA for", comp, "\n")
  anova_result <- propeller.anova(prop.list = props, design = design, coef = anova_comparisons[[comp]], robust = TRUE, trend = FALSE, sort = TRUE)
  write.csv(anova_result, file = paste0("Broad_ANOVA_", comp, ".csv"))
}

# Define t-test comparisons
ttest_comparisons <- list(
  als_vs_control_medulla = "conditionalsmed - conditionconmed",
  als_vs_control_oculomotor = "conditionalsocu - conditionconocu",
  als_vs_control_spinal_cord_cervical = "conditionalssc - conditionconsc",
  als_vs_control_spinal_cord_lumbar = "conditionalssl - conditionconsl",
  als_vs_control = "group_idals - group_idcon",
  med_vs_ocu = "regionmed - regionocu",
  sc_vs_ocu = "regionsc - regionocu",
  sl_vs_ocu = "regionsl - regionocu",
  med_vs_ocu_con = "conditionconmed - conditionconocu",
  sc_vs_ocu_con = "conditionconsc - conditionconocu",
  sl_vs_ocu_con = "conditionconsl - conditionconocu",
  med_vs_ocu_als = "conditionalsmed - conditionalsocu",
  sc_vs_ocu_als = "conditionalssc - conditionalsocu",
  sl_vs_ocu_als = "conditionalssl - conditionalsocu"
)

# Loop through each comparison and perform t-tests
for (comp in names(ttest_comparisons)) {
  cat("Performing t-test for", comp, "\n")
  contrast <- makeContrasts(contrasts = ttest_comparisons[[comp]], levels = design)
  ttest_result <- propeller.ttest(props, design, contrasts = contrast, robust = TRUE, trend = FALSE, sort = TRUE)
  write.csv(ttest_result, file = paste0("Broad_ttest_", comp, ".csv"))
}

sessionInfo()
