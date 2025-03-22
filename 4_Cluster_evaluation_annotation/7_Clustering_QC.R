library(Seurat)
library(future)
library(dplyr)
library(scclusteval)
library(purrr)
library(dittoSeq)
library(ggplot2)

# Load Data 
seurat_integrated <- readRDS("seurat_integrated_DB_removed.rds")

#Determine metrics to plot present in seurat_integrated@meta.data
metrics <- c("nCount_RNA", "nFeature_RNA")
pdf("4_feature_plot_4metrics.pdf", height = 10, width = 12)
plot <- FeaturePlot(
  seurat_integrated,
  reduction = "umap",
  features = metrics,
  pt.size = 0.5,
  order = TRUE,
  min.cutoff = 'q10',
  label = TRUE,
  ncol = 2  # Arranges the plots in 2 columns
)
print(plot)  # Ensures the plot is rendered in the PDF
dev.off()


# Generate the plot: nuclei proportion of different variables in different level of clustering

# Define the variables and grouping levels
variables <- c("group_id", "region", "batch")
group_levels <- c("main", "cluster_id")

# Iterate over each variable and grouping level
for (var in variables) {
  for (group in group_levels) {
    # Generate the dittoBarPlot
    plot <- dittoBarPlot(seurat_integrated, var, group.by = group)
    
    # Customize the plot's appearance
    plot <- plot + theme(
      axis.text = element_text(size = 6),   # Font size of axis labels
      axis.title = element_text(size = 6),  # Font size of axis titles
      plot.title = element_text(size = 8)   # Font size of plot title
    )
    
    # Define the filename based on the variable and grouping level
    filename <- paste0("dittoBarPlot_", var, "_", group, ".pdf")
    
    # Save the plot to a PDF file
    pdf(filename, height = 4, width = 6)
    print(plot)
    dev.off()
  }
}

sessionInfo()