# Loading Libraries 
library(Seurat)
library(future)
library(dplyr)
library(scclusteval)
library(purrr)

# Load Data 
fullsample_idents <- readRDS("gather_full_sample.rds")
subsample_idents <- readRDS("gather_subsample.rds")
als <- readRDS("seurat_integrated.rds")

## Explore Full Dataset 
# Display the full sample identities
print(fullsample_idents)

# Count the number of clusters for each combination of parameter set
fullsample_idents <- fullsample_idents %>%
  mutate(cluster_num = map_dbl(original_ident_full, ~ length(unique(.x))))

# Explore Subsampled Data 
# Assign stable clusters for all combinations of parameters
stable_clusters <- subsample_idents_list %>%
  mutate(stable_cluster = map(data, ~ AssignStableCluster(
    .x$original_ident,
    .x$recluster_ident,
    jaccard_cutoff = 0.8,
    method = "jaccard_percent", 
    percent_cutoff = 0.8
  )))

# Generate Plots 
# Create a scatter plot for parameter sets
pdf("ParameterSetScatterPlot.pdf", height = 20, width = 20)
ParameterSetScatterPlot(
  stable_clusters = stable_clusters,
  fullsample_idents = fullsample_idents,
  x_var = "k_param",
  y_var = "number",
  facet_rows = "resolution",
  facet_cols = "pc"
) +
  ggtitle("Parameter Set Scatter Plot")
dev.off()

# Create a plot for the percentage of cells in stable clusters
pdf("Percentage_of_Cells_in_Stable_Clusters.pdf", height = 20, width = 20)
ParameterSetScatterPlot(
  stable_clusters = stable_clusters,
  fullsample_idents = fullsample_idents,
  x_var = "k_param",
  y_var = "percentage",
  facet_rows = "resolution",
  facet_cols = "pc"
) +
  ggtitle("Percentage of Cells in Stable Clusters")
dev.off()

sessionInfo()
