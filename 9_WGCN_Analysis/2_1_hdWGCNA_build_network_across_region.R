# load librearies
library(WGCNA)
library(DESeq2)
library(GEOquery)
library(tidyverse)
library(CorLevelPlot)
library(gridExtra)
library(hdWGCNA)
library(Seurat)
library(tidyverse)
library(cowplot)
library(patchwork)
library(ggraph)
library(fgsea)
library(igraph)
library(farver)
library(ggforestplot)
library(RhpcBLASctl)
library(GeneOverlap)
library(corrplot)
omp_set_num_threads(1)
blas_set_num_threads(1)


# using the cowplot theme for ggplot
theme_set(theme_cowplot())

# set random seed for reproducibility
set.seed(12345)

# optionally enable multithreading
enableWGCNAThreads(nThreads = 8)

seurat_integrated <- readRDS("/lustre07/scratch/paria/als-project/analysis/ocu_med_sc/rpca/WGCNA/seurat_integrated_hdWGCNA.rds")
DefaultAssay(seurat_integrated) <- "RNA"


# Set up output directories
outDir <- "/lustre07/scratch/paria/als-project/analysis/ocu_med_sc/rpca/WGCNA_2/"


# Define the regions and cell types you want to analyze
#cell_types <- c("splatter", "Astrocyte", "Oligodendrocyte") 
cell_types <- c("MN")
params <- list(
  Microglia = list(k = 15, max_shared = 7),
  splatter = list(k = 15, max_shared = 7),
  Astrocyte = list(k = 15, max_shared = 7),
  Oligodendrocyte = list(k = 50, max_shared = 22),
  MN = list(k = 15, max_shared = 7)
)


for (cluster in cell_types) {
  # Subset the Seurat object by region
  seurat_obj <- subset(seurat_integrated, subset = main %in% cluster)
  seurat_obj$main <- as.factor(seurat_obj$main)
  seurat_obj$sample_id <- as.factor(seurat_obj$sample_id)
  Idents(seurat_obj) <- "main"
    # Create output directory for this specific region and cell type
  clusterDir <- file.path(outDir, paste0(cluster))
  dir.create(clusterDir, showWarnings = FALSE, recursive = TRUE)
    
  print(paste("Cell Type:", cluster))
  print("Setup For WGCNA")
  # we will select genes that are expressed in at least 5% of cells in this datase
  seurat_obj <- SetupForWGCNA(
  seurat_obj,
  gene_select = "fraction",
  fraction = 0.05,
  wgcna_name = cluster)

  # Extract the specific parameters for the current region
  k <- params[[cluster]]$k
  max_shared <- params[[cluster]]$max_shared
  
  print("construct metacells  in each group")
  # construct metacells  in each group
  seurat_obj <- MetacellsByGroups(
  seurat_obj = seurat_obj,
  group.by = c("main", "sample_id"),
  reduction = 'pca',
  k = k,
  max_shared = max_shared,
  ident.group = "main",
  min_cells = 100)


  print("normalize metacell expression matrix")
  seurat_obj <- NormalizeMetacells(seurat_obj)


  print("Set up the expression matrix")
  # Set up the expression matrix
  seurat_obj <- SetDatExpr(
  seurat_obj,
  group_name = cluster,
  group.by= 'main',
  assay = 'RNA',
  slot = 'data',
  wgcna_name = cluster)


  print("Select soft-power threshold")
  # Test different soft powers:
  seurat_obj <- TestSoftPowers(
  seurat_obj,
  networkType = "signed")


  power_table <- GetPowerTable(seurat_obj)
  head(power_table)


  print("construct co-expression network")
  # construct co-expression network:
  seurat_obj <- ConstructNetwork(
  seurat_obj,
  tom_name = cluster,
  overwrite_tom = TRUE)

  # Importantly, the “grey” module consists of genes that were not grouped into any co-expression module. 
  print("Optional: inspect the topoligcal overlap matrix (TOM)")
  TOM <- GetTOM(seurat_obj)


  print("ModuleEigengenes")
  # compute all MEs in the full single-cell dataset
  seurat_obj <- ModuleEigengenes(seurat_obj, vars.to.regress=c("PMI", "batchlib", "sex"))


  print("harmonized module eigengenes")
  # harmonized module eigengenes:
  hMEs <- GetMEs(seurat_obj)

  print("module eigengenes")
  # module eigengenes:
  MEs <- GetMEs(seurat_obj, harmonized=FALSE)

  print("Compute module connectivity")
  # Compute module connectivity
  # compute eigengene-based connectivity (kME):
  seurat_obj <- ModuleConnectivity(
  seurat_obj,
  group.by = 'main', group_name = cluster)


  print("Getting the module assignment table")
  # Getting the module assignment table
  # get the module assignment table:
  modules <- GetModules(seurat_obj) %>% subset(module != 'grey')

  # show the first 6 columns:
  head(modules[,1:6])

  print("save your output")
  # save your output
  rdsFilePath <- file.path(clusterDir, paste0("hdWGCNA_object_", cluster, ".rds"))
  saveRDS(seurat_obj, file = rdsFilePath)

}

sessionInfo()