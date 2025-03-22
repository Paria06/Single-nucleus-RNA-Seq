
# Load libraries
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

seurat_integrated <- readRDS("/lustre07/scratch/paria/als-project/analysis/ocu_med_sc/rpca/WGCNA/data/seurat_integrated_hdWGCNA.rds")
DefaultAssay(seurat_integrated) <- "RNA"


# Set up output directories
outDir <- "/lustre07/scratch/paria/als-project/analysis/ocu_med_sc/rpca/WGCNA/across_status/cluster_level"


regions <- c("ocu", "med", "sc", "sl")

#cell_types <- c("Microglia", "OPC", "splatter")

params <- list(
  "ocu" = list(k = 30, max_shared = 15, cell_type = "Mic1"),
  "med" = list(k = 30, max_shared = 15, cell_type = "OPC1"),
  "sc" = list(k = 25, max_shared = 10, cell_type = "ExInVentral1"),
  "sl" = list(k = 20, max_shared = 10, cell_type = "ExInVentral1"))


for (reg in regions) {

  print(paste("Processing region", reg))
  # Subset the Seurat object by region
  seurat_obj <- subset(seurat_integrated, subset = region == reg)
  DefaultAssay(seurat_obj) <- "RNA"
  
  # Extract the specific parameters for the current region
  k <- params[[reg]]$k
  max_shared <- params[[reg]]$max_shared
  cell_type <- params[[reg]]$cell_type

    print(paste("Processing Region:", reg, "-", "Cell Type:", cell_type))
  
    # Create output directory for this specific region and cell type
    clusterDir <- file.path(outDir, paste0(reg, "_", cell_type))
    dir.create(clusterDir, showWarnings = FALSE, recursive = TRUE)
    
 
    print("Setup For WGCNA")
    # we will select genes that are expressed in at least 5% of cells in this datase
    seurat_obj <- SetupForWGCNA(
    seurat_obj,
    gene_select = "fraction",
    fraction = 0.05,
    wgcna_name = paste0(reg, "_", cell_type))


    print("construct metacells  in each group")
    # construct metacells  in each group
    seurat_obj <- MetacellsByGroups(
    seurat_obj = seurat_obj,
    group.by = c("cluster_id", "sample_id"),
    reduction = 'pca',
    k = k,
    max_shared = max_shared,
    ident.group = "cluster_id",
    min_cells = 100)
    
    
    print("normalize metacell expression matrix")
    seurat_obj <- NormalizeMetacells(seurat_obj)


    print("Set up the expression matrix")
    # Set up the expression matrix
    seurat_obj <- SetDatExpr(
    seurat_obj,
    group_name = cell_type,
    group.by= 'cluster_id',
    assay = 'RNA',
    slot = 'data',
    wgcna_name = paste0(reg, "_", cell_type))


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
    tom_name = paste0(reg, "_", cell_type))

    # Importantly, the “grey” module consists of genes that were not grouped into any co-expression module. 
    # The grey module should be ignored for all downstream analysis and interpretation

    print("Optional: inspect the topoligcal overlap matrix (TOM)")
    TOM <- GetTOM(seurat_obj)


    print("Module Eigengenes and Connectivity")
    print("Compute harmonized module eigengenes")
    # need to run ScaleData first or else harmony throws an error:
    #DefaultAssay(seurat_obj) <- "integrated"
    #seurat_obj <- FindVariableFeatures(seurat_obj)
    #seurat_obj <- ScaleData(seurat_obj, features=VariableFeatures(seurat_obj))

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
    group.by = 'cluster_id', group_name = cell_type)



    print("save your output")
    # save your output
    rdsFilePath <- file.path(clusterDir, paste0("hdWGCNA_object_", reg, "_", cell_type, ".rds"))
    saveRDS(seurat_obj, file = rdsFilePath)

   }


sessionInfo()
