# single-cell analysis package
library(Seurat)

# plotting and data science packages
library(tidyverse)
library(cowplot)
library(patchwork)
library(MetBrewer)
# co-expression network analysis packages:
library(WGCNA)
library(hdWGCNA)

# gene enrichment packages
library(enrichR)
library(GeneOverlap)

# using the cowplot theme for ggplot
theme_set(theme_cowplot())

# set random seed for reproducibility
set.seed(12345)



# Set up output directories
outDir <- "/lustre07/scratch/paria/als-project/analysis/ocu_med_sc/rpca/WGCNA/across_status/cluster_level"


base_dirs <- c("ocu_Mic1" = "/ocu_Mic1/",
               "med_OPC1" = "/med_OPC1/", 
               "sc_ExInVentral1" = "/sc_ExInVentral1/",
               "sl_ExInVentral1" = "/sl_ExInVentral1/")

color_schemes <- list(
  "ocu_Mic1" = function(n) met.brewer("Hokusai3", n = n, type = 'continuous'),
  "med_OPC1" = function(n) met.brewer("Homer2", n = n, type = 'continuous'),
  "sc_ExInVentral1" = function(n) met.brewer("Benedictus", n = n, type = 'continuous'),
  "sl_ExInVentral1" = function(n) met.brewer("Redon", n = n, type = 'discrete'))



for (data_dir in names(base_dirs)) {
  #data_dir <- names(base_dirs)  
  # Load the Seurat object for the current dataset
  directory <- base_dirs[[data_dir]]
  seurat_obj <- readRDS(file = paste0(outDir, directory, "hdWGCNA_object_", data_dir, ".rds"))  
  
  # Set up an output directory for results
  output_dir <- file.path(paste0(outDir, directory, "results"))
  
  
  # reset module name and colors
  ### reset module name and colors
  seurat_obj <- ResetModuleNames(seurat_obj, new_name = paste0(data_dir, "-", "M"))

  modules <- GetModules(seurat_obj)
  mods <- levels(modules$module)
  mod_colors <- dplyr::select(modules, c(module, color)) %>%
    distinct %>% arrange(module) %>% .$color
  
  #n_colors <- length(mod_colors) -1
  # Define module names
  module_names <- mods[-1]  # Excluding the 'grey' module if present
  n_colors <- length(module_names)
  
  if (data_dir %in% names(color_schemes)) {
    new_colors <- color_schemes[[data_dir]](n_colors)
  } else {
    
  }
  new_colors <- sample(new_colors)
  # Map new colors to module names (excluding 'grey' if needed)
  module_color_map <- setNames(new_colors, module_names)
  seurat_obj <- ResetModuleColors(seurat_obj, new_colors)


  #dbs <- c('GO_Biological_Process_2021','GO_Cellular_Component_2021','GO_Molecular_Function_2021', 'REACTOME_2021', 'KEGG_2021')
  dbs <-c('GO_Biological_Process_2021','GO_Cellular_Component_2021','GO_Molecular_Function_2021', 'WikiPathway_2021_Human', 'KEGG_2021_Human')  
  # perform enrichment tests
  seurat_obj <- RunEnrichr(
    seurat_obj,
    dbs=dbs, # character vector of enrichr databases to test
    max_genes = 100 # number of genes per module to test. use max_genes = Inf to choose all genes!
  )

  # retrieve the output table
  enrich_df <- GetEnrichrTable(seurat_obj)
  write.csv(enrich_df, file = file.path(output_dir, "enrich_df.csv"))
  
  # make GO term plots:
 
  EnrichrBarPlot(
    seurat_obj,
    outdir = paste0(output_dir, "/enrichr_plots_", data_dir), # name of output directory
    n_terms = 10, # number of enriched terms to show (sometimes more show if there are ties!!!)
    plot_size = c(4,6), # width, height of the output .pdfs
    logscale=TRUE # do you want to show the enrichment as a log scale?
  )
 
  
  pdf(file.path(output_dir, "14.1_EnrichrDotPlot_GO_BP.pdf"), height=10, width=10)
  # enrichr dotplot
  EnrichrDotPlot(
    seurat_obj,
    mods = "all", # use all modules (this is the default behavior)
    database = "GO_Biological_Process_2021", # this has to be one of the lists we used above!!!
    n_terms=1 # number of terms for each module
  )
  dev.off()
  
  pdf(file.path(output_dir, "14.2_EnrichrDotPlot_GO_CC.pdf"), height=10, width=10)
  # enrichr dotplot
  EnrichrDotPlot(
    seurat_obj,
    mods = "all", # use all modules (this is the default behavior)
    database = "GO_Cellular_Component_2021", # this has to be one of the lists we used above!!!
    n_terms=1 # number of terms for each module
  )
  dev.off()


  pdf(file.path(output_dir, "14.3_EnrichrDotPlot_GO_MF.pdf"), height=10, width=10)
  # enrichr dotplot
  EnrichrDotPlot( 
    seurat_obj,
    mods = "all", # use all modules (this is the default behavior)
    database = "GO_Molecular_Function_2021", # this has to be one of the lists we used above!!!
    n_terms=1 # number of terms for each module
  )
  dev.off()

  

  pdf(file.path(output_dir, "14.4_EnrichrDotPlot_WikiPathway_2021_Human.pdf"), height=10, width=10)
  # enrichr dotplot
  EnrichrDotPlot(
    seurat_obj,
    mods = "all", # use all modules (this is the default behavior)
    database = "WikiPathway_2021_Human", # this has to be one of the lists we used above!!!
    n_terms=1 # number of terms for each module
  )
  dev.off()


  pdf(file.path(output_dir, "14.5_KEGG_2021_Human.pdf"), height=10, width=10)
  # enrichr dotplot
  EnrichrDotPlot(
    seurat_obj,
    mods = "all", # use all modules (this is the default behavior)
    database = "KEGG_2021_Human", # this has to be one of the lists we used above!!!
    n_terms=1 # number of terms for each module
  )
  dev.off()

 } 
  
  
 sessionInfo() 

