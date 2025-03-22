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
library(pheatmap)
library(MetBrewer)
omp_set_num_threads(1)
blas_set_num_threads(1)
theme_set(theme_cowplot())

# set random seed for reproducibility
set.seed(12345)

# optionally enable multithreading
enableWGCNAThreads(nThreads = 8)


# Set up output directories
outDir <- "/lustre07/scratch/paria/als-project/analysis/ocu_med_sc/rpca/WGCNA/across_region/main_level"


base_dirs <- c("Microglia" = "/Microglia/",
               "Astrocyte" = "/Astrocyte/", 
               "splatter" = "/splatter/",
               "Oligodendrocyte" = "/Oligodendrocyte/",
               "MN" = "/MN/")

color_schemes <- list(
  "Microglia" = function(n) met.brewer("Cross", n = n, type = 'discrete'),
  "Astrocyte" = function(n) met.brewer("Homer2", n = n, type = 'continuous'),
  "splatter" = function(n) met.brewer("Cassatt2", n = n, type = 'discrete'),
  "Oligodendrocyte" = function(n) met.brewer("Redon", n = n, type = 'discrete'),
  "MN" = function(n) met.brewer("Redon", n = n, type = 'discrete'))




cluster_sets <- list(
  "Microglia" = c("Mic1", "Mic2"),
  "Astrocyte" = c("WMAstrocyte", "GMAstrocyte"),
  "splatter" = c("splatter1", "splatter2", "MN", "InV"),
  "Oligodendrocyte" = c("OPC1", "OPC2", "OPC3", "COPC", "Olig1", "Olig2", "Olig3", "Olig4", "Olig5", "Olig6", "Olig7", "Olig8", "Olig9")
)


for (data_dir in names(base_dirs)) {

  # Load the Seurat object for the current dataset
  print("Load the Seurat object for the current dataset")
  directory <- base_dirs[[data_dir]]
  seurat_obj <- readRDS(file = paste0(outDir, directory, "hdWGCNA_object_", data_dir, ".rds"))  


  # Set up an output directory for results
  print("Set up an output directory for results")
  output_dir <- file.path(paste0(outDir, directory, "results"))
  dir.create(output_dir, showWarnings = FALSE) 

  print("Select soft-power threshold")
  # plot the results
  plot1 <- PlotSoftPowers(seurat_obj)
  # assemble with patchwork
  plot1 <- wrap_plots(plot1, ncol=2)


  pdf(file.path(output_dir, "1_PlotSoftPowers.pdf"), height=15, width=15)
  print(plot1)
  dev.off()
  
  
  power_table <- GetPowerTable(seurat_obj)

  print("PlotDendrogram")
  
  pdf(file.path(output_dir, "2_PlotDendrogram.pdf"), height=15, width=15)
  PlotDendrogram(seurat_obj, main = paste("hdWGCNA Dendrogram -", names(data_dir)))
  dev.off()
  

  # resset module name and colors
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

  
  print("Optional: inspect the topoligcal overlap matrix (TOM)")
  TOM <- GetTOM(seurat_obj)

 
  print("harmonized module eigengenes")
  # harmonized module eigengenes:
  hMEs <- GetMEs(seurat_obj)

  print("module eigengenes")
  # module eigengenes:
  MEs <- GetMEs(seurat_obj, harmonized=FALSE)


  print("plot genes ranked by kME for each module")
  # plot genes ranked by kME for each module
  plot3 <- PlotKMEs(seurat_obj)
  plot3 <- plot3 + theme(
  text = element_text(size = 8),   
  axis.text.x = element_text(size = 8),
  axis.text.y = element_text(size = 8),
  axis.title = element_text(size = 5))

  # Function to get module color
  get_module_color <- function(module_name, color_map) {
  return(color_map[module_name])
  }
  # Apply color to titles
  plot3 <- lapply(1:length(plot3), function(x) {
  module_name <- names(plot3)[x]
  module_color <- get_module_color(module_name, module_color_map)
  plot3[[x]] + NoLegend() +
    theme(
      plot.title = element_text(face='plain', color=module_color, size=6, vjust=0.25),
      plot.margin = margin(c(0, 0, 0, 0))
    )
  })

  pdf(file.path(output_dir, "3_PlotKMEs.pdf"), height=5, width=15)
  print(wrap_plots(plot3, ncol=3))
  dev.off()


  print("Getting the module assignment table")
  # Getting the module assignment table
  # get the module assignment table:
  modules <- GetModules(seurat_obj) %>% subset(module != 'grey')

  # show the first 6 columns:
  head(modules[,1:6])

  print("Compute hub gene signature scores")
  # Compute hub gene signature scores
  # compute gene scoring for the top 25 hub genes by kME for each module
  # with UCell method
  library(UCell)
  seurat_obj <- ModuleExprScore(
  seurat_obj,
  n_genes = 25,
  method='UCell')

  # plot Module Feature Plots
  print("plot Module Feature Plots")
  # Generate the feature plots
  plot4.1 <- ModuleFeaturePlot(seurat_obj, order=TRUE, raster=TRUE, raster_dpi=400, alpha=1, restrict_range=FALSE, raster_scale=0.25)
  
  # Ensure plot_list names match module_names order
  names(plot4.1) <- module_names
  
  # Function to get module color
  get_module_color <- function(module_name, color_map) {
    return(color_map[module_name])
  }
  
  # Apply color to titles
  plot4.1 <- lapply(1:length(plot4.1), function(x) {
    module_name <- names(plot4.1)[x]
    module_color <- get_module_color(module_name, module_color_map)
    plot4.1[[x]] + NoLegend() + 
      theme(
        plot.title = element_text(face='plain', color=module_color, size=5, vjust=0.25), 
        plot.margin = margin(c(0, 0, 0, 0))
      )
  })
  
  # Save the plots
  pdf(file.path(output_dir, "4.1_ModuleFeaturePlot.pdf"), height=5, width=5)
  print(wrap_plots(plot4.1, ncol=5))
  dev.off()

  # Hubgene circle plots:
  output_network_dir <- file.path(output_dir, "4.2_ModuleFeaturePlot_hub_genes")
  # Generate individual module networks
  ModuleNetworkPlot(
  seurat_obj,
  mods = "all",
  outdir = output_network_dir)

  # plot the hub gene signature score
  print("plot the hub gene signature score")
  # make a featureplot of hub scores for each module
  plot4.2 <- ModuleFeaturePlot(
  seurat_obj,
  features = 'scores',
  order_points = TRUE, 
  raster = TRUE, 
  raster_dpi = 400,
  alpha = 1,
  restrict_range = FALSE,
  raster_scale = 0.25)

  # Ensure plot_list names match module_names order
  names(plot4.2) <- module_names

  # Function to get module color
  get_module_color <- function(module_name, color_map) {
  return(color_map[module_name])
  }

  # Apply color to titles
  plot4.2 <- lapply(1:length(plot4.2), function(x) {
  module_name <- names(plot4.2)[x]
  module_color <- get_module_color(module_name, module_color_map)
  plot4.2[[x]] + NoLegend() +
    theme(
      plot.title = element_text(face='plain', color=module_color, size=5, vjust=0.25),
      plot.margin = margin(c(0, 0, 0, 0))
    )
  })

  # Save the plots
  pdf(file.path(output_dir, "4.2_ModuleFeaturePlot_hub_genes.pdf"), height=5, width=5)
  print(wrap_plots(plot4.2, ncol=5))
  dev.off()


  ##############plot module correlagram
  print("plot module correlagram")
  # plot module correlagram
  plot5 <- ModuleCorrelogram(seurat_obj)


  pdf(file.path(output_dir, "5_ModuleCorrelogram.pdf"), height=10, width=10)
  print(plot5)
  dev.off()

  # overlap with ALS_associated genes
  library(GeneOverlap)

  # ALS associated genes
  ALS.genes <- read.csv('/lustre07/scratch/paria/als-project/analysis/ocu_med_sc/rpca/WGCNA/data/ALS.csv')
  ALS.genes <- subset(ALS.genes, Gene %in% rownames(seurat_obj))

  # load modules
  modules <- GetModules(seurat_obj)
  saveRDS(modules, file.path(output_dir, "modules.rds"))

  mods <- levels(modules$module)
  genome.size <- nrow(modules)

  overlap_df <- do.call(rbind, lapply(mods, function(cur_mod){

  cur_genes <- modules %>% subset(module == cur_mod) %>% .$gene_name

  cur_overlap <- testGeneOverlap(newGeneOverlap(
      cur_genes,
      ALS.genes$Gene,
      genome.size=genome.size
  ))

  cur_overlap <- data.frame(
    'odds.ratio' = cur_overlap@odds.ratio,
    'pval' = cur_overlap@pval,
    'Jaccard' = cur_overlap@Jaccard,
    'size_intersection' = length(cur_overlap@intersection),
    'module' = cur_mod
  )

  cur_overlap

  })) %>% as.data.frame()

  overlap_df <- overlap_df %>% mutate(fdr=p.adjust(pval, method='BH'))
  overlap_df <- overlap_df %>% subset(module != 'grey')

  overlap_df$shape <- ifelse(overlap_df$fdr < 0.05, 21, 4)
  overlap_df <- overlap_df %>% arrange(odds.ratio, descending=TRUE)
  overlap_df$module <- factor(as.character(overlap_df$module), levels=as.character(overlap_df$module))

  mod_colors <- dplyr::select(modules, c(module, color)) %>%
  distinct
  cp <- mod_colors$color; names(cp) <- mod_colors$module

  plot6 <- overlap_df %>%
  ggplot(aes(y=module, x=odds.ratio, size= size_intersection, color=module)) +
  geom_segment(aes(y=module, yend=module, x=0, xend=odds.ratio), size=0.5, color='grey') +
  geom_point() +
  geom_point(shape=overlap_df$shape, color='black', fill=NA) +
  scale_color_manual(values=cp, guide='none') +
  ylab('') + xlab("Odds ratio") +
  scale_x_continuous(breaks = c(0, 1, 2,3)) +
  labs(size='Size\nintersection') +
  ggtitle('Overlap with ALS associated genes') +

  theme(
    panel.border = element_rect(size=1, color='black', fill=NA),
    axis.line.y = element_blank(),
    axis.line.x = element_blank(),
    plot.title = element_text(hjust=0.5, face='plain',size = 8),
    axis.title.x = element_text(size = 6),
    axis.title.y = element_text(size = 6),  
    axis.text.x = element_text(size = 6),    
    axis.text.y = element_text(size = 6),
    legend.title = element_text(size = 8), 
    legend.text = element_text(size = 8) )

  pdf(file.path(output_dir, "6_module_ALS_associated_genes_overlap.pdf"), width=4, height=3.5)
  print(plot6)
  dev.off()


  # UMAP visualization
  seurat_obj <- RunModuleUMAP(
  seurat_obj,
  n_hubs = 5,
  n_neighbors=15,
  min_dist=0.3,
  spread=1)

  
  # umap visualization no annotation
  # get the hub gene UMAP table from the seurat object
  umap_df <- GetModuleUMAP(seurat_obj)

  # plot with ggplot
  plot7.1 <- ggplot(umap_df, aes(x=UMAP1, y=UMAP2)) +
  geom_point(
   color=umap_df$color,
   size=umap_df$kME*2
  ) +
  umap_theme()

  pdf(file.path(output_dir, "7.1_umap_visualization.pdf"), height=3, width=4)
  print(plot7.1)
  dev.off()


  # umap visualization hub genes are labeled
  seurat_obj <- ResetModuleNames(
    seurat_obj,
    new_name = paste0("M"))

  # get the hub gene UMAP table from the seurat object
  umap_df <- GetModuleUMAP(seurat_obj)

  # plot with ggplot
  plot_df <- umap_df

  # compute coordinates for cluster labels
  centroid_df <- data.frame()
  for(cur_cluster in unique(plot_df[['module']])){
  cur_meta <- plot_df[plot_df[['module']] == cur_cluster,]
  df <- data.frame(
    cluster = cur_cluster,
    UMAP1 = mean(cur_meta$UMAP1),
    UMAP2 = mean(cur_meta$UMAP2)
  )
  centroid_df <- rbind(centroid_df, df)
  }
   # add annotation
  hub_genes <- GetHubGenes(seurat_obj, 3)
  anno_genes <- hub_genes$gene_name
  plot_df$anno <- ifelse(plot_df$gene %in% anno_genes, umap_df$gene, '')


  plot_df_anno <- subset(plot_df, anno != '')
  plot7.2 <-  plot_df %>%
  ggplot(aes(x=UMAP1, y=UMAP2, color=module)) +
  ggrastr::rasterise(
    geom_point(
      inherit.aes=FALSE,
      data=plot_df,
      aes(x=UMAP1, y=UMAP2, color=module),
      color=plot_df$color,
      size=plot_df$kME,
    ), dpi=800, dpi_scale=0.5) +
  geom_point(
    inherit.aes = FALSE,
    data = plot_df_anno,
    shape=21, color='black',
    fill=plot_df_anno$color,
    size=plot_df_anno$kME,
    aes(x=UMAP1, y=UMAP2, fill=module)
  ) +
  # add labels
  ggrepel::geom_text_repel(data = centroid_df, label=centroid_df$cluster, color='black', max.overlaps=Inf, size=2, fontface='bold') +
  geom_text_repel(label=plot_df$anno, max.overlaps=Inf, color='black', fontface='italic', size=1.5) +
  umap_theme() + NoLegend() +
  coord_equal() +
  theme(
    plot.margin = margin(0,0,0,0)
  )

  pdf(file.path(output_dir, "7.2_umap_visualization_hubgenes.pdf"), height=3, width=4)
  print(plot7.2)
  dev.off()


  # umap visualization ALS associated genes are labeled
  ALS.genes <- read.csv('/lustre07/scratch/paria/als-project/analysis/ocu_med_sc/rpca/WGCNA/data/ALS.csv')
  ALS.genes <- subset(ALS.genes, Gene %in% rownames(seurat_obj))

  hub_genes <- GetHubGenes(seurat_obj, 50)
  write.csv(hub_genes, file.path(output_dir, "hub_genes.csv"))
  label_genes <- intersect(hub_genes$gene_name, ALS.genes$Gene)

  ALS.hubs <- subset(GetHubGenes(seurat_obj, 50), gene_name %in% ALS.genes$Gene)
  write.csv(ALS.hubs, quote=FALSE, file.path(output_dir, 'ALS_hubgenes.csv'))

  # get the hub gene UMAP table from the seurat object
  umap_df <- GetModuleUMAP(seurat_obj)

  # get the hub gene UMAP table from the seurat object
  umap_df <- GetModuleUMAP(seurat_obj)

  # plot with ggplot
  plot_df <- umap_df

  # compute coordinates for cluster labels
  centroid_df <- data.frame()
  for(cur_cluster in unique(plot_df[['module']])){
  cur_meta <- plot_df[plot_df[['module']] == cur_cluster,]
  df <- data.frame(
    cluster = cur_cluster,
    UMAP1 = mean(cur_meta$UMAP1),
    UMAP2 = mean(cur_meta$UMAP2)
  )
  centroid_df <- rbind(centroid_df, df)
  }
 
  ### als hub genes
  hub_genes <- ALS.hubs
  anno_genes <- hub_genes$gene_name
  plot_df$anno <- ifelse(plot_df$gene %in% anno_genes, umap_df$gene, '')

  plot_df_anno <- subset(plot_df, anno != '')
  plot7.3 <-  plot_df %>%
  ggplot(aes(x=UMAP1, y=UMAP2, color=module)) +
  ggrastr::rasterise(
    geom_point(
      inherit.aes=FALSE,
      data=plot_df,
      aes(x=UMAP1, y=UMAP2, color=module),
      color=plot_df$color,
      size=plot_df$kME,
    ), dpi=800, dpi_scale=0.5) +
  geom_point(
    inherit.aes = FALSE,
    data = plot_df_anno,
    shape=21, color='black',
    fill=plot_df_anno$color,
    size=plot_df_anno$kME,
    aes(x=UMAP1, y=UMAP2, fill=module)
  ) +
  # add labels
  ggrepel::geom_text_repel(data = centroid_df, label=centroid_df$cluster, color='black', max.overlaps=Inf, size=2) +
  geom_text_repel(label=plot_df$anno, max.overlaps=Inf, color='black', fontface='italic', size=2.5, fontface='bold') +
  umap_theme() + NoLegend() +
  coord_equal() +
  theme(
    plot.margin = margin(0,0,0,0)
  )

  pdf(file.path(output_dir, "7.3_hubgene_umap_igraph.pdf"), height=3, width=4)
  print(plot7.3)
  dev.off()

  # save the risk genes
  # Convert gene names to uppercase for case insensitive matching
  genes_of_interest <- anno_genes
  
  # Create a data frame to store gene and module information
  gene_module_info <- data.frame(
    Gene = genes_of_interest,
    Module = sapply(genes_of_interest, function(gene) {
      module <- modules[which(modules$gene == gene), "module"]
      return(ifelse(length(module) == 0, NA, module))
    })
  )
  
  print(gene_module_info)
  write.csv(gene_module_info, file.path(output_dir, "gene_module_info.csv"))
  
  
  # refer to this to customize the network plot: https://smorabit.github.io/hdWGCNA/articles/network_visualizations.html
  # showing only the genes that are involved in different GO term based on Gene Ontology


  seurat_obj <- ResetModuleNames(
  seurat_obj,
  new_name = paste0(data_dir, "-", "M"))

  
  markers <- read.csv("/home/paria/scratch/als-project/analysis/ocu_med_sc/rpca/WGCNA/data/markers_main.csv")

  # compute marker gene overlaps
  overlap_df <- OverlapModulesDEGs(
  seurat_obj,
  deg_df = markers,
  fc_cutoff = 1)


  # overlap barplot, produces a plot for each cell type
  plot8.1 <- OverlapBarPlot(overlap_df)

  plot8.1 <- lapply(plot8.1, function(p) {
  p + theme(
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    plot.title = element_text(size = 14))
  })
  plot8.1 <- wrap_plots(plot8.1, ncol=4)
  
  pdf(file.path(output_dir, "8.1_OverlapBarPlot.pdf"), height=15, width=15)
  # stitch plots with patchwork
  print(plot8.1)
  dev.off()


  pdf(file.path(output_dir, "8.2_OverlapDotPlot.pdf"), height=10, width=10)
  # plot odds ratio of the overlap as a dot plot
  plot8.2 <- OverlapDotPlot(
  overlap_df,
  plot_var = 'odds_ratio') +
  ggtitle('Overlap of modules & cell-type markers')

  print(plot8.2)
  dev.off()

  #Differential module eigengene (DME) analysis
  seurat_obj <- ResetModuleNames(
  seurat_obj,
  new_name = paste0("M"))

  print("Differential module eigengene (DME) analysis")
  group1 <- seurat_obj@meta.data %>% subset(resistant == "V" & group_id == "als") %>% rownames
  group2 <- seurat_obj@meta.data %>% subset(resistant == "R" & group_id == "als") %>% rownames
  # region == "sl" & 
  # region == "ocu" & 
  head(group1)
  
  # Check if groups are empty and skip iteration if true
  if (length(group1) == 0 | length(group2) == 0) {
    print("One or more gene lists are empty. Skipping this iteration.")
    next
  }

  DMEs <- FindDMEs(
  seurat_obj,
  barcodes1 = group1,
  barcodes2 = group2,
  test.use='wilcox',
  wgcna_name= data_dir)

  head(DMEs)

  pdf(file.path(output_dir, "9_PlotDMEsLollipop_VvsR_ALS.pdf"), height=8, width=8)
  plot9.1 <- PlotDMEsLollipop(
  seurat_obj,
  DMEs,
  wgcna_name= data_dir,
  pvalue = "p_val_adj")
  print(plot9.1)


  plot9.2 <- PlotDMEsVolcano(
  seurat_obj,
  DMEs,
  wgcna_name = data_dir)
  print(plot9.2)

  plot9.3 <- PlotDMEsVolcano(
  seurat_obj,
  DMEs,
  wgcna_name = data_dir)
  print(plot9.3)
  
  dev.off()

  # Check if the value of data_dir exists in the cluster_sets
  #clusters <- c("Mic1", "Mic2")
  if (data_dir %in% names(cluster_sets)) {
    clusters <- cluster_sets[[data_dir]]
  } else {
    
  }
  
  
  DMEs <- data.frame()

  # loop through the clusters
  for(cur_cluster in clusters){
  # identify barcodes for group1 and group2 in each cluster
  # with the below code the als as group1 and con as group2 the con is considered the reference.
    group1 <- seurat_obj@meta.data %>% subset(resistant == "V" & group_id == "als") %>% rownames
    group2 <- seurat_obj@meta.data %>% subset(resistant == "R" & group_id == "als") %>% rownames 
  # region != "ocu"
  # region == "ocu" & 

  # run the DME test
  cur_DMEs <- FindDMEs(
   seurat_obj,
  barcodes1 = group1,
  barcodes2 = group2,
  test.use='wilcox',
  pseudocount.use=0.01, 
  wgcna_name = data_dir)
  
  # add the cluster info to the table
  cur_DMEs$cluster <- cur_cluster
  
  # append the table
  DMEs <- rbind(DMEs, cur_DMEs)

  }

  # get the modules table:
  modules <- GetModules(seurat_obj)
  mods <- levels(modules$module) 
  mods <- mods[mods != 'grey']

  # make a copy of the DME table for plotting
  plot_df <- DMEs

  # set the factor level for the modules so they plot in the right order:
  plot_df$module <- factor(as.character(plot_df$module), levels=mods)

  # set a min/max threshold for plotting
  maxval <- 3; minval <- -3
  plot_df$avg_log2FC <- ifelse(plot_df$avg_log2FC > maxval, maxval, plot_df$avg_log2FC)
  plot_df$avg_log2FC <- ifelse(plot_df$avg_log2FC < minval, minval, plot_df$avg_log2FC)

  # add significance levels
  plot_df$Significance <- gtools::stars.pval(plot_df$p_val_adj)

  # change the text color to make it easier to see 
  plot_df$textcolor <- ifelse(plot_df$avg_log2FC > 0.2, 'white', 'black')

  # make the heatmap with geom_tile
  plot10 <- plot_df %>% 
  ggplot(aes(y=cluster, x=module, fill=avg_log2FC)) +
  geom_tile() 

  # add the significance levels
  plot10 <- plot10 + geom_text(label=plot_df$Significance, color=plot_df$textcolor) 

  # customize the color and theme of the plot
  plot10 <- plot10 + 
  scale_fill_gradient2(low='white', mid='lightblue', high='midnightblue') +
  RotatedAxis() +
  theme(
    panel.border = element_rect(fill=NA, color='black', size=1),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    axis.text.x = element_text(size = 8),    # Customize x-axis text size
    axis.text.y = element_text(size = 8),    # Customize y-axis text size
     legend.text = element_text(size = 8),    # Customize legend text size
    legend.title = element_text(size = 8),
    plot.margin=margin(0,0,0,0)
  ) + xlab('') + ylab('') +
  coord_equal()


  pdf(file.path(output_dir, "10_DME_effect_sizes_heatmap_VvsR_ALS.pdf"), height=5, width=10)
  print(plot10)
  dev.off()
  
  ######## Module trait correlation
  print("Module trait correlation")
  # convert sex to factor
  seurat_obj$group_id <- as.factor(seurat_obj$group_id)
  seurat_obj$sex <- as.factor(seurat_obj$sex)
  

  plot11 <- PlotModuleTraitCorrelation(
  seurat_obj,
  label = 'fdr',
  label_symbol = 'stars',
  text_size = 2,
  text_digits = 2,
  text_color = 'white',
  high_color = 'midnightblue',
  mid_color = 'lightblue',
  low_color = 'white',
  plot_max = 0.75,
  combine=TRUE)

  plot11 <- plot11 +
  theme(
    axis.title.x = element_text(size = 10),
    axis.title.y = element_text(size = 10))

  pdf(file.path(output_dir, "11_PlotModuleTraitCorrelation.pdf"), height=18, width=8)
  print(plot11)
  dev.off()   

}

sessionInfo()

