
# help from this tutorial: http://htmlpreview.github.io/?https://github.com/RRHO2/RRHO2/blob/master/vignettes/RRHO2.html
closeAllConnections()
rm(list = ls())
library(ggpubr)
library(patchwork)
library(RRHO2)
library(tidyverse)
library(cowplot)
library(gplots)

setwd("/Users/paria/Library/CloudStorage/OneDrive-Personal/Rsession/analysis/")

res_ocu_cluster <- read_csv(file = "/Users/paria/Library/CloudStorage/OneDrive-Personal/Rsession/analysis/DEG/batch_regressed/ocu/cluster/1_raw_cluster.csv", show_col_types = FALSE) %>% 
  as.data.frame()

# Add an ID column to res_ocu_cluster
res_ocu_cluster$ID <- paste0(res_ocu_cluster$gene, "_", res_ocu_cluster$cluster_id)
ocu_cluster <- select(res_ocu_cluster, ID, pvalue, log2FoldChange)	
ocu_cluster$Score  <-  -log10(ocu_cluster$pvalue)
ocu_cluster$Score <- ocu_cluster$Score * (ocu_cluster$log2FoldChange)

gene_list1 <-
  data.frame(Genes = ocu_cluster$ID, DDE = ocu_cluster$Score)

res_sl_cluster <-
  read_csv(file = "/Users/paria/Library/CloudStorage/OneDrive-Personal/Rsession/analysis/DEG/batch_regressed/sl/cluster/1_raw_cluster.csv", show_col_types = FALSE) %>% 
  as.data.frame()
res_sl_cluster$ID <-
  paste0(res_sl_cluster$gene, "_", res_sl_cluster$cluster_id)
sl_cluster <- select(res_sl_cluster, ID, pvalue, log2FoldChange)
sl_cluster$Score  <-  -log10(sl_cluster$pvalue)
sl_cluster$Score <- sl_cluster$Score * (sl_cluster$log2FoldChange)

gene_list2 <-
  data.frame(Genes = sl_cluster$ID, DDE = sl_cluster$Score)

df <- inner_join(gene_list1, gene_list2,  by = c("Genes"), suffix = c(".ocu", ".sl"))
gene_list1 <- cbind.data.frame(Genes = df$Genes, DDE = df$DDE.ocu, stringsAsFactors = FALSE)
gene_list2 <- cbind.data.frame(Genes = df$Genes, DDE = df$DDE.sl, stringsAsFactors = FALSE)
gene_list1 <- na.omit(gene_list1)
gene_list2 <- na.omit(gene_list2)

#gene_list1 <- sample_n(gene_list1, 20000)
clusterID <- intersect(unique(res_ocu_cluster$cluster_id), unique(res_sl_cluster$cluster_id)) %>% sort()
for(method_type in c("hyper")) {
  #c("fisher", "hyper")) {
  for(correction in c("none", "BY")) {
    #c("none", "BY", "BH")) {
    for(scaling in c(TRUE)) {
      #c(TRUE, FALSE)) {
      pdf(paste0("/Users/paria/Library/CloudStorage/OneDrive-Personal/Rsession/analysis/RROH/Figure_",method_type, "_", correction, "_", scaling, "_RRHO2_correlation_cell_ocsl_cluster.pdf"),
          width = 6,
          height = 4,
          pointsize = 8)
      for (i in 1:length(clusterID)) {
        #i = 2
        print(paste0("Running clusterID = ", clusterID[i]))
        #i =1
        gene_list1_cluster <- gene_list1[grep(paste0("_", clusterID[i]),gene_list1$Genes),]
        gene_list2_cluster <- gene_list2[grep(paste0("_", clusterID[i]),gene_list2$Genes),]
        gene_list2_cluster <- gene_list2_cluster[match(gene_list1_cluster$Genes, gene_list2_cluster$Genes), ]
        
        RRHO_obj<- NULL
        RRHO_obj <-
          RRHO2_initialize(gene_list2_cluster,
                           gene_list1_cluster,
                           labels = c(paste0("sl_",clusterID[i]),
                                      #"_Batch_Age_PMI_pH"), 
                                      paste0("ocu_",clusterID[i])),
                           #"_Batch_Age_PMI_pH")),
                           log10.ind = TRUE, method = "hyper")
        
        print("Genes contributing to overlaps")
        head(RRHO_obj$genelist_dd$gene_list_overlap_dd) %>% print()
        head(RRHO_obj$genelist_uu$gene_list_overlap_uu) %>% print()
        head(RRHO_obj$genelist_du$gene_list_overlap_du) %>% print()
        head(RRHO_obj$genelist_ud$gene_list_overlap_ud) %>% print()
        max(RRHO_obj$hypermat, na.rm = TRUE) %>% print()	
        
        write.csv(RRHO_obj$hypermat, file = paste0(
          "/Users/paria/Library/CloudStorage/OneDrive-Personal/Rsession/analysis/RROH/RROH_", clusterID[i],
          method_type, "_", correction, "_", scaling, "_RRHO2_correlation_cell_ocusl_cluster.csv")) 
        
        summary_table <- rbind(summary_table, c(clusterID[i], method_type, correction, scaling, max(RRHO_obj$hypermat, na.rm = TRUE)))	
        
        try(RRHO2_heatmap(RRHO_obj))
        if(method_type == "hyper" & scaling == TRUE & correction == "none")  {
          try(RRHO2_heatmap(RRHO_obj, maximum = 50))
        }
        try(RRHO2_vennDiagram(RRHO_obj, "dd"))
        try(RRHO2_vennDiagram(RRHO_obj, "uu"))
        try(RRHO2_vennDiagram(RRHO_obj, "du"))
        try(RRHO2_vennDiagram(RRHO_obj, "ud"))
        
      }
      dev.off()
    }  
  }		
}	

##Broad

res_ocu_cluster <-
  read_csv(file = 
             "/Users/paria/Library/CloudStorage/OneDrive-Personal/Rsession/analysis/DEG/batch_regressed/ocu/main/1_raw_cluster.csv", show_col_types = FALSE) %>% 
  as.data.frame()
res_ocu_cluster$ID <-
  paste0(res_ocu_cluster$gene, "_", res_ocu_cluster$cluster_id)
ocu_cluster <- select(res_ocu_cluster, ID, pvalue, log2FoldChange)	
ocu_cluster$Score  <-  -log10(ocu_cluster$pvalue)
ocu_cluster$Score <- ocu_cluster$Score * (ocu_cluster$log2FoldChange)

gene_list1 <-
  data.frame(Genes = ocu_cluster$ID, DDE = ocu_cluster$Score)

res_sl_cluster <-
  read_csv(
    "/Users/paria/Library/CloudStorage/OneDrive-Personal/Rsession/analysis/DEG/batch_regressed/sl/main/1_raw_cluster.csv", show_col_types = FALSE) %>% 
  as.data.frame()
res_sl_cluster$ID <-
  paste0(res_sl_cluster$gene, "_", res_sl_cluster$cluster_id)
sl_cluster <- select(res_sl_cluster, ID, pvalue, log2FoldChange)
sl_cluster$Score  <-  -log10(sl_cluster$pvalue)
sl_cluster$Score <- sl_cluster$Score * (sl_cluster$log2FoldChange)

gene_list2 <-
  data.frame(Genes = sl_cluster$ID, DDE = sl_cluster$Score)

df <- inner_join(gene_list1, gene_list2,  by = c("Genes"), suffix = c(".ocu", ".sl"))
gene_list1 <- cbind.data.frame(Genes = df$Genes, DDE = df$DDE.ocu, stringsAsFactors = FALSE)
gene_list2 <- cbind.data.frame(Genes = df$Genes, DDE = df$DDE.sl, stringsAsFactors = FALSE)


#gene_list1 <- sample_n(gene_list1, 20000)
clusterID <- intersect(unique(res_ocu_cluster$cluster_id), unique(res_sl_cluster$cluster_id)) %>% sort()
for(method_type in c("hyper")) {
  #c("fisher", "hyper")) {
  for(correction in c("none", "BY")) {
    #c("none", "BY", "BH")) {
    for(scaling in c(TRUE)) {
      #c(TRUE, FALSE)) {
      pdf(paste0("/Users/paria/Library/CloudStorage/OneDrive-Personal/Rsession/analysis/RROH/Figure_",method_type, "_", correction, "_", scaling, "_RRHO2_correlation_cell_ocusl_main.pdf"),
          width = 6,
          height = 4)
      for (i in 1:length(clusterID)) {
        #i = 2
        print(paste0("Running clusterID = ", clusterID[i]))
        #i =1
        gene_list1_cluster <- gene_list1[grep(paste0("_", clusterID[i]),gene_list1$Genes),]
        gene_list2_cluster <- gene_list2[grep(paste0("_", clusterID[i]),gene_list2$Genes),]
        gene_list2_cluster <- gene_list2_cluster[match(gene_list1_cluster$Genes, gene_list2_cluster$Genes), ]
        
        RRHO_obj<- NULL
        try(RRHO_obj <-
              RRHO2_initialize(gene_list2_cluster,
                               gene_list1_cluster,
                               labels = c(paste0("sl_",clusterID[i]),
                                          #"_Batch_Age_PMI_pH"), 
                                          paste0("ocu_",clusterID[i])),
                               #"_Batch_Age_PMI_pH")),
                               log10.ind = scaling, method = method_type, multipleTesting = correction))
        
        print("Genes contributing to overlaps")
        head(RRHO_obj$genelist_dd$gene_list_overlap_dd) %>% print()
        head(RRHO_obj$genelist_uu$gene_list_overlap_uu) %>% print()
        head(RRHO_obj$genelist_du$gene_list_overlap_du) %>% print()
        head(RRHO_obj$genelist_ud$gene_list_overlap_ud) %>% print()
        max(RRHO_obj$hypermat, na.rm = TRUE) %>% print()	
        
        write.csv(RRHO_obj$hypermat, file = paste0(
          "/Users/paria/Library/CloudStorage/OneDrive-Personal/Rsession/analysis/RROH/RROH_", clusterID[i],
          method_type, "_", correction, "_", scaling, "_RRHO2_correlation_cell_ocusl_main.csv")) 
        
        summary_table <- rbind(summary_table, c(clusterID[i], method_type, correction, scaling, max(RRHO_obj$hypermat, na.rm = TRUE)))	
        
        try(RRHO2_heatmap(RRHO_obj))
        if(method_type == "hyper" & scaling == TRUE & correction == "none")  {
          try(RRHO2_heatmap(RRHO_obj, maximum = 50))
          #try(RRHO2_heatmap(RRHO_obj, maximum = 100))
        }
        try(RRHO2_vennDiagram(RRHO_obj, "dd"))
        try(RRHO2_vennDiagram(RRHO_obj, "uu"))
        try(RRHO2_vennDiagram(RRHO_obj, "du"))
        try(RRHO2_vennDiagram(RRHO_obj, "ud"))
        
      }
      dev.off()
    }  
  }		
}

colnames(summary_table) <- c("Cluster", "Method", "Correction", "Scaling", "Max_p")

write.csv(summary_table, "/Users/paria/Library/CloudStorage/OneDrive-Personal/Rsession/analysis/RROH/summary_table_RRHO.csv")

print(summary(warnings()))
sessionInfo()
