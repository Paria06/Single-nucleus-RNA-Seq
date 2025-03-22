
# help from this tutorial: http://htmlpreview.github.io/?https://github.com/RRHO2/RRHO2/blob/master/vignettes/RRHO2.html
closeAllConnections()
rm(list = ls())
library(ggpubr)
library(patchwork)
library(RRHO2)
library(tidyverse)
library(cowplot)
library(gplots)


setwd("/Users/paria/OneDrive/Rsession/")

res_ocu_cluster <-
  read_csv(file = 
             "/Users/paria/OneDrive/Rsession/raw_DEGs/1_raw_cluster_DEG_ocu.csv", show_col_types = FALSE) %>% 
  as.data.frame()
res_ocu_cluster$ID <-
  paste0(res_ocu_cluster$gene, "_", res_ocu_cluster$cluster_id)
ocu_cluster <- select(res_ocu_cluster, ID, pvalue, log2FoldChange)	
ocu_cluster$Score  <-  -log10(ocu_cluster$pvalue)
ocu_cluster$Score <- ocu_cluster$Score * (ocu_cluster$log2FoldChange)

gene_list1 <-
  data.frame(Genes = ocu_cluster$ID, DDE = ocu_cluster$Score)

res_sc_cluster <-
  read_csv(
    "/Users/paria/OneDrive/Rsession/raw_DEGs/1_raw_cluster_DEG_sc.csv", show_col_types = FALSE) %>% 
  as.data.frame()
res_sc_cluster$ID <-
  paste0(res_sc_cluster$gene, "_", res_sc_cluster$cluster_id)
sc_cluster <- select(res_sc_cluster, ID, pvalue, log2FoldChange)
sc_cluster$Score  <-  -log10(sc_cluster$pvalue)
sc_cluster$Score <- sc_cluster$Score * (sc_cluster$log2FoldChange)

gene_list2 <-
  data.frame(Genes = sc_cluster$ID, DDE = sc_cluster$Score)

df <- inner_join(gene_list1, gene_list2,  by = c("Genes"), suffix = c(".ocu", ".sc"))
gene_list1 <- cbind.data.frame(Genes = df$Genes, DDE = df$DDE.ocu, stringsAsFactors = FALSE)
gene_list2 <- cbind.data.frame(Genes = df$Genes, DDE = df$DDE.sc, stringsAsFactors = FALSE)
gene_list1 <- na.omit(gene_list1)
gene_list2 <- na.omit(gene_list2)

#gene_list1 <- sample_n(gene_list1, 20000)
clusterID <- intersect(unique(res_ocu_cluster$cluster_id), unique(res_sc_cluster$cluster_id)) %>% sort()
for(method_type in c("hyper")) {
  #c("fisher", "hyper")) {
  for(correction in c("none")) {
    #c("none", "BY", "BH")) {
    for(scaling in c(TRUE)) {
      #c(TRUE, FALSE)) {
      pdf(paste0("/Users/paria/OneDrive/Rsession/DEG/DE_correlation/Figure_hyper_none_TRUE_RRHO2_correlation_cell_ocusc_cluster.pdf"),
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
        
        #RRHO_obj<- NULL
        RRHO_obj <-
          RRHO2_initialize(gene_list2_cluster,
                           gene_list1_cluster,
                           labels = c(paste0("sc_",clusterID[i]),
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
          "/Users/paria/OneDrive/Rsession/DEG/DE_correlation/", clusterID[i],
          "hyper_none_TRUE_RRHO2_correlation_cell_ocusc_cluster.csv")) 
        
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

###############################Broad######################################

res_ocu_main <-
  read_csv(file = 
             "/Users/paria/OneDrive/Rsession/raw_DEGs/1_raw_main_DEG_ocu.csv", show_col_types = FALSE) %>% 
  as.data.frame()
res_ocu_main$ID <-
  paste0(res_ocu_main$gene, "_", res_ocu_main$cluster_id)
ocu_main <- select(res_ocu_main, ID, pvalue, log2FoldChange)	
ocu_main$Score  <-  -log10(ocu_main$pvalue)
ocu_main$Score <- ocu_main$Score * (ocu_main$log2FoldChange)

gene_list1 <-
  data.frame(Genes = ocu_main$ID, DDE = ocu_main$Score)

res_sc_main <-
  read_csv(
    "/Users/paria/OneDrive/Rsession/raw_DEGs/1_raw_main_DEG_ocu.csv", show_col_types = FALSE) %>% 
  as.data.frame()
res_sc_main$ID <-
  paste0(res_sc_main$gene, "_", res_sc_main$cluster_id)
sc_main <- select(res_sc_main, ID, pvalue, log2FoldChange)
sc_main$Score  <-  -log10(sc_main$pvalue)
sc_main$Score <- sc_main$Score * (sc_main$log2FoldChange)

gene_list2 <-
  data.frame(Genes = sc_main$ID, DDE = sc_main$Score)

df <- inner_join(gene_list1, gene_list2,  by = c("Genes"), suffix = c(".ocu", ".sc"))
gene_list1 <- cbind.data.frame(Genes = df$Genes, DDE = df$DDE.ocu, stringsAsFactors = FALSE)
gene_list2 <- cbind.data.frame(Genes = df$Genes, DDE = df$DDE.sc, stringsAsFactors = FALSE)


#gene_list1 <- sample_n(gene_list1, 20000)
mainID <- intersect(unique(res_ocu_main$cluster_id), unique(res_sc_main$cluster_id)) %>% sort()
for(method_type in c("hyper")) {
  #c("fisher", "hyper")) {
  for(correction in c("none")) {
    #c("none", "BY", "BH")) {
    for(scaling in c(TRUE)) {
      #c(TRUE, FALSE)) {
      pdf(paste0("/Users/paria/OneDrive/Rsession/DEG/DE_correlation/Figure_hyper_none_TRUE_RRHO2_correlation_cell_ocusc_main.pdf"),
          width = 6,
          height = 4,
          pointsize = 8)
      for (i in 1:length(mainID)) {
        #i = 2
        print(paste0("Running mainID = ", mainID[i]))
        #i =1
        gene_list1_main <- gene_list1[grep(paste0("_", mainID[i]),gene_list1$Genes),]
        gene_list2_main <- gene_list2[grep(paste0("_", mainID[i]),gene_list2$Genes),]
        gene_list2_main <- gene_list2_main[match(gene_list1_main$Genes, gene_list2_main$Genes), ]
        
        #RRHO_obj<- NULL
        RRHO_obj <-
          RRHO2_initialize(gene_list2_main,
                           gene_list1_main,
                           labels = c(paste0("sc_",mainID[i]),
                                      #"_Batch_Age_PMI_pH"), 
                                      paste0("ocu_",mainID[i])),
                           #"_Batch_Age_PMI_pH")),
                           log10.ind = TRUE, method = "hyper")
        
        print("Genes contributing to overlaps")
        head(RRHO_obj$genelist_dd$gene_list_overlap_dd) %>% print()
        head(RRHO_obj$genelist_uu$gene_list_overlap_uu) %>% print()
        head(RRHO_obj$genelist_du$gene_list_overlap_du) %>% print()
        head(RRHO_obj$genelist_ud$gene_list_overlap_ud) %>% print()
        max(RRHO_obj$hypermat, na.rm = TRUE) %>% print()	
        
        write.csv(RRHO_obj$hypermat, file = paste0(
          "/Users/paria/OneDrive/Rsession/DEG/DE_correlation/", mainID[i],
          "hyper_none_TRUE_RRHO2_correlation_cell_ocusc_main.csv")) 
        
        summary_table <- rbind(summary_table, c(mainID[i], method_type, correction, scaling, max(RRHO_obj$hypermat, na.rm = TRUE)))	
        
        try(RRHO2_heatmap(RRHO_obj))
      }
      try(RRHO2_vennDiagram(RRHO_obj, "dd"))
      try(RRHO2_vennDiagram(RRHO_obj, "uu"))
      try(RRHO2_vennDiagram(RRHO_obj, "du"))
      try(RRHO2_vennDiagram(RRHO_obj, "ud"))
      
    }
    dev.off()
  }  
}

sessionInfo()
