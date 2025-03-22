closeAllConnections()
rm(list = ls())
library(ggpubr)
library(patchwork)
library(RRHO2)
library(tidyverse)
library(cowplot)
library(gplots)

summary_table <- data.frame()

res_ocu_cluster <-
  read_csv(file = 
             "/Users/paria/OneDrive/Rsession/DEG/ocu/sig_cluster_DEG_ALS_vs_con_ocu.csv") %>% 
  as.data.frame()
res__cluster$ID <-
  paste0(res_ocu_cluster$gene, "_", res_ocu_cluster$cluster_id)
ocu_cluster <- select(res_ocu_cluster, ID, pvalue, log2FoldChange)	
ocu_cluster$Score  <-  -log10(ocu_cluster$p_val)
ocu_cluster$Score <- ocu_cluster$Score * (ocu_cluster$log2FoldChange)

gene_list1 <-
  data.frame(Genes = ocu_cluster$ID, DDE = ocu_cluster$Score)

res_med_cluster <-
  read_csv(
    "/Users/paria/OneDrive/Rsession/DEG/med/sig_cluster_DEG_ALS_vs_con_med.csv") %>% 
  as.data.frame()
res_med_cluster$ID <-
  paste0(res_med_cluster$gene, "_", res_med_cluster$cluster_id)
med_cluster <- select(res_med_cluster, ID, pvalue, log2FoldChange)
med_cluster$Score  <-  -log10(med_cluster$pvalue)
med_cluster$Score <- med_cluster$Score * (med_cluster$log2FoldChange)

gene_list2 <-
  data.frame(Genes = med_cluster$ID, DDE = med_cluster$Score)

df <- inner_join(gene_list1, gene_list2,  by = c("Genes"), suffix = c(".ocu", ".med"))
gene_list1 <- cbind.data.frame(Genes = df$Genes, DDE = df$DDE.ocu)
gene_list2 <- cbind.data.frame(Genes = df$Genes, DDE = df$DDE.med)


#gene_list1 <- sample_n(gene_list1, 20000)
clusterID <- intersect(unique(res_ocu_cluster$cluster_id), unique(res_ocu_cluster$cluster_id)) %>% sort()
for(method_type in c("hyper")) {
  #c("fisher", "hyper")) {
  for(correction in c("none", "BY")) {
    #c("none", "BY", "BH")) {
    for(scaling in c(TRUE)) {
      #c(TRUE, FALSE)) {
      pdf(paste0("/Users/paria/OneDrive/RsessionDEG/ocu_",method_type, "_", correction, "_", scaling, "_RRHO2_correlation_cell_ocumed_cluster.pdf"),
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
                               labels = c(paste0("med_",clusterID[i]),
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
          "/Users/paria/OneDrive/RsessionDEG/ocu", clusterID[i],
          method_type, "_", correction, "_", scaling, "_RRHO2_correlation_cell_ocumed_cluster.csv")) 
        
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


#Broad#####
res_ocu_cluster <-
  read_csv(file = 
             "/Users/paria/OneDrive/Rsession/DEG/ocu/sig_cluster_DEG_ALS_vs_con_ocu.csv") %>% 
  as.data.frame()
res_ocu_cluster$ID <-
  paste0(res_ocu_main$gene, "_", res_ocu_cluster$main)
ocu_main <- select(res_ocu_main, ID, pvalue, log2FoldChange)	
ocu_main$Score  <-  -log10(ocu_main$p_val)
ocu_main$Score <- ocu_main$Score * (ocu_main$log2FoldChange)

gene_list1 <-
  data.frame(Genes = ocu_main$ID, DDE = ocu_main$Score)

res_med_main <-
  read_csv(
    "/Users/paria/OneDrive/RsessionDEG/ocu/sig_cluster_DEG_ALS_vs_con_med.csv") %>% 
  as.data.frame()
res_med_main$ID <-
  paste0(res_med_cluster$gene, "_", res_med_main$main)
med_main <- select(res_med_main, ID, pvalue, log2FoldChange)
med_main$Score  <-  -log10(med_main$pvalue)
med_main$Score <- med_main$Score * (med_main$log2FoldChange)

gene_list2 <-
  data.frame(Genes = med_main$ID, DDE = med_main$Score)

df <- inner_join(gene_list1, gene_list2,  by = c("Genes"), suffix = c(".ocu", ".med"))
gene_list1 <- cbind.data.frame(Genes = df$Genes, DDE = df$DDE.ocu)
gene_list2 <- cbind.data.frame(Genes = df$Genes, DDE = df$DDE.med)


#gene_list1 <- sample_n(gene_list1, 20000)
clusterID <- intersect(unique(res_ocu_main$main), unique(res_ocu_main$main)) %>% sort()
for(method_type in c("hyper")) {
  #c("fisher", "hyper")) {
  for(correction in c("none", "BY")) {
    #c("none", "BY", "BH")) {
    for(scaling in c(TRUE)) {
      #c(TRUE, FALSE)) {
      pdf(paste0("/Users/paria/OneDrive/RsessionDEG/ocu_",method_type, "_", correction, "_", scaling, "_RRHO2_correlation_cell_ocumed_main.pdf"),
          width = 6,
          height = 4)
      for (i in 1:length(clusterID)) {
        #i = 2
        print(paste0("Running clusterID = ", clusterID[i]))
        #i =1
        gene_list1_main <- gene_list1[grep(paste0("_", clusterID[i]),gene_list1$Genes),]
        gene_list2_main <- gene_list2[grep(paste0("_", clusterID[i]),gene_list2$Genes),]
        gene_list2_main <- gene_list2_main[match(gene_list1_main$Genes, gene_list2_main$Genes), ]
        
        RRHO_obj<- NULL
        try(RRHO_obj <-
              RRHO2_initialize(gene_list2_main,
                               gene_list1_main,
                               labels = c(paste0("med_",clusterID[i]),
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
          "/Users/paria/OneDrive/RsessionDEG/ocu", clusterID[i],
          method_type, "_", correction, "_", scaling, "_RRHO2_correlation_cell_ocumed_main.csv")) 
        
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

write.csv(summary_table, "summary_table_RRHO.csv")

print(summary(warnings()))
sessionInfo()