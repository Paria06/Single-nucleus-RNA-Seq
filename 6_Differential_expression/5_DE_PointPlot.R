rm(list = ls(all = T))
cat("\014")
setwd("/Users/paria/OneDrive/Rsession/")
base_path <- "/Users/paria/OneDrive/Rsession/"

#### Packs ####
library(dplyr)
library(ggplot2)


#### Data & Plot ####
base_path <- "/Users/paria/OneDrive/Rsession/"

directories <- c("ocu_cluster" = "DEG/ocu/cluster/", 
                 "ocu_main" = "DEG/ocu/main/", 
                 "med_cluster" = "/DEG/med/cluster/", 
                 "med_main" = "DEG/med/main/",
                 "sc_cluster" = "DEG/sc/cluster/",
                 "sc_main" = "DEG/sc/main/",
                 "sl_cluster" = "DEG/sl/cluster/",
                 "sl_main" = "DEG/sl/main/")

for(directory_name in names(directories)) {
  #directory_name="sc_main"
  directory <- directories[[directory_name]]
  data <- read_csv(file = paste0(base_path, directory, "2_filtered_results.csv"))

plotdata <- data %>%
  filter(padj < 0.05) %>%
  select(log2FoldChange, padj, cluster_id) %>%
  mutate(pval_cat = cut(padj, breaks = 0.01*c(0:5), labels = 0.01*c(1:5)))

pdf(file = paste0(base_path, 'Finalized_outputs/', directory_name, 'log2foldchange_pointplot.pdf'), width = 8, height = 6)
plotdata %>%
  ggplot(aes(x = log2FoldChange, y = cluster_id, color = padj)) +
  geom_point(size = 2) +
  geom_vline(xintercept = 0) +
  theme_classic() +
  scale_colour_steps(n.breaks = 4, guide = 'legend') +
  ggtitle(label = "Title", subtitle = "Subtitle") +
  labs(color = "p-value")
dev.off()
}

sessionInfo()