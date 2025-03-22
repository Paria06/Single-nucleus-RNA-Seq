library(Seurat)
library(future)
library(dplyr)
library(scclusteval)
library(purrr)

# Load Data 
seurat_integrated <- readRDS("seurat_integrated_DB_removed")

## Dot plot depicting gene marker expression  and dendrogram

# Dendrogram
markers <- c("OLIG1", "PDGFRA", "VCAN", "PCDH15", "MOBP", "MYRF", "MOG", "MAG", "LRP2", "CNP", "CLDN11", "LYN", "LY86", "CX3CR1", "P2RY12", "CSF1R", "PTPRC", "GFAP", "AQP4", "SLC14A1", "CPAMD8", "FLT1", "CLDN5", "ERG", "NOSTRIN", "PDGFRB", "LAMC3", "NOTCH3", "DNAH11", "CFAP73", "CFAP157", "CD96", "CD247", "ITK", "RBFOX3", "GABRB2", "SYP", "SNAP25", "GATA3", "NEFH", "NEFL", "NEFM", "TUBA4A", "ACLY")
DefaultAssay(seurat_integrated) <- "RNA"

# Normalize RNA data for visualization purposes
seurat_integrated <- NormalizeData(seurat_integrated, verbose = FALSE)

seurat_integrated <- BuildClusterTree(seurat_integrated, features =markers, reduction ="pca")
ape::write.nexus(seurat_integrated@tools$BuildClusterTree, file = "10_dendrogram.nex")


# Dotplot
markers <- c("NEFH", "NEFL", "NEFM", "TUBA4A", "ACLY", "SYP", "SNAP25", "INA", "FLT1", "CLDN5", "PECAM1", "ERG", "NOSTRIN", "PDGFRA", "VCAN", "OLIG1" , "GFAP", "AQP4", "SLC14A1", "CPAMD8", "MOBP", "MYRF", "MOG", "MAG", "CD96", "CD247", "ITK", "DNAH11", "CFAP73", "CFAP157",  "FOXJ1", "RBFOX3", "GABRB2", "GATA3", "GAD1",  "GAD2",  "LYN", "LY86", "CX3CR1", "P2RY12", "CSF1R", "DCN", "LAMC3", "PDGFRB", "ACTA2", "NOTCH3")

DefaultAssay(seurat_integrated) <- "RNA"

# Normalize RNA data for visualization purposes
seurat_integrated <- NormalizeData(seurat_integrated, verbose = FALSE)
Idents(seurat_integrated) <- "seurat_clusters"
my_levels <- c(10,19,11,14,29,25,9,5,6,12,1,7,15,0,4,16,21,2,13,20,27,3,17,22,28,24,23,26)

seurat_integrated <- SetIdent(seurat_integrated, value = "seurat_clusters")
levels(seurat_integrated) <- my_levels

# Re-level object@ident
pdf("DotPlot.pdf", height=20, width=12)
DotPlot(object = seurat_integrated,
        features = markers, 
        dot.scale = 25, 
        col.min = 0.5,
        col.max = 15
) + scale_size_area(max_size = 14) + scale_colour_gradient2(low ="#F2F2F2",  mid ="#F0C9C0",  high ="#2DB600",  space ="Lab",  na.value ="white",  guide ="colourbar",  aesthetics ="colour")  + theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 12)) + coord_flip()
dev.off()

# Rename all identities
seurat_integrated_labelled <- RenameIdents(object = seurat_integrated,  "0" = "Olig6", "1" = "Olig3", "2" = "Mic1",  "3" = "ExInVentral1", "4" = "Olig5", "5" = "WMAstrocyte", "6" = "Olig4",
                                           "7" = "Olig2", "8" = "Olig9", "9" = "GMAstrocyte", "10" = "MN", "11" = "OPC1", "12" = "Olig1",  "13" = "Mic2", "14" = "OPC2", "15" = "Olig7",
                                           "16" = "EpendymalCell", "17" = "ExInVentral2", "18" = "Olig8", "19" = "EC", "20" = "ExDorsal1", "21" = "Lymphocyte", "22" = "InVentral", "23" = "ExDorsal2",
                                           "24" = "Meninge", "25" = "OPC3", "26" = "InDorsal", "27" = "ExDorsal3", "28" = "Pericyte", "29" = "COPC")


pdf("Cluster_labeled_UMAP_ocu_med_sc.pdf", height=30, width=40)
p0 <- DimPlot(seurat_integrated, reduction = "umap", group.by = "cluster_id", pt.size = 0.4, label = TRUE, label.size = 10, raster=FALSE) + theme(aspect.ratio=1)
p0
dev.off()

sessionInfo()