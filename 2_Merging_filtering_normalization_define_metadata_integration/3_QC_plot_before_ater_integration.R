library(Seurat)
library(gridExtra)

# Load your Seurat objects
seurat_filtered <- readRDS("seurat_filtered.rds")
seurat_integrated <- readRDS("seurat_integrated.rds")

# Generate UMAP plots for seurat_filtered
p1_before <- DimPlot(seurat_filtered, reduction = 'umap', group.by = 'status')
p2_before <- DimPlot(seurat_filtered, reduction = 'umap', group.by = 'individual')
p3_before <- DimPlot(seurat_filtered, reduction = 'umap', group.by = 'sex')
p4_before <- DimPlot(seurat_filtered, reduction = 'umap', group.by = 'batch')

# Generate UMAP plots for seurat_integrated
p1_after <- DimPlot(seurat_integrated, reduction = 'umap', group.by = 'status')
p2_after <- DimPlot(seurat_integrated, reduction = 'umap', group.by = 'individual')
p3_after <- DimPlot(seurat_integrated, reduction = 'umap', group.by = 'sex')
p4_after <- DimPlot(seurat_integrated, reduction = 'umap', group.by = 'batch')

# Arrange and save plots to PDF
pdf("UMAP_comparison.pdf", height = 20, width = 14)
grid.arrange(
  arrangeGrob(p1_before, p1_after, ncol = 1),
  arrangeGrob(p2_before, p2_after, ncol = 1),
  arrangeGrob(p3_before, p3_after, ncol = 1),
  arrangeGrob(p4_before, p4_after, ncol = 1),
  ncol = 1
)
dev.off()

sessionInfor()
