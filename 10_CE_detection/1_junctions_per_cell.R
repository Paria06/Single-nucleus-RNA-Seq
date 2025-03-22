
# Load necessary libraries
library(tidyverse)
library(data.table)  
library(optparse)
library(Seurat)

# Define sample list
samples <- c("T01499", "T01500", "T01487", "T01488", "T01491", "T01492", "T01479", "T01515", "T01516", "T01497", "T01489", "T01480", 
             "T01498", "T01486", "T01490", "T01477", "T01478", "T01513", "T01514")  # Add  sample IDs

# Loop through each sample
for (sample in samples) {
    prefix <- paste0("/home/paria/scratch/STARsolo/", sample)
    
    # Read in filtered GeneFull count matrix
    barcode_file <- paste0(prefix, "/", sample, "Solo.out/GeneFull/filtered/barcodes.tsv")
    barcodes <- readLines(barcode_file)
    
    # Read in raw count data assigned to junctions
    features <- read.csv(paste0(prefix, "/", sample, "SJ.out.tab"), sep = "\t", header = FALSE)
    features$new_column <- paste0(features$V1, ":", features$V2, "-", features$V3)
    write.table(features$new_column, file = paste0(prefix, "/", sample, "Solo.out/SJ/raw/features.tsv"),
                sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
    
    sj_folder <- file.path(prefix, paste0(sample, "Solo.out/SJ/raw"))
    sj <- ReadSTARsolo(data.dir = sj_folder, feature.column = 1)
    sj <- sj[, barcodes, drop = FALSE]
    sj <- sj[rowSums(sj) > 0, ]
    
    # Format and read the junction file
    junctions_file <- file.path(prefix, "output.sj") 
    junc <- vroom::vroom(junctions_file, 
                         col_names = c("chrom", "A", "B", "dot", "n_reads", "strand", "rA", "rb", "rgb", 
                                       "blockCount", "blockSize", "blockStarts"),
                         col_types = "cnncncnnnccc")
    
    # Get coordinates
    junc$Aoff <- as.numeric(gsub(",.*", "", junc$blockSize))
    junc$Boff <- as.numeric(gsub(".*,", "", junc$blockSize))
    junc$start <- junc$A + junc$Aoff
    junc$end <- junc$B - junc$Boff + 1
    junc$coord <- paste0(junc$chrom, ":", junc$start, "-", junc$end)
    junc <- junc %>% select(dot, coord)
    
    # Find detected junctions
    detected <- intersect(junc$coord, rownames(sj))
    print(paste0(" * ", length(detected), " cryptic junctions found in ", sample))
    
    # Process junctions for STMN2 and UNC13A
    coords <- rownames(sj)
    coords_df <- data.frame(
        chromosome = sapply(coords, function(x) strsplit(x, ":")[[1]][1]),
        start = as.integer(sapply(coords, function(x) strsplit(strsplit(x, ":")[[1]][2], "-")[[1]][1])),
        end = as.integer(sapply(coords, function(x) strsplit(strsplit(x, ":")[[1]][2], "-")[[1]][2]))
    )
    
    filtered_coords_sj_STMN <- coords_df[coords_df$chromosome == 'chr8' & 
                                             coords_df$start >= 79611210 & 
                                             coords_df$end <= 79636830, ]
    
    filtered_coords_sj_UNC13A <- coords_df[coords_df$chromosome == 'chr19' & 
                                               coords_df$start >= 17642845 & 
                                               coords_df$end <= 17641556, ]
    
    print(filtered_coords_sj_STMN)
    print(filtered_coords_sj_UNC13A)
    
    # Extract rows from sj based on filtered STMN2 coordinates
    filtered_coords_sj_STMN <- rownames(filtered_coords_sj_STMN)
    df <- sj[filtered_coords_sj_STMN, barcodes, drop = FALSE]
    
    # Transform data
    sj_df <- as.matrix(df) %>%
        as.data.frame() %>%
        tibble::rownames_to_column(var = "junction") %>%
        tidyr::pivot_longer(names_to = "barcode", values_to = "count", cols = !c(junction)) %>%
        dplyr::mutate(sample = sample) %>%
        filter(count > 0)
    
    # Load Seurat object for barcode-to-cell-type mapping
    seurat_object <- readRDS("/home/paria/scratch/als-project/analysis/ocu_med_sc/rpca/pseudo_bulk_DEG/seurat_integrated.rds")
    
    # Create a dataframe mapping barcodes to cell types
    barcode_to_celltype <- data.frame(
        barcode = rownames(seurat_object@meta.data),
        cell_type = seurat_object@meta.data$cluster_id
    )
    
    barcode_to_celltype$barcode <- sub("-1_1$", "", barcode_to_celltype$barcode)
    
    # Merge barcodes with cell types
    sj_df <- sj_df %>%
        left_join(barcode_to_celltype, by = "barcode")
    
    # Save output
    outFile <- paste0("/home/paria/scratch/STARsolo/", sample, "_sj_output.tsv")
    write_tsv(sj_df, outFile)
    
    print(paste0("Saved results for ", sample, " to ", outFile))
}
 sessionInfo()
 