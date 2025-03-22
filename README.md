# Project Overview

This repository contains various analysis methods and scripts for processing single-nucleus RNA sequencing (snRNA-seq) data. The following sections describe the steps, along with corresponding bash scripts and R scripts used in the project.

## 1. Preprocessing Methods
### Generation of Expression Matrix and Preprocessing of Single Nuclei Dataset
- **1_cellranger.sh**: FASTQ files are mapped to the GRCh38 reference genome to get the gene read count using the Cellranger count.

## 2. Merging, Filtering, Normalization, Metadata, and Integration
### Methods: Dimension reduction, integration, and clustering
- **bash_files**: Contains the bash scripts used to submit R scripts.
  - `Merging_filtering_normalization_metadata.sh`
  - `Integration_rPCA.sh`
- **Scripts**:
  - `1_Merging_filtering_normalization_metadata.R`: Create Seurat object, merge, and filter the dataset based on Satijalab vignettes.
  - `2_Integration_rPCA.R`: Normalization, finding variable features, integration, data scaling, and dimension reduction.
  - `3_QC_plot_before_after_integration.R`: Impact of different variables evaluated by plotting UMAP before and after integration.

## 3. ScClustEval Analysis
### Methods: Dimension reduction, integration, and clustering
- **pyflow_seurat_parameter_custom**: Files corresponding to the scClustEval Snakemake workflow.
  - `cluster.json`: Test the stability of clusters with different clustering parameters.
  - `config.yaml`: Configuration file for Snakemake workflow.
  - `pyflow-scBoot.sh`: Test cluster stability using bootstrapping.
- **scripts**:
  - `gather_fullsample.R`: Process full sample to test cluster stability.
  - `gather_subsample.R`: Process subsample to test cluster stability.
  - `preprocess.R`: Preprocess data for stability tests.
  - `subsample.R`: Handle subsampling for cluster stability testing.
- **Snakefile**: The Snakemake workflow file for managing the pipeline.

## 4. Cluster Evaluation and Annotation
### Methods: Quality control of clustering; Cluster annotation using automated and manual approaches
- **bash_files**: Contains bash scripts for running R scripts.
  - `Evaluation_of_scclusteval_output.sh`
  - `Doublet_removal.sh`
  - `Cluster_annotation_singleR_spinalcord.sh`
  - `Cluster_annotation_singleR_medulla.sh`
  - `Cluster_annotation_singleR_oculomotor.sh`
  - `Cluster_markers_manual.sh`
  - `Cluster_QC.sh`
- **Scripts**:
  - `1_Evaluation_of_scclusteval_output.R`: Visualize scClustEval output to determine the optimal combination of clustering parameters.
  - `2_Doublet_removal.R`: Remove doublets using the scDblFinder package.
  - `3_Cluster_annotation_singleR_spinalcord.R`: Automated approach using SingleR to annotate clusters with spinal cord dataset.
  - `4_Cluster_annotation_singleR_medulla.R`: Annotate clusters with midbrain dataset using SingleR.
  - `5_Cluster_annotation_singleR_oculomotor.R`: Annotate clusters with oculomotor dataset using SingleR.
  - `6_Cluster_markers_manual.R`: Generate cluster markers manually using Presto.
  - `7_Cluster_QC.R`: Evaluate QC parameters of clusters.

## 5. Cluster Visualization and Cell Type Proportion
### Methods: Cluster annotation using automated and manual approaches; Cell type proportions comparison
- **bash_files**:
  - `Clusters_visualization.sh`
  - `Cluster_celltype_proportion.sh`
- **Scripts**:
  - `1_Clusters_visualization.R`: Visualize cluster markers.
  - `2_Cluster_celltype_proportion.R`: Proportion of clusters and main cell types across different variables using Speckle package.

## 6. Differential Expression
### Methods: Differential Gene Expression Analysis using DESeq2
- **bash_files**:
  - `DEG_cluster_id_ALS_vs_con_oculomotor.sh`
  - `DEG_main_ALS_vs_con_oculomotor.sh`
  - `DEG_cluster_id_ALS_vs_con_medulla.sh`
  - `DEG_main_ALS_vs_con_medulla.sh`
  - `DEG_cluster_id_ALS_vs_con_cervicalspinalcord.sh`
  - `DEG_main_ALS_vs_con_cervicalspinalcord.sh`
  - `DE_PointPlot.sh`
- **Scripts**:
  - `1_1_DEG_cluster_id_ALS_vs_con_oculomotor.R`: DE analysis using DESeq2 comparing ALS vs control in oculomotor at cluster level.
  - `1_2_DEG_main_ALS_vs_con_oculomotor.R`: DE analysis using DESeq2 comparing ALS vs control in oculomotor at main level.
  - `2_1_DEG_cluster_id_ALS_vs_con_medulla.R`: DE analysis using DESeq2 comparing ALS vs control in medulla at cluster level.
  - `2_2_DEG_main_ALS_vs_con_medulla.R`: DE analysis using DESeq2 comparing ALS vs control in medulla at main level.
  - `3_1_DEG_cluster_id_ALS_vs_con_cervicalspinalcord.R`: DE analysis using DESeq2 comparing ALS vs control in cervical spinal cord at cluster level.
  - `3_2_DEG_main_ALS_vs_con_cervicalspinalcord.R`: DE analysis using DESeq2 comparing ALS vs control in cervical spinal cord at main level.
  - `4_1_DEG_cluster_id_ALS_vs_con_lumbarspinalcord.R`: DE analysis using DESeq2 comparing ALS vs control in lumbar spinal cord at cluster level.
  - `4_2_DEG_main_ALS_vs_con_lumbarspinalcord.R`: DE analysis using DESeq2 comparing ALS vs control in lumbar spinal cord at main level.
  - `5_DE_PointPlot.R`: Visualization of DEGs across different clusters.

## 7. RRHO Analysis
### Methods: Rank-Rank Hypergeometric Overlap without Thresholds
- **bash_files**:
  - `RROH_DE_correlation_ocumed.sh`
  - `RROH_DE_correlation_ocusc.sh`
  - `RROH_DE_correlation_ocusl.sh`
- **Scripts**:
  - `1_RROH_DE_correlation_ocumed.R`: Correlation of oculomotor DEG and medulla DEG without significance threshold.
  - `2_RROH_DE_correlation_ocusc.R`: Correlation of oculomotor DEG and cervical spinal cord DEG without significance threshold.
  - `3_RROH_DE_correlation_ocusl.R`: Correlation of oculomotor DEG and lumbar spinal cord DEG without significance threshold.

## 8. Differential Vulnerability Meta-analysis
### Methods: Differential Gene Expression Analysis using DESeq2
- **bash_files**:
  - `DEG_Differential_vulnerability_medocu_cluster.sh`
  - `DEG_Differential_vulnerability_medocu_main.sh`
  - `DEG_Differential_vulnerability_scocu_cluster.sh`
  - `DEG_Differential_vulnerability_scocu_main.sh`
  - `DEG_Differential_vulnerability_slocu_cluster.sh`
  - `DEG_Differential_vulnerability_slocu_main.sh`
  - `Differential_vulnerability_PointPlot.sh`
  - `meta_analysis.sh`
  - `meta_analysis_barplot.sh`
- **Scripts**:
  - `1_1_DEG_Differential_vulnerability_medocu_cluster.R`: Differential vulnerability analysis comparing medulla vs oculomotor in ALS at cluster level.
  - `1_2_DEG_Differential_vulnerability_medocu_main.R`: Differential vulnerability analysis comparing medulla vs oculomotor in ALS at main level.
  - `2_1_DEG_Differential_vulnerability_scocu_cluster.R`: Differential vulnerability analysis comparing cervical spinal cord vs oculomotor in ALS at cluster level.
  - `2_2_DEG_Differential_vulnerability_scocu_main.R`: Differential vulnerability analysis comparing cervical spinal cord vs oculomotor in ALS at main level.
  - `3_1_DEG_Differential_vulnerability_slocu_cluster.R`: Differential vulnerability analysis comparing lumbar spinal cord vs oculomotor in ALS at cluster level.
  - `3_2_DEG_Differential_vulnerability_slocu_main.R`: Differential vulnerability analysis comparing lumbar spinal cord vs oculomotor in ALS at main level.
  - `4_Differential_vulnerability_PointPlot.R`: Visualization of DEGs across different clusters.
  - `5_meta_analysis.R`: Combine p-values of DEG sets in comparisons 1, 2, and 3.
  - `6_meta_analysis_barplot.R`: Visualization of DEGs from meta-analysis.

## 9. WGCN Analysis
### Methods: Weighted Gene Co-expression Network Analysis
- **bash_files**:
  - `hdWGCNA_build_network_across_status.sh`
  - `hdWGCNA_downstream_across_status.sh`
  - `hdWGCNA_enrichment_analysis_across_status.sh`
  - `hdWGCNA_build_network_across_region.sh`
  - `hdWGCNA_downstream_across_region.sh`
  - `hdWGCNA_enrichment_analysis_across_region.sh`
- **Scripts**:
  - `1_1_hdWGCNA_build_network_across_status.R`: Build a WGCN incorporating cross-status data for each specific cluster.
  - `1_2_hdWGCNA_downstream_across_status.R`: Run downstream analysis on WGCN.
  - `1_3_hdWGCNA_enrichment_analysis_across_status.R`: Gene set enrichment analysis for different WGCNs.
  - `2_1_hdWGCNA_build_network_across_region.R`: Build a WGCN including MN dataset across regions in ALS.
  - `2_2_hdWGCNA_downstream_across_region.R`: Run downstream analysis on WGCN.
  - `2_3_hdWGCNA_enrichment_analysis_across_region.R`: Gene set enrichment analysis for WGCN.

## 10. Cryptic Exon Quantification
### Methods: Detection of Cryptic Exons from Junction Reads
- **bash_files**:
  - `Genome_generate.sh`: Generate the reference for alignment of FASTQ file using STARsolo.
  - `STARsolo_alignment.sh`: Align FASTQ file with created reference genome.
  - `Junction_extract_regtools.sh`: Extract junctions from BAM files.
  - `junctions_per_cell.sh`: Identify cryptic exons from junctions.
- **Scripts**:
  - `1_junctions_per_cell.R`: Identification of cryptic exons from junctions extracted using STAR
