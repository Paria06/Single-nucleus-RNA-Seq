


├── 1_Preprocessing								Methods sections: Generation of expression matrix and preprocessing of single nuclei dataset
│   ├── 1_cellranger.sh								FASTQ files are mapped to the GRCh38 reference genome to get the gene read count using the Cellranger count					             		                      
│   │
├── 2_Merging_filtering_normalization_define_metadata_integration		Methods sections: Dimension reduction, integration, and clustering 
│   ├── bash_files								Contains the bash scripts used to submit R scripts. 
│   │   ├── Merging_filtering_normalization_metadata.sh
│   │   ├── Integration_rPCA.sh
│   │   ├── Integration_rPCA.sh
│   └── Scripts
│       ├── 1_Merging_filtering_normalization_metadata.R			Create Seurat object, merge, and filtering dataset based on the Satijalab vignettes 
│       ├── 2_Integration_rPCA.R						Normalization, finding variable features, integration, data scaling, and dimension reduction  
│       └── 3_QC_plot_before_after_integration.R				The impact of different variables was evaluated plotting umap before and after integration 
│   
├── 3_Scclusteval_analysis							Methods sections: Dimension reduction, integration, and clustering 	
│   └── pyflow_seurat_parameter_custom						Files corresponding to the scclusteval Snakemake workflow
│       ├── cluster.json							Test the stability of clusters in different combinations of clustering parameters
│       ├── config.yaml								Test the stability of clusters in different combinations of clustering parameters
│       ├── pyflow-scBoot.sh							Test the stability of clusters in different combinations of clustering parameters
│       ├── scripts  								Includes scripts to process subsample and full sample
│       │   ├── gather_fullsample.R						Test the stability of clusters in different combinations of clustering parameters
│       │   ├── gather_subsample.R						Test the stability of clusters in different combinations of clustering parameters	
│       │   ├── preprocess.R							Test the stability of clusters in different combinations of clustering parameters
│       │   └── subsample.R							Test the stability of clusters in different combinations of clustering parameters
│       └── Snakefile								Test the stability of clusters in different combinations of clustering parameters
│
├── 4_Cluster_evaluation_annotation						Methods sections: Quality control of clustering;Cluster annotation using automated and manual approach
│   ├── bash_files								Contains the bash scripts used to submit R scripts.
│   │   ├── Evaluation_of_scclusteval_output.sh
│   │   ├── Doublet_removal.sh
│   │   ├── Cluster_annotation_singleR_spinalcord.sh
│   │   ├── Cluster_annotation_singleR_medulla.sh  
│   │   ├── Cluster_annotation_singleR_oculomotor.sh
│   │   ├── Cluster_markers_manual.sh
│   │   └── Cluster_QC.sh
│   ├── Scripts
│   │   ├── 1_Evaluation_of_scclusteval_output.R				Visualization of scclusteval output to determine the optimal combination
│   │   ├── 2_Doublet_removal.R                                         	Removal of doublets using scDblFinder packge
│   │   ├── 3_Cluster_annotation_singleR_spinalcord.R				Alignment of the dataset with a spinal cord dataset, automated approach using SingleR
│   │   ├── 4_Cluster_annotation_singleR_medulla.R				Alignment of the dataset with a midbrain dataset, automated approach using SingleR
│   │   ├── 5_Cluster_annotation_singleR_oculomotor.R				Alignment of the dataset with a midbrain dataset, automated approach using SingleR
│   │   ├── 6_Cluster_markers_manual.R						Generating clusters markers using Presto
│   │   └── 7_Cluster_QC.R							Evaluating the QC parameters of clusters
│   │
├── 5_Cluster_visualization 							Methods sections: Cluster annotation using automated and manual approach; Cell type proportions comparison
│   ├── bash_files								Contains the bash scripts used to submit R scripts.
│   │   ├── Clusters_visualization.sh  
│   │   └── Cluster_celltype_proportion.sh
│   └── Scripts
│       ├── 1_Clusters_visualization.R                                  	Cluster marker visualization
│       └── 2_Cluster_celltype_proportion.R					Proportion of clusters and main cell types across different variables using Speckle package
│   								
├── 6_Differential_expression							Methods sections: Differential Gene Expression analysis using DESeq2
│   ├── bash_files								Contains the bash scripts used to submit R scripts.
│   │   ├── DEG_cluster_id_ALS_vs_con_oculomotor.sh  
│   │   ├── DEG_main_ALS_vs_con_oculomotor.sh
│   │   ├── DEG_cluster_id_ALS_vs_con_medulla.sh 
│   │   ├── DEG_main_ALS_vs_con_medulla.sh
│   │   ├── DEG_cluster_id_ALS_vs_con_cervicalspinalcord.sh
│   │   ├── DEG_main_ALS_vs_con_cervicalspinalcord.sh 
│   │   └── DE_PointPlot.sh
│   └── Scripts
│       ├── 1_1_DEG_cluster_id_ALS_vs_con_oculomotor.R				DE analysis using DESeq2 comparing ALS vs control in oculomotor at cluster level
│       ├── 1_2_DEG_main_ALS_vs_con_oculomotor.R				DE analysis using DESeq2 comparing ALS vs control in oculomotor at main level
│       ├── 2_1_DEG_cluster_id_ALS_vs_con_medulla.R				DE analysis using DESeq2 comparing ALS vs control in medulla at cluster level
│       ├── 2_2_DEG_main_ALS_vs_con_medulla.R                           	DE analysis using DESeq2 comparing ALS vs control in medulla at main level
│       ├── 3_1_DEG_cluster_id_ALS_vs_con_cervicalspinalcord.R          	DE analysis using DESeq2 comparing ALS vs control in cervical spinal cord at cluster level
│       ├── 3_2_DEG_main_ALS_vs_con_cervicalspinalcord.R                	DE analysis using DESeq2 comparing ALS vs control in cervical spinal cord at main level
│       ├── 4_1_DEG_cluster_id_ALS_vs_con_lumbarspinalcord.R          		DE analysis using DESeq2 comparing ALS vs control in lumbar spinal cord at cluster level
│       ├── 4_2_DEG_main_ALS_vs_con_lumbarspinalcord.R                		DE analysis using DESeq2 comparing ALS vs control in lumbar spinal cord at main level
│       └── 5_DE_PointPlot.R				         		Visualization of DEGs across different clusters
│ 
├── 7_RRHO_analysis								Methods sections: Rank-Rank Hypergeometric Overlap without Thresholds
│   ├── bash_files								Contains the bash scripts used to submit R scripts.
│   │   ├── RROH_DE_correlation_ocumed.sh  
│   │   ├── RROH_DE_correlation_ocusc.sh 
│   │   └── RROH_DE_correlation_ocusl.sh
│   └── Scripts
│       ├── 1_RROH_DE_correlation_ocumed.R                              	Correlation of oculomotor DEG and medulla DEG without significance threshold
│       ├── 2_RROH_DE_correlation_ocusc.R                               	Correlation of oculomotor DEG and cervical spinal cord DEG without significance threshold
│       └── 3_RROH_DE_correlation_ocusl.R					Correlation of oculomotor DEG and lumbar spinal cord DEG without significance threshold
│   
├── 8_Differential_vulnerability_meta_analysis					Methods sections: Differential Gene Expression analysis using DESeq2
│   ├── bash_files								Contains the bash scripts used to submit R scripts.
│   │   ├── DEG_Differential_vulnerability_medocu_cluster.sh  
│   │   ├── DEG_Differential_vulnerability_medocu_main.sh
│   │   ├── DEG_Differential_vulnerability_scocu_cluster.sh 
│   │   ├── DEG_Differential_vulnerability_scocu_main.sh
│   │   ├── DEG_Differential_vulnerability_slocu_cluster.sh
│   │   ├── DEG_Differential_vulnerability_slocu_main.sh 
│   │   ├── Differential_vulnerability_PointPlot.sh  
│   │   ├── meta_analysis.sh 
│   │   └── meta_analysis_barplot.sh
│   └── Scripts
│       ├── 1_1_DEG_Differential_vulnerability_medocu_cluster.R			Differential vulnerability analysis using DESeq2 comparing medulla vs oculomotor in ALS at cluster level
│       ├── 1_2_DEG_Differential_vulnerability_medocu_main.R			Differential vulnerability analysis using DESeq2 comparing medulla vs oculomotor in ALS at main level
│       ├── 2_1_DEG_Differential_vulnerability_scocu_cluster.R			Differential vulnerability analysis using DESeq2 comparing cervical spinal cord vs oculomotor in ALS at cluster level
│       ├── 2_2_DEG_Differential_vulnerability_scocu_main.R             	Differential vulnerability analysis using DESeq2 comparing cervical spinal cord vs oculomotor in ALS at main level
│       ├── 3_1_DEG_Differential_vulnerability_slocu_cluster.R          	Differential vulnerability analysis using DESeq2 comparing lumbar spinal cord vs oculomotor in ALS at cluster level
│       ├── 3_2_DEG_Differential_vulnerability_slocu_main.R             	Differential vulnerability analysis using DESeq2 comparing lumbar spinal cord vs oculomotor in ALS at main level
│       ├── 4_Differential_vulnerability_PointPlot.R                    	Visualization of DEGs across different clusters, Figure 2d-g
│       ├── 5_meta_analysis.R							Combining the p values of DEG sets in 1, 2, and 3 comparisons 
│       └── 6_meta_analysis_barplot.R		                                Visualization of DEGs of meta-analysis
│
├── 9_WGCN_Analysis								Methods sections: Weighted Gene Co-expression Network Analysis	
│   ├── bash_files								Contains the bash scripts used to submit R scripts.
│   │   ├── hdWGCNA_build_network_across_status.sh
│   │   ├── hdWGCNA_downstream_across_status.sh
│   │   ├── hdWGCNA_enrichment_analysis_across_status.sh
│   │   ├── hdWGCNA_build_network_across_region.sh
│   │   ├── hdWGCNA_downstream_across_region.sh
│   │   └── hdWGCNA_enrichment_analysis_across_region.sh
│   └── Scripts
│      ├── 1_1_hdWGCNA_build_network_across_status.R				Build a WGCN incorporating cross-status data for each specific cluster
│      ├── 1_2_hdWGCNA_downstream_across_status.R				Run downstream analysis on WGCN 
│      ├── 1_3_hdWGCNA_enrichment_analysis_across_status.R			Gene set enrichment analysis for different WGCN
│      ├── 2_1_hdWGCNA_build_network_across_region.R				Build a WGCN including MN dataset across region in ALS
│      ├── 2_2_hdWGCNA_downstream_across_region.R				Run downstream analysis on WGCN												
│      └── 2_3_hdWGCNA_enrichment_analysis_across_region.R              	Gene set enrichment analysis for the WGCN
│		    																  		
├── 10_CrypticExon_quantification						Methods sections: Detection of Cryptic Exons (CE) from Junction Reads	
│   ├── bash_files								Contains the bash scripts.
│   │   ├── Genome_generate.sh                                          	Generating the reference for alignment of FASTQ file using STARsolo
│   │   ├── STARsoslo_alignment.sh                                      	Alignment of FASTQ file with the created reference genome
│   │   ├── Junction_extract_regtools.sh                                        Extraction of junctions from the BAM files
│   │   └── junctions_per_cell.sh
│   └── Script													
│      └── 1_junctions_per_cell.R                  				Identification of CEs from junction extracted in STARsolo and Regtools, formatting, and cell type mapping
│                                             
└── README

