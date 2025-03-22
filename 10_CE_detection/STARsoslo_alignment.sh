#!/bin/bash
#SBATCH --ntasks-per-node=4
#SBATCH --mem-per-cpu=10G
#SBATCH --time=7-00:00
#SBATCH --account=def-grouleau

# Load required modules
module load StdEnv/2023 star/2.7.11b samtools

# Base directories (adjust these paths as needed)
FASTQ_BASE="/home/paria/scratch/als-project/fastqs"
STAR_OUT_BASE="/home/paria/scratch/STARsolo"
GENOME_DIR="/home/paria/scratch/STARsolo/genomeDir"
WHITELIST="${STAR_OUT_BASE}/3M-february-2018_TRU.txt"

# Define an array of sample names
SAMPLES=("T01499" "T01500" "T01487" "T01488" "T01491" "T01492" "T01479" "T01515" "T01516" "T01497" "T01489" "T01480" "T01498" "T01486" "T01490" "T01477" "T01478" "T01513" "T01514")

for sample in "${SAMPLES[@]}"; do
    echo "Processing sample: ${sample}"

    # Define sample-specific directories and output prefix
    SAMPLE_FASTQ_DIR="${FASTQ_BASE}/${sample}"
    SAMPLE_OUT_DIR="${STAR_OUT_BASE}/${sample}"
    mkdir -p "${SAMPLE_OUT_DIR}"  # Create output directory if it doesn't exist
    OUT_PREFIX="${SAMPLE_OUT_DIR}/${sample}_"

    # Find all R2 and R1 files for the sample. The pattern assumes files contain _R2_ or _R1_ in their names.
    # Ensure all R2 files come first.
    R2_FILES=( $(ls ${SAMPLE_FASTQ_DIR}/*_R2_*.fastq.gz 2>/dev/null | sort) )
    R1_FILES=( $(ls ${SAMPLE_FASTQ_DIR}/*_R1_*.fastq.gz 2>/dev/null | sort) )

    # Check that we have at least one file in each group
    if [ ${#R2_FILES[@]} -eq 0 ] || [ ${#R1_FILES[@]} -eq 0 ]; then
        echo "Error: Could not find both R2 and R1 files for sample ${sample}."
        continue
    fi

    # Combine the file lists in the required order: first all R2, then all R1.
    # When expanding arrays, quoting them with "${ARRAY[@]}" preserves the separation.
    echo "Using the following R2 files:"
    printf '%s\n' "${R2_FILES[@]}"
    echo "Using the following R1 files:"
    printf '%s\n' "${R1_FILES[@]}"

    # Run STARsolo with the combined FASTQ files.
    STAR --soloType CB_UMI_Simple \
         --outSAMstrandField intronMotif \
         --outSAMattributes CR UR CB UB GX GN NH HI AS nM AS sM \
         --outSAMtype BAM SortedByCoordinate \
         --genomeDir "${GENOME_DIR}" \
         --soloCBwhitelist "${WHITELIST}" \
         --readFilesCommand zcat \
         --soloBarcodeReadLength 0 \
         --readFilesIn "${R2_FILES[@]}" "${R1_FILES[@]}" \
         --soloFeatures GeneFull SJ \
         --outFileNamePrefix "${OUT_PREFIX}" \
         --clipAdapterType CellRanger4 \
	     --outFilterScoreMin 30 \
	     --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
	     --soloUMIfiltering MultiGeneUMI_CR \
	     --soloUMIdedup 1MM_CR \
	     --soloCBstart 1 \
         --soloCBlen 16 \
         --soloUMIstart 17 \
         --soloUMIlen 12


done

