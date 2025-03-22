#!/bin/bash
#SBATCH --ntasks-per-node=4
#SBATCH --mem-per-cpu=10G
#SBATCH --time=7-00:00
#SBATCH --account=def-grouleau

# Load required modules
module load StdEnv/2023 star/2.7.11b samtools

for sample in "T01499" "T01500" "T01487" "T01488" "T01491" "T01492" "T01479" "T01515" "T01516" "T01497" "T01489" "T01480" "T01498" "T01486" "T01490" "T01477" "T01478" "T01513" "T01514"; do
    regtools junctions extract -b ${sample}_output.barcodes -o ${sample}_output.sj ${sample}Aligned.sortedByCoord.out.bam -s FR
done
