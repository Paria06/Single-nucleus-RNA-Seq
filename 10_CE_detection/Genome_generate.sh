#!/bin/bash


#SBATCH --mem-per-cpu=40G
#SBATCH --time=24:00:00
#SBATCH --account=def-grouleau

module load StdEnv/2023 star/2.7.11b

cd /home/paria/scratch/STARsolo

STAR --runMode genomeGenerate \
     --runThreadN 4 \
     --genomeDir /home/paria/scratch/STARsolo/genomeDir \
     --genomeFastaFiles /home/paria/scratch/STARsolo/genome/genome.fa \
     --sjdbGTFfile /home/paria/scratch/STARsolo/genome/genes.gtf \
     --sjdbOverhang 99
