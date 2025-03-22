#!/bin/bash
#SBATCH --time=9-00:00:00
#SBATCH --ntasks=2
#SBATCH --mem=8G
export PATH=/home/paria/cellranger-6.1.2:$PATH
source ~/.bashrc
for i in T01471 T01472 T01473 T01474 T01477 T01478 T01486 T01489 T01490 T01493 T01494 T01495 T01496 T01497 T01498 T01513 T01514 T01499 T01500 T01487 T01488 T01491 T01492 T01479 T01480 T01515 T01516 T01517 T01551 T01552 T01553 T01554 T01555 T01556 T01557 T01558 T01559 T01560
do
#Cellranger count version_6
      cellranger count --id=${i} \
      --transcriptome=/home/paria/refdata-gex-GRCh38-2020-A \
      --fastqs=/home/paria/runs/paria/als-project/raw-data/fastqs/${i} \
      --include-introns \
      --expect-cells=6000
Done

