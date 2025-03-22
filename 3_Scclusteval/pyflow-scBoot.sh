#! /bin/bash

module load r/4.1.2
source /home/paria/my_softwares/ENV/bin/activate
## make a folder for bsub_log if not exist
folder="sbatch_log"

if [ ! -d "$folder" ]; then
    mkdir -p "$folder"
fi

snakemake -j 99999 \
    -k \
    --rerun-incomplete \
    --latency-wait 10000 \
    --cluster-config cluster.json \
    --cluster "sbatch --mem-per-cpu {cluster.mem} --cpus-per-task {cluster.n} --time {cluster.time} --job-name {cluster.name} --output {cluster.output} --error {cluster.error} --account def-grouleau" \
    "$@" 


