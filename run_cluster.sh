#!/bin/bash
#SBATCH -p cclake
#SBATCH -A gottgens-sl2-cpu
#SBATCH -N 1
#SBATCH -n 6
#SBATCH --time 10:00:00
#SBATCH --job-name snakemake
#SBATCH --output snakemake-log-%J.txt
snakemake --cores 12 -j 99 --latency-wait 90 -p --cluster "sbatch -p cclake -A gottgens-sl2-cpu -n {threads} --time 10:00:00 --mem {resources.mem_mb}M"