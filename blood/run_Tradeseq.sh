#!/bin/bash
#SBATCH -p icelake #-himem #skylake-himem #icelake #skylake #-himem #cclake
#SBATCH -A gottgens-sl2-cpu
#SBATCH -N 1
#SBATCH -n 64
#SBATCH --time 10:00:00
#SBATCH --job-name Trajectory
#SBATCH --output logs/run_tradeseq.txt

Rscript Tradeseq.R