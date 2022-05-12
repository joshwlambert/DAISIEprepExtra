#!/bin/bash
#SBATCH --time=4-23:00:00
#SBATCH --partition=gelifes
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=performance
#SBATCH --output=/data/p287218/DAISIEmainland/logs/performance.log
#SBATCH --mem=5GB

ml R
Rscript /data/p287218/DAISIEprepExtra/scripts/performance.R
