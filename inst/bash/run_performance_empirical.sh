#!/bin/bash
#SBATCH --time=4-23:00:00
#SBATCH --partition=gelifes
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=performance_empirical
#SBATCH --output=/data/p287218/DAISIEprepExtra/logs/performance_empirical.log
#SBATCH --mem=5GB

ml R
Rscript /data/p287218/DAISIEprepExtra/inst/scripts/performance_empirical.R
