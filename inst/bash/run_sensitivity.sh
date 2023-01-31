#!/bin/bash
#SBATCH --time=1-23:00:00
#SBATCH --partition=gelifes
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=sensitivity
#SBATCH --output=/data/p287218/DAISIEprepExtra/logs/sensitivity.log
#SBATCH --mem=5GB

ml R
Rscript /data/p287218/DAISIEprepExtra/inst/scripts/sensitivity.R
