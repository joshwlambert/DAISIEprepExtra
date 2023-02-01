#!/bin/bash
#SBATCH --time=4-23:00:00
#SBATCH --partition=gelifes
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=performance
#SBATCH --output=/data/p287218/DAISIEprepExtra/logs/performance%a.log
#SBATCH --array=1-60
#SBATCH --mem=5GB

ml R
Rscript /data/p287218/DAISIEprepExtra/inst/scripts/performance_parallel.R ${SLURM_ARRAY_TASK_ID}
