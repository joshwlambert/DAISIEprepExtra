#!/bin/bash
#SBATCH --time=1-23:00:00
#SBATCH --partition=gelifes
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=test
#SBATCH --output=/data/p287218/DAISIEprepExtra/logs/test%a.log
#SBATCH --array=1-10
#SBATCH --mem=5GB

ml R
Rscript /data/p287218/DAISIEprepExtra/scripts/test.R ${SLURM_ARRAY_TASK_ID}
