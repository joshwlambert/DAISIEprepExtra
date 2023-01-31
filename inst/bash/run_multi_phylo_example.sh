#!/bin/bash
#SBATCH --time=4-23:00:00
#SBATCH --partition=gelifes
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=multi_phylo_example
#SBATCH --output=/data/p287218/DAISIEmainland/logs/multi_phylo_example.log
#SBATCH --mem=5GB

ml R
Rscript /data/p287218/DAISIEprepExtra/inst/scripts/single_phylo_example.R
