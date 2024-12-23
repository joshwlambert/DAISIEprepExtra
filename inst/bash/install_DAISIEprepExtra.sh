#!/bin/bash
#SBATCH --time=1-23:00:00
#SBATCH --partition=gelifes
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --job-name=install_DAISIEprepExtra
#SBATCH --output=/data/p287218/DAISIEprepExtra/install_DAISIEprepExtra.log
#SBATCH --mem=5GB

mkdir -p logs
mkdir -p results
ml R
Rscript -e "install.packages('renv')"
Rscript -e "renv::restore()"
Rscript -e "remotes::install_github('joshwlambert/DAISIEprepExtra')"
