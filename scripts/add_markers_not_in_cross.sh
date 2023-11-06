#!/bin/bash
#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=5
#SBATCH --mem=200gb
#SBATCH -J add_markers_not_in_cross
#SBATCH -o /home/hirschc1/della028/projects/genomic_prediction/simulation/analysis/projection_svs-snps/%x_%j.out
#SBATCH -e /home/hirschc1/della028/projects/genomic_prediction/simulation/analysis/projection_svs-snps/%x_%j.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=della028@umn.edu
#SBATCH --no-requeue

cd ~/projects/genomic_prediction/simulation

module load R/3.6.0

Rscript scripts/add_markers_not_in_cross.R analysis/projection_svs-snps
