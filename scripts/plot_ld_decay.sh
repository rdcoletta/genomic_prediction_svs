#!/bin/bash
#SBATCH --time=5:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=120gb
#SBATCH -J plot_ld_decay
#SBATCH -o /home/hirschc1/della028/projects/genomic_prediction/simulation/analysis/ld/%x_%j.out
#SBATCH -e /home/hirschc1/della028/projects/genomic_prediction/simulation/analysis/ld/%x_%j.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=della028@umn.edu
#SBATCH --no-requeue

module load R/3.6.0

# go to project folder
cd ~/projects/genomic_prediction/simulation

# plot ld decay
Rscript scripts/plot_ld_decay.R ${IN} ${OUT} --unequal-windows
