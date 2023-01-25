#!/bin/bash
#PBS -o /home/hirschc1/della028/projects/genomic_prediction/simulation/analysis/ld
#PBS -e /home/hirschc1/della028/projects/genomic_prediction/simulation/analysis/ld
#PBS -V
#PBS -N plot_ld_decay
#PBS -M della028@umn.edu
#PBS -m abe
#PBS -r n

module load R

# go to project folder
cd ~/projects/genomic_prediction/simulation

# plot ld decay
Rscript scripts/plot_ld_decay.R ${IN} ${OUT} ${OPT}
