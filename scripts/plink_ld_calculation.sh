#!/bin/bash
#PBS -o /home/hirschc1/della028/projects/genomic_prediction/simulation/analysis/ld
#PBS -e /home/hirschc1/della028/projects/genomic_prediction/simulation/analysis/ld
#PBS -V
#PBS -N plink_ld_calculation_${CHR}_${WINDOW}_${FILTER}
#PBS -M della028@umn.edu
#PBS -m abe
#PBS -r n

# go to project folder
cd /scratch.global/della028/hirsch_lab/genomic_prediction/ld

# calculate LD
plink --file usda_rils_projected-SVs-SNPs.chr${CHR}.plk --make-founders --r2 gz dprime with-freqs --ld-window-r2 0 --ld-window ${WINDOW}000 --ld-window-kb ${WINDOW} --geno ${FILTER} --out usda_rils_projected-SVs-SNPs.chr${CHR}.window-${WINDOW}kb.filter-${FILTER}
