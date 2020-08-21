#!/bin/bash
#PBS -l walltime=12:00:00,nodes=1:ppn=1,mem=150gb
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
plink --file usda_SNPs-SVs_rils.not-in-SVs.projected.chr${CHR}.reseq-SNPs.plk --make-founders --r2 gz --ld-window-r2 0 --ld-window ${WINDOW}000 --ld-window-kb ${WINDOW} --geno ${FILTER} --out usda_SNPs-SVs_rils.not-in-SVs.projected.chr${CHR}.reseq-SNPs.window-${WINDOW}kb.filter-${FILTER}
