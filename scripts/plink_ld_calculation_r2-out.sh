#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=150gb
#SBATCH -J plink_ld_calculation_r2-out_${CHR}
#SBATCH -o /home/hirschc1/della028/projects/genomic_prediction/simulation/analysis/ld/%x_%j.out
#SBATCH -e /home/hirschc1/della028/projects/genomic_prediction/simulation/analysis/ld/%x_%j.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=della028@umn.edu
#SBATCH --no-requeue

# go to project folder
cd /scratch.global/della028/hirsch_lab/genomic_prediction/ld

# calculate LD
plink --file usda_rils_projected-SVs-SNPs.chr${CHR}.plk --make-founders --r2 gz dprime with-freqs --ld-window-r2 ${R2OUT} --ld-window ${WINDOW}000 --ld-window-kb ${WINDOW} --geno ${FILTER} --out usda_rils_projected-SVs-SNPs.chr${CHR}.window-${WINDOW}kb.filter-${FILTER}.r2-${R2OUT}
