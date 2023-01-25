#!/bin/bash
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=10
#SBATCH --mem=150gb
#SBATCH -J plink_ld_calculation_random
#SBATCH -o /home/hirschc1/della028/projects/genomic_prediction/simulation/analysis/ld/%x_%j.out
#SBATCH -e /home/hirschc1/della028/projects/genomic_prediction/simulation/analysis/ld/%x_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=della028@umn.edu
#SBATCH --no-requeue

# go to project folder
cd /scratch.global/della028/hirsch_lab/genomic_prediction/ld

# calculate LD for each rep
for rep in {1..10}; do
  plink --file usda_rils_projected-SVs-SNPs.poly.chr${CHR}.rep${rep}.plk --make-founders --r2 gz dprime with-freqs --ld-window-r2 0 --ld-window ${WINDOW}000 --ld-window-kb ${WINDOW} --geno ${FILTER} --out usda_rils_projected-SVs-SNPs.poly.chr${CHR}.rep${rep}.window-${WINDOW}kb.filter-${FILTER} &
done
wait
