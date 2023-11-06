#!/bin/bash
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=10
#SBATCH --mem=100gb
#SBATCH -J add_mono_reseq-snps
#SBATCH -o /home/hirschc1/della028/projects/genomic_prediction/simulation/analysis/projection_svs-snps/%x_%j.out
#SBATCH -e /home/hirschc1/della028/projects/genomic_prediction/simulation/analysis/projection_svs-snps/%x_%j.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=della028@umn.edu
#SBATCH --no-requeue

# go to project folder
cd ~/projects/genomic_prediction/simulation

module load R/3.6.0

Rscript scripts/add_mono_reseq-snps.R ${CROSS} analysis/projection_svs-snps/usda_rils_SV-SNPchip-polySNPreseq.${CROSS}.projected.sliding-window.hmp.txt data/reseq_snps/widiv_snps_usda_parents.${CROSS}.not-in-SVs.hmp.txt analysis/projection_svs-snps/usda_rils_projected-SVs-SNPs.${CROSS}.hmp.txt
