#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --mem=120gb
#SBATCH -J project_svs-snps
#SBATCH -o /home/hirschc1/della028/projects/genomic_prediction/simulation/analysis/projection_svs-snps/%x_%j.out
#SBATCH -e /home/hirschc1/della028/projects/genomic_prediction/simulation/analysis/projection_svs-snps/%x_%j.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=della028@umn.edu
#SBATCH --no-requeue

# go to project folder
cd ~/projects/genomic_prediction/simulation

module load R/3.6.0

Rscript scripts/sliding_window_reseq-snps.R ${CROSS} analysis/projection_svs-snps/usda_rils_SV-SNPchip-polySNPreseq.${CROSS}.projected.hmp.txt data/hapmap_by_cross/usda_parents_SV-SNPchip-polySNPreseq.${CROSS}.hmp.txt --window_size=45 --window_step=1 --min_snps_per_window=15
