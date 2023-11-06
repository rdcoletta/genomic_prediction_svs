#!/bin/bash
#SBATCH --time=5:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=100gb
#SBATCH -J tassel_qc-summary
#SBATCH -o /home/hirschc1/della028/projects/genomic_prediction/simulation/analysis/projection_svs-snps/%x_%j.out
#SBATCH -e /home/hirschc1/della028/projects/genomic_prediction/simulation/analysis/projection_svs-snps/%x_%j.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=della028@umn.edu
#SBATCH --no-requeue

# go to project folder
cd ~/projects/genomic_prediction/simulation

run_pipeline.pl -Xmx100g -importGuess data/usda_rils_projected-SVs-SNPs.chr${CHR}.hmp.txt -GenotypeSummaryPlugin -endPlugin -export analysis/qc/tassel_usda-SNPs-SVs_chr${CHR}_summary
