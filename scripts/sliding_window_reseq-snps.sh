#!/bin/bash
#PBS -l walltime=12:00:00,nodes=1:ppn=20,mem=120gb
#PBS -o /home/hirschc1/della028/projects/genomic_prediction/simulation/analysis/projection_reseq-snps
#PBS -e /home/hirschc1/della028/projects/genomic_prediction/simulation/analysis/projection_reseq-snps
#PBS -V
#PBS -N sliding_window_reseq-snps_${CROSS}
#PBS -M della028@umn.edu
#PBS -m abe
#PBS -r n

# go to project folder
cd ~/projects/genomic_prediction/simulation

module load R

Rscript scripts/sliding_window_reseq-snps.R ${CROSS} analysis/projection_reseq-snps/biomAP_rils_SNPs-reseq_SNP-chip.${CROSS}.poly.projected.hmp.txt data/reseq_snps/${CROSS}/biomAP_parents_SNPs-reseq_SNP-chip.${CROSS}.poly.hmp.txt --window_size=45 --window_step=1 --min_snps_per_window=15
