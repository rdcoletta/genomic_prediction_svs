#!/bin/bash
#PBS -l walltime=3:00:00,nodes=1:ppn=1,mem=100gb
#PBS -o /home/hirschc1/della028/projects/genomic_prediction/simulation/data/reseq_snps
#PBS -e /home/hirschc1/della028/projects/genomic_prediction/simulation/data/reseq_snps
#PBS -V
#PBS -N merge_reseq-chip_projected_crosses_${CHR}
#PBS -M della028@umn.edu
#PBS -m abe
#PBS -r n

# go to project folder
cd ~/projects/genomic_prediction/simulation

module load R

Rscript scripts/merge_reseq-chip_projected_crosses.R data/usda_SNPs-SVs_rils.not-in-SVs.projected.chr${CHR}.hmp.txt data/reseq_snps/biomAP_rils_SNPs-reseq.projected.chr${CHR}.hmp.txt
