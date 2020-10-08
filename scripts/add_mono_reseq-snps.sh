#!/bin/bash
#PBS -l walltime=2:00:00,nodes=1:ppn=10,mem=100gb
#PBS -o /home/hirschc1/della028/projects/genomic_prediction/simulation/analysis/projection_svs-snps
#PBS -e /home/hirschc1/della028/projects/genomic_prediction/simulation/analysis/projection_svs-snps
#PBS -V
#PBS -N add_mono_reseq-snps_${CROSS}
#PBS -M della028@umn.edu
#PBS -m abe
#PBS -r n

# go to project folder
cd ~/projects/genomic_prediction/simulation

module load R

Rscript scripts/add_mono_reseq-snps.R ${CROSS} analysis/projection_svs-snps/usda_rils_SV-SNPchip-polySNPreseq.${CROSS}.projected.sliding-window.hmp.txt data/reseq_snps/biomAP_v_B73_SNPMatrix_parents.${CROSS}.not-in-SVs.hmp.txt analysis/projection_svs-snps/usda_rils_projected-SVs-SNPs.${CROSS}.hmp.txt
