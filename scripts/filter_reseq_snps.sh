#!/bin/bash
#PBS -l walltime=3:00:00,nodes=1:ppn=1,mem=50gb
#PBS -o /home/hirschc1/della028/projects/genomic_prediction/simulation/analysis
#PBS -e /home/hirschc1/della028/projects/genomic_prediction/simulation/analysis
#PBS -V
#PBS -N filter_reseq_snps_${CROSS}
#PBS -M della028@umn.edu
#PBS -m abe
#PBS -r n

module load R

# go to project folder
cd ~/projects/genomic_prediction/simulation

# remove SNPs inside deletions
Rscript scripts/remove_SNPs_within_dels.R data/reseq_snps/biomAP_v_B73_SNPMatrix_parents.${CROSS}.hmp.txt data/hapmap_by_cross/usda_SVs_parents.sorted.${CROSS}.hmp.txt

# keep only polymorphic
Rscript scripts/keep_only_poly_snps.R ${CROSS} data/reseq_snps/biomAP_v_B73_SNPMatrix_parents.${CROSS}.not-in-SVs.hmp.txt
