#!/bin/bash
#PBS -l walltime=3:00:00,nodes=1:ppn=1,mem=40gb
#PBS -o /home/hirschc1/della028/projects/genomic_prediction/simulation/analysis/projection_reseq-snps
#PBS -e /home/hirschc1/della028/projects/genomic_prediction/simulation/analysis/projection_reseq-snps
#PBS -V
#PBS -N project_reseq-snps_${CROSS}
#PBS -M della028@umn.edu
#PBS -m abe
#PBS -r n

# go to project folder
cd ~/projects/genomic_prediction/simulation/data/reseq_snps

# create donor files
run_pipeline.pl -Xmx10g -FILLINFindHaplotypesPlugin -hmp ${CROSS}/biomAP_parents_SNPs-reseq_SNP-chip.${CROSS}.poly.hmp.txt -o ~/projects/genomic_prediction/simulation/analysis/projection_reseq-snps/donors_${CROSS} -hapSize ${HAPSIZE} -minTaxa 1

# project reseq snps
run_pipeline.pl -Xmx10g -FILLINImputationPlugin -hmp ${CROSS}/biomAP_rils_SNPs-reseq_SNP-chip.${CROSS}.poly.hmp.txt -d ~/projects/genomic_prediction/simulation/analysis/projection_reseq-snps/donors_${CROSS} -o ~/projects/genomic_prediction/simulation/analysis/projection_reseq-snps/biomAP_rils_SNPs-reseq_SNP-chip.${CROSS}.poly.projected.hmp.txt -hapSize ${HAPSIZE} -hybNN false -accuracy
