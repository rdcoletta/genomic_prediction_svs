#!/bin/bash
#PBS -l walltime=5:00:00,nodes=1:ppn=10,mem=110gb
#PBS -o /home/hirschc1/della028/projects/genomic_prediction/simulation/analysis/projection_reseq-snps
#PBS -e /home/hirschc1/della028/projects/genomic_prediction/simulation/analysis/projection_reseq-snps
#PBS -V
#PBS -N add_snps_not_in_cross
#PBS -M della028@umn.edu
#PBS -m abe
#PBS -r n

cd ~/projects/genomic_prediction/simulation

module load R

Rscript scripts/add_snps_not_in_cross.R analysis/projection_reseq-snps
