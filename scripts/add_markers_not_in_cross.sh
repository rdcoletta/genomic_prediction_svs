#!/bin/bash
#PBS -l walltime=10:00:00,nodes=1:ppn=5,mem=200gb
#PBS -o /home/hirschc1/della028/projects/genomic_prediction/simulation/analysis/projection_svs-snps
#PBS -e /home/hirschc1/della028/projects/genomic_prediction/simulation/analysis/projection_svs-snps
#PBS -V
#PBS -N add_markers_not_in_cross
#PBS -M della028@umn.edu
#PBS -m abe
#PBS -r n

cd ~/projects/genomic_prediction/simulation

module load R

Rscript scripts/add_markers_not_in_cross.R analysis/projection_svs-snps
