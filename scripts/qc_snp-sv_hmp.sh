#!/bin/bash
#PBS -l walltime=2:00:00,nodes=1:ppn=1,mem=80gb
#PBS -o /home/hirschc1/della028/projects/genomic_prediction/simulation/analysis/trait_sim/MSI_dump
#PBS -e /home/hirschc1/della028/projects/genomic_prediction/simulation/analysis/trait_sim/MSI_dump
#PBS -V
#PBS -N qc_snp-sv_hmp
#PBS -M della028@umn.edu
#PBS -m abe
#PBS -r n

module load R/3.6.0

# go to project folder
cd ~/projects/genomic_prediction/simulation

# qc hapmap file
Rscript scripts/qc_snp-sv_hmp.R ${HMP} ${SVS} ${FOLDER}
