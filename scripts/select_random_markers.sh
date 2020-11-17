#!/bin/bash
#PBS -l walltime=1:00:00,nodes=1:ppn=1,mem=100gb
#PBS -o /home/hirschc1/della028/projects/genomic_prediction/simulation/analysis/trait_sim
#PBS -e /home/hirschc1/della028/projects/genomic_prediction/simulation/analysis/trait_sim
#PBS -V
#PBS -N select_random_markers
#PBS -M della028@umn.edu
#PBS -m abe
#PBS -r n

module load R/3.6.0

# go to project folder
cd ~/projects/genomic_prediction/simulation

# filter hapmap
OPT=$(echo $OPT | tr -d "'")
Rscript scripts/select_random_markers.R ${IN} ${OUT} ${NSAMPLE} ${OPT} --seed=${SEED}
