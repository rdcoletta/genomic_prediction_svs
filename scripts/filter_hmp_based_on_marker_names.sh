#!/bin/bash
#PBS -l walltime=4:00:00,nodes=1:ppn=1,mem=150gb
#PBS -o /home/hirschc1/della028/projects/genomic_prediction/simulation/analysis/trait_sim
#PBS -e /home/hirschc1/della028/projects/genomic_prediction/simulation/analysis/trait_sim
#PBS -V
#PBS -N filter_hmp_based_on_marker_names
#PBS -M della028@umn.edu
#PBS -m abe
#PBS -r n

module load R/3.6.0

# go to project folder
cd ~/projects/genomic_prediction/simulation

# filter hapmap
run_pipeline.pl -Xmx150g -importGuess ${IN} \
                -includeSiteNamesInFile ${LIST} \
                -export ${OUT} \
                -exportType HapmapDiploid
