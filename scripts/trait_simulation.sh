#!/bin/bash
#PBS -l walltime=2:00:00,nodes=1:ppn=1,mem=100gb
#PBS -o /home/hirschc1/della028/projects/genomic_prediction/simulation/analysis/trait_sim
#PBS -e /home/hirschc1/della028/projects/genomic_prediction/simulation/analysis/trait_sim
#PBS -V
#PBS -N trait_simulation
#PBS -M della028@umn.edu
#PBS -m abe
#PBS -r n

module load R/3.6.0

# go to project folder
cd ~/projects/genomic_prediction/simulation

# qc hapmap file
Rscript scripts/trait_simulation.R ${HMP} ${SVS} ${FOLDER} \
                                   --causal-variant=${VAR} \
                                   --rep=${REPS} \
                                   --ntraits=${ENVS} \
                                   --h2=${H2} \
                                   --model=${MODEL} \
                                   --add-QTN-num=${QTN} \
                                   --add-effect=${EFFECT} \
                                   --architecture=pleiotropic \
                                    --seed=${SEED}
