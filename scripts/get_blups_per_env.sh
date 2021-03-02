#!/bin/bash
#SBATCH --time=0:30:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=5gb
#SBATCH -J get_blups_per_env
#SBATCH -o /home/hirschc1/della028/projects/genomic_prediction/simulation/analysis/trait_sim/MSI_dump/%x_%j.out
#SBATCH -e /home/hirschc1/della028/projects/genomic_prediction/simulation/analysis/trait_sim/MSI_dump/%x_%j.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=della028@umn.edu
#SBATCH --no-requeue

module load R/3.6.0

# go to project folder
cd ~/projects/genomic_prediction/simulation

echo ${FOLDER}

for POP in $(seq 1 ${POPS}); do
  # get simulated traits file
  FILE=$(ls ${FOLDER}/pop${POP} | grep Simulated_Data)
  # get blups per env
  Rscript scripts/get_blups_per_env.R ${FOLDER}/pop${POP}/${FILE} ${FOLDER}/pop${POP}
done
