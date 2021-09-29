#!/bin/bash
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=5gb
#SBATCH -J anova_sim_traits
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
  # run anova
  Rscript scripts/anova_sim_traits.R ${FOLDER}/pop${POP}
  # plot pve
  Rscript scripts/plot_sim_pve_qtns.R ${FOLDER}/pop${POP} ${SVS}
done
