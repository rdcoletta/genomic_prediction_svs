#!/bin/bash
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=10
#SBATCH --mem=20gb
#SBATCH -J distribution_ld_snps-svs
#SBATCH -o /home/hirschc1/della028/projects/genomic_prediction/simulation/analysis/ld/%x_%j.out
#SBATCH -e /home/hirschc1/della028/projects/genomic_prediction/simulation/analysis/ld/%x_%j.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=della028@umn.edu
#SBATCH --no-requeue

module load R/3.6.0

# go to project folder
cd ~/projects/genomic_prediction/simulation

Rscript scripts/distribution_ld_snps-svs.R analysis/ld/ld_usda_rils_snp-sv_only.chr1.window-${WINDOW}kb.filter-${FILTER}.ld data/usda_SVs_parents.sorted.hmp.txt analysis/ld/window-${WINDOW}\kb_filter-${FILTER}
