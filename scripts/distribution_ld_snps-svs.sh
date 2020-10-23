#!/bin/bash
#PBS -l walltime=2:00:00,nodes=1:ppn=10,mem=10gb
#PBS -o /home/hirschc1/della028/projects/genomic_prediction/simulation/analysis/ld
#PBS -e /home/hirschc1/della028/projects/genomic_prediction/simulation/analysis/ld
#PBS -V
#PBS -N distribution_ld_snps-svs_${WINDOW}_${FILTER}
#PBS -M della028@umn.edu
#PBS -m abe
#PBS -r n

module load R

# go to project folder
cd ~/projects/genomic_prediction/simulation

Rscript scripts/distribution_ld_snps-svs.R analysis/ld/ld_usda_rils_snp-sv_only.chr1.window-${WINDOW}kb.filter-${FILTER}.ld data/usda_SVs_parents.sorted.hmp.txt analysis/ld/window-${WINDOW}\kb_filter-${FILTER}
