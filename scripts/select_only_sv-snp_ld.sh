#!/bin/bash
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=50gb
#SBATCH -J select_only_sv-snp_ld
#SBATCH -o /home/hirschc1/della028/projects/genomic_prediction/simulation/analysis/ld/%x_%j.out
#SBATCH -e /home/hirschc1/della028/projects/genomic_prediction/simulation/analysis/ld/%x_%j.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=della028@umn.edu
#SBATCH --no-requeue

# go to project folder
cd /scratch.global/della028/hirsch_lab/genomic_prediction/ld

# copy header for filtered file
zcat usda_rils_projected-SVs-SNPs.chr${CHR}.window-${WINDOW}kb.filter-${FILTER}.ld.gz | head -n 1 > ~/projects/genomic_prediction/simulation/analysis/ld/ld_usda_rils_snp-sv_only.chr${CHR}.window-${WINDOW}kb.filter-${FILTER}.ld

# keep only snp and sv r2 (excluding translocations)
zcat usda_rils_projected-SVs-SNPs.chr${CHR}.window-${WINDOW}kb.filter-${FILTER}.ld.gz | awk '$3 ~ /^del|^dup|^ins|^inv|^tra/ && $7 !~ /^del|^dup|^ins|^inv|^tra/ || $3 !~ /^del|^dup|^ins|^inv|^tra/ && $7 ~ /^del|^dup|^ins|^inv|^tra/' - >> ~/projects/genomic_prediction/simulation/analysis/ld/ld_usda_rils_snp-sv_only.chr${CHR}.window-${WINDOW}kb.filter-${FILTER}.ld
