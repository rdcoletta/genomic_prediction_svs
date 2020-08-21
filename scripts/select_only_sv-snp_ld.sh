#!/bin/bash
#PBS -l walltime=8:00:00,nodes=1:ppn=1,mem=100gb
#PBS -o /home/hirschc1/della028/projects/genomic_prediction/simulation/analysis/ld
#PBS -e /home/hirschc1/della028/projects/genomic_prediction/simulation/analysis/ld
#PBS -V
#PBS -N select_only_sv-snp_ld_${CHR}_${WINDOW}_${FILTER}
#PBS -M della028@umn.edu
#PBS -m abe
#PBS -r n

# go to project folder
cd /scratch.global/della028/hirsch_lab/genomic_prediction/ld

# copy header for filtered file
zcat usda_SNPs-SVs_rils.not-in-SVs.projected.chr${CHR}.reseq-SNPs.window-${WINDOW}kb.filter-${FILTER}.ld.gz | head -n 1 > ~/projects/genomic_prediction/simulation/analysis/ld/ld_usda_rils_snp-sv_only.chr${CHR}.window-${WINDOW}kb.filter-${FILTER}.ld

# keep only snp and sv r2 (excluding translocations)
zcat usda_SNPs-SVs_rils.not-in-SVs.projected.chr${CHR}.reseq-SNPs.window-${WINDOW}kb.filter-${FILTER}.ld.gz | awk '$3 ~ /^del|^dup|^ins|^inv|^tra/ && $6 ~ !/^del|^dup|^ins|^inv|^tra/ || $3 ~ !/^del|^dup|^ins|^inv|^tra/ && $6 ~ /^del|^dup|^ins|^inv|^tra/' - >> ~/projects/genomic_prediction/simulation/analysis/ld/ld_usda_rils_snp-sv_only.chr${CHR}.window-${WINDOW}kb.filter-${FILTER}.ld
