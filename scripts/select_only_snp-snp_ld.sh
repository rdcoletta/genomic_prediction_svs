#!/bin/bash
#PBS -l walltime=4:00:00,nodes=1:ppn=1,mem=50gb
#PBS -o /home/hirschc1/della028/projects/genomic_prediction/simulation/analysis/ld
#PBS -e /home/hirschc1/della028/projects/genomic_prediction/simulation/analysis/ld
#PBS -V
#PBS -N select_only_snp-snp_ld_${CHR}_${WINDOW}_${FILTER}
#PBS -M della028@umn.edu
#PBS -m abe
#PBS -r n

# go to project folder
cd /scratch.global/della028/hirsch_lab/genomic_prediction/ld

# keep only snp and snp LD (and exclude uncessary columns to reduce file size)
zcat usda_rils_projected-SVs-SNPs.chr${CHR}.window-${WINDOW}kb.filter-${FILTER}.ld.gz | awk '$3 !~ /^del|^dup|^ins|^inv|^tra/ && $7 !~ /^del|^dup|^ins|^inv|^tra/' - | sed 's/^ *//' | tr -s " " | tr " " "\t" | cut -f 1,2,6,9 >> ld_usda_rils_snp-snp_only.chr${CHR}.window-${WINDOW}kb.filter-${FILTER}.ld

# compress file
gzip ld_usda_rils_snp-snp_only.chr${CHR}.window-${WINDOW}kb.filter-${FILTER}.ld
