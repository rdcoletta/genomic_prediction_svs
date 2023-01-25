#!/bin/bash
#PBS -l walltime=5:00:00,nodes=1:ppn=1,mem=200gb
#PBS -o /home/hirschc1/della028/projects/genomic_prediction/simulation/analysis/ld
#PBS -e /home/hirschc1/della028/projects/genomic_prediction/simulation/analysis/ld
#PBS -V
#PBS -N calculate_ld_by_parent_${PARENT}_${CHR}
#PBS -M della028@umn.edu
#PBS -m abe
#PBS -r n

cd /scratch.global/della028/hirsch_lab/genomic_prediction/ld

# convert hapmap to plink format
run_pipeline.pl -Xmx200g -importGuess ~/projects/genomic_prediction/simulation/data/hapmap_by_common-parent/usda_rils_projected-SVs-SNPs.${PARENT}-parent.chr${CHR}.hmp.txt -export usda_rils_projected-SVs-SNPs.${PARENT}-parent.chr${CHR} -exportType Plink

# calculate LD
plink --file usda_rils_projected-SVs-SNPs.${PARENT}-parent.chr${CHR}.plk --make-founders --r2 gz dprime with-freqs --ld-window-r2 0 --ld-window 1000000 --ld-window-kb 1000 --geno 0.25 --out usda_rils_projected-SVs-SNPs.${PARENT}-parent.chr${CHR}.window-1000kb.filter-0.25

# copy header for filtered file
zcat usda_rils_projected-SVs-SNPs.${PARENT}-parent.chr${CHR}.window-1000kb.filter-0.25.ld.gz | head -n 1 > ld_usda_rils_snp-sv_only.${PARENT}-parent.chr${CHR}.window-1000kb.filter-0.25.ld
# keep only snp and sv r2 (excluding translocations)
zcat usda_rils_projected-SVs-SNPs.${PARENT}-parent.chr${CHR}.window-1000kb.filter-0.25.ld.gz | awk '$3 ~ /^del|^dup|^ins|^inv|^tra/ && $7 !~ /^del|^dup|^ins|^inv|^tra/ || $3 !~ /^del|^dup|^ins|^inv|^tra/ && $7 ~ /^del|^dup|^ins|^inv|^tra/' - >> ld_usda_rils_snp-sv_only.${PARENT}-parent.chr${CHR}.window-1000kb.filter-0.25.ld
# # compress ld file
# gzip ld_usda_rils_snp-sv_only.${PARENT}-parent.chr${CHR}.window-1000kb.filter-0.25.ld
