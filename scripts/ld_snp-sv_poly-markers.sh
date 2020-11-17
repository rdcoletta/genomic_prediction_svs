#!/bin/bash
#PBS -l walltime=1:00:00,nodes=1:ppn=1,mem=50gb
#PBS -o /home/hirschc1/della028/projects/genomic_prediction/simulation/analysis/ld
#PBS -e /home/hirschc1/della028/projects/genomic_prediction/simulation/analysis/ld
#PBS -V
#PBS -N ld_snp-sv_poly-markers_${CHR}
#PBS -M della028@umn.edu
#PBS -m abe
#PBS -r n

# go to project folder
cd ~/projects/genomic_prediction/simulation

# transform hmp into plink format
run_pipeline.pl -Xmx100g -importGuess data/usda_rils_projected-SVs-SNPs.chr${CHR}.poly.hmp.txt -export /scratch.global/della028/hirsch_lab/genomic_prediction/ld/usda_rils_projected-SVs-SNPs.chr${CHR}.poly -exportType Plink

# go to scratch folder
cd /scratch.global/della028/hirsch_lab/genomic_prediction/ld

# calculate LD
plink --file usda_rils_projected-SVs-SNPs.chr${CHR}.poly.plk --make-founders --r2 gz dprime with-freqs --ld-window-r2 0 --ld-window 1000 --ld-window-kb 1 --geno 0.25 --out usda_rils_projected-SVs-SNPs.chr${CHR}.poly

# keep only snp and sv r2 (excluding translocations)
zcat usda_rils_projected-SVs-SNPs.chr${CHR}.poly.ld.gz | head -n 1 > ~/projects/genomic_prediction/simulation/analysis/ld/ld_usda_rils_poly-snp-sv_only.chr${CHR}.ld
zcat usda_rils_projected-SVs-SNPs.chr${CHR}.poly.ld.gz | awk '$3 ~ /^del|^dup|^ins|^inv|^tra/ && $7 !~ /^del|^dup|^ins|^inv|^tra/ || $3 !~ /^del|^dup|^ins|^inv|^tra/ && $7 ~ /^del|^dup|^ins|^inv|^tra/' - >> ~/projects/genomic_prediction/simulation/analysis/ld/ld_usda_rils_poly-snp-sv_only.chr${CHR}.ld
