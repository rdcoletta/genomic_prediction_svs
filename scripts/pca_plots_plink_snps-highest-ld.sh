#!/bin/bash
#PBS -l walltime=4:00:00,nodes=1:ppn=1,mem=100gb
#PBS -o /home/hirschc1/della028/projects/genomic_prediction/simulation/analysis/ld/
#PBS -e /home/hirschc1/della028/projects/genomic_prediction/simulation/analysis/ld/
#PBS -V
#PBS -N pca_plots_plink_snps-highest-ld
#PBS -M della028@umn.edu
#PBS -m abe
#PBS -r n

cd ~/projects/genomic_prediction/simulation/

head -n 1 data/usda_rils_projected-SNPs-only.hmp.txt > analysis/ld/window-1kb_filter-0.25/usda_rils_projected-SNPs-only.highest-ld.hmp.txt
grep -f <(grep -v -P "^del|^dup|^ins|^inv|^tra" analysis/ld/window-1kb_filter-0.25/marker_info_highest-ld.txt | cut -f 1) data/usda_rils_projected-SNPs-only.hmp.txt >> analysis/ld/window-1kb_filter-0.25/usda_rils_projected-SNPs-only.highest-ld.hmp.txt

cd analysis/ld/window-1kb_filter-0.25/

# hmp to plink
run_pipeline.pl -Xmx10g -importGuess usda_rils_projected-SNPs-only.highest-ld.hmp.txt -export usda_rils_projected-SNPs-only.highest-ld -exportType Plink

# change first marker names (some are duplicated and plink doesn't like that)
module load R/3.6.0
Rscript ~/projects/genomic_prediction/simulation/scripts/change_marker_names_plink_map.R usda_rils_projected-SNPs-only.highest-ld.plk.map

# transform format to one that plink2 accepts
plink --file usda_rils_projected-SNPs-only.highest-ld.plk --mind 0.5 --make-bed --freq --out usda_rils_projected-SNPs-only.highest-ld --allow-extra-chr --make-founders

# create pca with plink2
plink2 --bfile usda_rils_projected-SNPs-only.highest-ld --make-grm-list --pca 10 --out pca_plink_usda_high-ld-snps --allow-extra-chr --read-freq usda_rils_projected-SNPs-only.highest-ld.frq

# plot pca
module load R
Rscript ~/projects/genomic_prediction/simulation/scripts/plot_PCA_plink.R pca_plink_usda_high-ld-snps.eigenvec ../pca_plink_usda_high-ld-snps_window-1kb_filter-0.25.pdf
