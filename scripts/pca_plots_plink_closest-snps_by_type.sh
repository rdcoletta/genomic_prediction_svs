#!/bin/bash
#PBS -l walltime=2:00:00,nodes=1:ppn=1,mem=50gb
#PBS -o /home/hirschc1/della028/projects/genomic_prediction/simulation/analysis/ld/
#PBS -e /home/hirschc1/della028/projects/genomic_prediction/simulation/analysis/ld/
#PBS -V
#PBS -N pca_plots_plink_closest-snps_by_${TYPE}
#PBS -M della028@umn.edu
#PBS -m abe
#PBS -r n

# go to project folder
cd ~/projects/genomic_prediction/simulation/analysis/ld/

run_pipeline.pl -Xmx10g -importGuess usda_SNPs-SVs_rils.not-in-SVs.projected.reseq-SNPs.closest-snps-to-svs.${TYPE}.hmp.txt -export usda_SNPs-SVs_rils.not-in-SVs.projected.reseq-SNPs.closest-snps-to-svs.${TYPE} -exportType Plink

# change first marker names (some are duplicated and plink doesn't like that)
module load R/3.6.0
Rscript ~/projects/genomic_prediction/simulation/scripts/change_marker_names_plink_map.R usda_SNPs-SVs_rils.not-in-SVs.projected.reseq-SNPs.closest-snps-to-svs.${TYPE}.plk.map

# transform format to one that plink2 accepts
plink --file usda_SNPs-SVs_rils.not-in-SVs.projected.reseq-SNPs.closest-snps-to-svs.${TYPE}.plk --mind 0.8 --make-bed --freq --out usda_SNPs-SVs_rils.not-in-SVs.projected.reseq-SNPs.closest-snps-to-svs.${TYPE} --allow-extra-chr --make-founders

# create pca with plink2
plink2 --bfile usda_SNPs-SVs_rils.not-in-SVs.projected.reseq-SNPs.closest-snps-to-svs.${TYPE} --make-grm-list --pca 10 --out pca_plink_usda_closest-snps_${TYPE} --allow-extra-chr --read-freq usda_SNPs-SVs_rils.not-in-SVs.projected.reseq-SNPs.closest-snps-to-svs.${TYPE}.frq

# plot pca
module load R
Rscript ~/projects/genomic_prediction/simulation/scripts/plot_PCA_plink.R pca_plink_usda_closest-snps_${TYPE}.eigenvec pca_plink_usda_closest-snps_${TYPE}.pdf
