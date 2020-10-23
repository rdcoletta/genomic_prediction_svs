#!/bin/bash
#PBS -l walltime=4:00:00,nodes=1:ppn=1,mem=100gb
#PBS -o /home/hirschc1/della028/projects/genomic_prediction/simulation/analysis/ld/
#PBS -e /home/hirschc1/della028/projects/genomic_prediction/simulation/analysis/ld/
#PBS -V
#PBS -N pca_plots_plink_svs
#PBS -M della028@umn.edu
#PBS -m abe
#PBS -r n

# go to project folder
cd ~/projects/genomic_prediction/simulation/

head -n 1 data/usda_rils_projected-SVs-SNPs.hmp.txt > data/usda_rils_projected-SVs-only.hmp.txt
grep -P "^del|^dup|^ins|^inv|^tra" data/usda_rils_projected-SVs-SNPs.hmp.txt >> data/usda_rils_projected-SVs-only.hmp.txt

run_pipeline.pl -Xmx10g -importGuess data/usda_rils_projected-SVs-only.hmp.txt -export analysis/ld/usda_rils_projected-SVs-only -exportType Plink

cd ~/projects/genomic_prediction/simulation/analysis/ld/

# change first marker names (some are duplicated and plink doesn't like that)
module load R/3.6.0
Rscript ~/projects/genomic_prediction/simulation/scripts/change_marker_names_plink_map.R usda_rils_projected-SVs-only.plk.map

# transform format to one that plink2 accepts
plink --file usda_rils_projected-SVs-only.plk --mind 0.5 --make-bed --freq --out usda_rils_projected-SVs-only --allow-extra-chr --make-founders

# create pca with plink2
plink2 --bfile usda_rils_projected-SVs-only --make-grm-list --pca 10 --out pca_plink_usda_svs --allow-extra-chr --read-freq usda_rils_projected-SVs-only.frq

# plot pca
module load R
Rscript ~/projects/genomic_prediction/simulation/scripts/plot_PCA_plink.R pca_plink_usda_svs.eigenvec pca_plink_usda_svs.pdf
