#!/bin/bash
#PBS -l walltime=5:00:00,nodes=1:ppn=1,mem=100gb
#PBS -o /home/hirschc1/della028/projects/genomic_prediction/simulation/analysis/qc
#PBS -e /home/hirschc1/della028/projects/genomic_prediction/simulation/analysis/qc
#PBS -V
#PBS -N tassel_qc-summary_${CHR}
#PBS -M della028@umn.edu
#PBS -m abe
#PBS -r n

# go to project folder
cd ~/projects/genomic_prediction/simulation

run_pipeline.pl -Xmx100g -importGuess data/usda_rils_projected-SVs-SNPs.chr${CHR}.hmp.txt -GenotypeSummaryPlugin -endPlugin -export analysis/qc/tassel_usda-SNPs-SVs_chr${CHR}_summary
