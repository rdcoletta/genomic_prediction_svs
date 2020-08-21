#!/bin/bash
#PBS -l walltime=2:00:00,nodes=1:ppn=1,mem=110gb
#PBS -o /home/hirschc1/della028/projects/genomic_prediction/simulation/analysis/ld
#PBS -e /home/hirschc1/della028/projects/genomic_prediction/simulation/analysis/ld
#PBS -V
#PBS -N hmp2plk_${CHR}
#PBS -M della028@umn.edu
#PBS -m abe
#PBS -r n

# go to project folder
cd ~/projects/genomic_prediction/simulation

# transform hmp into plink format
run_pipeline.pl -Xmx110g -importGuess data/usda_SNPs-SVs_rils.not-in-SVs.projected.chr${CHR}.reseq-SNPs.hmp.txt -export /scratch.global/della028/hirsch_lab/genomic_prediction/ld/usda_SNPs-SVs_rils.not-in-SVs.projected.chr${CHR}.reseq-SNPs -exportType Plink
