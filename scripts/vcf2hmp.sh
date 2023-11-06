#!/bin/bash
#SBATCH --time=6:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=110gb
#SBATCH -J vcf2hmp
#SBATCH -o /home/hirschc1/della028/projects/genomic_prediction/simulation/analysis/ld/%x_%j.out
#SBATCH -e /home/hirschc1/della028/projects/genomic_prediction/simulation/analysis/ld/%x_%j.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=della028@umn.edu
#SBATCH --no-requeue

# go to project folder
cd ~/projects/genomic_prediction/simulation

# transform hmp into plink format
run_pipeline.pl -Xmx100g -importGuess data/WiDiv508_B73v4_allchr_SNPS_maxmiss0.10.recode.vcf.gz \
                -includeTaxa LH82,PH207,PHG35,PHG39,PHG47,PHJ40 \
                -export data/widiv_snps.Qiu-et-al.no-B73.hmp.txt \
                -exportType HapmapDiploid
