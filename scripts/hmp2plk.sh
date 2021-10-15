#!/bin/bash
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=110gb
#SBATCH -J hmp2plk_${CHR}
#SBATCH -o /home/hirschc1/della028/projects/genomic_prediction/simulation/analysis/ld/%x_%j.out
#SBATCH -e /home/hirschc1/della028/projects/genomic_prediction/simulation/analysis/ld/%x_%j.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=della028@umn.edu
#SBATCH --no-requeue

# go to project folder
cd ~/projects/genomic_prediction/simulation

# transform hmp into plink format
run_pipeline.pl -Xmx110g -importGuess data/usda_rils_projected-SVs-SNPs.chr${CHR}.hmp.txt -export /scratch.global/della028/hirsch_lab/genomic_prediction/ld/usda_rils_projected-SVs-SNPs.chr${CHR} -exportType Plink
