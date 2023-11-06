#!/bin/bash
#SBATCH --time=4:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=50gb
#SBATCH -J sort_merged_projected_files
#SBATCH -o /home/hirschc1/della028/projects/genomic_prediction/simulation/analysis/projection_svs-snps/%x_%j.out
#SBATCH -e /home/hirschc1/della028/projects/genomic_prediction/simulation/analysis/projection_svs-snps/%x_%j.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=della028@umn.edu
#SBATCH --no-requeue

# go to project folder
cd ~/projects/genomic_prediction/simulation/analysis/projection_svs-snps

run_pipeline.pl -Xmx50g -SortGenotypeFilePlugin -inputFile usda_rils_projected-SVs-SNPs.${CROSS}.all-markers.hmp.txt -outputFile usda_rils_projected-SVs-SNPs.${CROSS}.all-markers.hmp.txt -fileType Hapmap

run_pipeline.pl -Xmx50g -importGuess usda_rils_projected-SVs-SNPs.${CROSS}.all-markers.hmp.txt -export usda_rils_projected-SVs-SNPs.${CROSS}.all-markers.hmp.txt -exportType HapmapDiploid
