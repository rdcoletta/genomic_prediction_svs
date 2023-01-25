#!/bin/bash
#PBS -l walltime=4:00:00,nodes=1:ppn=1,mem=50gb
#PBS -o /home/hirschc1/della028/projects/genomic_prediction/simulation/analysis/projection_svs-snps
#PBS -e /home/hirschc1/della028/projects/genomic_prediction/simulation/analysis/projection_svs-snps
#PBS -V
#PBS -N sort_merged_projected_files_${CROSS}
#PBS -M della028@umn.edu
#PBS -m abe
#PBS -r n

# go to project folder
cd ~/projects/genomic_prediction/simulation/analysis/projection_svs-snps

run_pipeline.pl -Xmx50g -SortGenotypeFilePlugin -inputFile usda_rils_projected-SVs-SNPs.${CROSS}.all-markers.hmp.txt -outputFile usda_rils_projected-SVs-SNPs.${CROSS}.all-markers.hmp.txt -fileType Hapmap

run_pipeline.pl -Xmx50g -importGuess usda_rils_projected-SVs-SNPs.${CROSS}.all-markers.hmp.txt -export usda_rils_projected-SVs-SNPs.${CROSS}.all-markers.hmp.txt -exportType HapmapDiploid
