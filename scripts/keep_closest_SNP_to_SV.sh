#!/bin/bash
#PBS -l walltime=5:00:00,nodes=1:ppn=1,mem=120gb
#PBS -o /home/hirschc1/della028/projects/genomic_prediction/simulation/analysis/ld
#PBS -e /home/hirschc1/della028/projects/genomic_prediction/simulation/analysis/ld
#PBS -V
#PBS -N keep_closest_SNP_to_SV_${CHR}_${TYPE}_${FILTER}
#PBS -M della028@umn.edu
#PBS -m abe
#PBS -r n

# go to project folder
cd ~/projects/genomic_prediction/simulation

# find closest low missing SNP to SV
module load R
Rscript scripts/get_closest-SNP_to_SV.R data/usda_SNPs-SVs_rils.not-in-SVs.projected.chr${CHR}.reseq-SNPs.hmp.txt analysis/ld/closest_low-missing-data-SNPs_to_SVs.chr-${CHR}.filter-${FILTER}.txt ${FILTER} ${TYPE}

# keep only SNPs with low missing data
run_pipeline.pl -Xmx120g -importGuess data/usda_SNPs-SVs_rils.not-in-SVs.projected.chr${CHR}.reseq-SNPs.hmp.txt \
                -includeSiteNamesInFile analysis/ld/closest_low-missing-data-SNPs_to_SVs.chr-${CHR}.filter-${FILTER}.${TYPE}.txt \
                -export analysis/ld/usda_SNPs-SVs_rils.not-in-SVs.projected.chr${CHR}.reseq-SNPs.closest-snps-to-svs.filter-${FILTER}.${TYPE}.hmp.txt \
                -exportType HapmapDiploid
