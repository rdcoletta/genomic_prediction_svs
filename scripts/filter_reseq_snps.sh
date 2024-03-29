#!/bin/bash
#SBATCH --time=3:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=50gb
#SBATCH -J vcf2hmp
#SBATCH -o /home/hirschc1/della028/projects/genomic_prediction/simulation/analysis/ld/%x_%j.out
#SBATCH -e /home/hirschc1/della028/projects/genomic_prediction/simulation/analysis/ld/%x_%j.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=della028@umn.edu
#SBATCH --no-requeue

module load R/3.6.0

# go to project folder
cd ~/projects/genomic_prediction/simulation

# remove SNPs inside deletions
Rscript scripts/remove_SNPs_within_dels.R data/reseq_snps/widiv_snps_usda_parents.${CROSS}.hmp.txt data/hapmap_by_cross/usda_SVs_parents.sorted.${CROSS}.hmp.txt

# keep only polymorphic
Rscript scripts/keep_only_poly_snps.R ${CROSS} data/reseq_snps/widiv_snps_usda_parents.${CROSS}.not-in-SVs.hmp.txt
