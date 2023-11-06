#!/bin/bash
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=50gb
#SBATCH -J ld_snp-sv_poly-markers_${CHR}
#SBATCH -o /home/hirschc1/della028/projects/genomic_prediction/simulation/analysis/trait_sim/MSI_dump/%x_%j.out
#SBATCH -e /home/hirschc1/della028/projects/genomic_prediction/simulation/analysis/trait_sim/MSI_dump/%x_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=della028@umn.edu
#SBATCH --no-requeue

# go to project folder
cd ~/projects/genomic_prediction/simulation

# transform hmp into plink format
run_pipeline.pl -Xmx100g -importGuess data/usda_rils_projected-SVs-SNPs.chr${CHR}.poly.hmp.txt -export /scratch.global/della028/hirsch_lab/genomic_prediction/ld/usda_rils_projected-SVs-SNPs.chr${CHR}.poly -exportType Plink

# go to scratch folder
cd /scratch.global/della028/hirsch_lab/genomic_prediction/ld

# calculate LD
plink --file usda_rils_projected-SVs-SNPs.chr${CHR}.poly.plk --make-founders --r2 gz dprime with-freqs --ld-window-r2 0 --ld-window ${WINDOW}000 --ld-window-kb ${WINDOW} --geno 0.25 --out usda_rils_projected-SVs-SNPs.chr${CHR}.poly

# keep only snp and sv r2 (excluding translocations)
zcat usda_rils_projected-SVs-SNPs.chr${CHR}.poly.ld.gz | head -n 1 > ~/projects/genomic_prediction/simulation/analysis/ld/ld_usda_rils_poly-snp-sv_only.chr${CHR}.ld
zcat usda_rils_projected-SVs-SNPs.chr${CHR}.poly.ld.gz | awk '$3 ~ /^del|^dup|^ins|^inv|^tra/ && $7 !~ /^del|^dup|^ins|^inv|^tra/ || $3 !~ /^del|^dup|^ins|^inv|^tra/ && $7 ~ /^del|^dup|^ins|^inv|^tra/' - >> ~/projects/genomic_prediction/simulation/analysis/ld/ld_usda_rils_poly-snp-sv_only.chr${CHR}.ld
