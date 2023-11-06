#!/bin/bash
#SBATCH --time=3:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=20gb
#SBATCH -J project_svs-snps
#SBATCH -o /home/hirschc1/della028/projects/genomic_prediction/simulation/analysis/projection_svs-snps/%x_%j.out
#SBATCH -e /home/hirschc1/della028/projects/genomic_prediction/simulation/analysis/projection_svs-snps/%x_%j.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=della028@umn.edu
#SBATCH --no-requeue

# go to project folder
cd ~/projects/genomic_prediction/simulation

# create donor files
run_pipeline.pl -Xmx20g -FILLINFindHaplotypesPlugin -hmp data/hapmap_by_cross/usda_parents_SV-SNPchip-polySNPreseq.${CROSS}.hmp.txt -o analysis/projection_svs-snps/donors_${CROSS} -hapSize ${HAPSIZE} -minTaxa 1

# project reseq snps
run_pipeline.pl -Xmx20g -FILLINImputationPlugin -hmp data/hapmap_by_cross/usda_rils_SV-SNPchip-polySNPreseq.${CROSS}.not-projected.hmp.txt -d analysis/projection_svs-snps/donors_${CROSS} -o analysis/projection_svs-snps/usda_rils_SV-SNPchip-polySNPreseq.${CROSS}.projected.hmp.txt -hapSize ${HAPSIZE} -hybNN false -accuracy
