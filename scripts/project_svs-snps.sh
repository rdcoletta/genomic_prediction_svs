#!/bin/bash
#PBS -l walltime=3:00:00,nodes=1:ppn=1,mem=20gb
#PBS -o /home/hirschc1/della028/projects/genomic_prediction/simulation/analysis/projection_svs-snps
#PBS -e /home/hirschc1/della028/projects/genomic_prediction/simulation/analysis/projection_svs-snps
#PBS -V
#PBS -N project_svs-snps_${CROSS}
#PBS -M della028@umn.edu
#PBS -m abe
#PBS -r n

# go to project folder
cd ~/projects/genomic_prediction/simulation

# create donor files
run_pipeline.pl -Xmx20g -FILLINFindHaplotypesPlugin -hmp data/hapmap_by_cross/usda_parents_SV-SNPchip-polySNPreseq.${CROSS}.hmp.txt -o analysis/projection_svs-snps/donors_${CROSS} -hapSize ${HAPSIZE} -minTaxa 1

# project reseq snps
run_pipeline.pl -Xmx20g -FILLINImputationPlugin -hmp data/hapmap_by_cross/usda_rils_SV-SNPchip-polySNPreseq.${CROSS}.not-projected.hmp.txt -d analysis/projection_svs-snps/donors_${CROSS} -o analysis/projection_svs-snps/usda_rils_SV-SNPchip-polySNPreseq.${CROSS}.projected.hmp.txt -hapSize ${HAPSIZE} -hybNN false -accuracy
