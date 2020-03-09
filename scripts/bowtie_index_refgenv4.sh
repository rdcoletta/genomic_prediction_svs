#!/bin/bash
#PBS -l walltime=2:00:00,nodes=1:ppn=5,mem=60gb
#PBS -o /home/hirschc1/della028/projects/genomic_prediction/simulation/data/refgen/B73v4
#PBS -e /home/hirschc1/della028/projects/genomic_prediction/simulation/data/refgen/B73v4
#PBS -V
#PBS -N bowtie_index_refgenv4
#PBS -M della028@umn.edu
#PBS -m abe
#PBS -r n

module load bowtie/1.1.2

# go to refgen v4 folder
cd ~/projects/genomic_prediction/simulation/data/refgen/B73v4
# build index
bowtie-build chr1.fa,chr2.fa,chr3.fa,chr4.fa,chr5.fa,chr6.fa,chr7.fa,chr8.fa,chr9.fa,chr10.fa index
