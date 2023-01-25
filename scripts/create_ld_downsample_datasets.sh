#!/bin/bash
#SBATCH --time=14:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem=100gb
#SBATCH -J create_ld_downsample_datasets
#SBATCH -o /home/hirschc1/della028/projects/genomic_prediction/simulation/analysis/trait_sim/MSI_dump/%x_%j.out
#SBATCH -e /home/hirschc1/della028/projects/genomic_prediction/simulation/analysis/trait_sim/MSI_dump/%x_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=della028@umn.edu
#SBATCH --no-requeue

module load R/3.6.0

WINDOW=5000
FILTER=0.25
SCRATCH=/scratch.global/della028/hirsch_lab/genomic_prediction/ld

# go to project folder
cd ~/projects/genomic_prediction/simulation

for chr in {1..10}; do
  # get all marker names in columns SNP_A and SNP_B
  zcat ${SCRATCH}/usda_rils_projected-SVs-SNPs.poly.chr${chr}.rep${REP}.window-${WINDOW}kb.filter-${FILTER}.ld.gz | sed 's/^ *//' | tr -s " " | tr " " "\t" | sed 1d | cut -f 3,7,9 > ${SCRATCH}/markers_with_ld_results.chr${chr}.rep${REP}.txt

  # get unique marker names for each chr -- will break files into
  # chunks first and then sort them to speed up process
  cd ${SCRATCH}
  cut -f 1,2 markers_with_ld_results.chr${chr}.rep${REP}.txt | sed 's/\t/\n/g' > marker_names_to_filter.chr${chr}.rep${REP}.txt
  split -l 1000000 marker_names_to_filter.chr${chr}.rep${REP}.txt small-chunk.chr${chr}.rep${REP}.
  for chunk in small-chunk.chr${chr}.rep${REP}.*; do
    sort -u $chunk > sorted-$chunk
  done
  sort -m -u sorted-small-chunk.chr${chr}.rep${REP}.* > ~/projects/genomic_prediction/simulation/analysis/ld_downsample/datasets/rep${REP}/marker_names_to_filter.no-dups.chr${chr}.rep${REP}.txt
  # remove tmp files and return to project folder
  rm small-chunk.chr${chr}.rep${REP}.* sorted-small-chunk.chr${chr}.rep${REP}.*
  cd ~/projects/genomic_prediction/simulation

  # filter hmp
  run_pipeline.pl -Xmx60g -importGuess analysis/ld_downsample/datasets/rep${REP}/usda_rils_projected-SVs-SNPs.poly.chr${chr}.rep${REP}.hmp.txt \
                  -includeSiteNamesInFile analysis/ld_downsample/datasets/rep${REP}/marker_names_to_filter.no-dups.chr${chr}.rep${REP}.txt \
                  -export analysis/ld_downsample/datasets/rep${REP}/usda_rils_projected-SVs-SNPs.poly.chr${chr}.rep${REP}.with-LD-info.hmp.txt \
                  -exportType HapmapDiploid
done

# merge chr
head -n 1 analysis/ld_downsample/datasets/rep${REP}/usda_rils_projected-SVs-SNPs.poly.chr1.rep${REP}.with-LD-info.hmp.txt > analysis/ld_downsample/datasets/rep${REP}/usda_rils_projected-SVs-SNPs.poly.rep${REP}.with-LD-info.hmp.txt
for chr in {1..10}; do
  sed 1d analysis/ld_downsample/datasets/rep${REP}/usda_rils_projected-SVs-SNPs.poly.chr${chr}.rep${REP}.with-LD-info.hmp.txt >> analysis/ld_downsample/datasets/rep${REP}/usda_rils_projected-SVs-SNPs.poly.rep${REP}.with-LD-info.hmp.txt
done

# remove chr files
for chr in {1..10}; do
  rm analysis/ld_downsample/datasets/rep${REP}/usda_rils_projected-SVs-SNPs.poly.chr${chr}.rep${REP}.with-LD-info.hmp.txt
done
