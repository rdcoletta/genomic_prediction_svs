#!/bin/bash
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=50gb
#SBATCH -J hmp2plk_random_${CHR}
#SBATCH -o /home/hirschc1/della028/projects/genomic_prediction/simulation/analysis/ld/%x_%j.out
#SBATCH -e /home/hirschc1/della028/projects/genomic_prediction/simulation/analysis/ld/%x_%j.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=della028@umn.edu
#SBATCH --no-requeue

# define scratch folder
SCRATCH=/scratch.global/della028/hirsch_lab/genomic_prediction/ld
# define seeds
get_seeded_random()
{
  seed="$1"
  openssl enc -aes-256-ctr -pass pass:"$seed" -nosalt < /dev/zero 2> /dev/null
}
SEEDS=($(shuf -i 1-100000 -n 10 --random-source=<(get_seeded_random 2021)))
# define starting seed index
SEED_IDX=0

# go to project folder
cd ~/projects/genomic_prediction/simulation

for rep in {1..10}; do
  # create output folder
  OUTFOLDER=analysis/ld_downsample/datasets/rep${rep}
  mkdir -p ${OUTFOLDER}
  # calculate number that represents 25% of data
  TOTAL=$(sed 1d data/usda_rils_projected-SVs-SNPs.chr${CHR}.poly.hmp.txt | wc -l)
  SAMPLE=$(( ${TOTAL} * 25 / 100 ))
  # randomly select markers to keep
  cut -f 1 data/usda_rils_projected-SVs-SNPs.chr${CHR}.poly.hmp.txt | sed 1d | shuf -n ${SAMPLE} --random-source=<(get_seeded_random ${SEEDS[${SEED_IDX}]}) > ${OUTFOLDER}/markers_to_keep.chr${CHR}.rep${rep}.txt
  # filter hmp to create dataset
  run_pipeline.pl -Xmx50g -importGuess data/usda_rils_projected-SVs-SNPs.chr${CHR}.poly.hmp.txt \
                  -includeSiteNamesInFile ${OUTFOLDER}/markers_to_keep.chr${CHR}.rep${rep}.txt \
                  -export ${OUTFOLDER}/usda_rils_projected-SVs-SNPs.poly.chr${CHR}.rep${rep}.hmp.txt \
                  -exportType HapmapDiploid
  # export to plink for LD calculation
  run_pipeline.pl -Xmx50g -importGuess ${OUTFOLDER}/usda_rils_projected-SVs-SNPs.poly.chr${CHR}.rep${rep}.hmp.txt \
                  -export ${SCRATCH}/usda_rils_projected-SVs-SNPs.poly.chr${CHR}.rep${rep} \
                  -exportType Plink
  # increment seed index variable
  SEED_IDX=$(( ${SEED_IDX} + 1 ))
done
