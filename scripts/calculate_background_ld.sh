#!/bin/bash
#SBATCH --time=5:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem=80gb
#SBATCH -J calculate_background_ld
#SBATCH -o /home/hirschc1/della028/projects/genomic_prediction/simulation/analysis/background_ld/%x_%j.out
#SBATCH -e /home/hirschc1/della028/projects/genomic_prediction/simulation/analysis/background_ld/%x_%j.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=della028@umn.edu
#SBATCH --no-requeue


# define function to sample different seeds
get_seeded_random()
{
  seed="$1"
  openssl enc -aes-256-ctr -pass pass:"$seed" -nosalt < /dev/zero 2> /dev/null
}

# bootstrap random sampling many times
for iter in {1..10}; do

  # create folder for iteration
  mkdir -p ${FOLDER}/iter${iter}

  # randomly select 50 markers per chr
  for chr in {1..10}; do
    # sample markers
    shuf -n 50 --random-source=<(get_seeded_random ${SEED}) <(awk -v chr="$chr" '$3 == chr' ${HMP} | cut -f 1) > ${FOLDER}/iter${iter}/markers_chr${chr}.txt
    # change seed
    SEED=($(shuf -i 1-100000 -n 1 --random-source=<(get_seeded_random ${SEED})))
  done

  # create file to store results
  BACKLDITER=${FOLDER}/iter${iter}/background_ld.txt
  echo -n "" > ${BACKLDITER}

  # calculate ld for all markers that are in different chromosomes
  for chr in {1..10}; do
    echo "iter${iter} - chr${chr}"
    # hmp2plk
    PLK=${FOLDER}/iter${iter}/inbreds_geno
    run_pipeline.pl -Xmx40g -importGuess ${HMP} \
                            -includeSiteNamesInFile <(cat ${FOLDER}/iter${iter}/markers_chr*.txt) \
                            -export ${PLK} \
                            -exportType Plink > /dev/null
    # calculate LD
    plink --file ${PLK}.plk \
          --make-founders \
          --r2 gz dprime with-freqs inter-chr \
          --ld-window-r2 0 \
          --geno 0.25 \
          --ld-snp-list ${FOLDER}/iter${iter}/markers_chr${chr}.txt \
          --out ${FOLDER}/iter${iter}/markers_chr${chr} > /dev/null
    # add R2 values into file
    zcat ${FOLDER}/iter${iter}/markers_chr${chr}.ld.gz | sed 's/^ *//' | tr -s " " | tr " " "\t" | awk -v chr=${chr} '$5 != chr' | cut -f 9 | sed 1d >> ${BACKLDITER}
  done

  # remove extra files
  rm ${FOLDER}/iter${iter}/*.log
  rm ${FOLDER}/iter${iter}/*.nosex

  # background ld as 95% of all these pairwise correlations

  # note the number of values in the sample (K)
  K=$(wc -l ${BACKLDITER} | cut -d " " -f 1)
  # calculate N. N = K x 0.95.
  N=$(printf "%.0f" $(echo "${K} * 0.95" | bc))
  # sort the samples in ascending order
  # the Nth value in the sorted list will be the 95th percentile value
  sort -g ${BACKLDITER} | head -n ${N} | tail -n 1 >> ${BACKLD}

  # https://www.manageengine.com/network-monitoring/faq/95th-percentile-calculation.html

done
