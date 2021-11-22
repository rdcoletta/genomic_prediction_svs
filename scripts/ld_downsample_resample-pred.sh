#!/bin/bash
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=10
#SBATCH --mem=100gb
#SBATCH -J ld_downsample_resample-pred
#SBATCH -o /home/hirschc1/della028/projects/genomic_prediction/simulation/analysis/ld_downsample/%x_%j.out
#SBATCH -e /home/hirschc1/della028/projects/genomic_prediction/simulation/analysis/ld_downsample/%x_%j.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=della028@umn.edu
#SBATCH --no-requeue

# function to get seed number
get_seeded_random()
{
  seed="$1"
  openssl enc -aes-256-ctr -pass pass:"$seed" -nosalt < /dev/zero 2> /dev/null
}
SEED=($(shuf -i 1-100000 -n 1 --random-source=<(get_seeded_random ${SEED})))

# go to project folder
cd ~/projects/genomic_prediction/simulation

echo -e "job started @ $(date)\n"

# get folder with original LD file and QTN list
LDFOLDER=analysis/ld_downsample/sim_traits/rep${REP}/${QTN}-QTNs_from_${VAR}/0.3-heritability/pop${POP}
# get folder with results for different LD levels
OUTFOLDER=analysis/ld_downsample/pred_${LD}/rep${REP}/${QTN}-QTNs_from_${VAR}/pop${POP}
echo ${OUTFOLDER}

# calculte number of markers to resample
NPREDS=$(wc -l ${OUTFOLDER}/predictors_${LD}-ld_to_QTNs.txt | cut -d " " -f 1)
# count number of iterations
i=1
# keep track of number of predictors from previous iteraction
# (exit while loop if number doesn't change for 5 times)
SAMEPREDS=0

while [[ ${NPREDS} < 500 && ${SAMEPREDS} < 5 ]]; do

  # keep initial number of predictors
  INITPRED=${NPREDS}
  # calculte number of markers to resample
  NRESAMPLE=$(( ${NMARKERS} - ${NPREDS} ))

  echo "  iteration $i"

  # randomly select markers to be predictors -- but exclude those already known
  # to be in the wrong LD category
  echo "  - resampling predictors"
  grep -v -F -f ${OUTFOLDER}/predictors_to_exclude.txt ${OUTFOLDER}/markers_${LD}-ld_to_QTNs.txt | shuf -n ${NRESAMPLE} --random-source=<(get_seeded_random ${SEED}) >> ${OUTFOLDER}/predictors_${LD}-ld_to_QTNs.txt

  # calculate LD between QTLs and predictors
  echo "  - calculating LD"
  for chr in {1..10}; do
    run_pipeline.pl -Xmx10g -importGuess data/usda_rils_projected-SVs-SNPs.chr${chr}.poly.hmp.txt \
                    -includeSiteNamesInFile <(cat ${LDFOLDER}/list_QTNs.chr${chr}.causal-pop${POP}.txt ${OUTFOLDER}/predictors_${LD}-ld_to_QTNs.txt) \
                    -export ${OUTFOLDER}/qtn-pred.chr${chr} \
                    -exportType Plink > /dev/null &
  done
  wait
  for chr in {1..10}; do
    if [[ -s ${OUTFOLDER}/qtn-pred.chr${chr}.plk.map ]]; then
      plink --file ${OUTFOLDER}/qtn-pred.chr${chr}.plk \
            --make-founders --r2 dprime with-freqs --ld-window-r2 0 \
            --ld-window 100000 --ld-window-kb 350000 --geno 0.25 \
            --out ${OUTFOLDER}/qtn-pred.chr${chr} > /dev/null &
    fi
  done
  wait
  # merge chr files
  cat ${OUTFOLDER}/qtn-pred.chr*.ld | head -n 1 > ${OUTFOLDER}/qtn-pred.ld
  for file in ${OUTFOLDER}/qtn-pred.chr*.ld; do
    sed 1d ${file} >> ${OUTFOLDER}/qtn-pred.ld
  done
  # remove unnecessary plink files
  rm ${OUTFOLDER}/qtn-pred.chr*

  # identify predictors above the LD threshold to any other QTL
  # and append to list of known predictors to exclude
  # (grep QTLs from ld file | grep predictors | filter by wrong LD >> predictor_to_exclude.txt)
  if [[ ${LD} == low ]]; then
    grep -F -f <(cat ${LDFOLDER}/list_QTNs.chr*.causal-pop${POP}.txt) ${OUTFOLDER}/qtn-pred.ld | grep -F -f ${OUTFOLDER}/predictors_${LD}-ld_to_QTNs.txt | awk '$9 > 0.5' | sed 's/^ *//' | tr -s " " | tr " " "\t" | cut -f 3,7 | sed 's/\t/\n/g' | sort | uniq | grep -v -F -f <(cat ${LDFOLDER}/list_QTNs.chr*.causal-pop${POP}.txt) >> ${OUTFOLDER}/predictors_to_exclude.txt
  elif [[ ${LD} == moderate ]]; then
    # identify and append to list of known predictors to exclude
    grep -F -f <(cat ${LDFOLDER}/list_QTNs.chr*.causal-pop${POP}.txt) ${OUTFOLDER}/qtn-pred.ld | grep -F -f ${OUTFOLDER}/predictors_${LD}-ld_to_QTNs.txt | awk '$9 > 0.9' | sed 's/^ *//' | tr -s " " | tr " " "\t" | cut -f 3,7 | sed 's/\t/\n/g' | sort | uniq | grep -v -F -f <(cat ${LDFOLDER}/list_QTNs.chr*.causal-pop${POP}.txt) >> ${OUTFOLDER}/predictors_to_exclude.txt
  fi
  # remove predictors
  echo "  - removing predictors above LD threshold"
  grep -v -f ${OUTFOLDER}/predictors_to_exclude.txt ${OUTFOLDER}/predictors_${LD}-ld_to_QTNs.txt > ${OUTFOLDER}/tmp_predictors_list
  mv ${OUTFOLDER}/tmp_predictors_list ${OUTFOLDER}/predictors_${LD}-ld_to_QTNs.txt

  # update number of predictors left
  NPREDS=$(wc -l ${OUTFOLDER}/predictors_${LD}-ld_to_QTNs.txt | cut -d " " -f 1)

  # update to number of times same number of predictors remain the same
  # otherwise, set variable back to 0
  if [[ ${INITPRED} == ${NPREDS} ]]; then
    SAMEPREDS=$((SAMEPREDS+1))
  else
    SAMEPREDS=0
  fi
  # update to number of iterations
  i=$((i+1))
  # use a different seed for next iteration
  SEED=($(shuf -i 1-100000 -n 1 --random-source=<(get_seeded_random ${SEED})))

done

# filter hmp to have only correct markers
run_pipeline.pl -Xmx10g -importGuess analysis/ld_downsample/datasets/rep${REP}/usda_rils_projected-SVs-SNPs.poly.rep${REP}.with-LD-info.hmp.txt \
                -includeSiteNamesInFile ${OUTFOLDER}/predictors_${LD}-ld_to_QTNs.txt \
                -export ${OUTFOLDER}/usda_rils_projected-SVs-SNPs.${LD}-ld_to_QTNs.hmp.txt \
                -exportType HapmapDiploid > /dev/null

echo -e "\njob finished @ $(date)"
