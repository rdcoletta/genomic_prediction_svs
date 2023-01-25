#!/bin/bash
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=10
#SBATCH --mem=80gb
#SBATCH -J ld_downsample_select-pred
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

echo "job started @ $(date)"

for QTN in 100; do
  for VAR in SNP SV both; do
    for LD in low moderate high; do

      # get folder with original LD file and QTN list
      LDFOLDER=analysis/ld_downsample/sim_traits/rep${REP}/${QTN}-QTNs_from_${VAR}/0.3-heritability/pop${POP}
      OUTFOLDER=analysis/ld_downsample/pred_${LD}/rep${REP}/${QTN}-QTNs_from_${VAR}/pop${POP}

      echo ${OUTFOLDER}

      # calculate LD between QTLs and predictors
      echo "    calculating LD"
      for chr in {1..10}; do
        run_pipeline.pl -Xmx10g -importGuess data/usda_rils_projected-SVs-SNPs.chr${chr}.poly.hmp.txt \
                        -includeSiteNamesInFile <(cat ${LDFOLDER}/list_QTNs.chr${chr}.causal-pop${POP}.txt ${OUTFOLDER}/markers_${LD}-ld_to_QTNs.txt) \
                        -export ${OUTFOLDER}/qtn-pred.chr${chr} \
                        -exportType Plink > /dev/null &
      done
      wait
      for chr in {1..10}; do
        if [[ -s ${OUTFOLDER}/qtn-pred.chr${chr}.plk.map ]]; then
          plink --file ${OUTFOLDER}/qtn-pred.chr${chr}.plk \
                --make-founders --r2 dprime with-freqs --ld-window-r2 0 \
                --ld-window 1000000 --ld-window-kb 350000 --geno 0.25 \
                --ld-snp-list ${LDFOLDER}/list_QTNs.chr${chr}.causal-pop${POP}.txt \
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


      # remove predictors in above the LD threshold to any other QTL
      # (grep QTLs from ld file | grep predictors | filter by wrong LD > predictor_to_exclude.txt)
      echo "    removing predictors above/below LD threshold"
      if [[ ${LD} == low ]]; then

        # identify
        grep -F -w -f <(cat ${LDFOLDER}/list_QTNs.chr*.causal-pop${POP}.txt) ${OUTFOLDER}/qtn-pred.ld | grep -F -w -f ${OUTFOLDER}/markers_${LD}-ld_to_QTNs.txt | awk '$9 > 0.5' | sed 's/^ *//' | tr -s " " | tr " " "\t" | cut -f 3,7 | sed 's/\t/\n/g' | sort | uniq | grep -v -F -w -f <(cat ${LDFOLDER}/list_QTNs.chr*.causal-pop${POP}.txt) > ${OUTFOLDER}/predictors_to_exclude.txt
        # remove
        comm -13 <(sort ${OUTFOLDER}/predictors_to_exclude.txt | uniq) <(sort ${OUTFOLDER}/markers_${LD}-ld_to_QTNs.txt | uniq) > ${OUTFOLDER}/predictors_${LD}-ld_to_QTNs.txt

      elif [[ ${LD} == moderate ]]; then

        # identify
        grep -F -w -f <(cat ${LDFOLDER}/list_QTNs.chr*.causal-pop${POP}.txt) ${OUTFOLDER}/qtn-pred.ld | grep -F -w -f ${OUTFOLDER}/markers_${LD}-ld_to_QTNs.txt | awk '$9 >= 0.9' | sed 's/^ *//' | tr -s " " | tr " " "\t" | cut -f 3,7 | sed 's/\t/\n/g' | sort | uniq | grep -v -F -w -f <(cat ${LDFOLDER}/list_QTNs.chr*.causal-pop${POP}.txt) > ${OUTFOLDER}/predictors_to_exclude.txt
        # remove
        comm -13 <(sort ${OUTFOLDER}/predictors_to_exclude.txt | uniq) <(sort ${OUTFOLDER}/markers_${LD}-ld_to_QTNs.txt | uniq) > ${OUTFOLDER}/predictors_${LD}-ld_to_QTNs.txt

      elif [[ ${LD} == high ]]; then

        # identify
        grep -F -w -f <(cat ${LDFOLDER}/list_QTNs.chr*.causal-pop${POP}.txt) ${OUTFOLDER}/qtn-pred.ld | grep -F -w -f ${OUTFOLDER}/markers_${LD}-ld_to_QTNs.txt | awk '$9 >= 0.9' | sed 's/^ *//' | tr -s " " | tr " " "\t" | cut -f 3,7 | sed 's/\t/\n/g' | sort | uniq | grep -v -F -w -f <(cat ${LDFOLDER}/list_QTNs.chr*.causal-pop${POP}.txt) > ${OUTFOLDER}/predictors_${LD}-ld_to_QTNs.txt

      fi

      # randomly select markers to be predictors
      echo "    randomly selecting ${NMARKERS} markers"
      if [[ ${LD} == low ]]; then

        shuf -n ${NMARKERS} --random-source=<(get_seeded_random ${SEED}) ${OUTFOLDER}/predictors_${LD}-ld_to_QTNs.txt > ${OUTFOLDER}/predictors_${LD}-ld_to_QTNs.${NMARKERS}-markers.txt

      elif [[ ${LD} == moderate ]]; then

        # for each qtn, select one predictor in moderate LD, but before:
        #   create empty file
        cat /dev/null > ${OUTFOLDER}/predictors_${LD}-ld_to_QTNs.${NMARKERS}-markers.txt
        #   define seeds for each qtn before
        QTNSEED=($(shuf -i 1-100000 -n 1 --random-source=<(get_seeded_random ${SEED})))
        for qtn in $(cat ${LDFOLDER}/list_QTNs.chr*.causal-pop${POP}.txt); do
          # find all predictors in ld to that qtn
          grep -F -w ${qtn} ${OUTFOLDER}/qtn-pred.ld | grep -F -w -f ${OUTFOLDER}/predictors_${LD}-ld_to_QTNs.txt | awk '$9 > 0.5 && $9 < 0.9' | sed 's/^ *//' | tr -s " " | tr " " "\t" | sed 1d | cut -f 3,7 | sed 's/\t/\n/g' | sort | uniq | grep -v -F -w ${qtn} | shuf -n 1 --random-source=<(get_seeded_random ${QTNSEED}) >> ${OUTFOLDER}/predictors_${LD}-ld_to_QTNs.${NMARKERS}-markers.txt
          # change seed number
          QTNSEED=($(shuf -i 1-100000 -n 1 --random-source=<(get_seeded_random ${QTNSEED})))
        done
        # then, randomly select the rest of markers that are not in high LD to QTL
        SAMPLED=$(wc -l ${OUTFOLDER}/predictors_${LD}-ld_to_QTNs.${NMARKERS}-markers.txt | cut -d " " -f 1)
        REMAINING=$(( ${NMARKERS} - ${SAMPLED} ))
        shuf -n ${REMAINING} --random-source=<(get_seeded_random ${SEED}) analysis/ld_downsample/pred_low/rep${REP}/${QTN}-QTNs_from_${VAR}/pop${POP}/predictors_low-ld_to_QTNs.txt >> ${OUTFOLDER}/predictors_${LD}-ld_to_QTNs.${NMARKERS}-markers.txt


      elif [[ ${LD} == high ]]; then

        # for each qtn, select one predictor in high LD, but before:
        #   create empty file
        cat /dev/null > ${OUTFOLDER}/predictors_${LD}-ld_to_QTNs.${NMARKERS}-markers.txt
        #   define seeds for each qtn before
        QTNSEED=($(shuf -i 1-100000 -n 1 --random-source=<(get_seeded_random ${SEED})))
        for qtn in $(cat ${LDFOLDER}/list_QTNs.chr*.causal-pop${POP}.txt); do
          # find all predictors in ld to that qtn
          grep -F -w ${qtn} ${OUTFOLDER}/qtn-pred.ld | grep -F -w -f ${OUTFOLDER}/predictors_${LD}-ld_to_QTNs.txt | awk '$9 >= 0.9' | sed 's/^ *//' | tr -s " " | tr " " "\t" | sed 1d | cut -f 3,7 | sed 's/\t/\n/g' | sort | uniq | grep -v -F -w ${qtn} | shuf -n 1 --random-source=<(get_seeded_random ${QTNSEED}) >> ${OUTFOLDER}/predictors_${LD}-ld_to_QTNs.${NMARKERS}-markers.txt
          # change seed number
          QTNSEED=($(shuf -i 1-100000 -n 1 --random-source=<(get_seeded_random ${QTNSEED})))
        done
        # then, randomly select the rest of markers between low/moderate ld
        SAMPLED=$(wc -l ${OUTFOLDER}/predictors_${LD}-ld_to_QTNs.${NMARKERS}-markers.txt | cut -d " " -f 1)
        REMAINING=$(( ${NMARKERS} - ${SAMPLED} ))
        shuf -n ${REMAINING} --random-source=<(get_seeded_random ${SEED}) <(cat analysis/ld_downsample/pred_low/rep${REP}/${QTN}-QTNs_from_${VAR}/pop${POP}/predictors_low-ld_to_QTNs.txt analysis/ld_downsample/pred_moderate/rep${REP}/${QTN}-QTNs_from_${VAR}/pop${POP}/predictors_moderate-ld_to_QTNs.txt) >> ${OUTFOLDER}/predictors_${LD}-ld_to_QTNs.${NMARKERS}-markers.txt

      fi


      # filter hmp to have only these markers
      echo "    creating hapmap file with selected predictors"
      run_pipeline.pl -Xmx10g -importGuess analysis/ld_downsample/datasets/rep${REP}/usda_rils_projected-SVs-SNPs.poly.rep${REP}.with-LD-info.hmp.txt \
                      -includeSiteNamesInFile ${OUTFOLDER}/predictors_${LD}-ld_to_QTNs.${NMARKERS}-markers.txt \
                      -export ${OUTFOLDER}/usda_rils_projected-SVs-SNPs.${LD}-ld_to_QTNs.hmp.txt \
                      -exportType HapmapDiploid > /dev/null
      # use a different seed for next iteration
      SEED=($(shuf -i 1-100000 -n 1 --random-source=<(get_seeded_random ${SEED})))

    done
  done
done

echo "job finished @ $(date)"
