#!/bin/bash
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=10
#SBATCH --mem=100gb
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

for QTN in 10 100; do
  for VAR in SNP SV; do
    for LD in low moderate high; do

      # get folder with original LD file and QTN list
      LDFOLDER=analysis/ld_downsample/sim_traits/rep${REP}/${QTN}-QTNs_from_${VAR}/0.3-heritability/pop${POP}
      # get folder with results for different LD levels
      OUTLOW=analysis/ld_downsample/pred_low/rep${REP}/${QTN}-QTNs_from_${VAR}/pop${POP}
      OUTMOD=analysis/ld_downsample/pred_moderate/rep${REP}/${QTN}-QTNs_from_${VAR}/pop${POP}
      OUTHIGH=analysis/ld_downsample/pred_high/rep${REP}/${QTN}-QTNs_from_${VAR}/pop${POP}

      if [[ ${LD} == low ]]; then

        echo "  ${OUTLOW}"

        # randomly select markers to be predictors
        shuf -n ${NMARKERS} --random-source=<(get_seeded_random ${SEED}) ${OUTLOW}/markers_${LD}-ld_to_QTNs.txt > ${OUTLOW}/predictors_${LD}-ld_to_QTNs.txt

      elif [[ ${LD} == moderate ]]; then

        echo "  ${OUTMOD}"

        # randomly select markers to be predictors
        shuf -n ${NMARKERS} --random-source=<(get_seeded_random ${SEED}) ${OUTMOD}/markers_${LD}-ld_to_QTNs.txt > ${OUTMOD}/predictors_${LD}-ld_to_QTNs.txt

      elif [[ ${LD} == high ]]; then

        echo "  ${OUTHIGH}"

        # for each qtn, select one predictor in high LD, but before:
        #   create empty file
        cat /dev/null > ${OUTHIGH}/predictors_${LD}-ld_to_QTNs.txt
        #   define seeds for each qtn before
        QTNSEED=($(shuf -i 1-100000 -n 1 --random-source=<(get_seeded_random ${SEED})))
        for qtn in $(cat ${LDFOLDER}/list_QTNs.chr*.causal-pop${POP}.txt); do
          # find all predictors in ld to that qtn
          sed 's/^ *//' ${LDFOLDER}/ld_info.QTNs-rep${REP}-pop${POP}.${LD}.ld | tr -s " " | tr " " "\t" | sed 1d | cut -f 3,7 | grep ${qtn} | sed 's/\t/\n/g' | sort | uniq | grep -v -F ${qtn} | shuf -n 1 --random-source=<(get_seeded_random ${QTNSEED}) >> ${OUTHIGH}/predictors_${LD}-ld_to_QTNs.txt
          # change seed number
          QTNSEED=($(shuf -i 1-100000 -n 1 --random-source=<(get_seeded_random ${QTNSEED})))
        done
        # then, randomly select the rest of markers between low/moderate ld
        REMAINING=$(( ${NMARKERS} - ${QTN} ))
        shuf -n ${REMAINING} --random-source=<(get_seeded_random ${SEED}) <(cat ${OUTLOW}/markers_low-ld_to_QTNs.txt ${OUTMOD}/markers_moderate-ld_to_QTNs.txt) >> ${OUTHIGH}/predictors_${LD}-ld_to_QTNs.txt

      fi

      # correct output folder name to save plink and hmp files
      OUTFOLDER=analysis/ld_downsample/pred_${LD}/rep${REP}/${QTN}-QTNs_from_${VAR}/pop${POP}

      # calculate LD between QTLs and predictors
      echo "    calculating LD"
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

      # remove predictors in above the LD threshold to any other QTL
      # (grep QTLs from ld file | grep predictors | filter by wrong LD > predictor_to_exclude.txt)
      echo "    removing predictors above LD threshold"
      if [[ ${LD} == low ]]; then
        # identify
        grep -F -f <(cat ${LDFOLDER}/list_QTNs.chr*.causal-pop${POP}.txt) ${OUTFOLDER}/qtn-pred.ld | grep -F -f ${OUTFOLDER}/predictors_${LD}-ld_to_QTNs.txt | awk '$9 > 0.5' | sed 's/^ *//' | tr -s " " | tr " " "\t" | cut -f 3,7 | sed 's/\t/\n/g' | sort | uniq | grep -v -F -f <(cat ${LDFOLDER}/list_QTNs.chr*.causal-pop${POP}.txt) > ${OUTFOLDER}/predictors_to_exclude.txt
        # remove
        grep -v -f ${OUTFOLDER}/predictors_to_exclude.txt ${OUTFOLDER}/predictors_${LD}-ld_to_QTNs.txt > ${OUTFOLDER}/tmp_predictors_list
        mv ${OUTFOLDER}/tmp_predictors_list ${OUTFOLDER}/predictors_${LD}-ld_to_QTNs.txt
      elif [[ ${LD} == moderate ]]; then
        # identify
        grep -F -f <(cat ${LDFOLDER}/list_QTNs.chr*.causal-pop${POP}.txt) ${OUTFOLDER}/qtn-pred.ld | grep -F -f ${OUTFOLDER}/predictors_${LD}-ld_to_QTNs.txt | awk '$9 > 0.9' | sed 's/^ *//' | tr -s " " | tr " " "\t" | cut -f 3,7 | sed 's/\t/\n/g' | sort | uniq | grep -v -F -f <(cat ${LDFOLDER}/list_QTNs.chr*.causal-pop${POP}.txt) > ${OUTFOLDER}/predictors_to_exclude.txt
        # remove
        grep -v -f ${OUTFOLDER}/predictors_to_exclude.txt ${OUTFOLDER}/predictors_${LD}-ld_to_QTNs.txt > ${OUTFOLDER}/tmp_predictors_list
        mv ${OUTFOLDER}/tmp_predictors_list ${OUTFOLDER}/predictors_${LD}-ld_to_QTNs.txt
      fi

      # filter hmp to have only these markers
      run_pipeline.pl -Xmx10g -importGuess analysis/ld_downsample/datasets/rep${REP}/usda_rils_projected-SVs-SNPs.poly.rep${REP}.with-LD-info.hmp.txt \
                      -includeSiteNamesInFile ${OUTFOLDER}/predictors_${LD}-ld_to_QTNs.txt \
                      -export ${OUTFOLDER}/usda_rils_projected-SVs-SNPs.${LD}-ld_to_QTNs.hmp.txt \
                      -exportType HapmapDiploid > /dev/null
      # use a different seed for next iteration
      SEED=($(shuf -i 1-100000 -n 1 --random-source=<(get_seeded_random ${SEED})))

    done
  done
done

echo "job finished @ $(date)"
