#!/bin/bash

# go to project folder
cd /home/candy/rafa/genomic_prediction/simulation

echo "job started"

echo "--- Part 1 @ $(date) ---"
for H2 in 0.3 0.7; do
  for QTN in 10 100; do
    for VAR in SNP SV; do
      for MODEL in A AD D; do
        # determine correct number of add/dom QTNs
        if [[ ${MODEL} == A ]]; then
          ADDQTN=${QTN}
          DOMQTN=0
        elif [[ ${MODEL} == AD ]]; then
          ADDQTN=$(( ${QTN} / 2 ))
          DOMQTN=$(( ${QTN} / 2 ))
        elif [[ ${MODEL} == D ]]; then
          ADDQTN=0
          DOMQTN=${QTN}
        fi
        for POP in $(seq 1 20); do
          # get folder name
          FOLDER=analysis/trait_sim_hybrids/multi_env/with_gxe/${MODEL}_model/${ADDQTN}-add-QTNs_${DOMQTN}-dom-QTNs_from_${VAR}/${H2}-heritability/pop${POP}
          INFILE=$(ls ${FOLDER}/Simulated_Data_*)
          echo "  ${FOLDER}"
          # get BLUPs -- without adding environmental weights
          Rscript scripts/get_blups_per_env.R ${INFILE} ${FOLDER} > ${FOLDER}/get_blups_per_env.log 2>&1 &
          # get BLUPs -- adding environmental weights
          Rscript scripts/get_blups_per_env.R ${INFILE} ${FOLDER} --envs-weight > ${FOLDER}/get_blups_per_env_weighted.log 2>&1 &
          wait
        done
      done
    done
  done
done


echo "--- Part 2 @ $(date) ---"
VAR=both
RATIO=0.5
ADDEFFECT=0.1
SVEFFECT=0.1
for H2 in 0.3 0.7; do
  for QTN in 10 100; do
    for MODEL in A AD D; do
      # determine correct number of add/dom QTNs
      if [[ ${MODEL} == A ]]; then
        ADDQTN=${QTN}
        DOMQTN=0
      elif [[ ${MODEL} == AD ]]; then
        ADDQTN=$(( ${QTN} / 2 ))
        DOMQTN=$(( ${QTN} / 2 ))
      elif [[ ${MODEL} == D ]]; then
        ADDQTN=0
        DOMQTN=${QTN}
      fi
      for POP in $(seq 1 20); do
        # get folder name
        FOLDER=analysis/trait_sim_hybrids/multi_env/with_gxe/${MODEL}_model/${ADDQTN}-add-QTNs_${DOMQTN}-dom-QTNs_from_${VAR}/SNP-SV-ratio_${RATIO}/effects_SNP-${ADDEFFECT}_SV-${SVEFFECT}/${H2}-heritability/pop${POP}
        INFILE=$(ls ${FOLDER}/Simulated_Data_*)
        echo "  ${FOLDER}"
        # get BLUPs -- without adding environmental weights
        Rscript scripts/get_blups_per_env.R ${INFILE} ${FOLDER} > ${FOLDER}/get_blups_per_env.log 2>&1 &
        # get BLUPs -- adding environmental weights
        Rscript scripts/get_blups_per_env.R ${INFILE} ${FOLDER} --envs-weight > ${FOLDER}/get_blups_per_env_weighted.log 2>&1 &
        wait
      done
    done
  done
done

echo "job finished @ $(date)"
