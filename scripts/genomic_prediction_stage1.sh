#!/bin/bash

# go to project folder
cd /home/candy/rafa/genomic_prediction/simulation

echo "job started"

echo "--- Part 1 @ $(date) ---"
for H2 in 0.3 0.7; do
  for QTN in 10 100; do
    for VAR in SNP SV; do
      for GXE in no with; do
        for POP in $(seq 1 20); do
          # get folder name
          FOLDER=analysis/trait_sim/multi_env/${GXE}_gxe/additive_model/equal_effects/${QTN}-QTNs_from_${VAR}/${H2}-heritability/pop${POP}
          INFILE=$(ls ${FOLDER}/Simulated_Data_3_Reps_Herit_*)
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
EFFECTSIZE=0.1
for H2 in 0.3 0.7; do
  for QTN in 10 100; do
    for GXE in no with; do
      for RATIO in 0.5 0.8; do
        for SVEFFECT in 0.1 0.2 0.5; do
          for POP in $(seq 1 20); do
            # get folder name
            FOLDER=analysis/trait_sim/multi_env/${GXE}_gxe/additive_model/equal_effects/${QTN}-QTNs_from_${VAR}/SNP-SV-ratio_${RATIO}/effects_SNP-${EFFECTSIZE}_SV-${SVEFFECT}/${H2}-heritability/pop${POP}
            INFILE=$(ls ${FOLDER}/Simulated_Data_3_Reps_Herit_*)
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
done

echo "--- Part 3 @ $(date) ---"
VAR=both
EFFECTSIZE=0.1
for H2 in 0.3 0.7; do
  for QTN in 10 100; do
    for GXE in with; do
      for RATIO in 0.5 0.8; do
        for SVEFFECT in 0.2 0.5; do
          for POP in $(seq 1 20); do
            # get folder name
            FOLDER=analysis/trait_sim/multi_env/${GXE}_gxe/additive_model/equal_effects/${QTN}-QTNs_from_${VAR}/SNP-SV-ratio_${RATIO}/effects_SNP-${EFFECTSIZE}_SV-${SVEFFECT}_diff-dist-gxe/${H2}-heritability/pop${POP}
            INFILE=$(ls ${FOLDER}/Simulated_Data_3_Reps_Herit_*)
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
done

echo "job finished @ $(date)"
