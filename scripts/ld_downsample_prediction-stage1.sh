#!/bin/bash

# go to project folder
cd /home/candy/rafa/genomic_prediction/simulation

echo "job started @ $(date)"

for H2 in 0.3 0.7; do
  for QTN in 100; do
    for VAR in SNP SV both; do
      for REP in $(seq 1 10); do
        for POP in {1..3}; do
          # get folder name
          FOLDER=analysis/ld_downsample/sim_traits/rep${REP}/${QTN}-QTNs_from_${VAR}/${H2}-heritability/pop${POP}
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

echo "job finished @ $(date)"
