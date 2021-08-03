#!/bin/bash

function usage {
  echo
  echo "Usage: $(basename $0) [-d INT] [-e INT] [-f INT] [-i INT] [-s INT] [-h] [-w]" 2>&1
  echo
  echo '    -d      choose the iteration of dataset with markers for prediction'
  echo '    -e      total number of environments to predict'
  echo '    -f      number of folds for cross-validation'
  echo '    -h      print help'
  echo '    -i      number of cross-validation iterations'
  echo '    -s      seed number'
  echo '    -w      add this flag to use weighted BLUPs from 1st stage'
  echo
  exit 1
}

function warning {
  echo
  echo "Argument for -$1 should be an integer" 2>&1
  echo
  exit 1
}

# set default options
NDATASET=1
NENVS=5
NFOLDS=5
NITER=3
SEED=2021
WEIGHT=

# get command-line optional arguments
optstring=":d:e:f:i:hs:w"

while getopts ${optstring} arg; do
  case ${arg} in
    d) [[ $OPTARG =~ ^[0-9]+$ ]] && NDATASET=${OPTARG} || warning ${arg} ;;
    e) [[ $OPTARG =~ ^[0-9]+$ ]] && NENVS=${OPTARG} || warning ${arg} ;;
    f) [[ $OPTARG =~ ^[0-9]+$ ]] && NFOLDS=${OPTARG} || warning ${arg} ;;
    h) usage ;;
    i) [[ $OPTARG =~ ^[0-9]+$ ]] && NITER=${OPTARG} || warning ${arg} ;;
    s) [[ $OPTARG =~ ^[0-9]+$ ]] && SEED=${OPTARG} || warning ${arg} ;;
    w) WEIGHT=--envs-weight ;;

    :) echo "Error: -${OPTARG} requires an argument"; usage ;;
    *) echo "Unknown option: -${OPTARG}"; usage ;;

  esac
done

# shift argument index so positional arguments can start at $1
shift $(($OPTIND - 1))

# assert positional arguments
if [ $# -ge 1 ]; then
    echo -e "Ignoring positional arguments (none needed)\n"
fi

# # debug
# echo -e "positional args:\n  ${@}"
# echo "optional args:"
# for opt in NDATASET NENVS NFOLDS NITER SEED WEIGHT; do
#   echo "  ${opt}=${!opt}"
# done


# run predictions (stage 2)

# go to project folder
cd /home/candy/rafa/genomic_prediction/simulation

echo "job started"

echo "--- Part 1 @ $(date) (all_markers, h2 0.2, qtn 10 100, gxe no with already completed)---"
for PREDICTOR in all_markers; do
  for H2 in 0.2; do
    for QTN in 200; do
      for VAR in SNP SV; do
        for GXE in no with; do
          for POP in $(seq 1 3); do
            # get folder name
            FOLDER=analysis/trait_sim/multi_env/${GXE}_gxe/additive_model/equal_effects/${QTN}-QTNs_from_${VAR}/${H2}-heritability/pop${POP}
            # get genotypic data
            GDATA=analysis/trait_sim/datasets/iter${NDATASET}/usda_rils.${PREDICTOR}.adjusted-n-markers.hmp.txt
            # get pheno data and set output folder
            if [[ -z ${WEIGHT} ]]; then
              # if using blups without weights
              PDATA=${FOLDER}/blups_1st_stage.txt
              OUT=${FOLDER}/prediction_${PREDICTOR}
              LOG=${FOLDER}/genomic_prediction_from_blups.log
            else
              # if using blups with weights
              PDATA=${FOLDER}/blups_1st_stage_weighted.txt
              OUT=${FOLDER}/prediction_${PREDICTOR}_weighted
              LOG=${FOLDER}/genomic_prediction_from_blups_weighted.log
            fi
            echo "  ${FOLDER}"
            # run prediction with cv1
            Rscript scripts/genomic_prediction_from_blups.R ${GDATA} ${PDATA} ${OUT} --cv-type=CV1 --n-folds=${NFOLDS} --cv-iter=${NITER} --total-envs=${NENVS} --seed=${SEED} ${WEIGHT} > ${LOG%.log}.CV1.log 2>&1 &
            # run prediction with cv2
            Rscript scripts/genomic_prediction_from_blups.R ${GDATA} ${PDATA} ${OUT} --cv-type=CV2 --n-folds=${NFOLDS} --cv-iter=${NITER} --total-envs=${NENVS} --seed=${SEED} ${WEIGHT} > ${LOG%.log}.CV2.log 2>&1 &
            wait
          done
        done
      done
    done
  done
done

# echo "--- Part 1 @ $(date) ---"
# for PREDICTOR in all_markers snp_ld_markers sv_markers snp_markers snp_not_ld_markers; do
#   for H2 in 0.2 0.5 0.9; do
#     for QTN in 10 100 200; do
#       for VAR in SNP SV; do
#         for GXE in no with; do
#           for POP in $(seq 1 3); do
#             # get folder name
#             FOLDER=analysis/trait_sim/multi_env/${GXE}_gxe/additive_model/equal_effects/${QTN}-QTNs_from_${VAR}/${H2}-heritability/pop${POP}
#             # get genotypic data
#             GDATA=analysis/trait_sim/datasets/iter${NDATASET}/usda_rils.${PREDICTOR}.adjusted-n-markers.hmp.txt
#             # get pheno data and set output folder
#             if [[ -z ${WEIGHT} ]]; then
#               # if using blups without weights
#               PDATA=${FOLDER}/blups_1st_stage.txt
#               OUT=${FOLDER}/prediction_${PREDICTOR}
#               LOG=${FOLDER}/genomic_prediction_from_blups.log
#             else
#               # if using blups with weights
#               PDATA=${FOLDER}/blups_1st_stage_weighted.txt
#               OUT=${FOLDER}/prediction_${PREDICTOR}_weighted
#               LOG=${FOLDER}/genomic_prediction_from_blups_weighted.log
#             fi
#             echo "  ${FOLDER}"
#             # run prediction with cv1
#             Rscript scripts/genomic_prediction_from_blups.R ${GDATA} ${PDATA} ${OUT} --cv-type=CV1 --n-folds=${NFOLDS} --cv-iter=${NITER} --total-envs=${NENVS} --seed=${SEED} ${WEIGHT} > ${LOG%.log}.CV1.log 2>&1 &
#             # run prediction with cv2
#             Rscript scripts/genomic_prediction_from_blups.R ${GDATA} ${PDATA} ${OUT} --cv-type=CV2 --n-folds=${NFOLDS} --cv-iter=${NITER} --total-envs=${NENVS} --seed=${SEED} ${WEIGHT} > ${LOG%.log}.CV2.log 2>&1 &
#             wait
#           done
#         done
#       done
#     done
#   done
# done


# echo "--- Part 2 @ $(date) ---"
# VAR=both
# EFFECTSIZE=0.1
# for H2 in 0.2 0.5 0.9; do
#   for QTN in 10 100 200; do
#     for GXE in no with; do
#       for RATIO in 0.5 0.8; do
#         for SVEFFECT in 0.1 0.2 0.5; do
#           for POP in $(seq 1 10); do
#             # get folder name
#             FOLDER=analysis/trait_sim/multi_env/${GXE}_gxe/additive_model/equal_effects/${QTN}-QTNs_from_${VAR}/SNP-SV-ratio_${RATIO}/effects_SNP-${EFFECTSIZE}_SV-${SVEFFECT}/${H2}-heritability/pop${POP}
#             INFILE=$(ls ${FOLDER}/Simulated_Data_3_Reps_Herit_*)
#             echo "  ${FOLDER}"
#             # get BLUPs -- without adding environmental weights
#             Rscript scripts/get_blups_per_env.R ${INFILE} ${FOLDER} > ${FOLDER}/get_blups_per_env.log 2>&1 &
#             # get BLUPs -- adding environmental weights
#             Rscript scripts/get_blups_per_env.R ${INFILE} ${FOLDER} --envs-weight > ${FOLDER}/get_blups_per_env_weighted.log 2>&1 &
#             wait
#           done
#         done
#       done
#     done
#   done
# done
#
# echo "--- Part 3 @ $(date) ---"
# VAR=both
# EFFECTSIZE=0.1
# for H2 in 0.2 0.5 0.9; do
#   for QTN in 10 100 200; do
#     for GXE in with; do
#       for RATIO in 0.5 0.8; do
#         for SVEFFECT in 0.2 0.5; do
#           for POP in $(seq 1 10); do
#             # get folder name
#             FOLDER=analysis/trait_sim/multi_env/${GXE}_gxe/additive_model/equal_effects/${QTN}-QTNs_from_${VAR}/SNP-SV-ratio_${RATIO}/effects_SNP-${EFFECTSIZE}_SV-${SVEFFECT}_diff-dist-gxe/${H2}-heritability/pop${POP}
#             INFILE=$(ls ${FOLDER}/Simulated_Data_3_Reps_Herit_*)
#             echo "  ${FOLDER}"
#             # get BLUPs -- without adding environmental weights
#             Rscript scripts/get_blups_per_env.R ${INFILE} ${FOLDER} > ${FOLDER}/get_blups_per_env.log 2>&1 &
#             # get BLUPs -- adding environmental weights
#             Rscript scripts/get_blups_per_env.R ${INFILE} ${FOLDER} --envs-weight > ${FOLDER}/get_blups_per_env_weighted.log 2>&1 &
#             wait
#           done
#         done
#       done
#     done
#   done
# done

echo "job finished @ $(date)"
