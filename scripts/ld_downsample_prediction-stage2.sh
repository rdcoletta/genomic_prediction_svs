#!/bin/bash

function usage {
  echo
  echo "Usage: $(basename $0) [-d INT] [-e INT] [-f INT] [-i INT] [-p INT] [-P INT] [-s INT] [-h] [-w]" 2>&1
  echo
  echo '    -d      choose the iteration of dataset with markers for prediction'
  echo '    -e      total number of environments to predict'
  echo '    -f      number of folds for cross-validation'
  echo '    -h      print help'
  echo '    -i      number of cross-validation iterations'
  echo '    -p      start predictions at this population number'
  echo '    -P      end predictions at this population number'
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
POPSTART=1
POPEND=20
SEED=2021
WEIGHT=

# get command-line optional arguments
optstring=":d:e:f:hi:p:P:s:w"

while getopts ${optstring} arg; do
  case ${arg} in
    d) [[ $OPTARG =~ ^[0-9]+$ ]] && NDATASET=${OPTARG} || warning ${arg} ;;
    e) [[ $OPTARG =~ ^[0-9]+$ ]] && NENVS=${OPTARG} || warning ${arg} ;;
    f) [[ $OPTARG =~ ^[0-9]+$ ]] && NFOLDS=${OPTARG} || warning ${arg} ;;
    h) usage ;;
    i) [[ $OPTARG =~ ^[0-9]+$ ]] && NITER=${OPTARG} || warning ${arg} ;;
    p) [[ $OPTARG =~ ^[0-9]+$ ]] && POPSTART=${OPTARG} || warning ${arg} ;;
    P) [[ $OPTARG =~ ^[0-9]+$ ]] && POPEND=${OPTARG} || warning ${arg} ;;
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
# for opt in NDATASET NENVS NFOLDS NITER POPSTART POPEND SEED WEIGHT; do
#   echo "  ${opt}=${!opt}"
# done
# echo "pops to evaluate:"
# echo "$(seq ${POPSTART} ${POPEND})"


# run predictions (stage 2)

# go to project folder
cd /home/candy/rafa/genomic_prediction/simulation

for PREDICTOR in low moderate high; do
  for H2 in 0.3 0.7; do
    for QTN in 100; do
      for VAR in SNP SV both; do
        for POP in $(seq ${POPSTART} ${POPEND}); do
          # get folder name
          FOLDER=analysis/ld_downsample/sim_traits/rep${NDATASET}/${QTN}-QTNs_from_${VAR}/${H2}-heritability/pop${POP}
          # get genotypic data
          GDATA=analysis/ld_downsample/pred_${PREDICTOR}/rep${NDATASET}/${QTN}-QTNs_from_${VAR}/pop${POP}/usda_rils_projected-SVs-SNPs.${PREDICTOR}-ld_to_QTNs.hmp.txt
          # get pheno data and set output folder
          if [[ -z ${WEIGHT} ]]; then
            # if using blups without weights
            PDATA=${FOLDER}/blups_1st_stage.txt
            OUT=${FOLDER}/pred_${PREDICTOR}_ld
            LOG=${FOLDER}/genomic_prediction_from_blups.${PREDICTOR}-ld.log
            mkdir -p ${OUT}
          else
            # if using blups with weights
            PDATA=${FOLDER}/blups_1st_stage_weighted.txt
            OUT=${FOLDER}/pred_${PREDICTOR}_ld_weighted
            LOG=${FOLDER}/genomic_prediction_from_blups_weighted.${PREDICTOR}-ld.log
            mkdir -p ${OUT}
          fi
          # echo "    ${FOLDER}"
          # run prediction with cv1
          echo "Rscript scripts/genomic_prediction_from_blups.R ${GDATA} ${PDATA} ${OUT} --cv-type=CV1 --n-folds=${NFOLDS} --cv-iter=${NITER} --total-envs=${NENVS} --seed=${SEED} ${WEIGHT} > ${LOG%.log}.CV1.log" >> scripts/commands_ld_downsample_prediction-stage2.txt
          # run prediction with cv2
          echo "Rscript scripts/genomic_prediction_from_blups.R ${GDATA} ${PDATA} ${OUT} --cv-type=CV2 --n-folds=${NFOLDS} --cv-iter=${NITER} --total-envs=${NENVS} --seed=${SEED} ${WEIGHT} > ${LOG%.log}.CV2.log" >> scripts/commands_ld_downsample_prediction-stage2.txt
        done
      done
    done
  done
done
