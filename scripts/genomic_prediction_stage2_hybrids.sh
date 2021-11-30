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

for PREDICTOR in all_markers snp_ld_markers sv_markers snp_markers snp_not_ld_markers; do
  for H2 in 0.3 0.7; do
    for QTN in 10 100; do
      for VAR in SNP SV; do
        for MODEL in A AD D; do
          # determine correct number of add/dom QTNs
          if [[ ${MODEL} == A ]]; then
            ADDQTN=${QTN}
            DOMQTN=0
            IMPUTEEFFECT=Add
          elif [[ ${MODEL} == AD ]]; then
            ADDQTN=$(( ${QTN} / 2 ))
            DOMQTN=$(( ${QTN} / 2 ))
            IMPUTEEFFECT=Dom
          elif [[ ${MODEL} == D ]]; then
            ADDQTN=0
            DOMQTN=${QTN}
            IMPUTEEFFECT=Dom
          fi
          for POP in $(seq ${POPSTART} ${POPEND}); do
            # get folder name
            FOLDER=analysis/trait_sim_hybrids/multi_env/with_gxe/${MODEL}_model/${ADDQTN}-add-QTNs_${DOMQTN}-dom-QTNs_from_${VAR}/${H2}-heritability/pop${POP}
            # get genotypic data
            GDATA=analysis/trait_sim_hybrids/datasets/iter${NDATASET}/usda_hybrids.${PREDICTOR}.adjusted-n-markers.hmp.txt
            # get pheno data and set output folder
            if [[ -z ${WEIGHT} ]]; then
              # if using blups without weights
              PDATA=${FOLDER}/blups_1st_stage.txt
              OUT=${FOLDER}/prediction_iter${NDATASET}/${PREDICTOR}
              LOG=${FOLDER}/prediction_iter${NDATASET}/genomic_prediction_from_blups.${PREDICTOR}.log
              mkdir -p ${FOLDER}/prediction_iter${NDATASET}
            else
              # if using blups with weights
              PDATA=${FOLDER}/blups_1st_stage_weighted.txt
              OUT=${FOLDER}/prediction_iter${NDATASET}_weighted/${PREDICTOR}
              LOG=${FOLDER}/prediction_iter${NDATASET}_weighted/genomic_prediction_from_blups_weighted.${PREDICTOR}.log
              mkdir -p ${FOLDER}/prediction_iter${NDATASET}_weighted
            fi
            # echo "    ${FOLDER}"
            # run prediction with cv1
            echo "Rscript scripts/genomic_prediction_from_blups_hybrids.R ${GDATA} ${PDATA} ${OUT} --cv-type=CV1 --n-folds=${NFOLDS} --cv-iter=${NITER} --total-envs=${NENVS} --impute-effect=${IMPUTEEFFECT} --impute-type=Middle --seed=${SEED} ${WEIGHT} > ${LOG%.log}.CV1.log" >> scripts/commands_genomic_prediction_stage2_hybrids.txt
            # run prediction with cv2
            echo "Rscript scripts/genomic_prediction_from_blups_hybrids.R ${GDATA} ${PDATA} ${OUT} --cv-type=CV2 --n-folds=${NFOLDS} --cv-iter=${NITER} --total-envs=${NENVS} --impute-effect=${IMPUTEEFFECT} --impute-type=Middle --seed=${SEED} ${WEIGHT} > ${LOG%.log}.CV2.log" >> scripts/commands_genomic_prediction_stage2_hybrids.txt
          done
        done
      done
    done
  done
done

VAR=both
RATIO=0.5
ADDEFFECT=0.1
SVEFFECT=0.1
for PREDICTOR in all_markers snp_ld_markers sv_markers snp_markers snp_not_ld_markers; do
  for H2 in 0.3 0.7; do
    for QTN in 10 100; do
      for MODEL in A AD D; do
        # determine correct number of add/dom QTNs
        if [[ ${MODEL} == A ]]; then
          ADDQTN=${QTN}
          DOMQTN=0
          IMPUTEEFFECT=Add
        elif [[ ${MODEL} == AD ]]; then
          ADDQTN=$(( ${QTN} / 2 ))
          DOMQTN=$(( ${QTN} / 2 ))
          IMPUTEEFFECT=Dom
        elif [[ ${MODEL} == D ]]; then
          ADDQTN=0
          DOMQTN=${QTN}
          IMPUTEEFFECT=Dom
        fi
        for POP in $(seq ${POPSTART} ${POPEND}); do
          # get folder name
          FOLDER=analysis/trait_sim_hybrids/multi_env/with_gxe/${MODEL}_model/${ADDQTN}-add-QTNs_${DOMQTN}-dom-QTNs_from_${VAR}/SNP-SV-ratio_${RATIO}/effects_SNP-${ADDEFFECT}_SV-${SVEFFECT}/${H2}-heritability/pop${POP}
          # get genotypic data
          GDATA=analysis/trait_sim_hybrids/datasets/iter${NDATASET}/usda_hybrids.${PREDICTOR}.adjusted-n-markers.hmp.txt
          # get pheno data and set output folder
          if [[ -z ${WEIGHT} ]]; then
            # if using blups without weights
            PDATA=${FOLDER}/blups_1st_stage.txt
            OUT=${FOLDER}/prediction_iter${NDATASET}/${PREDICTOR}
            LOG=${FOLDER}/prediction_iter${NDATASET}/genomic_prediction_from_blups.${PREDICTOR}.log
            mkdir -p ${FOLDER}/prediction_iter${NDATASET}
          else
            # if using blups with weights
            PDATA=${FOLDER}/blups_1st_stage_weighted.txt
            OUT=${FOLDER}/prediction_iter${NDATASET}_weighted/${PREDICTOR}
            LOG=${FOLDER}/prediction_iter${NDATASET}_weighted/genomic_prediction_from_blups_weighted.${PREDICTOR}.log
            mkdir -p ${FOLDER}/prediction_iter${NDATASET}_weighted
          fi
          # echo "    ${FOLDER}"
          # run prediction with cv1
          echo "Rscript scripts/genomic_prediction_from_blups_hybrids.R ${GDATA} ${PDATA} ${OUT} --cv-type=CV1 --n-folds=${NFOLDS} --cv-iter=${NITER} --total-envs=${NENVS} --impute-effect=${IMPUTEEFFECT} --impute-type=Middle --seed=${SEED} ${WEIGHT} > ${LOG%.log}.CV1.log" >> scripts/commands_genomic_prediction_stage2_hybrids.txt
          # run prediction with cv2
          echo "Rscript scripts/genomic_prediction_from_blups_hybrids.R ${GDATA} ${PDATA} ${OUT} --cv-type=CV2 --n-folds=${NFOLDS} --cv-iter=${NITER} --total-envs=${NENVS} --impute-effect=${IMPUTEEFFECT} --impute-type=Middle --seed=${SEED} ${WEIGHT} > ${LOG%.log}.CV2.log" >> scripts/commands_genomic_prediction_stage2_hybrids.txt
        done
      done
    done
  done
done
