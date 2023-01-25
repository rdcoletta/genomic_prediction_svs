#!/bin/bash
#PBS -l walltime=24:00:00,nodes=1:ppn=5,mem=100gb
#PBS -o /home/hirschc1/della028/projects/genomic_prediction/simulation/analysis/trait_sim
#PBS -e /home/hirschc1/della028/projects/genomic_prediction/simulation/analysis/trait_sim
#PBS -V
#PBS -N predict_sim_trait_${EFFECT}_${MARKER}
#PBS -M della028@umn.edu
#PBS -m abe
#PBS -r n

module load R/3.6.0

# go to project folder
cd ~/projects/genomic_prediction/simulation

# check if single or multiple environments
[[ ${ENVS} == "single" ]] && ENVS_FOLDER='trait_sim' && OPT=''
[[ ${ENVS} == "multiple" ]] && ENVS_FOLDER='trait_sim_mult-env' && OPT='--multiple-envs'

echo "Options selected: ${EFFECT} ${MARKER} ${REPS} ${ENVS} ${ENVS_FOLDER} ${OPT}"
echo ""

# predict simulated traits
for var_source in SNP SV both; do
  for qtn in 3 25 75; do
    for h2 in 0.2 0.5 0.9; do
      for (( rep=1; rep<=${REPS}; rep++ )); do

        echo "--- Causative variant: ${var_source} | QTNs: ${qtn} | Heritability: ${h2} | Rep: ${rep} ---"
        # set up arguments
        GENO=analysis/trait_sim/datasets/usda_rils.${MARKER}_markers.adjusted-n-markers.hmp.txt
        TRAIT=analysis/${ENVS_FOLDER}/additive_model/${EFFECT}/${qtn}-QTNs_from_${var_source}/${h2}-heritability/rep${rep}/Simulated_Data_1_Reps_Herit_${h2}*.txt
        FOLDER=analysis/${ENVS_FOLDER}/additive_model/${EFFECT}/${qtn}-QTNs_from_${var_source}/gs_with_max-number_markers/${MARKER}_markers/${h2}-heritability/rep${rep}
        SEED=677
        FOLDS=5
        # run script
        Rscript scripts/predict_sim_trait.R ${GENO} ${TRAIT} ${FOLDER} --number-folds=${FOLDS} --seed=${SEED} ${OPT}
        echo ""

      done
    done
  done
done
