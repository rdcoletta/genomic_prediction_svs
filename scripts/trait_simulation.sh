#!/bin/bash
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=100gb
#SBATCH -J trait_simulation
#SBATCH -o /home/hirschc1/della028/projects/genomic_prediction/simulation/analysis/trait_sim/MSI_dump/%x_%j.out
#SBATCH -e /home/hirschc1/della028/projects/genomic_prediction/simulation/analysis/trait_sim/MSI_dump/%x_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=della028@umn.edu
#SBATCH --no-requeue

module load R/3.6.0

# check for optional variables -- if not requested, set it blank

# SV options
[[ -z ${RATIO} ]] && RATIO="" || RATIO="--SNP-SV-ratio=${RATIO}"
[[ -z ${SVEFFECT} ]] && SVEFFECT="" || SVEFFECT="--SV-effect=${SVEFFECT}"
# mult env options
[[ -z ${GENCOR} ]] && GENCOR="" || GENCOR="--gen-cor-matrix=${GENCOR}"
[[ -z ${RESCOR} ]] && RESCOR="" || RESCOR="--res-cor-matrix=${RESCOR}"
[[ -z ${GXE} ]] && GXE="" || GXE="--${GXE}"
[[ -z ${DIFFDIST} ]] && DIFFDIST="" || DIFFDIST="--${DIFFDIST}"
[[ -z ${DIFFMEAN} ]] && DIFFMEAN="" || DIFFMEAN="--${DIFFMEAN}"
# report option
[[ -z ${QTNVAR} ]] && QTNVAR="" || QTNVAR="--${QTNVAR}"

# go to project folder
cd ~/projects/genomic_prediction/simulation

# simulate traits
Rscript scripts/trait_simulation.R ${HMP} ${SVS} ${FOLDER} \
                                   --causal-variant=${VAR} ${RATIO} \
                                   --pops=${POPS} \
                                   --reps=${REPS} \
                                   --envs=${ENVS} \
                                   --h2=${H2} \
                                   --model=${MODEL} \
                                   --add-QTN-num=${QTN} \
                                   --add-effect-type=${EFFECTTYPE} \
                                   --marker-effect=${EFFECTSIZE} ${SVEFFECT} \
                                   --architecture=pleiotropic \
                                   --seed=${SEED} \
                                   ${GENCOR} ${RESCOR} ${GXE} ${DIFFDIST} ${DIFFMEAN} ${QTNVAR}
