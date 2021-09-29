#!/bin/bash
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=10
#SBATCH --mem=100gb
#SBATCH -J ld_predictor_causal-vars
#SBATCH -o /home/hirschc1/della028/projects/genomic_prediction/simulation/analysis/trait_sim/MSI_dump/%x_%j.out
#SBATCH -e /home/hirschc1/della028/projects/genomic_prediction/simulation/analysis/trait_sim/MSI_dump/%x_%j.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=della028@umn.edu
#SBATCH --no-requeue

module load R/3.6.0

# go to project folder
cd ~/projects/genomic_prediction/simulation

echo "job started"

echo "--- Part 1 @ $(date) ---"
H2=0.3
GXE=with
for QTN in 10 100; do
  for VAR in SNP SV; do
    # get folder name
    FOLDER=analysis/trait_sim/multi_env/${GXE}_gxe/additive_model/equal_effects/${QTN}-QTNs_from_${VAR}/${H2}-heritability/pop${POP}
    echo "  ${FOLDER}"
    # get hapmap file with predictors
    PREDFILE=analysis/trait_sim/datasets/iter${NDATASET}/usda_rils.${PREDICTOR}.adjusted-n-markers.hmp.txt
    # get file with causative variants
    CAUSALFILE=${FOLDER}/Additive_Selected_QTNs.txt
    # create folder to save files
    OUTFOLDER=analysis/trait_sim/multi_env/${GXE}_gxe/additive_model/equal_effects/${QTN}-QTNs_from_${VAR}/ld_causative-vars_predictors/pop${POP}
    mkdir -p ${OUTFOLDER}
    # get markers to calculate ld
    cat <(cut -f 1 ${PREDFILE} | sed 1d) <(cut -f 2 ${CAUSALFILE} | sed 1d) > ${OUTFOLDER}/markers_for_ld.${PREDICTOR}.pred-iter${NDATASET}.causal-pop${POP}.txt
    # filter hmp file
    echo "    filtering hmp file"
    for chr in $(seq 1 10); do
      run_pipeline.pl -Xmx10g -importGuess data/usda_rils_projected-SVs-SNPs.chr${chr}.poly.hmp.txt \
                      -includeSiteNamesInFile ${OUTFOLDER}/markers_for_ld.${PREDICTOR}.pred-iter${NDATASET}.causal-pop${POP}.txt \
                      -export ${OUTFOLDER}/${PREDICTOR}.pred-iter${NDATASET}.causal-pop${POP}.chr${chr}.hmp.txt \
                      -exportType HapmapDiploid > /dev/null &
    done
    wait
    # transform hmp into plink format
    echo "    transforming hmp to plink format"
    for chr in $(seq 1 10); do
      run_pipeline.pl -Xmx10g -importGuess ${OUTFOLDER}/${PREDICTOR}.pred-iter${NDATASET}.causal-pop${POP}.chr${chr}.hmp.txt \
                      -export ${OUTFOLDER}/${PREDICTOR}.pred-iter${NDATASET}.causal-pop${POP}.chr${chr} -exportType Plink > /dev/null &
    done
    wait
    # calculate LD
    echo "    calculating LD"
    for chr in $(seq 1 10); do
      plink --file ${OUTFOLDER}/${PREDICTOR}.pred-iter${NDATASET}.causal-pop${POP}.chr${chr}.plk \
            --make-founders --r2 gz dprime with-freqs --ld-window-r2 0 \
            --ld-window 10000 --ld-window-kb 350000 --geno 0.25 \
            --out ${OUTFOLDER}/${PREDICTOR}.pred-iter${NDATASET}.causal-pop${POP}.chr${chr} > /dev/null &
    done
    wait
    # remove unnecessary plink files
    for chr in $(seq 1 10); do
      rm ${OUTFOLDER}/${PREDICTOR}.pred-iter${NDATASET}.causal-pop${POP}.chr${chr}.log
      rm ${OUTFOLDER}/${PREDICTOR}.pred-iter${NDATASET}.causal-pop${POP}.chr${chr}.plk.*
      rm ${OUTFOLDER}/${PREDICTOR}.pred-iter${NDATASET}.causal-pop${POP}.chr${chr}.nosex
    done
    # define LD and outfile names
    LDFILE=${OUTFOLDER}/${PREDICTOR}.pred-iter${NDATASET}.causal-pop${POP}.chr${chr}.ld.gz
    OUTFILE=${OUTFOLDER}/ld_summary.${PREDICTOR}.pred-iter${NDATASET}.causal-pop${POP}.txt
    # get LD
    echo "    summarizing LD"
    Rscript scripts/ld_predictor_causal-vars.R ${PREDFILE} ${CAUSALFILE} ${LDFILE} ${OUTFILE}
  done
done
echo ""

echo "--- Part 2 @ $(date) ---"
VAR=both
EFFECTSIZE=0.1
SVEFFECT=0.1
for QTN in 10 100; do
  for RATIO in 0.5 0.8; do
    # get folder name
    FOLDER=analysis/trait_sim/multi_env/${GXE}_gxe/additive_model/equal_effects/${QTN}-QTNs_from_${VAR}/SNP-SV-ratio_${RATIO}/effects_SNP-${EFFECTSIZE}_SV-${SVEFFECT}/${H2}-heritability/pop${POP}
    echo "  ${FOLDER}"
    # get hapmap file with predictors
    PREDFILE=analysis/trait_sim/datasets/iter${NDATASET}/usda_rils.${PREDICTOR}.adjusted-n-markers.hmp.txt
    # get file with causative variants
    CAUSALFILE=${FOLDER}/Additive_Selected_QTNs.txt
    # create folder to save files
    OUTFOLDER=analysis/trait_sim/multi_env/${GXE}_gxe/additive_model/equal_effects/${QTN}-QTNs_from_${VAR}/SNP-SV-ratio_${RATIO}/ld_causative-vars_predictors/pop${POP}
    mkdir -p ${OUTFOLDER}
    # get markers to calculate ld
    cat <(cut -f 1 ${PREDFILE} | sed 1d) <(cut -f 2 ${CAUSALFILE} | sed 1d) > ${OUTFOLDER}/markers_for_ld.${PREDICTOR}.pred-iter${NDATASET}.causal-pop${POP}.txt
    # filter hmp file
    echo "    filtering hmp file"
    for chr in $(seq 1 10); do
      run_pipeline.pl -Xmx10g -importGuess data/usda_rils_projected-SVs-SNPs.chr${chr}.poly.hmp.txt \
                      -includeSiteNamesInFile ${OUTFOLDER}/markers_for_ld.${PREDICTOR}.pred-iter${NDATASET}.causal-pop${POP}.txt \
                      -export ${OUTFOLDER}/${PREDICTOR}.pred-iter${NDATASET}.causal-pop${POP}.chr${chr}.hmp.txt \
                      -exportType HapmapDiploid > /dev/null &
    done
    wait
    ## transform hmp into plink format
    echo "    transforming hmp to plink format"
    for chr in $(seq 1 10); do
      run_pipeline.pl -Xmx10g -importGuess ${OUTFOLDER}/${PREDICTOR}.pred-iter${NDATASET}.causal-pop${POP}.chr${chr}.hmp.txt \
                      -export ${OUTFOLDER}/${PREDICTOR}.pred-iter${NDATASET}.causal-pop${POP}.chr${chr} -exportType Plink > /dev/null &
    done
    wait
    # calculate LD
    echo "    calculating LD"
    for chr in $(seq 1 10); do
      plink --file ${OUTFOLDER}/${PREDICTOR}.pred-iter${NDATASET}.causal-pop${POP}.chr${chr}.plk \
            --make-founders --r2 gz dprime with-freqs --ld-window-r2 0 \
            --ld-window 10000 --ld-window-kb 350000 --geno 0.25 \
            --out ${OUTFOLDER}/${PREDICTOR}.pred-iter${NDATASET}.causal-pop${POP}.chr${chr} > /dev/null &
    done
    wait
    # remove unnecessary plink files
    for chr in $(seq 1 10); do
      rm ${OUTFOLDER}/${PREDICTOR}.pred-iter${NDATASET}.causal-pop${POP}.chr${chr}.log
      rm ${OUTFOLDER}/${PREDICTOR}.pred-iter${NDATASET}.causal-pop${POP}.chr${chr}.plk.*
      rm ${OUTFOLDER}/${PREDICTOR}.pred-iter${NDATASET}.causal-pop${POP}.chr${chr}.nosex
    done
    # define LD and outfile names
    LDFILE=${OUTFOLDER}/${PREDICTOR}.pred-iter${NDATASET}.causal-pop${POP}.chr${chr}.ld.gz
    OUTFILE=${OUTFOLDER}/ld_summary.${PREDICTOR}.pred-iter${NDATASET}.causal-pop${POP}.txt
    # get LD
    echo "    summarizing LD"
    Rscript scripts/ld_predictor_causal-vars.R ${PREDFILE} ${CAUSALFILE} ${LDFILE} ${OUTFILE}
  done
done

echo "job finished @ $(date)"
