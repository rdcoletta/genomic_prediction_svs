#!/bin/bash
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=80gb
#SBATCH -J qc_snp-sv_hmp
#SBATCH -o /home/hirschc1/della028/projects/genomic_prediction/simulation/analysis/trait_sim/MSI_dump/%x_%j.out
#SBATCH -e /home/hirschc1/della028/projects/genomic_prediction/simulation/analysis/trait_sim/MSI_dump/%x_%j.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=della028@umn.edu
#SBATCH --no-requeue


module load R/3.6.0

# go to project folder
cd ~/projects/genomic_prediction/simulation

# create folder to save results
mkdir -p ${FOLDER}
# get base for output name
outbase=$(echo ${HMP} | cut -d "." -f 1 | tr "/" "\n" | tail -n 1)
# get summary of missing data per marker
run_pipeline.pl -Xmx80g -importGuess ${HMP} \
                -GenotypeSummaryPlugin -endPlugin \
                -export ${FOLDER}/${outbase}\_OverallSummary,${FOLDER}/${outbase}\_AlleleSummary,${FOLDER}/${outbase}\_SiteSummary,${FOLDER}/${outbase}\_TaxaSummary

# qc hapmap file
Rscript scripts/qc_snp-sv_hmp.R ${FOLDER}/${outbase}\_SiteSummary.txt ${SVS} ${OUT} ${MISS}
