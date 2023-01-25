#!/bin/bash
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=50gb
#SBATCH -J hmp2plk_random_${CHR}
#SBATCH -o /home/hirschc1/della028/projects/genomic_prediction/simulation/analysis/ld/%x_%j.out
#SBATCH -e /home/hirschc1/della028/projects/genomic_prediction/simulation/analysis/ld/%x_%j.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=della028@umn.edu
#SBATCH --no-requeue

module load R/3.6.0

# define scratch folder
SCRATCH=/scratch.global/della028/hirsch_lab/genomic_prediction/ld

# go to project folder
cd ~/projects/genomic_prediction/simulation

# summarize selected marker data
run_pipeline.pl -Xmx50g -importGuess data/usda_rils_projected-SVs-SNPs.chr${CHR}.poly.hmp.txt \
                -GenotypeSummaryPlugin -endPlugin \
                -export analysis/ld_downsample/datasets/markers_chr${CHR}_OverallSummary,analysis/ld_downsample/datasets/markers_chr${CHR}_AlleleSummary,analysis/ld_downsample/datasets/markers_chr${CHR}_SiteSummary,analysis/ld_downsample/datasets/markers_chr${CHR}_TaxaSummary

# downsample number of SNPs while matching MAF distribution of SVs
Rscript scripts/downsample_markers_by_MAF.R \
        analysis/ld_downsample/datasets/markers_chr${CHR}_SiteSummary.txt \
        data/SVs_IDs_poly.txt \
        analysis/ld_downsample/datasets/markers_chr${CHR}_MAF-downsampled.txt \
        --prop-downsample=0.25 --bin-size=0.01 --reps=${REPS} --seed=${SEED}

for rep in $(seq 1 ${REPS}); do
  # create output folder
  OUTFOLDER=analysis/ld_downsample/datasets/rep${rep}
  mkdir -p ${OUTFOLDER}
  # filter hmp to create dataset
  run_pipeline.pl -Xmx50g -importGuess data/usda_rils_projected-SVs-SNPs.chr${CHR}.poly.hmp.txt \
                  -includeSiteNamesInFile analysis/ld_downsample/datasets/markers_chr${CHR}_MAF-downsampled.rep${rep}.txt \
                  -export ${OUTFOLDER}/usda_rils_projected-SVs-SNPs.poly.chr${CHR}.rep${rep}.hmp.txt \
                  -exportType HapmapDiploid
  # export to plink for LD calculation
  run_pipeline.pl -Xmx50g -importGuess ${OUTFOLDER}/usda_rils_projected-SVs-SNPs.poly.chr${CHR}.rep${rep}.hmp.txt \
                  -export ${SCRATCH}/usda_rils_projected-SVs-SNPs.poly.chr${CHR}.rep${rep} \
                  -exportType Plink
done
