#!/bin/bash
#SBATCH --time=2:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=10
#SBATCH --mem=100gb
#SBATCH -J ld_filter_selected-qtns
#SBATCH -o /home/hirschc1/della028/projects/genomic_prediction/simulation/analysis/ld_downsample/%x_%j.out
#SBATCH -e /home/hirschc1/della028/projects/genomic_prediction/simulation/analysis/ld_downsample/%x_%j.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=della028@umn.edu
#SBATCH --no-requeue

# go to project folder
cd ~/projects/genomic_prediction/simulation

echo "job started @ $(date)"

H2=0.3
WINDOW=5000
FILTER=0.25
SCRATCH=/scratch.global/della028/hirsch_lab/genomic_prediction/ld

for QTN in 100; do
  for VAR in SNP SV both; do
    # get folder name
    FOLDER=analysis/ld_downsample/sim_traits/rep${REP}/${QTN}-QTNs_from_${VAR}/${H2}-heritability/pop${POP}
    echo "  ${FOLDER}"
    # select QTNs per chr
    for chr in {1..10}; do
      # get ld file
      LDFILE=${SCRATCH}/usda_rils_projected-SVs-SNPs.poly.chr${chr}.rep${REP}.window-${WINDOW}kb.filter-${FILTER}.ld.gz
      # get file with causative variants
      awk -v chr="$chr" '$9 == chr' ${FOLDER}/Additive_QTNs.txt | cut -f 7 > ${FOLDER}/list_QTNs.chr${chr}.causal-pop${POP}.txt
      # keep only QTNs
      zcat ${LDFILE} | head -n 1 > ${FOLDER}/ld_info.QTNs-rep${REP}-pop${POP}.chr${chr}.ld
      zcat ${LDFILE} | grep -F -w -f ${FOLDER}/list_QTNs.chr${chr}.causal-pop${POP}.txt >> ${FOLDER}/ld_info.QTNs-rep${REP}-pop${POP}.chr${chr}.ld &
    done
    wait
    # merge results from chrs
    cp ${FOLDER}/ld_info.QTNs-rep${REP}-pop${POP}.chr1.ld ${FOLDER}/ld_info.QTNs-rep${REP}-pop${POP}.ld
    for chr in {2..10}; do
      sed 1d ${FOLDER}/ld_info.QTNs-rep${REP}-pop${POP}.chr${chr}.ld >> ${FOLDER}/ld_info.QTNs-rep${REP}-pop${POP}.ld
    done
    # remove chr files
    for chr in {1..10}; do
      rm ${FOLDER}/ld_info.QTNs-rep${REP}-pop${POP}.chr${chr}.ld
    done
    # # copy results to higher h2 folder
    # cp ${FOLDER}/ld_info.QTNs-rep${REP}-pop${POP}.ld analysis/ld_downsample/sim_traits/rep${REP}/${QTN}-QTNs_from_${VAR}/0.7-heritability/pop${POP}/
  done
done

echo "job finished @ $(date)"
