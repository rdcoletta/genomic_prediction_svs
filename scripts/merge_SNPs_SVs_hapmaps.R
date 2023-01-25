#### arguments for command line ----

args <- commandArgs(trailingOnly = TRUE)

# help
if (all(length(args) == 1 & args == "-h" | args == "--help")) {
  cat("
      Description: this script merges SNPs and SVs hapmap files from usda parents to be used in
                   Tassel 5 when projecting SVs into RILs

      Usage: Rscript merge_SNPs_SVs_hapmaps.R [...]")
  quit()
}

# make sure the correct number of arguments are used
if (length(args) != 5) {
  stop("incorrect number of arguments provided.

       Usage: Rscript merge_SNPs_SVs_hapmaps.R [...]
       ")
}

# assign arguments to variables
sv.file <- args[1]
chip.parents.file <- args[2]
chip.rils.file <- args[3]
output.parents <- args[4]
output.rils <- args[5]

# setwd("~/projects/genomic_prediction/simulation")
# sv.file <- "data/hapmap_by_cross/usda_SVs_parents.sorted.B73xLH82.hmp.txt"
# chip.parents.file <- "data/hapmap_by_cross/usda_22kSNPs_parents.sorted.diploid.PHG35-reconstructed.filtered.v4.B73xLH82.not-in-SVs.hmp.txt"
# chip.rils.file <- "data/hapmap_by_cross/usda_22kSNPs_rils.sorted.diploid.filtered.v4.B73xLH82.not-in-SVs.sliding-window.hmp.txt"
# output.parents <- "data/hapmap_by_cross/usda_parents_SV-SNPchip.B73xLH82.hmp.txt"
# output.rils <- "data/hapmap_by_cross/usda_rils_SV-SNPchip.B73xLH82.not-projected.hmp.txt"

# setwd("~/projects/genomic_prediction/simulation")
# sv.file <- "data/reseq_snps/biomAP_v_B73_SNPMatrix_parents.B73xLH82.poly.not-in-SVs.hmp.txt"
# chip.parents.file <- "data/hapmap_by_cross/usda_parents_SV-SNPchip.B73xLH82.hmp.txt"
# chip.rils.file <- "data/hapmap_by_cross/usda_rils_SV-SNPchip.B73xLH82.not-projected.hmp.txt"
# output.parents <- "data/hapmap_by_cross/usda_parents_SV-SNPchip-SNPreseq.B73xLH82.hmp.txt"
# output.rils <- "data/hapmap_by_cross/usda_rils_SV-SNPchip-SNPreseq.B73xLH82.not-projected.hmp.txt"


#### libraries used ----

if(!require("data.table")) install.packages("data.table")




#### function ----

MergeHapmaps <- function(snp_file, parental_sv_file, merge_RILs = FALSE) {

  # load SNP and SV datasets
  snp.hmp <- fread(snp_file, header = TRUE, data.table = FALSE)
  sv.hmp <- fread(parental_sv_file, header = TRUE, data.table = FALSE)

  # filter out SVs that are missing in both parents
  sv.hmp <- sv.hmp[which(sv.hmp[, 12] != "NN" & sv.hmp[, 13] != "NN"), ]

  # if want to merge SNPs and SVs from parents...
  if (merge_RILs == FALSE) {

    # remove anything after "_" in SNP colnames
    for (i in 12:NCOL(snp.hmp)) {
      colnames(snp.hmp)[i] <- unlist(strsplit(colnames(snp.hmp)[i], split = "_"))[1]
    }

    # make sure columns have the same names and order in both files
    sv.hmp <- cbind(sv.hmp[, 1:11],
                    sv.hmp[, colnames(snp.hmp)[12:NCOL(snp.hmp)]])
    colnames(sv.hmp) <- colnames(snp.hmp)

    # merge hapmaps
    merged.hmp <- rbind(snp.hmp, sv.hmp)
  }

  # if want to create dataset for projection of SVs from parents to RILs...
  if (merge_RILs == TRUE) {

    # since i want to impute SVs from parents to RILs, I need to also include the SV marker names in
    # the RIL hapmap, but the actual genotypes of these markers are unknown (so far) and need to be
    # set to NA

    # create empty data frame with nrow = number of SVs, and ncol = number of RILs to be imputed
    rils.genos <- data.frame(matrix("NN", nrow = NROW(sv.hmp), ncol = length(12:NCOL(snp.hmp))))
    colnames(rils.genos) <- colnames(snp.hmp)[12:NCOL(snp.hmp)]

    # make sure values for alleles and genotypes are NAs before merging with RIL hmp
    rils.sv.hmp <- cbind(sv.hmp[, 1:11], rils.genos)
    rils.sv.hmp$alleles <- "NN"
    # also make sure the hapmap columns are the same
    colnames(rils.sv.hmp)[1:11] <- colnames(snp.hmp)[1:11]

    # merge hapmaps
    merged.hmp <- rbind(snp.hmp, rils.sv.hmp)
  }

  # sort by chromosome and position
  merged.hmp <- merged.hmp[order(merged.hmp$chrom, merged.hmp$pos),]

  return(merged.hmp)

}



#### merge hapmaps from RILs ----

cat("Merging SNPs and SVs in parents...\n")
parents.merged <- MergeHapmaps(snp_file = chip.parents.file, parental_sv_file = sv.file,
                               merge_RILs = FALSE)

fwrite(parents.merged, output.parents, quote = FALSE, sep = "\t", na = NA, row.names = FALSE)
cat("Done!\n")

#### merge hapmaps from RILs ----

cat("Merging SNPs and SVs in RILs\n")
rils.merged <- MergeHapmaps(snp_file = chip.rils.file, parental_sv_file = sv.file,
                            merge_RILs = TRUE)

fwrite(rils.merged, output.rils, quote = FALSE, sep = "\t", na = NA, row.names = FALSE)
cat("Done!\n")



#### make sure parents and rils have the same markers ----

if (NROW(parents.merged) > NROW(rils.merged)) {
  cat("Making sure parents and rils have the same SNPs\n")
  parents.merged.filtered <- parents.merged[which(parents.merged[, 1] %in% rils.merged[, 1]), ]
  fwrite(parents.merged.filtered, output.parents, quote = FALSE, sep = "\t", na = NA, row.names = FALSE)
  cat("Done!\n\n")
}
