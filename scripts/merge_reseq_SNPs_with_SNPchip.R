#### arguments for command line ----

args <- commandArgs(trailingOnly = TRUE)

# help
if (all(length(args) == 1 & args == "-h" | args == "--help")) {
  cat("
      Description: this script removes resequencing SNPs from deletions, subset SNPs by population,
                   get only polymorphic SNPs for each population, and finally merge these SNPs with
                   SNPs from chip data

      Usage: Rscript merge_reseq_SNPs_with_SNPchip.R [...]")
  quit()
}

# make sure the correct number of arguments are used
if (length(args) != 7) {
  stop("incorrect number of arguments provided.

       Usage: Rscript merge_reseq_SNPs_with_SNPchip.R [...]
       ")
}

# assign arguments to variables
file.snps.parents <- args[1]
file.svs.parents <- args[2]
cross.info <- args[3]
sv.proj.folder <- args[4]
sv.donors.folder <- args[5]
outfile.parents <- args[6]
outfile.rils <- args[7]

# file.snps.parents <- "data/reseq_snps/biomAP_v_B73_SNPMatrix_7parents.hmp.txt"
# file.svs.parents <- "data/usda_SVs_7parents.sorted.hmp.txt"
# cross.info <- "data/usda_biparental-crosses.txt"
# sv.proj.folder <- "analysis/projection"
# sv.donors.folder <- "data/merged_hapmaps_by_cross"
# outfile.parents <- "data/reseq_snps/biomAP_parents_SNPs-reseq_SNP-chip_SVs-proj.hmp.txt"
# outfile.rils <- "data/reseq_snps/biomAP_rils_SNPs-reseq_SNP-chip_SVs-proj.hmp.txt"



#### libraries used ----

if(!require("data.table")) install.packages("data.table")
if(!require("doParallel")) install.packages("doParallel")



#### function ----

RemoveSNPsInPAVs <- function(snp_file, parental_sv_file) {

  # load data
  snp.hmp <- fread(snp_file, header = TRUE, data.table = FALSE)
  sv.hmp <- fread(parental_sv_file, header = TRUE, data.table = FALSE)

  # create list to store SNPs within SVs
  SNPs.in.SVs <- list()

  # parse by chromosome
  for (chr in 1:length(unique(snp.hmp[, "chrom"]))) {

    cat("  Analyzing chromosome ", chr, ":  ", sep = "")

    # get each SV ID in a chromosome
    sv.hmp.subset <- sv.hmp[which(sv.hmp[, "chrom"] == chr), 1]
    # get each SNP name and position in a chromosome
    snp.hmp.subset <- snp.hmp[which(snp.hmp[, "chrom"] == chr), c(1, 4)]

    for (SV in sv.hmp.subset) {

      # check SNPs that fall in deletions only (not dups)
      SV.type <- unlist(strsplit(SV, split = ".", fixed = TRUE))[1]

      if (SV.type == "del") {

        # get start and end positions of SV
        SV.start <- as.numeric(unlist(strsplit(SV, split = ".", fixed = TRUE))[3])
        SV.end <- as.numeric(unlist(strsplit(SV, split = ".", fixed = TRUE))[4])
        SV.size <- SV.end - SV.start

        # check if a SNP is in a SV boundary
        SNP.within.del <- snp.hmp.subset[which(snp.hmp.subset[, 2] %in% seq(SV.start, SV.end, by = 1)), 1]

        if (length(SNP.within.del) > 0 & SV.size <= 100000) {
          # append SNPs to SV list of chromosome
          SNPs.in.SVs[[paste0("chr", chr)]] <- append(SNPs.in.SVs[[paste0("chr", chr)]], SNP.within.del)
        }
      }

    }

    # remove redundancy
    SNPs.in.SVs[[paste0("chr", chr)]] <- unique(SNPs.in.SVs[[paste0("chr", chr)]])

    cat(length(SNPs.in.SVs[[paste0("chr", chr)]]), "out of", NROW(snp.hmp.subset), "SNPs are within a PAV\n")

  }

  # filter hapmap so it doesn't have the SNPs in the SNPs.in.SVs list
  snp.hmp.filtered <- snp.hmp[which(!snp.hmp[, 1] %in% as.character(unlist(SNPs.in.SVs))), ]

  return(snp.hmp.filtered)

}


MergeHapmaps <- function(snp_file, reseq_file, merge_RILs = FALSE) {

  # load SNP and SV datasets
  snp.hmp <- fread(snp_file, header = TRUE, data.table = FALSE)
  reseq.hmp <- fread(reseq_file, header = TRUE, data.table = FALSE)

  # if want to merge SNPs and SVs from parents...
  if (merge_RILs == FALSE) {

    # make sure columns have the same names and order in both files
    reseq.hmp <- cbind(reseq.hmp[, 1:11],
                    reseq.hmp[, colnames(snp.hmp)[12:NCOL(snp.hmp)]])
    colnames(reseq.hmp) <- colnames(snp.hmp)

    # merge hapmaps
    merged.hmp <- rbind(snp.hmp, reseq.hmp)
  }

  # if want to create dataset for projection of SVs from parents to RILs...
  if (merge_RILs == TRUE) {

    # since i want to impute SVs from parents to RILs, I need to also include the SV marker names in
    # the RIL hapmap. The values for these markers will be NA since we don't know yet if there is
    # a SV there

    # create empty data frame with nrow = number of SVs, and ncol = number of RILs to be imputed
    rils.genos <- data.frame(matrix("NN", nrow = NROW(reseq.hmp), ncol = length(12:NCOL(snp.hmp))))
    colnames(rils.genos) <- colnames(snp.hmp)[12:NCOL(snp.hmp)]

    # make sure values for alleles and genotypes are NAs before merging with RIL hmp
    rils.reseq.hmp <- cbind(reseq.hmp[, 1:11], rils.genos)
    rils.reseq.hmp$alleles <- NA
    # also make sure the hapmap columns are the same
    colnames(rils.reseq.hmp)[1:11] <- colnames(snp.hmp)[1:11]

    # merge hapmaps
    merged.hmp <- rbind(snp.hmp, rils.reseq.hmp)
  }

  # sort by chromosome and position
  merged.hmp <- merged.hmp[order(merged.hmp$chrom, merged.hmp$pos),]

  return(merged.hmp)

}


#### filter SNPs within SVs boundaries ----

cat("Removing SNPs within 100kb deletions\n")

# snp.hmp.parents <- RemoveSNPsInPAVs(snp_file = file.snps.parents, parental_sv_file = file.svs.parents)

outfile.parents.not.in.dels <- unlist(strsplit(file.snps.parents, split = ".", fixed = TRUE))[1]
outfile.parents.not.in.dels <- paste0(outfile.parents.not.in.dels, ".not-in-PAVs.hmp.txt")
# fwrite(snp.hmp.parents, outfile.parents.not.in.dels, sep = "\t", quote = FALSE, na = NA)



#### subset data by cross ----

# For each population, merge these SNPs with best markers + projected SVs
# - Parents: best markers + projected SVs + SNPs in reseq
# - RILs: best markers + projected SVs + SNPs in reseq (but all converted to NN)

snp.hmp.parents <- fread("data/reseq_snps/biomAP_v_B73_SNPMatrix_7parents.not-in-PAVs.hmp.txt",
                          header = TRUE, data.table = FALSE)


# load table with cross information
df.crosses <- fread(cross.info, header = TRUE, data.table = FALSE)

rows.to.keep <- c()
# filter cross information to have only projected crosses
for (row in 1:NROW(df.crosses)) {

  # get cross name, and change * by x in the name
  cross <- df.crosses[row, "cross"]
  cross <- gsub(pattern = "*", replacement = "x", x = cross, fixed = TRUE)

  filename.after.proj <- list.files(path = sv.proj.folder,
                                    pattern = paste0(cross, "_RILs.projected.hmp.txt"),
                                    full.names = TRUE)

  if (length(filename.after.proj) > 0) {
    rows.to.keep <- append(rows.to.keep, row)
  }

}
df.crosses <- df.crosses[rows.to.keep, ]

# write merged RILs hapmap by cross
for (row in 1:NROW(df.crosses)) {

  # get cross name, and change * by x in the name
  cross <- df.crosses[row, "cross"]
  cross <- gsub(pattern = "*", replacement = "x", x = cross, fixed = TRUE)
  # get RIL names for that cross
  RILs <- unlist(strsplit(df.crosses[row, "RILs"], split = ","))

  # subset merged hapmaps from parents
  parents <- unlist(strsplit(cross, split = "x"))
  subset.parents <- cbind(snp.hmp.parents[, 1:11],
                          snp.hmp.parents[, which(colnames(snp.hmp.parents) %in% parents)])
  # since i'm subsetting, the alleles column might be wrong...
  subset.parents$alleles <- NA

  cat("Analyzing cross ", cross, ":\n", sep = "")

  # get parents column numbers in resequencing data
  p1.col <- grep(parents[1], colnames(subset.parents))
  p2.col <- grep(parents[2], colnames(subset.parents))

  # remove SNPs that are missing in at least one parent
  subset.parents <- subset.parents[which(subset.parents[, p1.col] != "NN" & subset.parents[, p2.col] != "NN"), ]
  
  # create folder to store results for that cross
  cross.folder <- unlist(strsplit(file.snps.parents, split = "/"))
  cross.folder <- paste0(cross.folder[-length(cross.folder)], collapse = "/")
  cross.folder <- paste0(cross.folder, "/", cross)
  
  if (!dir.exists(cross.folder)) dir.create(cross.folder)
  
  
  #### keep only polymorphic SNPs ----

  cat("  Removing monomorphic SNPs\n")

  registerDoParallel(cores = 10)
  subset.parents.poly <- foreach(chr=1:10, .combine = rbind) %dopar% {

    # separate into chromosomes
    subset.parents.chr <- subset(subset.parents, chrom == chr)
    # get type of each marker
    marker.type.chr <- apply(X = subset.parents.chr[, c(p1.col, p2.col)],
                             MARGIN = 1, FUN = function(snp) {
                               # get unique genotypes between parents
                               genotypes <- unique(snp)
                               genotypes <- genotypes[genotypes != "NN"]

                               if (length(genotypes) == 0) {

                                 # if there is no genotype, snp is missing
                                 return("missing")

                               } else if (length(genotypes) == 1) {

                                 # if there is one genotype, it's monomorphic
                                 # but distinguish if SNP is het
                                 alleles <- unlist(strsplit(genotypes, split = ""))
                                 if (alleles[1] == alleles[2]) {
                                   return("mono")
                                 } else {
                                   return("het")
                                 }

                               } else {

                                 # if there are two genotypes, it's polymorphic
                                 # but distiguish if one of the genotypes is het
                                 p1.alleles <- unlist(strsplit(genotypes[1], split = ""))
                                 p2.alleles <- unlist(strsplit(genotypes[2], split = ""))
                                 if (p1.alleles[1] == p1.alleles[2] & p2.alleles[1] == p2.alleles[2]) {
                                   return("poly")
                                 } else {
                                   return("het")
                                 }

                               }
                             })
    # keep only homozygous polymorphic markers between parents
    subset.parents.poly.chr <- subset.parents.chr[which(marker.type.chr == "poly"), ]

    # return filter chomosome set
    subset.parents.poly.chr

  }
  stopImplicitCluster()

  # write results
  subset.parents.out <- gsub(pattern = "not-in-PAVs.hmp.txt", replacement = paste0(cross, ".poly.not-in-PAVs.hmp.txt"),
                             x = outfile.parents.not.in.dels)
  subset.parents.out <- rev(unlist(strsplit(subset.parents.out, split = "/")))[1]
  subset.parents.out <- paste0(cross.folder, "/", subset.parents.out)
  fwrite(subset.parents.poly, subset.parents.out, quote = FALSE, sep = "\t", na = NA, row.names = FALSE)


  # open parental files with 22kSNPs + SVs
  SNPs_SVs_parents_file <- paste0(sv.donors.folder, "/usda_SNPs-SVs_", cross, "_parents.sorted.hmp.txt")
  # SNPs_SVs_parents <- fread(SNPs_SVs_parents_file, header = TRUE, data.table = FALSE)
  # open RIL files with 22kSNPs + projected SVs
  SNPs_SVs_rils_file <- paste0(sv.proj.folder, "/usda_SNPs-SVs_", cross, "_RILs.projected.hmp.txt")
  # SNPs_SVs_rils <- fread(SNPs_SVs_rils_file, header = TRUE, data.table = FALSE)


  #### merge hapmaps from parents ----

  cat("  Merging resequencing SNPs of parents with SNP chip data and SVs\n")
  parents.merged <- MergeHapmaps(snp_file = SNPs_SVs_parents_file,
                                 reseq_file = subset.parents.out,
                                 merge_RILs = FALSE)


  #### merge hapmaps from RILs ----

  cat("  Merging resequencing SNPs of RILs with SNP chip data and projected SVs\n")
  rils.merged <- MergeHapmaps(snp_file = SNPs_SVs_rils_file,
                              reseq_file = subset.parents.out,
                              merge_RILs = TRUE)



  #### remove any duplicated SNP between reseq and SNPchip ----

  # find duplicated snps from reseq and chip
  dup.rows.to.remove <- duplicated(parents.merged[, c("chrom", "pos")])

  # since i merged files with rbind(snp.hmp, rils.reseq.hmp), the first duplicate will be
  # the SNP from chip, and the second duplicate will be from the resequencing data
  # also, I will remove the duplicates from the resequencing data to maintain compatibility
  # with previous SV projections
  parents.merged.filtered <- parents.merged[!dup.rows.to.remove, ]
  rils.merged.filtered <- rils.merged[!dup.rows.to.remove, ]


  cat("  Making sure parental and RIL data have the same markers\n")
  if (NROW(parents.merged.filtered) == NROW(rils.merged.filtered)) {
    all(parents.merged.filtered[, 1] == rils.merged.filtered[, 1])
  } else {
    stop("Parents and RILs don't have the same markers\n")
  }


  #### write results ----

  cat("  Writing files\n")
  
  parents.merged.outfile <- unlist(strsplit(outfile.parents, split = ".", fixed = TRUE))[1]
  parents.merged.outfile <- paste0(parents.merged.outfile, ".", cross, ".poly.hmp.txt")
  parents.merged.outfile <- rev(unlist(strsplit(parents.merged.outfile, split = "/")))[1]
  parents.merged.outfile <- paste0(cross.folder, "/", parents.merged.outfile)
  
  fwrite(parents.merged.filtered, parents.merged.outfile, quote = FALSE, sep = "\t", na = NA, row.names = FALSE)

  rils.merged.outfile <- unlist(strsplit(outfile.rils, split = ".", fixed = TRUE))[1]
  rils.merged.outfile <- paste0(rils.merged.outfile, ".", cross, ".poly.hmp.txt")
  rils.merged.outfile <- rev(unlist(strsplit(rils.merged.outfile, split = "/")))[1]
  rils.merged.outfile <- paste0(cross.folder, "/", rils.merged.outfile)
  
  fwrite(rils.merged.filtered, rils.merged.outfile, quote = FALSE, sep = "\t", na = NA, row.names = FALSE)

  cat("  Done!\n")
}
