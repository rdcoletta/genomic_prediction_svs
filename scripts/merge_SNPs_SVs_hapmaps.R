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
if (length(args) != 7) {
  stop("incorrect number of arguments provided.

       Usage: Rscript merge_SNPs_SVs_hapmaps.R [...]
       ")
}

# assign arguments to variables
file.snps.parents <- args[1]
file.snps.rils <- args[2]
file.svs.parents <- args[3]
outfile.parents <- args[4]
outfile.rils <- args[5]
cross.info <- args[6]
dir.hmp.crosses <- args[7]

# file.snps.parents <- "data/usda_22kSNPs_7parents.sorted.diploid.v4.PHG35-reconstructed.hmp.txt"
# file.snps.rils <- "data/usda_22kSNPs_325rils.sorted.diploid.v4.hmp.txt"
# file.svs.parents <- "data/usda_SVs_7parents.sorted.hmp.txt"
# outfile.parents <- "data/usda_SNPs-SVs_7parents.not-in-PAVs.hmp.txt"
# outfile.rils <- "data/usda_SNPs-SVs_325rils.not-in-SVs.hmp.txt"
# cross.info <- "data/usda_biparental-crosses.txt"
# dir.hmp.crosses <- "data/merged_hapmaps_by_cross"


#### libraries used ----

if(!require("data.table")) install.packages("data.table")




#### function ----

RemoveSNPsInPAVs <- function(snp_file, parental_sv_file) {

  # load data
  snp.hmp <- fread(snp_file, header = TRUE, data.table = FALSE)
  sv.hmp <- fread(parental_sv_file, header = TRUE, data.table = FALSE)

  # create list to store SNPs within SVs
  SNPs.in.SVs <- list()

  # parse by chromosome
  for (chr in 1:length(unique(snp.hmp[, "chrom"]))) {

    cat("Analyzing chromosome ", chr, "...\n", sep = "")

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
        # print(SV.start)

        # check if a SNP is in a SV boundary
        SNP.within.del <- snp.hmp.subset[which(snp.hmp.subset[, 2] %in% seq(SV.start, SV.end, by = 1)), 1]

        if (length(SNP.within.del) > 0 & SV.size <= 100000) {
          # append SNPs to SV list of chromosome (only if deletion is smaller than 100kb)
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


MergeHapmaps <- function(snp_file, parental_sv_file, merge_RILs = FALSE) {

  # load SNP and SV datasets
  snp.hmp <- fread(snp_file, header = TRUE, data.table = FALSE)
  sv.hmp <- fread(parental_sv_file, header = TRUE, data.table = FALSE)

  # if want to merge SNPs and SVs from parents...
  if (merge_RILs == FALSE) {

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
    # the RIL hapmap. The values for these markers will be NA since we don't know yet if there is
    # a SV there

    # create empty data frame with nrow = number of SVs, and ncol = number of RILs to be imputed
    rils.genos <- data.frame(matrix("NN", nrow = NROW(sv.hmp), ncol = length(12:NCOL(snp.hmp))))
    colnames(rils.genos) <- colnames(snp.hmp)[12:NCOL(snp.hmp)]

    # make sure values for alleles and genotypes are NAs before merging with RIL hmp
    rils.sv.hmp <- cbind(sv.hmp[, 1:11], rils.genos)
    rils.sv.hmp$alleles <- NA
    # also make sure the hapmap columns are the same
    colnames(rils.sv.hmp)[1:11] <- colnames(snp.hmp)[1:11]


    # merge hapmaps
    merged.hmp <- rbind(snp.hmp, rils.sv.hmp)
  }

  # sort by chromosome and position
  merged.hmp <- merged.hmp[order(merged.hmp$chrom, merged.hmp$pos),]

  return(merged.hmp)

}



#### filter SNPs within SVs boundaries ----

# from parents
snp.hmp.parents <- RemoveSNPsInPAVs(snp_file = file.snps.parents, parental_sv_file = file.svs.parents)

outfile.parents.not.in.dels <- unlist(strsplit(file.snps.parents, split = ".", fixed = TRUE))[1]
outfile.parents.not.in.dels <- paste0(outfile.parents.not.in.dels, ".not-in-SVs.hmp.txt")
fwrite(snp.hmp.parents, outfile.parents.not.in.dels, sep = "\t", quote = FALSE, na = NA)


# from RILs
snp.hmp.rils <- RemoveSNPsInPAVs(snp_file = file.snps.rils, parental_sv_file = file.svs.parents)

outfile.rils.not.in.dels <- unlist(strsplit(file.snps.rils, split = ".", fixed = TRUE))[1]
outfile.rils.not.in.dels <- paste0(outfile.rils.not.in.dels, ".not-in-SVs.hmp.txt")
fwrite(snp.hmp.rils, outfile.rils.not.in.dels, sep = "\t", quote = FALSE, na = NA)




#### merge hapmaps from parents ----

parents.merged <- MergeHapmaps(snp_file = outfile.parents.not.in.dels,
                               parental_sv_file = file.svs.parents,
                               merge_RILs = FALSE)

fwrite(parents.merged, outfile.parents,  sep = "\t", quote = FALSE, na = NA)



#### merge hapmaps from RILs ----

rils.merged <- MergeHapmaps(snp_file = outfile.rils.not.in.dels,
                            parental_sv_file = file.svs.parents,
                            merge_RILs = TRUE)

fwrite(rils.merged, outfile.rils, sep = "\t", quote = FALSE, na = NA)




#### subset by cross ----

# load table with cross information
df.crosses <- fread(cross.info, header = TRUE, data.table = FALSE)

# make dir to store results
if (!dir.exists(dir.hmp.crosses)) dir.create(dir.hmp.crosses)

# write merged RILs hapmap by cross
for (row in 1:NROW(df.crosses)) {

  # get cross name, and change * by x in the name
  cross <- df.crosses[row, "cross"]
  cross <- gsub(pattern = "*", replacement = "x", x = cross, fixed = TRUE)
  # get RIL names for that cross
  RILs <- unlist(strsplit(df.crosses[row, "RILs"], split = ","))

  # subset merged hapmaps from parents
  parents <- unlist(strsplit(cross, split = "x"))
  subset.parents.merged <- cbind(parents.merged[, 1:11],
                                 parents.merged[, which(colnames(parents.merged) %in% parents)])
  # since i'm subsetting, the alleles column might be wrong...
  subset.parents.merged$alleles <- NA

  # subset merged hapmaps from RILs
  subset.rils.merged <- cbind(rils.merged[, 1:11],
                              rils.merged[, which(colnames(rils.merged) %in% RILs)])
  # since i'm subsetting, the alleles column might be wrong...
  subset.rils.merged$alleles <- NA

  # write results only if RILs were actually planted (i.e., if RILs are in the merged file)
  if (NCOL(subset.rils.merged) > 11) {
    # parents
    out_parents <- paste0(dir.hmp.crosses, "/usda_SNPs-SVs_", cross, "_parents.sorted.hmp.txt")
    fwrite(subset.parents.merged, file = out_parents, sep = "\t", na = "NN", quote = FALSE)

    # RILs
    out_rils <- paste0(dir.hmp.crosses, "/usda_SNPs-SVs_", cross, "_RILs.sorted.hmp.txt")
    fwrite(subset.rils.merged, file = out_rils, sep = "\t", na = "NN", quote = FALSE)
  }
}
