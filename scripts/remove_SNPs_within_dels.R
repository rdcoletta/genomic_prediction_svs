#### arguments for command line ----

args <- commandArgs(trailingOnly = TRUE)

# help
if (all(length(args) == 1 & args == "-h" | args == "--help")) {
  cat("
      Description: this script removes SNPs within deletions up to 100kb

      Usage: Rscript remove_SNPs_within_dels.R [...]")
  quit()
}

# make sure the correct number of arguments are used
if (length(args) != 2) {
  stop("incorrect number of arguments provided.

       Usage: Rscript remove_SNPs_within_dels.R [...]
       ")
}

# assign arguments to variables
file.snps <- args[1]
file.svs <- args[2]

# file.snps <- "data/hapmap_by_cross/usda_22kSNPs_parents.sorted.diploid.PHG35-reconstructed.filtered.v4.B73xLH82.hmp.txt"
# file.snps <- "data/hapmap_by_cross/usda_22kSNPs_rils.sorted.diploid.filtered.v4.B73xLH82.hmp.txt"
# file.svs <- "data/hapmap_by_cross/usda_SVs_parents.sorted.B73xLH82.hmp.txt"



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
    sv.hmp.subset <- sv.hmp[which(sv.hmp[, "chrom"] == chr), c(1, 12, 13)]
    # get each SNP name and position in a chromosome
    snp.hmp.subset <- snp.hmp[which(snp.hmp[, "chrom"] == chr), c(1, 4)]

    for (row in 1:NROW(sv.hmp.subset)) {

      # get sv name and parental calls
      SV <- sv.hmp.subset[row, 1]
      p1.call <- sv.hmp.subset[row, 2]
      p2.call <- sv.hmp.subset[row, 3]

      # check SNPs that fall in deletions only (not dups)
      SV.type <- unlist(strsplit(SV, split = ".", fixed = TRUE))[1]

      if (SV.type == "del") {

        if (p1.call == "TT" | p2.call == "TT") {

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
    }

    # remove redundancy
    SNPs.in.SVs[[paste0("chr", chr)]] <- unique(SNPs.in.SVs[[paste0("chr", chr)]])

    cat(length(SNPs.in.SVs[[paste0("chr", chr)]]), "out of", NROW(snp.hmp.subset), "SNPs are within a PAV\n")

  }

  # filter hapmap so it doesn't have the SNPs in the SNPs.in.SVs list
  snp.hmp.filtered <- snp.hmp[which(!snp.hmp[, 1] %in% as.character(unlist(SNPs.in.SVs))), ]

  return(snp.hmp.filtered)

}



#### filter SNPs within SVs boundaries ----

snp.hmp <- RemoveSNPsInPAVs(snp_file = file.snps, parental_sv_file = file.svs)

outfile <- gsub(".hmp.txt", ".not-in-SVs.hmp.txt", file.snps)
fwrite(snp.hmp, outfile, sep = "\t", quote = FALSE, na = NA)
