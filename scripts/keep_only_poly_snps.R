#### arguments for command line ----

args <- commandArgs(trailingOnly = TRUE)

# help
if (all(length(args) == 1 & args == "-h" | args == "--help")) {
  cat("
      Description: this script selects only polymorphic SNPs for each population

      Usage: Rscript keep_only_poly_snps.R [...]")
  quit()
}

# make sure the correct number of arguments are used
if (length(args) != 2) {
  stop("incorrect number of arguments provided.

       Usage: Rscript keep_only_poly_snps.R [...]
       ")
}

# assign arguments to variables
cross <- args[1]
reseq.snps.file <- args[2]

# cross <- "B73xLH82"
# reseq.snps.file <- "data/reseq_snps/biomAP_v_B73_SNPMatrix_parents.B73xLH82.not-in-SVs.hmp.txt"



#### libraries ----

if(!require("data.table")) install.packages("data.table")



#### keeping only polymorphic ----

# load data
reseq.snps <- fread(reseq.snps.file, header = TRUE, data.table = FALSE)

# get parents names
cross <- gsub(pattern = "*", replacement = "x", x = cross, fixed = TRUE)
cross <- unlist(strsplit(cross, split = "x"))
parent1 <- cross[1]
parent2 <- cross[2]

# get parents column numbers in resequencing data
p1.col <- grep(parent1, colnames(reseq.snps))
p2.col <- grep(parent2, colnames(reseq.snps))

# get type of each marker
marker.type <- apply(X = reseq.snps[, c(p1.col, p2.col)],
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
reseq.snps.poly <- reseq.snps[which(marker.type == "poly"), ]

# write results
reseq.snps.out <- gsub("not-in-SVs.hmp.txt", "poly.not-in-SVs.hmp.txt", reseq.snps.file)
fwrite(reseq.snps.poly, reseq.snps.out, quote = FALSE, sep = "\t",
       na = NA, row.names = FALSE)
