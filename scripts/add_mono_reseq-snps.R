#### arguments for command line ----

args <- commandArgs(trailingOnly = TRUE)

# help
if (all(length(args) == 1 & args == "-h" | args == "--help")) {
  cat("
Description: this script adds monomorphic resequencing SNPs back to projected SNP data.

Usage: Rscript add_mono-reseq-SNPs_after_projection.R [...]")
  quit()
}

# make sure the correct number of arguments are used
# you should provide 3 arguments
if (length(args) != 3) {
  stop("incorrect number of arguments provided.

Usage: Rscript add_mono-reseq-SNPs_after_projection.R [...]
       ")
}

# assign arguments to variables
cross <- args[1]
hmp.proj.filename <- args[2]
hmp.parents.filename <- args[3]


# cross <- "B73xLH82"
# hmp.proj.filename <- "analysis/projection_reseq-snps/biomAP_rils_SNPs-reseq_SNP-chip.B73xLH82.poly.projected.sliding-window.hmp.txt"
# hmp.parents.filename <- "data/reseq_snps/biomAP_v_B73_SNPMatrix_7parents.not-in-PAVs.hmp.txt"



#### libraries ----

if(!require("data.table")) install.packages("data.table")
if(!require("doParallel")) install.packages("doParallel")

if (detectCores() < 10) {
  num.cores <- detectCores()
} else {
  num.cores <- 10
}



#### add monomorphic snps by chromosome ----

# get parents names
parent1 <- unlist(strsplit(cross, "x"))[1]
parent2 <- unlist(strsplit(cross, "x"))[2]

# open file with reseq SNPs positions
hmp.parents <- fread(hmp.parents.filename, header = TRUE, data.table = FALSE)
hmp.parents <- cbind(hmp.parents[, 1:11], hmp.parents[, c(parent1, parent2)])

# after projection
hmp.proj <- fread(hmp.proj.filename, header = TRUE, data.table = FALSE)

registerDoParallel(cores = num.cores)
merged.hmp <- foreach(chr=1:10, .combine = rbind) %dopar% {

  # get non-missing monomorphic markers from chromosome
  hmp.parents.chr <- subset(hmp.parents, chrom == chr)
  hmp.parents.chr.mono <- hmp.parents.chr[which(hmp.parents.chr[, parent1] == hmp.parents.chr[, parent2]), ]
  hmp.parents.chr.mono <- hmp.parents.chr.mono[which(hmp.parents.chr.mono[, parent1] != "NN"), ]
  
  # remove SNP chip data from projected file
  hmp.proj.chr <- subset(hmp.proj, chrom == chr)
  hmp.proj.chr <- hmp.proj.chr[grep("^snp", hmp.proj.chr[, 1], perl = TRUE), ]
  
  # replicate columns of monomorphic parental SNPs to the have the number of RILs after projection
  reps <- length(12:NCOL(hmp.proj.chr)) - 2
  hmp.proj.chr.mono <- cbind(hmp.parents.chr.mono, replicate(reps, hmp.parents.chr.mono[, parent1]))
  # make sure columns have the same name as rils
  colnames(hmp.proj.chr.mono) <- colnames(hmp.proj.chr)
  
  # merge mono and polymorphic SNPs
  merged.hmp.chr <- rbind(hmp.proj.chr.mono, hmp.proj.chr)
  # sort by chromosome and position
  merged.hmp.chr <- merged.hmp.chr[order(merged.hmp.chr$pos), ]
  
  # return merged hmp for chromosome
  merged.hmp.chr
  
}
stopImplicitCluster()

# write merged files
outfile <- gsub("poly.projected.sliding-window.hmp.txt", "projected.sliding-window.hmp.txt", hmp.proj.filename)
fwrite(merged.hmp, outfile, quote = FALSE, sep = "\t", na = NA, row.names = FALSE)



