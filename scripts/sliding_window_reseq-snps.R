#### arguments for command line ----

args <- commandArgs(trailingOnly = TRUE)

# help
if (all(length(args) == 1 & args == "-h" | args == "--help")) {
  cat("
Description: This script applies a sliding window approach to correct possible miscalls after projection.

Usage: Rscript sliding_window_reseq-snps.R [...]")
  quit()
}

# make sure the correct number of arguments are used
# you should provide 6 arguments
if (length(args) != 6) {
  stop("incorrect number of arguments provided.
       
       Usage: Rscript sliding_window_reseq-snps.R [...]
       ")
}
cross <- args[1]
data.filename <- args[2]
reseq.parents.filename <- args[3]

if (grepl("--window_size=", args[4])) {
  window.size <- as.numeric(unlist(strsplit(args[4], split = "="))[2])
} else {
  stop("Invalid argument 4")
}

if (grepl("--window_step=", args[5])) {
  window.step <- as.numeric(unlist(strsplit(args[5], split = "="))[2])
} else {
  stop("Invalid argument 5")
}

if (grepl("--min_snps_per_window=", args[6])) {
  min.snps.per.window <- as.numeric(unlist(strsplit(args[6], split = "="))[2])
} else {
  stop("Invalid argument 6")
}

# cross <- "B73xLH82"
# data.filename <- "analysis/projection_reseq-snps/biomAP_rils_SNPs-reseq_SNP-chip.B73xLH82.poly.projected.hmp.txt"
# reseq.parents.filename <- "data/reseq_snps/B73xLH82/biomAP_parents_SNPs-reseq_SNP-chip.B73xLH82.poly.hmp.txt"
# window.size <- 45
# window.step <- 1
# min.snps.per.window <- 15




#### libraries ----

if(!require("data.table")) install.packages("data.table")
if(!require("foreach")) install.packages("foreach")
if(!require("doParallel")) install.packages("doParallel")



#### sliding window approach ----

cat("\nLoading cross ", cross, "...\n", sep = "")


# get available cores for paralellizing
num.cores <- detectCores()



# load data
geno.data.infile <- fread(data.filename, header = TRUE, data.table = FALSE)
reseq.parents.infile <- fread(reseq.parents.filename, header = TRUE, data.table = FALSE)

# get parents names
parent1 <- unlist(strsplit(cross, "x"))[1]
parent2 <- unlist(strsplit(cross, "x"))[2]


# get parents column numbers in resequencing data
p1.col.reseq <- grep(parent1, colnames(reseq.parents.infile))
p2.col.reseq <- grep(parent2, colnames(reseq.parents.infile))


# make sure that gbs data has the snps as reseq
if (all(geno.data.infile[,1] != reseq.parents.infile[,1])) stop("Data have different length")

geno.data.outfile <- data.frame(stringsAsFactors = FALSE)
for (chr in 1:10) {
  
  cat("Analyze chromosome ", chr, "\n", sep = "")
  # subset by chr
  reseq.parents.infile.chr <- subset(reseq.parents.infile, chrom == chr)
  geno.data.infile.chr <- subset(geno.data.infile, chrom == chr)
  
  registerDoParallel(cores = num.cores)
  geno.data.infile.window <- foreach(ril.col=12:NCOL(geno.data.infile.chr), .combine = cbind) %dopar% {
    
    # set up first window
    window.start <- 1
    window.stop <- window.start + (window.size - 1)
    
    # create a vector to store consenus genotype for each window
    window.consensus <- c()
    
    # use slide window approach until end of the window reaches the last SNP
    while (window.stop <= NROW(geno.data.infile.chr)) {
      
      # get genotypes from parents and ril for that window
      window <- cbind(reseq.parents.infile.chr[window.start:window.stop, p1.col.reseq],
                      reseq.parents.infile.chr[window.start:window.stop, p2.col.reseq],
                      geno.data.infile.chr[window.start:window.stop, ril.col])
      
      # define from which parents the SNP in ril comes from
      window.calls <- apply(window, MARGIN = 1, FUN = function(genotypes) {
        if (genotypes[3] == "NN") {
          # if ril snp is NN
          return("missing")
        } else if (genotypes[1] == genotypes[3]) {
          # if ril snp is the same as parent1
          return("p1")
        } else if (genotypes[2] == genotypes[3]) {
          # if ril snp is the same as parent2
          return("p2")
        } else {
          # check if ril snp is a het or if it has a different allele from parents
          p1.alleles <- unlist(strsplit(genotypes[1], split = ""))
          p2.alleles <- unlist(strsplit(genotypes[2], split = ""))
          ril.alleles <- unlist(strsplit(genotypes[3], split = ""))
          if (all(unique(p1.alleles) %in% ril.alleles) & all(unique(p2.alleles) %in% ril.alleles) & ril.alleles[1] != ril.alleles[2]) {
            return("het")
          } else {
            return("missing")
          }
        }
      })
      
      # check if there is enough ril snps genotyped
      if (sum(window.calls != "missing") >= min.snps.per.window) {
        
        # get number of alleles for each parent (multiply by 2 because a homozygous has 2 alleles)
        n.p1.alleles <- (sum(window.calls == "p1") * 2) + (sum(window.calls == "het"))
        n.p2.alleles <- (sum(window.calls == "p2") * 2) + (sum(window.calls == "het"))
        total.alleles <- sum(window.calls != "missing") * 2
        # from those not missing, what proportion is p1, p2 and het?
        prop.p1 <- n.p1.alleles/total.alleles
        
        # define consensus based on threshold (p1: p1>0.7, p2: p1<0.3, het: 0.3<p1<0.7)
        if (prop.p1 >= 0.7) {
          # assign parent1 genotype of first snp on window to consensus
          window.consensus <- append(window.consensus, window[1:window.step, 1])
        } else if (prop.p1 <= 0.3) {
          # assign parent2 genotype of first snp on window to consensus
          window.consensus <- append(window.consensus, window[1:window.step, 2])
        } else {
          # assign het genotype to consensus
          p1.allele <- sapply(window[1:window.step, 1], function(x) return(unlist(strsplit(x, split = ""))[1]))
          p2.allele <- sapply(window[1:window.step, 2], function(x) return(unlist(strsplit(x, split = ""))[1]))
          window.consensus <- append(window.consensus, paste0(p1.allele, p2.allele))
        }
      } else {
        # if there is very few or none ril SNPs genotyped, consider consensus as missing data
        window.consensus <- append(window.consensus, rep("NN", times = window.step))
      }
      
      # set up the start of next window
      window.start <- window.start + window.step
      window.stop <- window.start + (window.size - 1)
    }
    
    # after that, add (window.size - 1) NNs in the consensus
    window.consensus <- append(window.consensus, rep("NN", times = NROW(geno.data.infile.chr) - length(window.consensus)))
    
    window.consensus
    
  }
  stopImplicitCluster()
  
  # correct column names
  colnames(geno.data.infile.window) <- colnames(geno.data.infile.chr)[12:NCOL(geno.data.infile.chr)]
  
  # create hapmap again
  geno.data.outfile.chr <- cbind(geno.data.infile.chr[, 1:11], geno.data.infile.window, stringsAsFactors = FALSE)
  
  # bind to final data frame
  geno.data.outfile <- rbind(geno.data.outfile, geno.data.outfile.chr)
}



# write file
outfile <- gsub(".projected.hmp.txt", ".projected.sliding-window.hmp.txt", data.filename, fixed = TRUE)
fwrite(geno.data.outfile, outfile, quote = FALSE, sep = "\t", na = "NA", row.names = FALSE)



