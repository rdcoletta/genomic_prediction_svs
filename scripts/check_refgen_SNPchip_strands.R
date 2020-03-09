#### arguments for command line ----

args <- commandArgs(trailingOnly = TRUE)

# help
if (all(length(args) == 1 & args == "-h" | args == "--help")) {
  cat("
      Description: check which reference genome the SNP chip data belongs, while accounting for probes in different strands.
      
      Usage: Rscript check_refgen_SNPchip_strands.R [file_ref-markers_chr1]")
  quit()
}

# make sure the correct number of arguments are used
if (length(args) != 1) {
  stop("incorrect number of arguments provided.
       
       Usage: Rscript check_refgen_SNPchip_strands.R [file_ref-markers_chr1]
       ")
}

# assign arguments to variables
infile <- args[1]

# infile <- "data/check_refgen_SNPchip/ref-markers_chr1.txt"




#### check strand ----

# read data
chr1.refs <- read.table(infile, header = TRUE, stringsAsFactors = FALSE)

# do the following for each version of ref gen
for (version in c("v1", "v2", "v3", "v4")) {
  # filter dataset to have only the genotyped marker and the ref gen version
  df <- chr1.refs[, c("marker", version)]
  # keep only markers that differ
  df.diff.alleles <- df[which(df[, "marker"] != df[, version]), ]
  
  # get proportion of markers that match between genotyped and ref gen version
  n.diff <- NROW(df.diff.alleles)
  n.total <- NROW(df)
  prop.same <- 1 - (n.diff / n.total)
  
  cat("Same alleles between marker and", version, "is", prop.same, "\n")
  
  # transform one column to its complement
  df.diff.alleles[, "marker"] <- sapply(df.diff.alleles[, "marker"], FUN = function(base) {
    if (base == "A") {
      base <- "T"
    } else if (base == "T") {
      base <- "A"
    } else if (base == "C") {
      base <- "G"
    } else if (base == "G") {
      base <- "C"
    } else {
      base <- NA
    }
  })
  
  # remove any missing data
  n.missing <- sum(is.na(df.diff.alleles[, "marker"]))
  cat("Missing alleles removed:", n.missing, "\n")
  df.diff.alleles <- df.diff.alleles[which(!is.na(df.diff.alleles[, "marker"])), ]
  
  # # get proportion of markers that match between complement genotyped and ref gen version
  n.diff.revcomp <- NROW(df[which(df.diff.alleles[, "marker"] != df.diff.alleles[, version]), ])  # 0
  n.total.revcomp <- NROW(df.diff.alleles)
  prop.same.revcomp <- 1 - (n.diff.revcomp / n.total.revcomp)
  
  cat("Reverse complement alleles between marker and", version, "is", prop.same.revcomp, "\n\n")
}
