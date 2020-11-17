library(data.table)

usage <- function() {
  cat("
description: remove markers that are monomorphic across all individuals in a hapmap file.

usage: Rscript remove_mono_markers.R [infile_name] [outfile_name]

positional arguments:
  infile_name             hapmap file containing SNPs and SVs
  outfile_name            name of output file

optional argument:
  --help                  show this helpful message

"
  )
}

removeMonoMarkers <- function(hmp) {
  
  hmp_poly <- apply(hmp, MARGIN = 1, function(marker) {
    
    geno <- unique(marker[12:length(marker)])
    geno <- geno[geno != "NN"]
    if (length(geno) > 1)  return(marker)
    
  })
  hmp_poly <- data.frame(do.call(rbind, hmp_poly), stringsAsFactors = FALSE, check.names = FALSE)
  # fix chromsome and position columns (remove whitespace added by do.call function)
  hmp_poly[, "chrom"] <- as.integer(gsub("^ +", "", hmp_poly[, "chrom"], perl = TRUE))
  hmp_poly[, "pos"] <- as.integer(gsub("^ +", "", hmp_poly[, "pos"], perl = TRUE))
  
  return(hmp_poly)
  
}




#### command line options ----

args <- commandArgs(trailingOnly = TRUE)
if ("--help" %in% args) usage() & q(save = "no")

# assert to have the correct optional arguments
if (length(args) < 2) stop(usage(), "missing positional argument(s)")

# get arguments
infile_name <- args[1]
outfile_name <- args[2]
# infile_name <- "data/test_usda_rils_projected-SVs-only.hmp.txt"
# out_folder <- "data/test_usda_rils_projected-SVs-only.poly.hmp.txt"



#### remove monomorphic markers ----

# load data
geno_data <- fread(infile_name, header = TRUE, data.table = FALSE)

# remove monomorphic markers (across all individuals)
geno_data <- removeMonoMarkers(geno_data)

# write file
fwrite(x = geno_data, file = outfile_name, quote = FALSE, sep = "\t", na = NA, row.names = FALSE)
