library(data.table)

usage <- function() {
  cat("
description: filter a hapmap file to have only polymorphic markers from a list of marker IDs.

usage: Rscript create_dataset_for_prediction.R [infile_name] [outfile_name] [...]

positional arguments:
  infile_name             hapmap file containing SNPs and SVs
  outfile_name            name of filtered dataset
  --markers-list=FILE     single-column file containing marker IDs to keep.
                          Use '--markers-list=ALL' to keep all markers from 'infile_name'
  
optional argument:
  --help                  show this helpful message
  
    
"
  )
}

getArgValue <- function(arg) {
  
  # get value from command-line arguments
  arg <- unlist(strsplit(arg, "="))
  if (length(arg) == 1) return(TRUE)
  if (length(arg) > 1) return(arg[2])
  
}


#### command line options ----

args <- commandArgs(trailingOnly = TRUE)
if ("--help" %in% args) usage() & q(save = "no")

# assert to have the correct number of arguments
if (length(args) < 3) stop(usage(), "missing positional argument(s)")

# get positional arguments
infile_name <- args[1]
outfile_name <- args[2]
markers_list <- getArgValue(args[3])

# infile_name <- "data/test_usda_rils_projected-SVs-only.hmp.txt"
# outfile_name <- "analysis/datasets/test_svs-only.hmp.txt"
# markers_list <- "analysis/ld/test_snp-names_highest-ld.no-duplicates.txt"



#### create dataset ----

# load hapmap file
geno_data <- fread(infile_name, header = TRUE, data.table = FALSE)

if (markers_list != "ALL") {
  # get markers to keep
  markers_to_keep <- fread(markers_list, header = TRUE, data.table = FALSE)
  markers_to_keep <- markers_to_keep[, 1]
  
  # filter hapmap
  geno_data <- geno_data[which(geno_data[, 1] %in% markers_to_keep), ]
}

# write file
fwrite(x = geno_data, file = outfile_name, quote = FALSE, sep = "\t", na = NA, row.names = FALSE)
