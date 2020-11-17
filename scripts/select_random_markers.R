library(data.table)

usage <- function() {
  cat("
description: select random markers from a hapmap file.

usage: Rscript select_random_markers.R [infile_name] [outfile_name] [...]

positional arguments:
  infile_name                   hapmap file containing SNPs and/or SVs
  outfile_name                  name of filtered dataset
  markers_to_sample             number of markers to randomly sample
  
optional argument:
  --help                        show this helpful message
  --proportion-SV-SNP=[VALUE]   if hapmap contains both SNPs and SVs, choose the proportion of 
                                random markers to be sampled for each marker type. For example,
                                a value of 0.5 means that half of the random markers will be SNPs
                                and the other half will be SVs.
  --SVs-list=FILE               single-column file containing SV IDs. Providing SV IDs is necessary
                                if '--proportion-SV-SNP' option was requested.
  --seed=VALUE                  value for set.seed (default: NULL; random number is selected)
  
    
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

# set default
proportion_SV_SNP <- NULL
seed <- NULL
SVs_list <- NULL

args <- commandArgs(trailingOnly = TRUE)
if ("--help" %in% args) usage() & q(save = "no")

# assert to have the correct optional arguments
if (length(args) < 3) stop(usage(), "missing positional argument(s)")

if (length(args) > 3) {
  
  opt_args <- args[-1:-3]
  opt_args_allowed <- c("--proportion-SV-SNP", "--SVs-list", "--seed")
  opt_args_requested <- as.character(sapply(opt_args, function(x) unlist(strsplit(x, split = "="))[1]))
  if (any(!opt_args_requested %in% opt_args_allowed)) stop(usage(), "wrong optional argument(s)")
  
  # change default based on the argument provided
  for (argument in opt_args_allowed) {
    if (any(grepl(argument, opt_args_requested))) {
      arg_name <- gsub("-", "_", gsub("--", "", argument))
      arg_value <- getArgValue(opt_args[grep(argument, opt_args)])
      assign(arg_name, arg_value)
    }
  }
  
}

if (!is.null(proportion_SV_SNP)) {
  if (suppressWarnings(!any(is.na(as.numeric(proportion_SV_SNP))))) {
    proportion_SV_SNP <- as.numeric(proportion_SV_SNP)
    if (proportion_SV_SNP < 0 | proportion_SV_SNP > 1) {
      stop("Optional argument '--proportion-SV-SNP' should be a number between 0 and 1")
    }
    if (is.null(SVs_list)) {
      stop("Provide a file for '--SVs-list' option. This is required when calling '--proportion-SV-SNP'")
    }
  } else {
    stop("Optional argument '--proportion-SV-SNP' should be a number")
  }
}

if (is.null(seed)) {
  seed <- ceiling(runif(1, 0, 1000000))
} else {
  if (suppressWarnings(!any(is.na(as.numeric(seed))))) {
    seed <- as.numeric(seed) 
  } else {
    stop("Optional argument '--seed' should be a number")
  }
}


# get positional arguments
infile_name <- args[1]
outfile_name <- args[2]
markers_to_sample <- as.numeric(args[3])

if (suppressWarnings(!any(is.na(as.numeric(markers_to_sample))))) {
  seed <- as.numeric(markers_to_sample) 
} else {
  stop("Positional argument 'markers_to_sample' should be a number")
}


# infile_name <- "analysis/trait_sim/datasets/test_usda_rils.all_markers.hmp.txt"
# outfile_name <- "analysis/trait_sim/datasets/test_usda_rils.all_markers.adjusted-n-markers.hmp.txt"
# markers_to_sample <- 334
# proportion_SV_SNP <- 0.5
# SVs_list <- "data/test_SVs_IDs.txt"
# seed <- 1990



#### select random markers ----

# load hapmap
hmp_file <- fread(infile_name, header = TRUE, data.table = FALSE)

if (is.null(proportion_SV_SNP)) {

  # sample random markers
  set.seed(seed)
  rows_to_keep <- sort(sample(1:NROW(hmp_file), size = markers_to_sample, replace = FALSE))

} else {

  # load SV IDs
  SVs <- fread(SVs_list, header = FALSE, data.table = FALSE)
  SVs <- SVs[, 1]

  # get row numbers for SVs and SNPs
  SV_rows <- which(hmp_file[, 1] %in% SVs)
  SNP_rows <- which(!hmp_file[, 1] %in% SVs)

  # sample proportion of SVs and SNPs
  set.seed(seed)
  SVs_to_keep <- sample(SV_rows, size = floor(markers_to_sample * proportion_SV_SNP), replace = FALSE)
  set.seed(seed)
  SNPs_to_keep <- sample(SNP_rows, size = ceiling(markers_to_sample * proportion_SV_SNP), replace = FALSE)

  # merge rows to keep
  rows_to_keep <- sort(c(SVs_to_keep, SNPs_to_keep))

}

# sample dataset
hmp_file <- hmp_file[rows_to_keep, ]

# write file
fwrite(x = hmp_file, file = outfile_name, quote = FALSE, sep = "\t", na = NA, row.names = FALSE)
