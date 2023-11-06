library(data.table)


usage <- function() {
  cat("
description: add B73 alleles to hapmap file.

usage: Rscript add_B73_alleles_to_hmp.R [hmp_file] [b73_pos_file] [outfile] [...]

positional arguments:
  hmp_file                path to folder containing simulated traits from simplePHENOTYPES
  b73_pos_file            single-column file containing only SV IDs
  outfile                 name of output

optional argument:
  --help                  show this helpful message

"
  )
}



#### command line options ----

args <- commandArgs(trailingOnly = TRUE)
if ("--help" %in% args) usage() & q(save = "no")

# assert to have the correct optional arguments
if (length(args) != 3) stop(usage(), "missing positional argument(s)")

# get positional arguments
hmp_file <- args[1]
b73_pos_file <- args[2]
outfile <- args[3]
# hmp_file <- "data/widiv_snps.Qiu-et-al.no-B73.hmp.txt"
# b73_pos_file <- "data/widiv_snps.B73-alleles-only.txt"
# outfile <- "data/widiv_snps.Qiu-et-al.hmp.txt"


#### add alleles to hmp ----

# load files
hmp <- fread(hmp_file, header = TRUE, data.table = FALSE)
b73_pos <- fread(b73_pos_file, header = FALSE, data.table = FALSE)

# transform ref alleles into homozygous state
b73_pos[, 3] <- paste0(b73_pos[, 3], b73_pos[, 3])

# order files by chrom and pos
hmp <- hmp[order(hmp$chrom, hmp$pos), ]
b73_pos <- b73_pos[order(b73_pos[, 1], b73_pos[, 2]), ]

# confirm order is the same
if (any(all(hmp$chrom != b73_pos[, 1]) & all(hmp$pos != b73_pos[, 2]))) {
  stop("marker chrom and pos don't match")
}

# merge B73 to hmp
hmp <- data.frame(hmp[, 1:11], B73 = b73_pos[, 3], hmp[, 12:ncol(hmp)])

# write file
fwrite(hmp, file = outfile, sep = "\t", quote = FALSE, na = NA, row.names = FALSE)