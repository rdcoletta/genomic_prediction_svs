library(data.table)

# assign arguments to variables
args <- commandArgs(trailingOnly = TRUE)
hmp_file <- args[1]
# hmp_file <- "data/hapmap_by_cross/usda_parents_SV-SNPchip-polySNPreseq.B73xLH82.hmp.txt"

# remove duplicates
hmp <- fread(hmp_file, header = TRUE, data.table = FALSE)

hmp_filtered <- data.frame(stringsAsFactors = FALSE)
for (chr in 1:10) {
  
  hmp_chr <- subset(hmp, chrom == chr)
  
  # find snps with duplicated positions (i.e. one snp that is present in SNP chip and resequencing)
  duplicates <- hmp_chr[which(duplicated(hmp_chr[, 4]) | duplicated(hmp_chr[, 4], fromLast = TRUE)), 1]
  # keep only name of snps to be removed (i.e. from resequencing data)
  duplicates <- duplicates[grep(paste0("^snp.", chr), duplicates, perl = TRUE)]
  
  # filter df
  hmp_chr <- hmp_chr[which(!hmp_chr[, 1] %in% duplicates), ]
  hmp_filtered <- rbind(hmp_filtered, hmp_chr)
  
}

# write filtered file
fwrite(hmp_filtered, hmp_file, quote = FALSE, sep = "\t", na = NA, row.names = FALSE)
