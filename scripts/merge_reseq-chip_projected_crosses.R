#### arguments for command line ----

args <- commandArgs(trailingOnly = TRUE)

# help
if (all(length(args) == 1 & args == "-h" | args == "--help")) {
  cat("
      Description: this script merges the results of resequencing SNP and SNP chip projections of RILs into one file

      
      Usage: Rscript merge_reseq-chip_projected_crosses.R [...]")
  quit()
}

# make sure the correct number of arguments are used
if (length(args) != 2) {
  stop("incorrect number of arguments provided.
       
       Usage: Rscript merge_reseq-chip_projected_crosses.R [...]
       ")
}

# assign arguments to variables
hmp.chip.file <- args[1]
hmp.reseq.file <- args[2]


# hmp.chip.file <- "data/usda_SNPs-SVs_325rils.not-in-SVs.projected.chr1.hmp.txt"
# hmp.reseq.file <- "data/reseq_snps/biomAP_rils_SNPs-reseq.projected.chr1.hmp.txt"



#### libraries ----

if(!require("data.table")) install.packages("data.table")



#### merge files ----


# load files
cat("Loading data\n")
hmp.chip <- fread(hmp.chip.file, header = TRUE, data.table = FALSE)
hmp.reseq <- fread(hmp.reseq.file, header = TRUE, data.table = FALSE)

# rearrange SNP chip column names based on reseq data
hmp.chip <- hmp.chip[, colnames(hmp.reseq)]

# merge files
cat("Merging data\n")
hmp.merged <- rbind(hmp.reseq, hmp.chip, stringsAsFactors = FALSE)


# write results
outfile <- gsub(".hmp.txt", ".reseq-SNPs.hmp.txt", hmp.chip.file, fixed = TRUE)
fwrite(hmp.merged, outfile, sep = "\t", quote = FALSE, na = NA)


# fix allele columns using TASSEL
cat("Sorting file\n")
commands.tassel.sort <- paste0("/home/hirschc1/della028/software/tassel-5-standalone/run_pipeline.pl",
                               " -Xmx100g -SortGenotypeFilePlugin -inputFile ", outfile,
                               " -outputFile ", outfile,
                               " -fileType Hapmap")
system(commands.tassel.sort)

cat("Converting to hapmap\n")
commands.tassel.hapdip <- paste0("/home/hirschc1/della028/software/tassel-5-standalone/run_pipeline.pl",
                                 " -Xmx100g -importGuess ", outfile,
                                 " -export ", outfile,
                                 " -exportType HapmapDiploid")
system(commands.tassel.hapdip)
