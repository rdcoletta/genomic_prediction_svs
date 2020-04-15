#### arguments for command line ----

args <- commandArgs(trailingOnly = TRUE)

# help
if (all(length(args) == 1 & args == "-h" | args == "--help")) {
  cat("
Description: this script identifies duplicated individuals in a hapmap file and collapses
             their genotypic information.

Usage: Rscript collapse_geno-info_dup-RILs.R hapmap_file")
  quit()
}

# make sure the correct number of arguments are used
if (length(args) != 1) {
  stop("incorrect number of arguments provided.

Usage: Rscript collapse_geno-info_dup-RILs.R hapmap_file
       ")
}

# assign arguments to variables
hmp.file <- args[1]

# hmp.file <- "data/Final_22kSNPs_DAS_UIUC_RILsGenotypeData_122318.hmp.txt"




#### libraries ----

if(!require("data.table")) install.packages("data.table")



#### collapse genotypes ----

# load data
hmp <- fread(hmp.file, header = TRUE, data.table = FALSE)

# identify duplicated RILs
dup.rils <- colnames(hmp)[duplicated(colnames(hmp))]

# create data frame to store results
collapsed.genos <- data.frame(matrix(nrow = NROW(hmp), ncol = 0), stringsAsFactors = FALSE)
for (ril in dup.rils) {
  # collapse genotypes for that RIL (set to missing if calls disagree between duplicates)
  collapsed.ril <- apply(hmp[, which(colnames(hmp) == ril)], MARGIN = 1, function(x) {
    if (length(unique(x)) == 1) return(unique(x))
    else return("NN")
  })
  # append to data frame
  collapsed.genos <- cbind(collapsed.genos, collapsed.ril)
}
# change column names
colnames(collapsed.genos) <- dup.rils

# create final final with collapsed genotypes
hmp.final <- hmp[, unique(colnames(hmp))]
for (ril in dup.rils) {
  hmp.final[, ril] <- collapsed.genos[, ril]
}

# overwrite hapmap file
fwrite(hmp.final, file = hmp.file, sep = "\t", quote = FALSE, na = NA)
