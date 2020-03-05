#### arguments for command line ----

args <- commandArgs(trailingOnly = TRUE)

# help
if (all(length(args) == 1 & args == "-h" | args == "--help")) {
  cat("
      Description: remove extra genotypic data from hapmap files that will not be used in the USDA project

      Usage: Rscript remove_extra_geno-data.R [...]")
  quit()
}

# make sure the correct number of arguments are used
if (length(args) != 5) {
  stop("incorrect number of arguments provided.

       Usage: Rscript remove_extra_geno-data.R [...]
       ")
}

# assign arguments to variables
infile.parents <- args[1]
outfile.parents <- args[2]
infile.rils <- args[3]
outfile.rils <- args[4]
rils.to.keep <- args[5]


# infile.parents <- "data/Final_22kSNPs_DAS_UIUC_ParentalGenotypeData_122318.hmp.txt"
# outfile.parents <- "data/usda_22kSNPs_7parents.hmp.txt"
# infile.rils <- "data/Final_22kSNPs_DAS_UIUC_RILsGenotypeData_122318.hmp.txt"
# outfile.rils <- "data/usda_22kSNPs_325rils.hmp.txt"
# rils.to.keep <- "data/usda_RILs_2018.txt"




#### libraries used ----

if(!require("data.table")) install.packages("data.table")


#### parental data ----

# read in hmp file
hmp_parents <- fread(infile.parents, header = TRUE, data.table = FALSE)

# parents to keep
keep_parents <- c("B73", "PHJ40", "PHG39", "PHG47", "PH207", "PHG35", "LH82")

# remove extra parents
hmp_parents_filtered <- hmp_parents[,1:11]
hmp_parents_filtered <- cbind(hmp_parents_filtered, hmp_parents[,keep_parents])

# write filtered version
fwrite(hmp_parents_filtered, file = outfile.parents, sep = "\t", na = NA, quote = FALSE)



#### RIL data ----

# read in hmp file
hmp_rils <- fread(infile.rils, header = TRUE, data.table = FALSE)

# read id table with RIL names and IDs
name2id_table <- list.files(path = "data", pattern = "id_table", full.names = TRUE, recursive = TRUE)
name2id_table <- name2id_table[grep("RIL", name2id_table)]
name2id_table <- fread(name2id_table, header = TRUE, data.table = FALSE)

# read in ril names to keep
keep_rils <- fread(rils.to.keep, header = FALSE, data.table = FALSE)
keep_rils_name <- unique(keep_rils$V1)

# select ril IDs to keep based on ril names
keep_rils_id <- as.character(name2id_table[which(name2id_table[,"genotype_name"] %in% keep_rils_name), "genotype_id"])

# remove extra rils
hmp_rils_filtered <- hmp_rils[,1:11]
hmp_rils_filtered <- cbind(hmp_rils_filtered, hmp_rils[,keep_rils_id])

# write filtered version
fwrite(hmp_rils_filtered, file = outfile.rils, sep = "\t", na = NA, quote = FALSE)
