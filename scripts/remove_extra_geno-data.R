#### description ----

# remove extra genotypic data from hapmap files that will not be used in the USDA project



#### libraries used ----

library(data.table)



#### parental data ----

# read in hmp file
hmp_parents <- fread("data/Final_22kSNPs_DAS_UIUC_ParentalGenotypeData_122318.hmp.txt",
                     header = TRUE, data.table = FALSE)

# parents to keep
keep_parents <- c("B73", "PHJ40", "PHG39", "PHG47", "PH207", "PHG35", "LH82")

# remove extra parents
hmp_parents_filtered <- hmp_parents[,1:11]
hmp_parents_filtered <- cbind(hmp_parents_filtered, hmp_parents[,keep_parents])

# write filtered version
fwrite(hmp_parents_filtered, file = "data/usda_22kSNPs_7parents.hmp.txt", sep = "\t", na = NA, quote = FALSE)



#### RIL data ----

# read in hmp file
hmp_rils <- fread("data/Final_22kSNPs_DAS_UIUC_RILsGenotypeData_122318.hmp.txt", header = TRUE,
                  data.table = FALSE)

# read id table with RIL names and IDs
name2id_table <- fread("data/id_table_22kSNPs_DAS_UIUC_RILsGenotypeData.txt", header = TRUE,
                       data.table = FALSE)

# read in ril names to keep
keep_rils <- fread("data/usda_RILs_2018.txt", header = FALSE, data.table = FALSE)
keep_rils_name <- unique(keep_rils$V1)

# select ril IDs to keep based on ril names
keep_rils_id <- as.character(name2id_table[which(name2id_table[,"genotype_name"] %in% keep_rils_name), "genotype_id"])

# remove extra rils
hmp_rils_filtered <- hmp_rils[,1:11]
hmp_rils_filtered <- cbind(hmp_rils_filtered, hmp_rils[,keep_rils_id])

# write filtered version
fwrite(hmp_rils_filtered, file = "data/usda_22kSNPs_325rils.hmp.txt", sep = "\t", na = NA, quote = FALSE)
