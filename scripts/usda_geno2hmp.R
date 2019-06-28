#### description ----

# tranform genotypic data of 22k SNPs from both parental lines and RILs from USDA project



#### libraries used ----

library(data.table)



#### functions ----

geno2hmp <- function(geno_file, outfile_path, id_table_outfile_path, names = TRUE) {
  # read in table
  geno_data <- fread(geno_file, header = FALSE, data.table = F)
  
  # create hapmap structrue
  hmp_data <- data.frame("rs" = geno_data[-(1:2),1],
                         "alleles" = NA,
                         "chrom" = geno_data[-(1:2),2],
                         "pos" = geno_data[-(1:2),3],
                         "strand" = NA,
                         "assembly" = NA,
                         "center" = NA,
                         "protLSID" = NA,
                         "assayLSID" = NA,
                         "panel" = NA,
                         "QCcode" = NA)
  
  # add genotypes
  hmp_data <- cbind(hmp_data, geno_data[-(1:2),4:NCOL(geno_data)])
  
  if (names == TRUE) {
    # include genotype names as column names
    colnames(hmp_data)[12:NCOL(hmp_data)] <- as.character(geno_data[1,4:NCOL(geno_data)])
  }
  if (names == FALSE) {
    # include genotype IDs as column names
    colnames(hmp_data)[12:NCOL(hmp_data)] <- as.character(geno_data[2,4:NCOL(geno_data)])
  }
  
  # print hapmap converted genotypic data
  fwrite(hmp_data, file = outfile_path, sep = "\t", na = NA, quote = FALSE)
  
  
  # also print a table relating parent name and its ID when genotyping
  id_table <- data.frame("genotype_name" = as.character(geno_data[1,4:NCOL(geno_data)]),
                         "genotype_id" = as.character(geno_data[2,4:NCOL(geno_data)]))
  fwrite(id_table, file = id_table_outfile_path, sep = "\t")
}



#### parental data ----

geno2hmp(geno_file = "data/SNP_chip/Final_22kSNPs_DAS_UIUC_ParentalGenotypeData_122318.csv",
         outfile_path = "data/Final_22kSNPs_DAS_UIUC_ParentalGenotypeData_122318.hmp.txt",
         id_table_outfile_path = "data/id_table_22kSNPs_DAS_UIUC_ParentalGenotypeData.txt",
         names = TRUE)



#### RIL data ----

# this dataset is divided in two separate files. First, i will use function above to convert part 1
# to hapmap, and then will append only genotypic information of part 2 into part 1

geno2hmp(geno_file = "data/SNP_chip/Final_22kSNPs_DAS_UIUC_RILsGenotypeData-Part1_122318.csv",
         outfile_path = "data/Final_22kSNPs_DAS_UIUC_RILsGenotypeData_122318.hmp.txt",
         id_table_outfile_path = "data/id_table_22kSNPs_DAS_UIUC_RILsGenotypeData.txt",
         names = FALSE)


# read ril data again
ril_data1 <- fread("data/Final_22kSNPs_DAS_UIUC_RILsGenotypeData_122318.hmp.txt",
                   header = TRUE, data.table = F)

ril_data2 <- fread("data/SNP_chip/Final_22kSNPs_DAS_UIUC_RILsGenotypeData-Part2_122318.csv",
                   header = FALSE, data.table = F)

# filter part 2 so it has only the genotypic information
ril_data2_filter <- ril_data2[-(1:2),4:NCOL(ril_data2)]
colnames(ril_data2_filter) <- ril_data2[2,4:NCOL(ril_data2)]

# append genotypic info of part 2 to 
ril_data_combined <- cbind(ril_data1, ril_data2_filter)
View(ril_data_combined)

# overwrite RIL data so it has both parts now
fwrite(ril_data_combined, file = "data/Final_22kSNPs_DAS_UIUC_RILsGenotypeData_122318.hmp.txt",
       sep = "\t", na = NA, quote = FALSE)


# append the rest of ids of part 2 on table that already has the part 1
id_table_rils <- fread("data/id_table_22kSNPs_DAS_UIUC_RILsGenotypeData.txt",
                       header = TRUE, data.table = F)

id_table_ril_part2 <- data.frame("genotype_name" = as.character(ril_data2[1,4:NCOL(ril_data2)]),
                                 "genotype_id" = as.character(ril_data2[2,4:NCOL(ril_data2)]))

id_table_rils_combined <- rbind(id_table_rils, id_table_ril_part2)

# overwrite previous table file
fwrite(id_table_rils_combined, file = "data/id_table_22kSNPs_DAS_UIUC_RILsGenotypeData.txt", sep = "\t")
