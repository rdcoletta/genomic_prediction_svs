#### description ----

# this script merges SNPs and SVs hapmap files from usda parents to be used in Tassel 5 when
# projecting SVs into RILs



#### libraries ----

library(data.table)



#### function ----

MergeHapmaps <- function(snp_file, parental_sv_file, merge_RILs = FALSE) {
  
  # load SNP and SV datasets
  snp.hmp <- fread(snp_file, header = TRUE, data.table = FALSE)
  sv.hmp <- fread(parental_sv_file, header = TRUE, data.table = FALSE)
  
  # if want to merge SNPs and SVs from parents...
  if (merge_RILs == FALSE) {
    
    # make sure columns have the same names and order in both files
    sv.hmp <- cbind(sv.hmp[, 1:11],
                    sv.hmp[, colnames(snp.hmp)[12:NCOL(snp.hmp)]])
    colnames(sv.hmp) <- colnames(snp.hmp)
    
    # merge hapmaps
    merged.hmp <- rbind(snp.hmp, sv.hmp)
  }
  
  # if want to create dataset for projection of SVs from parents to RILs...
  if (merge_RILs == TRUE) {
    
    # since i want to impute SVs from parents to RILs, I need to also include the SV marker names in
    # the RIL hapmap. The values for these markers will be NA since we don't know yet if there is
    # a SV there
    
    # create empty data frame with nrow = number of SVs, and ncol = number of RILs to be imputed
    rils.genos <- data.frame(matrix(NA, nrow = NROW(sv.hmp), ncol = length(12:NCOL(snp.hmp))))
    colnames(rils.genos) <- colnames(snp.hmp)[12:NCOL(snp.hmp)]
    
    # make sure values for alleles and genotypes are NAs before merging with RIL hmp
    rils.sv.hmp <- cbind(sv.hmp[, 1:11], rils.genos)
    rils.sv.hmp$alleles <- NA
    # also make sure the hapmap columns are the same
    colnames(rils.sv.hmp)[1:11] <- colnames(snp.hmp)[1:11]

    
    # merge hapmaps
    merged.hmp <- rbind(snp.hmp, rils.sv.hmp)
  }
  
  # sort by chromosome and position
  merged.hmp <- merged.hmp[order(merged.hmp$chrom, merged.hmp$pos),]
  
  return(merged.hmp)
  
}



#### merge hapmaps from parents ----

parents.merged <- MergeHapmaps(snp_file = "tests/data/usda_22kSNPs_7parents.sorted.diploid.hmp.txt",
                               parental_sv_file = "tests/data/usda_SVs_7parents.sorted.hmp.txt",
                               merge_RILs = FALSE)



#### merge hapmaps from RILs ----

rils.merged <- MergeHapmaps(snp_file = "tests/data/usda_22kSNPs_325rils.sorted.diploid.hmp.txt",
                            parental_sv_file = "tests/data/usda_SVs_7parents.sorted.hmp.txt",
                            merge_RILs = TRUE)



#### subset by cross ----

# load table with cross information
df.crosses <- fread("tests/data/usda_biparental-crosses.txt",
                           header = TRUE, data.table = FALSE)

# make dir to store results
dir.create("tests/data/merged_hapmaps_by_cross")

# write merged RILs hapmap by cross
for (row in 1:NROW(df.crosses)) {
  
  # get cross name, and change * by x in the name
  cross <- df.crosses[row, "cross"]
  cross <- gsub(pattern = "*", replacement = "x", x = cross, fixed = TRUE)
  # get RIL names for that cross
  RILs <- unlist(strsplit(df.crosses[row, "RILs"], split = ","))
  
  # subset merged hapmaps from parents
  parents <- unlist(strsplit(cross, split = "x"))
  subset.parents.merged <- cbind(parents.merged[, 1:11],
                                 parents.merged[, which(colnames(parents.merged) %in% parents)])
  # since i'm subsetting, the alleles column might be wrong...
  subset.parents.merged$alleles <- NA
  
  # subset merged hapmaps from RILs
  subset.rils.merged <- cbind(rils.merged[, 1:11],
                              rils.merged[, which(colnames(rils.merged) %in% RILs)])
  # since i'm subsetting, the alleles column might be wrong...
  subset.rils.merged$alleles <- NA
  
  # write results only if RILs were actually planted (i.e., if RILs are in the merged file)
  if (NCOL(subset.rils.merged) > 11) {
    # parents
    out_parents <- paste0("tests/data/merged_hapmaps_by_cross/usda_SNPs-SVs_", cross, "_parents.sorted.hmp.txt")
    fwrite(subset.parents.merged, file = out_parents, sep = "\t", na = "NN", quote = FALSE)
    
    # RILs
    out_rils <- paste0("tests/data/merged_hapmaps_by_cross/usda_SNPs-SVs_", cross, "_RILs.sorted.hmp.txt")
    fwrite(subset.rils.merged, file = out_rils, sep = "\t", na = "NN", quote = FALSE)
  }
}
