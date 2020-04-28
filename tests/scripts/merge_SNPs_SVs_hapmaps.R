#### description ----

# this script merges SNPs and SVs hapmap files from usda parents to be used in Tassel 5 when
# projecting SVs into RILs



#### libraries ----

library(data.table)



#### function ----

RemoveSNPsInPAVs <- function(snp_file, parental_sv_file) {
  
  # load data
  snp.hmp <- fread(snp_file, header = TRUE, data.table = FALSE)
  sv.hmp <- fread(parental_sv_file, header = TRUE, data.table = FALSE)
  
  # create list to store SNPs within SVs
  SNPs.in.SVs <- list()
  
  # parse by chromosome
  for (chr in 1:length(unique(snp.hmp[, "chrom"]))) {
    
    cat("Analyzing chromosome ", chr, "...\n", sep = "")
    
    # get each SV ID in a chromosome
    sv.hmp.subset <- sv.hmp[which(sv.hmp[, "chrom"] == chr), 1]
    # get each SNP name and position in a chromosome
    snp.hmp.subset <- snp.hmp[which(snp.hmp[, "chrom"] == chr), c(1, 4)]
    
    for (SV in sv.hmp.subset) {
      
      # check SNPs that fall in deletions only (not dups)
      SV.type <- unlist(strsplit(SV, split = ".", fixed = TRUE))[1]
      
      if (SV.type == "del") {
        # get start and end positions of SV
        SV.start <- unlist(strsplit(SV, split = ".", fixed = TRUE))[3]
        SV.end <- unlist(strsplit(SV, split = ".", fixed = TRUE))[4]
        # check if a SNP is in a SV boundary
        SNP.within.del <- snp.hmp.subset[which(snp.hmp.subset[, 2] %in% seq(SV.start, SV.end, by = 1)), 1]
        
        if (length(SNP.within.del) > 0) {
          # append SNPs to SV list of chromosome
          SNPs.in.SVs[[paste0("chr", chr)]] <- append(SNPs.in.SVs[[paste0("chr", chr)]], SNP.within.del)
        }
      }
      
    }
    
    # remove redundancy
    SNPs.in.SVs[[paste0("chr", chr)]] <- unique(SNPs.in.SVs[[paste0("chr", chr)]])
    
    cat(length(SNPs.in.SVs[[paste0("chr", chr)]]), "out of", NROW(snp.hmp.subset), "SNPs are within a PAV\n")
    
  }
  
  # filter hapmap so it doesn't have the SNPs in the SNPs.in.SVs list
  snp.hmp.filtered <- snp.hmp[which(!snp.hmp[, 1] %in% as.character(unlist(SNPs.in.SVs))), ]
  
  return(snp.hmp.filtered)

}


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



#### filter SNPs within SVs boundaries ----

# from parents
snp.hmp.parents <- RemoveSNPsInPAVs(snp_file = "data/usda_22kSNPs_7parents.sorted.diploid.v4.PHG35-corrected.hmp.txt",
                                     parental_sv_file = "tests/data/usda_SVs_7parents.sorted.hmp.txt")

fwrite(snp.hmp.parents, "tests/data/usda_18kSNPs_7parents.not-in-PAVs.hmp.txt",
       sep = "\t", quote = FALSE, na = NA)


# from RILs
snp.hmp.rils <- RemoveSNPsInPAVs(snp_file = "data/usda_22kSNPs_325rils.sorted.diploid.v4.hmp.txt",
                                 parental_sv_file = "tests/data/usda_SVs_7parents.sorted.hmp.txt")

fwrite(snp.hmp.rils, "tests/data/usda_18kSNPs_325rils.not-in-SVs.hmp.txt",
       sep = "\t", quote = FALSE, na = NA)




#### merge hapmaps from parents ----

parents.merged <- MergeHapmaps(snp_file = "tests/data/usda_18kSNPs_7parents.not-in-PAVs.hmp.txt",
                               parental_sv_file = "tests/data/usda_SVs_7parents.sorted.hmp.txt",
                               merge_RILs = FALSE)

fwrite(parents.merged, "tests/data/usda_SNPs-SVs_7parents.not-in-PAVs.hmp.txt",
       sep = "\t", quote = FALSE, na = NA)



#### merge hapmaps from RILs ----

rils.merged <- MergeHapmaps(snp_file = "tests/data/usda_18kSNPs_325rils.not-in-SVs.hmp.txt",
                            parental_sv_file = "tests/data/usda_SVs_7parents.sorted.hmp.txt",
                            merge_RILs = TRUE)

fwrite(rils.merged, "tests/data/usda_SNPs-SVs_325rils.not-in-SVs.hmp.txt",
       sep = "\t", quote = FALSE, na = NA)



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






#### USING ALL PARENTS FROM RESEQ ----

# filter SNPs within SVs boundaries
# from parents...
snp.hmp.parents <- RemoveSNPsInPAVs(snp_file = "data/usda_22kSNPs_7parents.sorted.diploid.v4.all-parents-corrected.hmp.txt",
                                    parental_sv_file = "tests/data/usda_SVs_7parents.sorted.hmp.txt")

fwrite(snp.hmp.parents, "tests/data/usda_18kSNPs_7parents.reseq.not-in-PAVs.hmp.txt",
       sep = "\t", quote = FALSE, na = NA)

# from RILs...
snp.hmp.rils <- RemoveSNPsInPAVs(snp_file = "data/usda_22kSNPs_325rils.sorted.diploid.v4.hmp.txt",
                                 parental_sv_file = "tests/data/usda_SVs_7parents.sorted.hmp.txt")

fwrite(snp.hmp.rils, "tests/data/usda_18kSNPs_325rils.not-in-SVs.hmp.txt",
       sep = "\t", quote = FALSE, na = NA)



# merge hapmaps from parents
parents.merged <- MergeHapmaps(snp_file = "tests/data/usda_18kSNPs_7parents.reseq.not-in-PAVs.hmp.txt",
                               parental_sv_file = "tests/data/usda_SVs_7parents.sorted.hmp.txt",
                               merge_RILs = FALSE)
# merge hapmaps from rils
rils.merged <- MergeHapmaps(snp_file = "tests/data/usda_18kSNPs_325rils.not-in-SVs.hmp.txt",
                            parental_sv_file = "tests/data/usda_SVs_7parents.sorted.hmp.txt",
                            merge_RILs = TRUE)

# subset by cross
df.crosses <- fread("tests/data/usda_biparental-crosses.txt",
                    header = TRUE, data.table = FALSE)

# make dir to store results
dir.create("tests/data/merged_hapmaps_by_cross_reseq-parents/")

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
    out_parents <- paste0("tests/data/merged_hapmaps_by_cross_reseq-parents/usda_SNPs-SVs_", cross, "_parents.sorted.hmp.txt")
    fwrite(subset.parents.merged, file = out_parents, sep = "\t", na = "NN", quote = FALSE)
    
    # RILs
    out_rils <- paste0("tests/data/merged_hapmaps_by_cross_reseq-parents/usda_SNPs-SVs_", cross, "_RILs.sorted.hmp.txt")
    fwrite(subset.rils.merged, file = out_rils, sep = "\t", na = "NN", quote = FALSE)
  }
}





#### USING DATA WITH 15K ONLY ----

# filter SNPs within SVs boundaries

# parents...
snp.hmp.parents <- RemoveSNPsInPAVs(snp_file = "data/usda_15kSNPs_7parents.sorted.diploid.v4.PHG35-corrected.hmp.txt",
                                    parental_sv_file = "tests/data/usda_SVs_7parents.sorted.hmp.txt")

fwrite(snp.hmp.parents, "tests/data/usda_14kSNPs_7parents.not-in-PAVs.hmp.txt",
       sep = "\t", quote = FALSE, na = NA)

# RILs...
snp.hmp.rils <- RemoveSNPsInPAVs(snp_file = "data/usda_15kSNPs_325rils.sorted.diploid.v4.PHG35-corrected.hmp.txt",
                                 parental_sv_file = "tests/data/usda_SVs_7parents.sorted.hmp.txt")

fwrite(snp.hmp.rils, "tests/data/usda_14kSNPs_325rils.not-in-SVs.hmp.txt",
       sep = "\t", quote = FALSE, na = NA)


# merge hapmaps from parents
parents.merged <- MergeHapmaps(snp_file = "tests/data/usda_14kSNPs_7parents.not-in-PAVs.hmp.txt",
                               parental_sv_file = "tests/data/usda_SVs_7parents.sorted.hmp.txt",
                               merge_RILs = FALSE)

# merge hapmaps from RILs
rils.merged <- MergeHapmaps(snp_file = "tests/data/usda_14kSNPs_325rils.not-in-SVs.hmp.txt",
                            parental_sv_file = "tests/data/usda_SVs_7parents.sorted.hmp.txt",
                            merge_RILs = TRUE)



# subset by cross
# load table with cross information
df.crosses <- fread("tests/data/usda_biparental-crosses.txt",
                    header = TRUE, data.table = FALSE)

# make dir to store results
dir.create("tests/data/merged_hapmaps_by_cross_15k")

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
    out_parents <- paste0("tests/data/merged_hapmaps_by_cross_15k/usda_SNPs-SVs_", cross, "_parents.sorted.hmp.txt")
    fwrite(subset.parents.merged, file = out_parents, sep = "\t", na = "NN", quote = FALSE)
    
    # RILs
    out_rils <- paste0("tests/data/merged_hapmaps_by_cross_15k/usda_SNPs-SVs_", cross, "_RILs.sorted.hmp.txt")
    fwrite(subset.rils.merged, file = out_rils, sep = "\t", na = "NN", quote = FALSE)
  }
}






#### USING ALL PARENTS FROM RESEQ WITH 15k SNPs ----

# filter SNPs within SVs boundaries
# from parents...
snp.hmp.parents <- RemoveSNPsInPAVs(snp_file = "data/usda_15kSNPs_7parents.sorted.diploid.v4.all-parents-corrected.hmp.txt",
                                    parental_sv_file = "tests/data/usda_SVs_7parents.sorted.hmp.txt")

fwrite(snp.hmp.parents, "tests/data/usda_14kSNPs_7parents.reseq.not-in-PAVs.hmp.txt",
       sep = "\t", quote = FALSE, na = NA)

# from RILs...
snp.hmp.rils <- RemoveSNPsInPAVs(snp_file = "data/usda_15kSNPs_325rils.sorted.diploid.v4.all-parents-corrected.hmp.txt",
                                 parental_sv_file = "tests/data/usda_SVs_7parents.sorted.hmp.txt")

fwrite(snp.hmp.rils, "tests/data/usda_14kSNPs_325rils.reseq.not-in-SVs.hmp.txt",
       sep = "\t", quote = FALSE, na = NA)



# merge hapmaps from parents
parents.merged <- MergeHapmaps(snp_file = "tests/data/usda_14kSNPs_7parents.reseq.not-in-PAVs.hmp.txt",
                               parental_sv_file = "tests/data/usda_SVs_7parents.sorted.hmp.txt",
                               merge_RILs = FALSE)
# merge hapmaps from rils
rils.merged <- MergeHapmaps(snp_file = "tests/data/usda_14kSNPs_325rils.reseq.not-in-SVs.hmp.txt",
                            parental_sv_file = "tests/data/usda_SVs_7parents.sorted.hmp.txt",
                            merge_RILs = TRUE)

# subset by cross
df.crosses <- fread("tests/data/usda_biparental-crosses.txt",
                    header = TRUE, data.table = FALSE)

# make dir to store results
dir.create("tests/data/merged_hapmaps_by_cross_reseq-parents_15k/")

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
    out_parents <- paste0("tests/data/merged_hapmaps_by_cross_reseq-parents_15k/usda_SNPs-SVs_", cross, "_parents.sorted.hmp.txt")
    fwrite(subset.parents.merged, file = out_parents, sep = "\t", na = "NN", quote = FALSE)
    
    # RILs
    out_rils <- paste0("tests/data/merged_hapmaps_by_cross_reseq-parents_15k/usda_SNPs-SVs_", cross, "_RILs.sorted.hmp.txt")
    fwrite(subset.rils.merged, file = out_rils, sep = "\t", na = "NN", quote = FALSE)
  }
}
