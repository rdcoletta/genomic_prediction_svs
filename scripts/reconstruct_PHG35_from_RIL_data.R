

crosses.file <- "tests/data/usda_biparental-crosses.txt"
data.per.cross <- "tests/data/merged_hapmaps_by_cross"
output.folder <- "tests/data/reconstructed_PHG35/"





if(!require("data.table")) install.packages("data.table")
if(!require("ggplot2")) install.packages("ggplot2")




# open file with crosses information and select only crosses with PHG35
crosses <- fread(crosses.file, header = TRUE, data.table = FALSE)
crosses <- crosses[grep("PHG35", crosses[, "cross"]), ]
crosses[, "cross"] <- gsub("*", "x", crosses[, "cross"], fixed = TRUE)

# create empty list to store PHG35 reconstruction for each cross
reconstructed.PHG35.all.crosses <- list()

for (cross in crosses[, "cross"]) {

  # only parse PHG35 crosses that will be used for projection
  donors.list <- list.files(path = data.per.cross, pattern = "parents.sorted.hmp.txt",
                            recursive = TRUE, full.names = FALSE)
  donors.list <- gsub("usda_SNPs-SVs_", "", donors.list)
  donors.list <- gsub("_parents.sorted.hmp.txt", "", donors.list)
  donors.list <- donors.list[grep("PHG35", donors.list)]

  if (cross %in% donors.list) {

    cat("Analyzing cross", cross, "\n")

    # load parental data
    parents.hmp <- list.files(path = data.per.cross, pattern = paste0(cross, "_parents.sorted.hmp.txt"),
                              recursive = TRUE, full.names = TRUE)
    parents.hmp <- fread(parents.hmp, header = TRUE, data.table = FALSE)

    # get column number of non-PHG35 parent
    parent2 <- unlist(strsplit(cross, split = "x"))
    parent2 <- parent2[grep("PHG35", parent2, invert = TRUE)]
    col.reseq.parent2 <- which(colnames(parents.hmp) == parent2)

    # load ril data
    rils.hmp <- list.files(path = data.per.cross, pattern = paste0(cross, "_RILs.sorted.hmp.txt"),
                           recursive = TRUE, full.names = TRUE)
    rils.hmp <- fread(rils.hmp, header = TRUE, data.table = FALSE)

    # create empty vector to store reconstructed PHG35 alleles
    reconstructed.PHG35 <- rep("NN", times = NROW(parents.hmp))

    for (row in 1:NROW(parents.hmp)) {

      # skip SVs
      if (!grepl("^del|^dup|^inv|^tra", parents.hmp[row, 1])) {

        # get allele on resequencing data from non-PHG35 parent
        alleles.parent2 <- unlist(strsplit(parents.hmp[row, col.reseq.parent2], split = ""))

        # only proceed if allele in non-PHG35 parent is not missing and is homozygous
        if (alleles.parent2[1] == alleles.parent2[2] & all(unique(alleles.parent2) != "N")) {

          # get alleles on ril data
          alleles.rils <- as.character(rils.hmp[row, 12:NCOL(rils.hmp)])
          alleles.rils <- unlist(strsplit(paste0(alleles.rils, collapse = ""), split = ""))
          # count allele frequency
          df.allele.count <- data.frame(table(alleles.rils), stringsAsFactors = FALSE)
          df.allele.count <- cbind(df.allele.count, percent = NA)
          # but first exclude missing allele ("N") if there are other two alleles
          if (NROW(df.allele.count) > 2) {
            df.allele.count <- df.allele.count[which(df.allele.count[, "alleles.rils"] != "N"), ]
          }
          # get percent allele frequency
          for (allele in 1:NROW(df.allele.count)) {
            allele.freq <- df.allele.count[allele, "Freq"] / sum(df.allele.count[, "Freq"])
            df.allele.count[allele, "percent"] <- round(allele.freq, digits = 2)
          }

          # keep only alleles with allele frequency > 0.05 (since it can be an artifact from SNP calling)
          alleles.rils <- df.allele.count[which(df.allele.count[, "percent"] > 0.05), "alleles.rils"]
          alleles.rils <- as.character(alleles.rils)

          # exclude allele from non-PHG35 parent to get PHG35 allele
          allele.PHG35 <- alleles.rils[!alleles.rils %in% alleles.parent2]

          # if there's only one allele in "alleles.rils" it means it is monomorphic between the two
          # parents, thus "allele.PHG35" should have the same allele as the other parent
          if (length(allele.PHG35) == 0) {
            allele.PHG35 <- unique(alleles.parent2)
          }
          # if there are two alleles remaining in "allele.PHG35" it means that there is either a
          # third allele for this locus or there is missing data on rils and I can't confirm whether
          # this missing allele is from PHG35 or parent2
          if (length(allele.PHG35) > 1) {
            allele.PHG35 <- "N"
          }

          # debug
          if (length(reconstructed.PHG35[row]) != length(paste0(allele.PHG35, allele.PHG35))) {
            print(row)
          }


          reconstructed.PHG35[row] <- paste0(allele.PHG35, allele.PHG35)
          # note: there won't be any heterozygote for reconstructed PHG35

        }
      }
    }

    # append reconstructed PHG35 for this cross into list
    reconstructed.PHG35.all.crosses[[cross]] <- reconstructed.PHG35

  }
}

# get final reconstructed PHG35 by merging results from all crosses and transforming disagreements
# in missing data
reconstructed.PHG35.df <- data.frame(do.call(cbind, reconstructed.PHG35.all.crosses), stringsAsFactors = FALSE)
reconstructed.PHG35.final <- apply(reconstructed.PHG35.df[, ], MARGIN = 1, function(genotypes) {
  # get unique genotypes among all crosses for that SNP
  genotypes <- unique(genotypes)
  # if there are more than one genotype
  if (length(genotypes) > 1) {
    # exclude missing data
    genotypes <- genotypes[genotypes != "NN"]
    # even if after excluding missing data, there is more than one unique genotype, make SNP as NN
    if (length(genotypes) > 1) {
      genotypes <- "NN"
    }
  }
  return(genotypes)
})


# sum(reconstructed.PHG35.final != "NN")


# fix all parental data with the new reconstructed PHG35 data
# make sure to not replace the SV calls though

for (cross in crosses[, "cross"]) {

  # only parse PHG35 crosses that will be used for projection
  donors.list <- list.files(path = data.per.cross, pattern = "parents.sorted.hmp.txt",
                            recursive = TRUE, full.names = FALSE)
  donors.list <- gsub("usda_SNPs-SVs_", "", donors.list)
  donors.list <- gsub("_parents.sorted.hmp.txt", "", donors.list)
  donors.list <- donors.list[grep("PHG35", donors.list)]

  if (cross %in% donors.list) {

    cat("Analyzing cross", cross, "\n")

    # load parental data
    parents.hmp <- list.files(path = data.per.cross, pattern = paste0(cross, "_parents.sorted.hmp.txt"),
                              recursive = TRUE, full.names = TRUE)
    parents.hmp <- fread(parents.hmp, header = TRUE, data.table = FALSE)

    for (row in 1:NROW(parents.hmp)) {
      # skip SVs
      if (!grepl("^del|^dup|^inv|^tra", parents.hmp[row, 1])) {
        # change resequencing data
        parents.hmp[row, "PHG35"] <- reconstructed.PHG35.final[row]
      }
    }

    # write fixed hmp file
    fixed.outfile <- list.files(path = data.per.cross, pattern = paste0(cross, "_parents.sorted.hmp.txt"),
                                recursive = TRUE, full.names = TRUE)
    fixed.outfile <- gsub("sorted.hmp.txt", "sorted.PHG35-reconstructed.hmp.txt", fixed.outfile)
    fwrite(parents.hmp, file = fixed.outfile, quote = FALSE, sep = "\t", na = "NA", row.names = FALSE)

  }
}
