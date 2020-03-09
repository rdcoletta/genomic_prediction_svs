

crosses.file <- "tests/data/usda_biparental-crosses.txt"
parents.filename <- "data/usda_22kSNPs_7parents.sorted.diploid.v4.hmp.txt"
rils.filename <- "data/usda_22kSNPs_325rils.sorted.diploid.v4.hmp.txt"
output.folder <- "tests/data/reconstructed_PHG35/"





if(!require("data.table")) install.packages("data.table")
if(!require("ggplot2")) install.packages("ggplot2")



# load hapmap files
parents.hmp <- fread(parents.filename, header = TRUE, data.table = FALSE)
rils.hmp <- fread(rils.filename, header = TRUE, data.table = FALSE)

# open file with crosses information and select only crosses with PHG35
crosses <- fread(crosses.file, header = TRUE, data.table = FALSE)
crosses <- crosses[grep("PHG35", crosses[, "cross"]), ]
crosses[, "cross"] <- gsub("*", "x", crosses[, "cross"], fixed = TRUE)

# create empty list to store PHG35 reconstruction for each cross
reconstructed.PHG35.all.crosses <- list()


for (cross in crosses[, "cross"]) {

  # select the RILs that were actually genotyped for that cross
  cross.rils <- unlist(strsplit(crosses[which(crosses[, "cross"] == cross), "RILs"], split = ","))
  cross.rils <- cross.rils[cross.rils %in% colnames(rils.hmp)]

  # only parse PHG35 crosses that have genotyped rils
  if (length(cross.rils) > 0) {

    cat("Analyzing cross", cross, "\n")

    # filter parental data to have only parents of the cross
    parent2 <- unlist(strsplit(cross, split = "x"))
    parent2 <- parent2[grep("PHG35", parent2, invert = TRUE)]
    parents.hmp.cross <- cbind(parents.hmp[, 1:11], parents.hmp[, c("PHG35", parent2)])

    # filter ril data to have only genotyped rils for that cross
    rils.hmp.cross <- cbind(rils.hmp[, 1:11], rils.hmp[, cross.rils])

    # create empty vector to store reconstructed PHG35 alleles
    reconstructed.PHG35 <- rep("NN", times = NROW(parents.hmp.cross))

    for (row in 1:NROW(parents.hmp.cross)) {

      # get allele on resequencing data from non-PHG35 parent
      alleles.parent2 <- unlist(strsplit(parents.hmp.cross[row, parent2], split = ""))

      # only proceed if allele in non-PHG35 parent is not missing and is homozygous
      if (alleles.parent2[1] == alleles.parent2[2] & all(unique(alleles.parent2) != "N")) {

        # get alleles on ril data
        alleles.rils <- as.character(rils.hmp.cross[row, 12:NCOL(rils.hmp.cross)])
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


# sum(parents.hmp[, "PHG35"] != reconstructed.PHG35.final)


# fix parental data with the new reconstructed PHG35 data
parents.hmp[, "PHG35"] <- reconstructed.PHG35.final

# write fixed hmp file
fixed.outfile <- gsub(pattern = "7parents.sorted.diploid.v4.hmp.txt",
                      replacement = "7parents.sorted.diploid.v4.PHG35-reconstructed.hmp.txt",
                      parents.filename)
fwrite(parents.hmp, file = fixed.outfile, quote = FALSE, sep = "\t", na = "NA", row.names = FALSE)




#### compare reconstructed with previous PHG35 from SNPchip and resequencing ----

# add reconstructed versions in the same hapmap
hmp.reconstructed <- cbind(parents.hmp[, 1:11],
                           reconstructed.PHG35.df,
                           final_reconstruction = reconstructed.PHG35.final)

# load other versions
hmp.reseq <- fread("data/usda_22kSNPs_7parents.sorted.diploid.v4.PHG35-corrected.hmp.txt",
                   header = TRUE, data.table = FALSE)
hmp.chip <- fread("data/usda_22kSNPs_7parents.sorted.diploid.hmp.txt",
                  header = TRUE, data.table = FALSE)

# resequencing data has the same SNPs as the reconstructed version, but SNP chip have a few more
all(hmp.reconstructed[,1] == hmp.reseq[,1])
# so need to exclude these extra snps from chip
snps.to.keep <- as.character(hmp.reconstructed[, 1])
hmp.chip <- hmp.chip[which(hmp.chip[, 1] %in% snps.to.keep), ]
hmp.chip <- hmp.chip[match(hmp.reconstructed[,1], hmp.chip[,1]), ]
all(hmp.reconstructed[,1] == hmp.chip[,1])

# merge all PHG35 versions into one dataframe
PHG35.versions <- cbind(hmp.reconstructed[, 12:NCOL(hmp.reconstructed)],
                        resequencing = hmp.reseq[, "PHG35"],
                        snp_chip = hmp.chip[, "PHG35"])

# create directory to save matrices
dir.create("analysis/qc/PHG35_reconstruction")


# make pairwise comparisons of similarities (don't count if missing in one of the genomes)
CalcSimMatrix <- function(PHG35.versions, report = c("percent", "total")) {
  
  # create empty matrix
  sim.matrix <- matrix(NA, nrow = NCOL(PHG35.versions), ncol = NCOL(PHG35.versions))
  rownames(sim.matrix) <- colnames(PHG35.versions)
  colnames(sim.matrix) <- colnames(PHG35.versions)
  
  for (v1 in 1:NCOL(PHG35.versions)) {
    for (v2 in 1:NCOL(PHG35.versions)) {
      
      # select only two versions
      comparison <- cbind(as.character(PHG35.versions[, v1]), as.character(PHG35.versions[, v2]))
      # remove missing data
      comparison <- comparison[which(comparison[, 1] != "NN" & comparison[, 2] != "NN"), ]
      # calculate similarity
      if (report == "percent")  similarity <- sum(comparison[, 1] == comparison[, 2]) / NROW(comparison)
      if (report == "total")  similarity <- sum(comparison[, 1] == comparison[, 2])
      # add results to correct place at matrix
      name.row <- colnames(PHG35.versions)[v1]
      name.col <- colnames(PHG35.versions)[v2]
      sim.matrix[name.row, name.col] <- round(similarity, digits = 2)
      
    }
  }
  
  return(sim.matrix)
}

sim.matrix.percent <- CalcSimMatrix(PHG35.versions, report = "percent")
sim.matrix <- CalcSimMatrix(PHG35.versions, report = "total")
sim.matrix[upper.tri(sim.matrix, diag = FALSE)] <- sim.matrix.percent[upper.tri(sim.matrix.percent, diag = FALSE)]

fwrite(data.frame(sim.matrix), "analysis/qc/PHG35_reconstruction/similarity_matrix.txt",
       quote = FALSE, sep = "\t", row.names = TRUE)



# make pairwise comparisons of disagreements (don't count if missing in one of the genomes)
CalcDisagreeMatrix <- function(PHG35.versions, report = c("percent", "total"), type = c("reverse", "het", "diff_allele")) {
  
  # create empty matrix
  matrix.disagree <- matrix(NA, nrow = NCOL(PHG35.versions), ncol = NCOL(PHG35.versions))
  rownames(matrix.disagree) <- colnames(PHG35.versions)
  colnames(matrix.disagree) <- colnames(PHG35.versions)
  
  for (v1 in 1:NCOL(PHG35.versions)) {
    for (v2 in 1:NCOL(PHG35.versions)) {
      
      # select only two versions of PHG35
      comparison <- cbind(as.character(PHG35.versions[, v1]), as.character(PHG35.versions[, v2]))
      # remove missing data
      comparison <- comparison[which(comparison[, 1] != "NN" & comparison[, 2] != "NN"), ]
      # select only disagreements
      comparison <- comparison[which(comparison[, 1] != comparison[, 2]), ]
      # transform second column to reverse column
      comparison[, 2] <- sapply(comparison[, 2], FUN = function(geno) {
        alleles <- unlist(strsplit(geno, split = ""))
        comp.alleles <- sapply(alleles, function(allele) {
          if (allele == "A") return("T")
          else if (allele == "T") return("A")
          else if (allele == "C") return("G")
          else if (allele == "G") return("C")
        })
        comp.alleles <- paste0(comp.alleles, collapse = "")
        return(comp.alleles)
      })
      
      # only proceed if there are disagreements
      if (NROW(comparison) > 0) {
        
        if (type == "reverse") {
          # calculate how much they disagree after phasing (i.e. controlling for reverse strand)
          if (report == "percent") disagree <- sum(comparison[, 1] == comparison[, 2]) / NROW(comparison)
          if (report == "total") disagree <- sum(comparison[, 1] == comparison[, 2])
        }
        
        if (type == "het" | type == "diff_allele") {
          # select disagreements after phasing
          comparison <- comparison[which(comparison[, 1] != comparison[, 2]), ]
          # check if homo or het
          # note: i'm not checking if one of the alleles in the het the same as the one in the other genotype
          comparison.type <- apply(comparison, MARGIN = 1, FUN = function(geno) {
            alleles.v1 <- unlist(strsplit(geno[1], split = ""))
            alleles.v2 <- unlist(strsplit(geno[2], split = ""))
            if (alleles.v1[1] == alleles.v1[2]) type.v1 <- "homo"
            if (alleles.v1[1] != alleles.v1[2]) type.v1 <- "het"
            if (alleles.v2[1] == alleles.v2[2]) type.v2 <- "homo"
            if (alleles.v2[1] != alleles.v2[2]) type.v2 <- "het"
            return(c(type.v1, type.v2))
          })
          comparison.type <- t(comparison.type)
          
          if (type == "het") {
            # calculate how much they disagree due to hets in one of the genomes
            if (report == "percent") disagree <- sum(comparison.type[, 1] == "het" | comparison.type[, 2] == "het") / NROW(comparison.type)
            if (report == "total") disagree <- sum(comparison.type[, 1] == "het" | comparison.type[, 2] == "het")
          }
          if (type == "diff_allele") {
            # calculate how much they disagree due to hets in one of the genomes
            if (report == "percent") disagree <- sum(comparison.type[, 1] != "het" & comparison.type[, 2] != "het") / NROW(comparison.type)
            if (report == "total") disagree <- sum(comparison.type[, 1] != "het" & comparison.type[, 2] != "het")
          }
        }
        
      } else {
        # if there's no disagreement...
        disagree <- 0
      }
      
      # add results to correct place at matrix
      name.row <- colnames(PHG35.versions)[v1]
      name.col <- colnames(PHG35.versions)[v2]
      matrix.disagree[name.row, name.col] <- round(disagree, digits = 2)
      
    }
  }
  
  return(matrix.disagree)
  
}

# of all that disagree, how many were actually reverse strand
disagree.rev.matrix.percent <- CalcDisagreeMatrix(PHG35.versions, report = "percent", type = "reverse")
disagree.rev.matrix <- CalcDisagreeMatrix(PHG35.versions, report = "total", type = "reverse")
disagree.rev.matrix[upper.tri(disagree.rev.matrix, diag = FALSE)] <- disagree.rev.matrix.percent[upper.tri(disagree.rev.matrix.percent, diag = FALSE)]

fwrite(data.frame(disagree.rev.matrix), "analysis/qc/PHG35_reconstruction/disagreements_rev-comp.txt",
       quote = FALSE, sep = "\t", row.names = TRUE)


# of all that disagree after controlling for phasing, how many were due to hets
disagree.het.matrix.percent <- CalcDisagreeMatrix(PHG35.versions, report = "percent", type = "het")
disagree.het.matrix <- CalcDisagreeMatrix(PHG35.versions, report = "total", type = "het")
disagree.het.matrix[upper.tri(disagree.het.matrix, diag = FALSE)] <- disagree.het.matrix.percent[upper.tri(disagree.het.matrix.percent, diag = FALSE)]

fwrite(data.frame(disagree.het.matrix), "analysis/qc/PHG35_reconstruction/disagreements_hets.txt",
       quote = FALSE, sep = "\t", row.names = TRUE)

# of all that disagree after controlling for phasing, how many were due to hets
disagree.diff.matrix.percent <- CalcDisagreeMatrix(PHG35.versions, report = "percent", type = "diff_allele")
disagree.diff.matrix <- CalcDisagreeMatrix(PHG35.versions, report = "total", type = "diff_allele")
disagree.diff.matrix[upper.tri(disagree.diff.matrix, diag = FALSE)] <- disagree.diff.matrix.percent[upper.tri(disagree.diff.matrix.percent, diag = FALSE)]

fwrite(data.frame(disagree.diff.matrix), "analysis/qc/PHG35_reconstruction/disagreements_diff-allele.txt",
       quote = FALSE, sep = "\t", row.names = TRUE)
