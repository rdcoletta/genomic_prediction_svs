

crosses.file <- "tests/data/usda_biparental-crosses.txt"
data.per.cross <- "tests/data/merged_hapmaps_by_cross"
output.folder <- "tests/data/reconstructed_PHG35/"





if(!require("data.table")) install.packages("data.table")
if(!require("ggplot2")) install.packages("ggplot2")



CorrectAlleleColumn <- function(genotypes) {
  
  genotypes <- unique(genotypes)
  alleles <- c()
  for (geno in genotypes) {
    allele <- unlist(strsplit(geno, split = ""))
    alleles <- append(alleles,unique(allele))
  }
  alleles <- paste0(unique(alleles), collapse = "/")
  return(alleles)
  
}



CountGenotypes <- function(genotypes, allele1, allele2) {
  
  alleles <- unlist(strsplit(genotypes, split = ""))
  if ("N" %in% alleles) {
    return("missing")
  } else {
    if (alleles[1] == alleles[2] & alleles[1] == allele1) {
      return("homo1")
    } else if (alleles[1] == alleles[2] & alleles[1] == allele2) {
      return("homo2")
    } else {
      return("het")
    }
  }
  
}


# open file with crosses information and select only crosses with PHG35
crosses <- fread(crosses.file, header = TRUE, data.table = FALSE)
crosses <- crosses[grep("PHG35", crosses[, "cross"]), ]
crosses[, "cross"] <- gsub("*", "x", crosses[, "cross"], fixed = TRUE)


for (cross in crosses[, "cross"]) {
  
  # only parse PHG35 crosses that will be used for projection
  donors.list <- list.files(path = data.per.cross, pattern = "parents.sorted.hmp.txt",
                            recursive = TRUE, full.names = FALSE)
  donors.list <- gsub("usda_SNPs-SVs_", "", donors.list)
  donors.list <- gsub("_parents.sorted.hmp.txt", "", donors.list)
  donors.list <- donors.list[grep("PHG35", donors.list)]
  
  if (cross %in% donors.list) {
    
    cat("Analyzing cross", cross)
    
    # load parental data
    parents.hmp <- list.files(path = data.per.cross, pattern = paste0(cross, "_parents.sorted.hmp.txt"),
                             recursive = TRUE, full.names = TRUE)
    parents.hmp <- fread(parents.hmp, header = TRUE, data.table = FALSE)
    
    # load ril data
    rils.hmp <- list.files(path = data.per.cross, pattern = paste0(cross, "_RILs.sorted.hmp.txt"),
                           recursive = TRUE, full.names = TRUE)
    rils.hmp <- fread(rils.hmp, header = TRUE, data.table = FALSE)
    
    # correct allelle columns on parental and ril data
    parents.hmp[, "alleles"] <- apply(parents.hmp[, 12:NCOL(parents.hmp)], MARGIN = 1, FUN = CorrectAlleleColumn)
    rils.hmp[, "alleles"] <- apply(rils.hmp[, 12:NCOL(rils.hmp)], MARGIN = 1, FUN = CorrectAlleleColumn)
    
    # check which SNP calls disagree
    disagree <- c()
    for (row in 1:NROW(parents.hmp)) {
      
      parents.alleles <- unlist(strsplit(parents.hmp[row, "alleles"], split = "/"))
      rils.alleles <- unlist(strsplit(rils.hmp[row, "alleles"], split = "/"))

      # ignore missing data
      rils.alleles <- rils.alleles[grep("N", rils.alleles, invert = TRUE)]
      
      if (length(rils.alleles) > 0 & !all(rils.alleles %in% parents.alleles)) {
        disagree <- append(disagree, row)
      }
      
    }
    
    # for each snp that disagrees:
    #   - see if one of the parents call is missing
    #   - if parents alleles are not missing, see what's the proportion of allele1 and allele2 and how many homo1, het and homo2
    
    
    df.disagree <- data.frame(matrix(nrow = 0, ncol = 9))
    colnames(df.disagree) <- list("SNP", "alleles_parents", "alleles_rils", "missing_in_PHG35_reseq",
                                  "missing_in_other_parent_reseq", "homo1", "het", "homo2", "missing")
    
    parents <- unlist(strsplit(cross, split = "x"))
    not.PHG35.parent <- parents[parents != "PHG35"]
    
    for (snp in disagree) {
      
      if ("NN" %in% parents.hmp[snp, "PHG35"] & "NN" %in% parents.hmp[snp, not.PHG35.parent]) {
        df.disagree <- rbind(df.disagree, list(SNP = parents.hmp[snp, 1],
                                               alleles_parents = parents.hmp[snp, 2],
                                               alleles_rils = rils.hmp[snp, 2],
                                               missing_in_PHG35_reseq = TRUE,
                                               missing_in_other_parent_reseq = TRUE,
                                               homo1 = NA,
                                               het = NA,
                                               homo2 = NA,
                                               missing = NA),
                             stringsAsFactors = FALSE)
        
      } else if ("NN" %in% parents.hmp[snp, "PHG35"]) {
        df.disagree <- rbind(df.disagree, list(SNP = parents.hmp[snp, 1],
                                               alleles_parents = parents.hmp[snp, 2],
                                               alleles_rils = rils.hmp[snp, 2],
                                               missing_in_PHG35_reseq = TRUE,
                                               missing_in_other_parent_reseq = FALSE,
                                               homo1 = NA,
                                               het = NA,
                                               homo2 = NA,
                                               missing = NA),
                             stringsAsFactors = FALSE)
        
      } else if ("NN" %in% parents.hmp[snp, not.PHG35.parent]) {
        df.disagree <- rbind(df.disagree, list(SNP = parents.hmp[snp, 1],
                                               alleles_parents = parents.hmp[snp, 2],
                                               alleles_rils = rils.hmp[snp, 2],
                                               missing_in_PHG35_reseq = FALSE,
                                               missing_in_other_parent_reseq = TRUE,
                                               homo1 = NA,
                                               het = NA,
                                               homo2 = NA,
                                               missing = NA),
                             stringsAsFactors = FALSE)
        
      } else {
        
        ril.alleles <- unlist(strsplit(rils.hmp[snp, 2], split = "/"))
        ril.alleles <- ril.alleles[ril.alleles != "N"]
        allele1 <- ril.alleles[1]
        allele2 <- ril.alleles[2]
        
        # count genotypes
        geno.count <- sapply(rils.hmp[snp, 12:NCOL(rils.hmp)], FUN = CountGenotypes, allele1, allele2)
        # tally up genotypes
        homo1 <- round(sum(geno.count == "homo1") / length(geno.count), digits = 2)
        het <- round(sum(geno.count == "het") / length(geno.count), digits = 2)
        homo2 <- round(sum(geno.count == "homo2") / length(geno.count), digits = 2)
        missing <- round(sum(geno.count == "missing") / length(geno.count), digits = 2)
        
        df.disagree <- rbind(df.disagree, list(SNP = parents.hmp[snp, 1],
                                               alleles_parents = parents.hmp[snp, 2],
                                               alleles_rils = rils.hmp[snp, 2],
                                               missing_in_PHG35_reseq = FALSE,
                                               missing_in_other_parent_reseq = FALSE,
                                               homo1 = homo1,
                                               het = het,
                                               homo2 = homo2,
                                               missing = missing),
                             stringsAsFactors = FALSE)
        
      }
    }
    
    total.missing.PHG35 <- sum(df.disagree[, "missing_in_PHG35_reseq"] == TRUE)
    percent.missing.PHG35 <- (total.missing.PHG35 / NROW(df.disagree)) * 100
    percent.missing.PHG35 <- round(percent.missing.PHG35, digits = 2)
    
    # cat("  ", total.missing.PHG35, " (", percent.missing.PHG35, "%) of disagreements are due to missing data on resequencing of PHG35 parent\n", sep = "")
    cat(" (", NROW(df.disagree), " disagreements)\n", sep = "")
    
    
    filter.alleles.both.parents <- sapply(df.disagree[, "alleles_parents"], FUN = function(alleles) {
      alleles <- unlist(strsplit(alleles, split = "/"))
      return(length(alleles))
    })
    
    
    # missing data in both parents
    summary.missing.both <- NROW(df.disagree[which(filter.alleles.both.parents == 1 & df.disagree[, "missing_in_PHG35_reseq"] == TRUE), ])
    cat("  ", summary.missing.both, " (", round((summary.missing.both / NROW(df.disagree)) * 100, digits = 2),
        "%) of SNPs that disagree between parents and RILs are due to missing data both parents\n", sep = "")
    
    # missing data only in PHG35 parent
    summary.missing.PHG35 <- NROW(df.disagree[which(filter.alleles.both.parents > 1 & df.disagree[, "missing_in_PHG35_reseq"] == TRUE), ])
    cat("  ", summary.missing.PHG35, " (", round((summary.missing.PHG35 / NROW(df.disagree)) * 100, digits = 2),
        "%) of SNPs that disagree between parents and RILs are due to missing data in PHG35 parent\n", sep = "")
    
    # missing data only in other parent
    summary.missing.notPHG35 <- NROW(df.disagree[which(filter.alleles.both.parents > 1 & df.disagree[, "missing_in_PHG35_reseq"] == FALSE & df.disagree[, "missing_in_other_parent_reseq"] == TRUE), ])
    
    cat("  ", summary.missing.notPHG35, " (", round((summary.missing.notPHG35 / NROW(df.disagree)) * 100, digits = 2),
        "%) of SNPs that disagree between parents and RILs are due to missing data in the non-PHG35 parent\n", sep = "")
    
    # what if I have monomorphic call for a SNP in RIL data mixed with missing data, and the parents
    # data have the allele called in RILs? I can't know if the SNP is actually monomorphic, or if
    # the missing data would actually be a different allele.
    summary.undefined <- NROW(df.disagree[which(filter.alleles.both.parents == 1 & df.disagree[, "missing_in_PHG35_reseq"] == FALSE), ])
    # View(df.disagree[which(filter.alleles.both.parents == 1 & df.disagree[, "missing_in_PHG35_reseq"] == FALSE), ])
    
    cat("  ", summary.undefined, " (", round((summary.undefined / NROW(df.disagree)) * 100, digits = 2),
        "%) of SNPs that disagree between parents and RILs are due to uncertainty in SNP being mono or polymorphic\n", sep = "")
    
    # actual disagreement 
    summary.actual.disagree <- NROW(df.disagree[which(filter.alleles.both.parents > 1 & df.disagree[, "missing_in_PHG35_reseq"] == FALSE & df.disagree[, "missing_in_other_parent_reseq"] == FALSE), ])
    # View(df.disagree[which(filter.alleles.both.parents > 1 & df.disagree[, "missing_in_PHG35_reseq"] == FALSE & df.disagree[, "missing_in_other_parent_reseq"] == FALSE), ])
    
    cat("  ", summary.actual.disagree, " (", round((summary.actual.disagree / NROW(df.disagree)) * 100, digits = 2),
        "%) of SNPs that disagree between parents and RILs are due to PHG35 allele being different between resequencing and SNP data\n\n", sep = "")
    
    
    summary.df <- data.frame(reason_disagreement = factor(c("missing both\nparents", "missing PHG35", "missing other\nparent",
                                       "uncertainty", "actual disagree"), levels = (c("missing both\nparents", "missing PHG35", "missing other\nparent",
                                                                                      "uncertainty", "actual disagree"))),
                             counts = c(summary.missing.both, summary.missing.PHG35, summary.missing.notPHG35,
                                        summary.undefined, summary.actual.disagree))
    
    
    if (!dir.exists("tests/analysis/PHG35_disagreements")) {
      dir.create("tests/analysis/PHG35_disagreements")
    }
    
    summary.plot <- ggplot(summary.df, aes(x = reason_disagreement, y = counts)) +
      geom_col()
    ggsave(filename = paste0("tests/analysis/PHG35_disagreements/", cross, "_summary.png"),
           plot = summary.plot, device = "png")
    
  }
  
}
