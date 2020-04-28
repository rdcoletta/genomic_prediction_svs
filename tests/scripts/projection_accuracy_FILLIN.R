library(data.table)


ProjectionAccuracy <- function(parents_data, ril_data, crosses_info, output_dir,
                               markers_to_mask = 0.1, only_homo_markers = TRUE,
                               remove_tmp_files = TRUE) {
  
  cat("Projection accuracy:\n\n")
  cat("- FILLIN plugin from TASSEL 5\n")
  cat("- SNPs only\n")
  cat("- Masked", paste0(markers_to_mask*100, "% of markers"), "\n")
  
  # load genotypic datasets
  geno.parents <- fread(parents_data, header = TRUE, data.table = FALSE)
  geno.rils <- fread(ril_data, header = TRUE, data.table = FALSE)
  
  # exclude hets and missing data from parents?
  if (only_homo_markers == TRUE) {
    
    # create empty vector to store homozygous markers for all parents
    markers.to.keep <- vector(mode = "numeric")
    
    # check each marker
    for (row in 1:NROW(geno.parents)) {
      # get marker types
      marker.types <- sapply(X = geno.parents[row, 12:NCOL(geno.parents)],
                             FUN = function(markers) {
                               # check if marker has allele(s) missing
                               if (grepl("N", markers)) {
                                 type <- "missing"
                               } else {
                                 # if not missing, check if it's homo or het
                                 markers <- unlist(strsplit(markers, ""))
                                 if (markers[1] == markers[2]) {
                                   type <- "homo"
                                 } else {
                                   type <- "het"
                                 }
                               }
                               return(type)
                             })
      # keep marker (row) if it is homozygous for all parents
      if (all(marker.types == "homo")) {
        markers.to.keep <- append(markers.to.keep, row)
      }
    }
    
    # just confirm that all markers from parents and rils are the same
    # all(geno.parents[, "rs#"] == geno.rils[, "rs#"])
    
    # filter parents and rils
    geno.parents <- geno.parents[markers.to.keep, ]
    geno.rils <- geno.rils[markers.to.keep, ]
    
    cat("- Removed heterozygous and missing markers among all parents\n")
  }

  # next:
  # - subset by cross
  # - randomly select SNPs to mask
  # - call FILLIN from R using system()
  # - load results and calculate accuracy
  # - go to next cross
  
  # load table with cross information
  df.crosses <- fread(crosses_info, header = TRUE, data.table = FALSE)
  
  # make dir to store results
  out.dir <- output_dir
  if (!dir.exists(out.dir)) {
    dir.create(out.dir, recursive = TRUE)
  }
  
  # create a df to store results for each cross
  df.accuracy.per.cross <- data.frame(matrix(nrow = 0, ncol = 16), stringsAsFactors = FALSE)
  colnames(df.accuracy.per.cross) <- c("cross", "lines_per_cross", "total_markers", "markers_masked", 
                                       "projected_genotypes", "matches", "major_to_minor",
                                       "minor_to_major", "het_to_homo", "homo_to_het", "not_imputed",
                                       "missing_to_missing", "missing_to_major", "missing_to_minor",
                                       "missing_to_het", "predicted")
  
  # calculate projection accuracy by cross
  for (row in 1:NROW(df.crosses)) {
    
    # get cross name, and change * by x in the name
    cross <- df.crosses[row, "cross"]
    cross <- gsub(pattern = "*", replacement = "x", x = cross, fixed = TRUE)
    
    # get RIL names for that cross
    RILs <- unlist(strsplit(df.crosses[row, "RILs"], split = ","))
    
    # subset parents for that cross
    parents <- unlist(strsplit(cross, split = "x"))
    geno.parents.cross <- cbind(geno.parents[, 1:11],
                                geno.parents[, which(colnames(geno.parents) %in% parents)])
    # since there are only two parents,i can't know which allele is major or minor
    geno.parents.cross$alleles <- NA
    
    # subset rils for that cross
    geno.rils.cross <- cbind(geno.rils[, 1:11],
                             geno.rils[, which(colnames(geno.rils) %in% RILs)])
    
    
    # make sure only RILs that were actually planted and genotyped are analyzed
    if (NCOL(geno.rils.cross) > 11) {
      
      cat("\n------------------\n")
      cat(cross, "summary")
      cat("\n------------------\n\n")
      
      # fix major and minor alleles for that cross
      for (i in 1:NROW(geno.rils.cross)) {
        # get all genotypes from ea
        alleles.all <- as.character(geno.rils.cross[i, 12:NCOL(geno.rils.cross)])
        alleles.types <- table(unlist(strsplit(alleles.all, "")))
        # remove Ns
        if ("N" %in% names(alleles.types)) {
          alleles.types <- alleles.types[names(alleles.types) != "N"]
        }
        # sort by allele count
        alleles.types <- sort(alleles.types, decreasing = TRUE)
        # get major allele
        major <- names(alleles.types[1])
        # get minor allele (if present) and correct alleles in the geno.rils.cross
        if (length(alleles.types) > 1) {
          minor <- names(alleles.types[2])
          geno.rils.cross[i, "alleles"] <- paste0(major, "/", minor)
        } else {
          geno.rils.cross[i, "alleles"] <- major
        }
      }
      
      # randomly select SNPs to mask in the RILs
      set.seed(4613)
      markers.to.mask <- sort(sample(x = 1:NROW(geno.rils.cross),
                                     size = round(NROW(geno.rils.cross) * markers_to_mask),  # mask 10% of markers
                                     replace = FALSE))
      geno.rils.masked <- geno.rils.cross
      geno.rils.masked[markers.to.mask, 12:NCOL(geno.rils.masked)] <- "NN"
      # keep track of which markers were masked
      markers.masked <- geno.rils.masked[markers.to.mask, 1]
      
      # write datasets from parents and rils so FILLIN can read it
      filename.parents <- paste0(out.dir, "/usda_SNPs_", cross, "_parents.sorted.hmp.txt")
      filename.rils <- paste0(out.dir, "/usda_SNPs_", cross, "_RILs.sorted.hmp.txt")
      filename.masked.rils <- paste0(out.dir, "/usda_SNPs_", cross, "_RILs.masked.hmp.txt")
      
      fwrite(geno.parents.cross, file = filename.parents, sep = "\t", quote = FALSE, na = "NA")
      fwrite(geno.rils.cross, file = filename.rils, sep = "\t", quote = FALSE, na = "NA")
      fwrite(geno.rils.masked, file = filename.masked.rils, sep = "\t", quote = FALSE, na = "NA")
      
      # write TASSEL 5 commands
      commands.find.hap <- paste0("/Users/rafael/Documents/tassel-5-standalone/run_pipeline.pl",
                                  " -FILLINFindHaplotypesPlugin",
                                  " -hmp ", filename.parents,
                                  " -o ", out.dir, "/donors_", cross,
                                  " -hapSize 1000",
                                  " -minTaxa 1")
      
      commands.impute <- paste0("/Users/rafael/Documents/tassel-5-standalone/run_pipeline.pl",
                                " -FILLINImputationPlugin",
                                " -hmp ", filename.masked.rils,
                                " -d ", out.dir, "/donors_", cross,
                                " -o ", out.dir, "/usda_SNPs_", cross, "_RILs.projected.hmp.txt",
                                " -hapSize 1000",
                                " -accuracy")
      
      commands.hapdip <- paste0("/Users/rafael/Documents/tassel-5-standalone/run_pipeline.pl",
                                " -importGuess ", out.dir, "/usda_SNPs_", cross, "_RILs.projected.hmp.txt",
                                " -export ", out.dir, "/usda_SNPs_", cross, "_RILs.projected.hmp.txt",
                                " -exportType HapmapDiploid")
      
      # create haplotypes from parents
      system(commands.find.hap)
      cat(paste0("Haplotypes created: ", out.dir, "/donors_", cross), "\n")
      # impute masked ril genotypes based on parents 
      system(commands.impute)
      # transform projected rils to hapmap diploid
      system(commands.hapdip)
      cat(paste0("Projection complete: ", out.dir, "/usda_SNPs_", cross, "_RILs.projected.hmp.txt"), "\n")
      
      if (remove_tmp_files == TRUE) {
        # remove temporay files
        system(paste0("rm ", filename.parents))
        system(paste0("rm ", filename.rils))
        system(paste0("rm ", filename.masked.rils))
      }
      
      # load projected RILs
      geno.rils.projected <- fread(paste0(out.dir, "/usda_SNPs_", cross, "_RILs.projected.hmp.txt"),
                                   header = TRUE, data.table = FALSE)
      
      # calculate accuracy of imputation based on how many lines were  correctly genotyped for each
      # masked marker, and also calculate different mistakes from imputation
      
      df.projection.counts <- data.frame(matrix(nrow = length(markers.masked), ncol = 13),
                                         stringsAsFactors = FALSE)
      colnames(df.projection.counts) <- c("marker", "lines", "matches", "major_to_minor",
                                          "minor_to_major", "het_to_homo", "homo_to_het",
                                          "not_imputed", "missing_to_missing", "missing_to_major",
                                          "missing_to_minor", "missing_to_het", "predicted")
      
      for (marker in 1:length(markers.masked)) {
        
        # get expected and observed genotypes
        geno.exp <- geno.rils.cross[marker, 12:NCOL(geno.rils.cross)]
        geno.obs <- geno.rils.projected[marker, 12:NCOL(geno.rils.projected)]
        
        # add name and number of lines genotyped in this population
        df.projection.counts[marker, "marker"] <- markers.masked[marker]
        df.projection.counts[marker, "lines"] <- length(geno.exp)
        
        # # count perfect matches and add to df
        # df.projection.counts[marker, "matches"] <- sum(geno.exp == geno.obs)
        
        # get major and minor alleles
        alleles <- geno.rils.cross[marker, "alleles"]
        major.allele <- unlist(strsplit(alleles, "/"))[1]
        minor.allele <- unlist(strsplit(alleles, "/"))[2]
        # but if all alleles are missing from the masked marker, use alleles from projection
        if (is.na(alleles)) {
          alleles <- geno.rils.projected[marker, "alleles"]
          major.allele <- unlist(strsplit(alleles, "/"))[1]
          minor.allele <- unlist(strsplit(alleles, "/"))[2]
          # and if all alleles are also missing in projection, make sure to transform NA into NN
          if (is.na(alleles)) {
            major.allele <- "N"
            minor.allele <- "N"
          }
        }
        
        # transform genotypes into different allele types
        alleles.exp <- sapply(X = geno.exp,
                              FUN = function(markers, major, minor) {
                                if (grepl("N", markers)) {
                                  type <- "missing"
                                } else {
                                  if (markers == paste0(major, major)) { type <- "major" }
                                  if (markers == paste0(minor, minor)) { type <- "minor" }
                                  if (markers == paste0(major, minor)) { type <- "het" }
                                  if (markers == paste0(minor, major)) { type <- "het" }
                                }
                                return(type)
                              }, major.allele, minor.allele)
        
        alleles.obs <- sapply(X = geno.obs,
                              FUN = function(markers, major, minor) {
                                if (grepl("N", markers)) {
                                  type <- "missing"
                                } else {
                                  if (markers == paste0(major, major)) { type <- "major" }
                                  if (markers == paste0(minor, minor)) { type <- "minor" }
                                  if (markers == paste0(major, minor)) { type <- "het" }
                                  if (markers == paste0(minor, major)) { type <- "het" }
                                  exp.alleles <- c(major, minor)
                                  obs.alleles <- unlist(strsplit(markers, ""))
                                  # if observed allele from projection is not one of the alleles
                                  # in the parents, call them "predicted"
                                  if (!obs.alleles[1] %in% exp.alleles | !obs.alleles[2] %in% exp.alleles) {
                                    type <- "predicted"
                                  }
                                }
                                return(type)
                              }, major.allele, minor.allele)
        
        # count how many perfect matches, how many were not imputed, how many majors were projected
        # as minor alleles, minor to major, het to homozygous, and homo to het.
        
        count.miss2miss <- 0
        count.miss2major <- 0
        count.miss2minor <- 0
        count.miss2het <- 0
        count.not.imp <- 0
        count.matches <- 0
        count.ma2mi <- 0
        count.mi2ma <- 0
        count.het2hom <- 0
        count.homo2het <- 0
        count.pred <- 0
        
        
        for (i in 1:length(alleles.exp)) {
          if (alleles.exp[i] == "missing" & alleles.obs[i] == "missing") {
            count.miss2miss <- count.miss2miss + 1
          } else if (alleles.exp[i] == "missing" & alleles.obs[i] == "major") {
            count.miss2major <- count.miss2major + 1
          } else if (alleles.exp[i] == "missing" & alleles.obs[i] == "minor") {
            count.miss2minor <- count.miss2minor + 1
          } else if (alleles.exp[i] == "missing" & alleles.obs[i] == "het") {
            count.miss2het <- count.miss2het + 1
          } else {
            if (alleles.obs[i] == "missing") {
              count.not.imp <- count.not.imp + 1
            } else if (alleles.obs[i] == "predicted") {
              count.pred <- count.pred + 1
            } else {
              if (alleles.exp[i] == alleles.obs[i]) { count.matches <- count.matches + 1 }
              if (alleles.exp[i] == "major" & alleles.obs[i] == "minor") { count.ma2mi <- count.ma2mi + 1 }
              if (alleles.exp[i] == "minor" & alleles.obs[i] == "major") { count.mi2ma <- count.mi2ma + 1 }
              if (alleles.exp[i] == "het" & alleles.obs[i] != "het") { count.het2hom <- count.het2hom + 1 }
              if (alleles.exp[i] != "het" & alleles.obs[i] == "het") { count.homo2het <- count.homo2het + 1 }
            }
          }
        }
        
        # add results to df
        df.projection.counts[marker, 3:13] <- c(count.matches, count.ma2mi, count.mi2ma, count.het2hom,
                                                count.homo2het, count.not.imp, count.miss2miss,
                                                count.miss2major, count.miss2minor, count.miss2het,
                                                count.pred)
        
      }
      
      
      # check if the counts don't sum to the total number of lines genotyped
      for (j in 1:NROW(df.projection.counts)) {
        if (df.projection.counts[j, "lines"] != sum(df.projection.counts[j, 3:13])) {
          cat(paste0("\nError with marker ", df.projection.counts[j, "marker"], ". Counts don't sum to ",
                     df.projection.counts[j, "lines"], "(the number of lines genotyped)\n\n"))
        }
      }
      
      # add summary to main df
      df.accuracy.per.cross <- rbind(df.accuracy.per.cross,
                                     list(cross = cross,
                                          lines_per_cross = length(12:NCOL(geno.rils.cross)),
                                          total_markers = NROW(geno.rils.cross),
                                          markers_masked = length(markers.masked),
                                          projected_genotypes = sum(df.projection.counts$lines),
                                          matches = sum(df.projection.counts$matches),
                                          major_to_minor = sum(df.projection.counts$major_to_minor),
                                          minor_to_major = sum(df.projection.counts$minor_to_major),
                                          het_to_homo = sum(df.projection.counts$het_to_homo),
                                          homo_to_het = sum(df.projection.counts$homo_to_het),
                                          not_imputed = sum(df.projection.counts$not_imputed),
                                          missing_to_missing = sum(df.projection.counts$missing_to_missing),
                                          missing_to_major = sum(df.projection.counts$missing_to_major),
                                          missing_to_minor = sum(df.projection.counts$missing_to_minor),
                                          missing_to_het = sum(df.projection.counts$missing_to_het),
                                          predicted = sum(df.projection.counts$predicted)),
                                     stringsAsFactors = FALSE)
      
      
      # print accuracy for this cross
      cat("\nMarkers masked:", paste0(round((length(markers.masked) / NROW(geno.rils.cross)) * 100, digits = 0), "%"), "\n\n")
      
      cat("Projection accuracy:", paste0(round((sum(df.projection.counts$matches) / sum(df.projection.counts$lines)) * 100, digits = 2), "%"), "\n")
      cat("Not imputed:", paste0(round((sum(df.projection.counts$not_imputed) / sum(df.projection.counts$lines)) * 100, digits = 2), "%"), "\n")
      cat("Predicted:", paste0(round((sum(df.projection.counts$predicted) / sum(df.projection.counts$lines)) * 100, digits = 2), "%"), "\n")
      cat("Missing in both:", paste0(round((sum(df.projection.counts$missing_to_missing) / sum(df.projection.counts$lines))*100, digits = 2), "%"), "\n")
      cat("Missing to major, minor or het:", paste0(round(((sum(df.projection.counts$missing_to_major) + sum(df.projection.counts$missing_to_minor) + sum(df.projection.counts$missing_to_het)) / sum(df.projection.counts$lines))*100, digits = 2), "%"), "\n")
      cat("Major to minor:", paste0(round((sum(df.projection.counts$major_to_minor) / sum(df.projection.counts$lines))*100, digits = 2), "%"), "\n")
      cat("Minor to major:", paste0(round((sum(df.projection.counts$minor_to_major) / sum(df.projection.counts$lines))*100, digits = 2), "%"), "\n")
      cat("Het to homo:", paste0(round((sum(df.projection.counts$het_to_homo) / sum(df.projection.counts$lines))*100, digits = 2), "%"), "\n")
      cat("Homo to het:", paste0(round((sum(df.projection.counts$homo_to_het) / sum(df.projection.counts$lines))*100, digits = 2), "%"), "\n\n")
      
    }
    
  }
  
  
  # check if counts sum up to the number of projected genotypes
  for (j in 1:NROW(df.accuracy.per.cross)) {
    if(df.accuracy.per.cross[j, "projected_genotypes"] != sum(df.accuracy.per.cross[j, 6:16])) {
      cat(paste0("\nError with marker ", df.accuracy.per.cross[j, "cross"], ". Counts don't sum to ",
                 df.accuracy.per.cross[j, "projected_genotypes"], " (the total number of projected genotypes)\n\n"))
    }
  }
  
  # write final results
  fwrite(df.accuracy.per.cross, file = paste0(out.dir, "/accuracy_projection_all_crosses.txt"),
         sep = "\t", quote = FALSE, na = "NA")
  
  # print overall accuracy
  cat("\n-------------------------\n")
  cat("Summary among all crosses")
  cat("\n-------------------------\n")
  
  attach(df.accuracy.per.cross)
  cat("Markers masked:", paste0(round((sum(markers_masked) / sum(total_markers)) * 100, digits = 0), "%"), "\n")
  
  cat("Projection accuracy:", paste0(round((sum(matches) / sum(projected_genotypes)) * 100, digits = 2), "%"), "\n")
  cat("Not imputed:", paste0(round((sum(not_imputed) / sum(projected_genotypes)) * 100, digits = 2), "%"), "\n")
  cat("Predicted:", paste0(round((sum(predicted) / sum(projected_genotypes)) * 100, digits = 2), "%"), "\n")
  cat("Missing in both:", paste0(round((sum(missing_to_missing) / sum(projected_genotypes)) * 100, digits = 2), "%"), "\n")
  cat("Missing to major, minor or het:", paste0(round(((sum(missing_to_major) + sum(missing_to_minor) + sum(missing_to_het)) / sum(projected_genotypes)) * 100, digits = 2), "%"), "\n")
  cat("Major to minor:", paste0(round((sum(major_to_minor) / sum(projected_genotypes))*100, digits = 2), "%"), "\n")
  cat("Minor to major:", paste0(round((sum(minor_to_major) / sum(projected_genotypes))*100, digits = 2), "%"), "\n")
  cat("Het to homo:", paste0(round((sum(het_to_homo) / sum(projected_genotypes))*100, digits = 2), "%"), "\n")
  cat("Homo to het:", paste0(round((sum(homo_to_het) / sum(projected_genotypes))*100, digits = 2), "%"), "\n\n")
  detach(df.accuracy.per.cross)
  
}





# calculate accuracy of projection using different parameter combinations

for (percent_mask in c(0.01, 0.1, 0.5)) {
  for (homo_only in c(TRUE, FALSE)) {
    if (homo_only == TRUE) {
      output.dir <- paste0("tests/analysis/imputation_reseq-parents_15k/accuracy/only-homo_", percent_mask*100, "-percent-masked")
    } else {
      output.dir <- paste0("tests/analysis/imputation_reseq-parents_15k/accuracy/all-markers_", percent_mask*100, "-percent-masked")
    }
    
    # create directory
    dir.create(output.dir, recursive = TRUE)
    
    # write log of this program
    sink(paste0(output.dir, "/accuracy_projection_log.txt"))
    
    ProjectionAccuracy(parents_data = "data/usda_15kSNPs_7parents.sorted.diploid.v4.all-parents-corrected.hmp.txt",
                       ril_data = "data/usda_15kSNPs_325rils.sorted.diploid.v4.all-parents-corrected.hmp.txt",
                       crosses_info = "tests/data/usda_biparental-crosses.txt",
                       markers_to_mask = percent_mask,
                       output_dir = output.dir,
                       only_homo_markers = homo_only,
                       remove_tmp_files = TRUE)
    
    # close sink connection
    sink()
    
  }
}
closeAllConnections()



# analyze results
library(ggplot2)
library(dplyr)

results.dir <- "tests/analysis/imputation_reseq-parents_15k/accuracy"
results.files <- list.files(path = results.dir,
                           pattern = "accuracy_projection_all_crosses.txt",
                           recursive = TRUE)

master.summary <- data.frame()

for (results in results.files) {
  
  parameters <- unlist(strsplit(results, split = "/"))[1]
  
  summary <- fread(paste0(results.dir, "/", results), header = TRUE, data.table = FALSE)
  summary <- cbind(parameters, summary)
  
  master.summary <- rbind(master.summary, summary)
  
}

# plot summary
master.summary %>% 
  group_by(parameters, cross) %>% 
  summarize(accuracy = matches/projected_genotypes) %>% 
  ggplot(aes(x = cross, y = accuracy, color = parameters)) +
  geom_point() +
  ylim(0,1) +
  labs(x = "cross", y = "Projection accuracy") +
  scale_color_manual(values = c("#fdcc8a", "#fc8d59", "#d7301f", "#bdc9e1", "#74a9cf", "#0570b0"),
                     labels = c("all markers (1% masked)", "all markers (10% masked)",
                                "all markers (50% masked)", "only homo (1% masked)",
                                "only homo (10% masked)", "only homo (50% masked)"),
                     name = "Parameters") +
  scale_x_discrete(labels = function(cross) gsub("x", "\nx\n", cross))

plot.name <- paste0(results.dir, "/accuracy_summary_plot.png")
ggsave(filename = plot.name, device = "png", height = 8, width = 10)
