#### arguments for command line ----

args <- commandArgs(trailingOnly = TRUE)

# help
if (all(length(args) == 1 & args == "-h" | args == "--help")) {
  cat("
      Description: This script creates hybrid genotypes based on RIL genotypes used for the cross.

      Usage: Rscript create_hybrid_genotypes.R [...]")
  quit()
}

# make sure the correct number of arguments are used
# you should provide 6 arguments
if (length(args) != 3) {
  stop("incorrect number of arguments provided.

       Usage: Rscript create_hybrid_genotypes.R [pheno_hybrids_excel] [hmp_rils_filename] [outfile]
       ")
}

pheno_hybrids_excel <- args[1]
hmp_rils_filename <- args[2]
outfile <- args[3]



#### libraries ----

if(!require("data.table")) install.packages("data.table")
if(!require("readxl")) install.packages("readxl")
if(!require("gtools")) install.packages("gtools")


#### get hybrid information ----

# load pheno data
hybrid_info <- read_excel(pheno_hybrids_excel)

# keep only pedigree information
hybrid_info <- hybrid_info[, c("Hybrid", "ParentA", "ParentB")]
hybrid_info <- hybrid_info[!duplicated(hybrid_info$Hybrid), ]
# remove checks
hybrid_info <- hybrid_info[grep("check", hybrid_info$Hybrid, invert = TRUE, ignore.case = TRUE), ]
# adjust format
hybrid_info <- as.data.frame(hybrid_info[mixedorder(hybrid_info$Hybrid), ])



#### create genotypes ----

# load ril data
hmp_rils <- fread(hmp_rils_filename, header = TRUE, data.table = FALSE)

# create empty hmp to store hybrid genotypes
hmp_hybrids <- hmp_rils[, 1:11]

cat("Creating hybrid genotypes:\n")

for (hybrid in hybrid_info$Hybrid) {
  
  cat("  hybrid ", hybrid, "\n")
  
  # get RIL names used in for single cross name
  p1 <- hybrid_info[which(hybrid_info$Hybrid == hybrid), "ParentA"]
  p2 <- hybrid_info[which(hybrid_info$Hybrid == hybrid), "ParentB"]
  
  # make sure there is genotypic data for both parents
  if (p1 %in% colnames(hmp_rils) & p2 %in% colnames(hmp_rils)) {
    
    cat("  ", p1, " x ", p2, "\n", sep = "")
    
    # check parental marker type
    marker_type <- apply(X = hmp_rils[, c(p1, p2)], MARGIN = 1, FUN = function(marker) {
      
      # get unique genotypes between parents
      genotypes <- unique(marker)
      
      if (any(grepl("NN", genotypes))) {
        
        # if NN, marker is missing
        return("missing")
        
      } else if (length(genotypes) == 1) {
        
        # if there is one genotype, it's monomorphic
        # but distinguish if marker is het
        alleles <- unlist(strsplit(genotypes, split = ""))
        if (alleles[1] == alleles[2]) {
          return("mono")
        } else {
          return("het")
        }
        
      } else {
        
        # if there are two genotypes, it's polymorphic
        # but distiguish if one of the genotypes is het
        p1_alleles <- unlist(strsplit(genotypes[1], split = ""))
        p2_alleles <- unlist(strsplit(genotypes[2], split = ""))
        if (p1_alleles[1] == p1_alleles[2] & p2_alleles[1] == p2_alleles[2]) {
          return("poly")
        } else {
          return("het")
        }
        
      }
    })
    
    # add marker type to single cross info
    geno_parents <- cbind(hmp_rils[, c(p1, p2)], marker_type)
    
    # create hybrid genotype
    geno_hybrid <- apply(geno_parents, MARGIN = 1, function(marker) {
      
      if (marker[3] == "mono" | marker[3] == "poly") {
        # for markers that are poly or monomorphic, get one allele from each parent
        p1_allele <- unlist(strsplit(marker[1], split = ""))[1]
        p2_allele <- unlist(strsplit(marker[2], split = ""))[1]
        geno_hybrid <- paste0(p1_allele, p2_allele)
      } else {
        # for markers that are missing or het in at least one parent, set marker to NN
        geno_hybrid <- "NN"
      }
      
      return(geno_hybrid)
      
    })
    
    # View(cbind(geno_parents, geno_hybrid))
    
    # append hybrid genotype to final hmp
    hmp_hybrids <- cbind(hmp_hybrids, geno_hybrid)
    # correct column name for that hybrid
    colnames(hmp_hybrids)[NCOL(hmp_hybrids)] <- hybrid
    
  }
  else {
    # debug
    cat("  missing genotypic data in at least one of the parents\n")
  }
  
}

# write file
fwrite(hmp_hybrids, outfile, quote = FALSE, sep = "\t", na = "NA", row.names = FALSE)

# fix allele columns using TASSEL
commands.hapdip <- paste0("/home/hirschc1/della028/software/tassel-5-standalone/run_pipeline.pl",
                          " -importGuess ", outfile,
                          " -export ", outfile,
                          " -exportType HapmapDiploid")
system(commands.hapdip)




##### debug ----

# pheno_hybrids_excel <- "data/NIFA2020_CombinedData.xlsx"
# hmp_rils_filename <- "../single_cross_prediction/data/usda_22kSNPs_rils.sorted.diploid.filtered.v4.sliding-window.hmp.txt"
# outfile <- "data/usda_22kmarkers_hybrids.hmp.txt"
