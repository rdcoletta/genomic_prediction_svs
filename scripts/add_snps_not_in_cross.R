#### arguments for command line ----

args <- commandArgs(trailingOnly = TRUE)

# help
if (all(length(args) == 1 & args == "-h" | args == "--help")) {
  cat("
Description: this script add back all SNPs not in present in a cross that were projected in others

Usage: ")
  quit()
}

# make sure the correct number of arguments are used
# you should provide 2 arguments
if (length(args) != 1) {
  stop("incorrect number of arguments provided.
       
       Usage:
       ")
}

# assign arguments to variables
proj.folder <- args[1]


# proj.folder <- "analysis/projection_reseq-snps"




#### libraries ----

if(!require("data.table")) install.packages("data.table")
if(!require("doParallel")) install.packages("doParallel")

if (detectCores() < 10) {
  num.cores <- detectCores()
} else {
  num.cores <- 10
}


# get list with all NAM families
cross.list <- system("ls -d ~/projects/genomic_prediction/simulation/data/reseq_snps/[BLP]* | xargs -n 1 basename",
                     intern = TRUE)



#### create empty df with all snps ----

# get a list with all sv positions per chromosome for all crosses
all.proj.snp.pos <- list()

cat("Getting list with SV positions for all populations...\n")

for (cross in cross.list) {
  
  cat("  ", cross, "\n")
  
  # get parents names
  parent1 <- unlist(strsplit(cross, "x"))[1]
  parent2 <- unlist(strsplit(cross, "x"))[2]
  
  # load hapmap after projection
  filename.after.proj <- list.files(path = proj.folder,
                                    pattern = paste0(cross, ".projected.sliding-window.hmp.txt"),
                                    full.names = TRUE)
  hmp.after <- fread(filename.after.proj, header = TRUE, data.table = FALSE)
  # filter hmp by svs and select only chromosome and positions
  hmp.after <- hmp.after[, c(1, 3, 4)]
  
  for (chr in 1:10) {
    
    hmp.after.chr <- subset(hmp.after, chrom == chr)
    
    if (!as.character(chr) %in% names(all.proj.snp.pos)) {
      all.proj.snp.pos[[as.character(chr)]] <- hmp.after.chr[, "pos"]
    } else {
      all.proj.snp.pos[[as.character(chr)]] <- append(all.proj.snp.pos[[as.character(chr)]], hmp.after.chr[, "pos"])
      all.proj.snp.pos[[as.character(chr)]] <- sort(unique(all.proj.snp.pos[[as.character(chr)]]))
    }
    
  }
}

cat("Done!\n\n")


# prepare snp name column
vector.snp.names <- lapply(seq_along(all.proj.snp.pos), function(chr) {
  return(paste0("snp.", names(all.proj.snp.pos)[[chr]], ".", unlist(all.proj.snp.pos[[chr]])))
})
vector.snp.names <- unlist(vector.snp.names)

# prepare chrom column
vector.chr <- lapply(seq_along(all.proj.snp.pos), function(chr) {
  return(rep(names(all.proj.snp.pos)[[chr]], times = length(unlist(all.proj.snp.pos[[chr]]))))
})
vector.chr <- unlist(vector.chr)

# create an empty hmp file with hmp columns
all.snps.hmp <- data.frame(`rs#` = vector.snp.names,
                           alleles = "N",
                           chrom = vector.chr,
                           pos = unlist(all.proj.snp.pos),
                           strand = "+",
                           `assembly#` = "NA",
                           center = "NA", 
                           protLSID = "NA",
                           assayLSID = "NA",
                           panelLSID = "NA",
                           QCcode = "NA",
                           check.names = FALSE, stringsAsFactors = FALSE)




#### add back snps not in cross ----


cat("Creating file with SVs only for each cross...\n")

for (cross in cross.list) {
  
  cat("  ", cross, "\n")
  
  # load hapmap after projection
  filename.after.proj <- list.files(path = proj.folder,
                                    pattern = paste0(cross, ".projected.sliding-window.hmp.txt"),
                                    full.names = TRUE)
  hmp.after <- fread(filename.after.proj, header = TRUE, data.table = FALSE)
  
  
  # subset by chr and run in parallel
  registerDoParallel(cores = num.cores)
  all.snps.hmp.data <- foreach(chr=1:10, .combine = rbind) %dopar% {
    
    all.snps.hmp.chr <- subset(all.snps.hmp, chrom == chr)
    hmp.after.chr <- subset(hmp.after, chrom == chr)
    
    # cbind NN matrix with column number = number of genotypes in a cross
    all.snps.hmp.chr <- cbind(all.snps.hmp.chr,
                          data.frame(matrix("NN", nrow = NROW(all.snps.hmp.chr),
                                            ncol = NCOL(hmp.after.chr) - 11),
                                     stringsAsFactors = FALSE))
    colnames(all.snps.hmp.chr) <- colnames(hmp.after.chr)
    
    # filter dataframe with NN to have only SNPs that were not projected for that particular cross
    all.snps.hmp.chr.NN <- all.snps.hmp.chr[which(!all.snps.hmp.chr[, 1] %in% hmp.after.chr[, 1]), ]
    all.snps.hmp.chr.not.NN <- all.snps.hmp.chr[which(all.snps.hmp.chr[, 1] %in% hmp.after.chr[, 1]), ]
    
    # remove any duplicated row
    hmp.after.chr <- hmp.after.chr[!duplicated(hmp.after.chr[, 1]), ]
    
    # make sure the data has the same length and same positions 
    if (all(hmp.after.chr[, 1] == all.snps.hmp.chr.not.NN[, 1])) {
      # merge missing data
      merged.hmp <- rbind(all.snps.hmp.chr.NN, hmp.after.chr)
      merged.hmp <- merged.hmp[order(merged.hmp$pos), ]
    } else {
      stop("Data have different length")
    }
    
    # return merged dataset for chromosome
    merged.hmp
    
  }
  stopImplicitCluster()
  
  NROW(all.snps.hmp) == NROW(all.snps.hmp.data)
  
  # write results for cross
  out.filename <- gsub(".projected.sliding-window.hmp.txt", ".projected.snps-all-crosses.hmp.txt", filename.after.proj, fixed = TRUE)
  fwrite(all.snps.hmp.data, file = out.filename, quote = FALSE, sep = "\t", na = "NA", row.names = FALSE)
  
}

cat("Done!\n\n")
