#### arguments for command line ----

args <- commandArgs(trailingOnly = TRUE)

# help
if (all(length(args) == 1 & args == "-h" | args == "--help")) {
  cat("
Description: this script add back all markers not in present in a cross that were projected in others

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


# proj.folder <- "analysis/projection_svs-snps"




#### libraries ----

if(!require("data.table")) install.packages("data.table")
if(!require("dplyr")) install.packages("dplyr")
if(!require("doParallel")) install.packages("doParallel")

if (detectCores() < 10) {
  num.cores <- detectCores()
} else {
  num.cores <- 10
}


# get list with all families
cross.list <- system("ls data/hapmap_by_cross/usda_22kSNPs_rils.sorted.diploid.filtered.v4.*.hmp.txt | xargs -n 1 basename |  cut -d '.' -f 6 | uniq",
                     intern = TRUE)



#### create empty df with all snps ----

# get a list with all sv positions per chromosome for all crosses
all.proj.markers <- data.frame(stringsAsFactors = FALSE)

cat("Getting list with SV positions for all populations...\n")

for (cross in cross.list) {

  cat("  ", cross, "\n")

  # get parents names
  parent1 <- unlist(strsplit(cross, "x"))[1]
  parent2 <- unlist(strsplit(cross, "x"))[2]

  # load hapmap after projection
  filename.after.proj <- list.files(path = proj.folder,
                                    pattern = paste0("usda_rils_projected-SVs-SNPs.", cross, ".hmp.txt"),
                                    full.names = TRUE)
  hmp.after <- fread(filename.after.proj, header = TRUE, data.table = FALSE)
  # filter hmp by svs and select only chromosome and positions
  hmp.after <- hmp.after[, c(1, 3, 4)]

  for (chr in 1:10) {

    # subset by chr
    hmp.after.chr <- subset(hmp.after, chrom == chr)
    # get marker info and keep only unique info among crosses
    all.proj.markers <- rbind(all.proj.markers, hmp.after.chr)
    all.proj.markers <- distinct(all.proj.markers)

  }

  # sort df by chr and pos
  all.proj.markers <- all.proj.markers[order(all.proj.markers$chrom, all.proj.markers$pos), ]

}

cat("Done!\n\n")


# create an empty hmp file with hmp columns
all.proj.markers <- data.frame(`rs#` = all.proj.markers[, 1],
                               alleles = "N",
                               chrom = all.proj.markers[, 2],
                               pos = all.proj.markers[, 3],
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
                                    pattern = paste0("usda_rils_projected-SVs-SNPs.", cross, ".hmp.txt"),
                                    full.names = TRUE)
  hmp.after <- fread(filename.after.proj, header = TRUE, data.table = FALSE)


  # subset by chr and run in parallel
  registerDoParallel(cores = num.cores)
  all.markers.data <- foreach(chr=1:10, .combine = rbind) %dopar% {

    all.proj.markers.chr <- subset(all.proj.markers, chrom == chr)
    hmp.after.chr <- subset(hmp.after, chrom == chr)

    # cbind NN matrix with column number = number of genotypes in a cross
    all.proj.markers.chr <- cbind(all.proj.markers.chr,
                          data.frame(matrix("NN", nrow = NROW(all.proj.markers.chr),
                                            ncol = NCOL(hmp.after.chr) - 11),
                                     stringsAsFactors = FALSE))
    colnames(all.proj.markers.chr) <- colnames(hmp.after.chr)

    # filter dataframe with NN to have only SNPs that were not projected for that particular cross
    all.proj.markers.chr.NN <- all.proj.markers.chr[which(!all.proj.markers.chr[, 1] %in% hmp.after.chr[, 1]), ]

    # merge missing data
    merged.hmp <- rbind(all.proj.markers.chr.NN, hmp.after.chr)
    merged.hmp <- merged.hmp[order(merged.hmp$pos), ]

    # return merged dataset for chromosome
    merged.hmp

  }
  stopImplicitCluster()

  if (NROW(all.proj.markers) != NROW(all.markers.data)) cat(cross, "has wrong number of rows\n")

  # write results for cross
  out.filename <- gsub(".hmp.txt", ".all-markers.hmp.txt", filename.after.proj, fixed = TRUE)
  fwrite(all.markers.data, file = out.filename, quote = FALSE, sep = "\t", na = "NA", row.names = FALSE)

}

cat("Done!\n\n")
