library(data.table)
library(ggplot2)

usage <- function() {
  cat("
description: plot how many SNPs and SVs are missing per marker and per RIL, and also
             plot the distribution of these markers along the chromosome

usage: Rscript qc_snp-sv_hmp.R [infile_name] [sv_list] [out_folder]

positional arguments:
  infile_name             hapmap file containing SNPs and SVs
  sv_list                 single-column file containing only SV IDs
  out_folder              path to folder to save plots

optional argument:
  --help                  show this helpful message

"
  )
}

missingPerMarker <- function(geno_data = NULL,
                             list_of_SV_IDs = NULL,
                             marker_type = c("SNPs", "SVs")) {
  
  # make sure list of SVs is provided
  if (is.null(list_of_SV_IDs)) stop("No list of SV IDs provided!")
  
  # separate data to have only SNPs or SVs
  if (marker_type == "SNPs") geno_data <- geno_data[which(!geno_data[,1] %in% list_of_SV_IDs), ]
  if (marker_type == "SVs") geno_data <- geno_data[which(geno_data[,1] %in% list_of_SV_IDs), ]
  
  list_missing <- list()
  for (margin in c(1, 2)) {
    
    # count missing data per marker (MARGIN = 1) or per RIL (MARGIN = 2)
    amount_missing <- apply(geno_data[, 12:NCOL(geno_data)], MARGIN = margin, function(genotypes) {
      return(sum(genotypes == "NN"))
    })
    # amount_missing <- amount_missing[amount_missing > 0]
    if (margin == 1) list_missing[["per_marker"]] <- amount_missing
    if (margin == 2) list_missing[["per_ril"]] <- amount_missing
    
  }
  
  return(list_missing)
  
}


getMarkersPositions <- function(geno_data = NULL,
                                list_of_SV_IDs = NULL) {
  
  # make sure list of SVs is provided
  if (is.null(list_of_SV_IDs)) stop("No list of SV IDs provided!")
  
  # create list of SVs
  SV_pos <- geno_data[which(geno_data[, 1] %in% list_of_SV_IDs), c(1, 3, 4)]
  SV_pos[, 1] <- as.character(sapply(SV_pos[, 1], function(x) {
    x <- unlist(strsplit(x, split = ".", fixed = TRUE))
    x <- toupper(x)
    return(x[1])
  }))
  SV_pos <- cbind(SV_pos, marker = "SV")
  
  # create list of SNPs
  SNP_pos <- geno_data[which(!geno_data[,1] %in% list_of_SV_IDs), c(1, 3, 4)]
  SNP_pos[, 1] <- "SNP"
  SNP_pos <- cbind(SNP_pos, marker = "SNP")
  
  markers_pos <- rbind(SNP_pos, SV_pos)
  colnames(markers_pos)[1] <- "specific_marker_type"
  
  return(markers_pos)
  
}



#### command line options ----

args <- commandArgs(trailingOnly = TRUE)
if ("--help" %in% args) usage() & q(save = "no")

# assert to have the correct optional arguments
if (length(args) < 3) stop(usage(), "missing positional argument(s)")

# get arguments
infile_name <- args[1]
sv_list <- args[2]
out_folder <- args[3]
# infile_name <- "data/usda_rils_projected-SVs-SNPs.chr1.poly.hmp.txt"
# sv_list <- "data/SVs_IDs.txt"
# out_folder <- "analysis/trait_sim/qc"



#### data qc ----

# load SV information
SVs <- fread(sv_list, header = FALSE, data.table = FALSE)
SVs <- SVs[, 1]

# analyze snps and svs separately
for (marker in c("SNPs", "SVs")) {
  
  # create empty list to store results
  missing <- list()
  
  for (chr in 1:10) {
    
    # load genotypic data from one chromosome at a time
    infile_name_chr <- gsub("chr[0-9]+", paste0("chr", chr), infile_name)
    geno_data <- fread(infile_name_chr, header = TRUE, data.table = FALSE)
    
    # calculate how many markers are missing per ril and per marker
    missing_chr <- missingPerMarker(geno_data = geno_data,
                                    list_of_SV_IDs = SVs,
                                    marker_type = marker)
    
    # for missing data per marker, add them to list with all chromosomes
    missing[["per_marker"]] <- append(missing[["per_marker"]], missing_chr[["per_marker"]])
    # for missing data per ril, sum their numbers with all chromosomes
    if (length(missing[["per_ril"]]) == 0) missing[["per_ril"]] <- rep(0, times = length(missing_chr[["per_ril"]]))
    missing[["per_ril"]] <- missing[["per_ril"]] + missing_chr[["per_ril"]]
    
  }
  
  # total number of missing markers / (number of RILs x number of markers)
  total_missing <- sum(missing[["per_marker"]]) / (length(missing[["per_ril"]]) * length(missing[["per_marker"]]))
  total_missing <- round(total_missing * 100, digits = 2)
  
  cat("Total amount of ", gsub("s", "", marker), " calls missing: ", total_missing, "%\n", sep = "")
  
  # transform missing data into percentage
  missing[["per_marker"]] <- as.numeric(round(missing[["per_marker"]] / length(missing[["per_ril"]]), digits = 2))
  missing[["per_ril"]] <- as.numeric(round(missing[["per_ril"]] / length(missing[["per_marker"]]), digits = 2))
  
  for (missing_type in c("per_marker", "per_ril")) {
    
    # create a df for ploting
    amount_missing <- data.frame(missing = as.numeric(missing[[missing_type]]))
    # plot
    plot_missing <- ggplot(data = amount_missing, aes(x = missing)) +
      geom_histogram(binwidth = 0.01) +
      coord_cartesian(xlim = c(0, 1)) +
      ylab(paste0("Number of ",
                  ifelse(missing_type == "per_marker", marker, "RILs"))) +
      xlab(paste0("Percent missing data",
                  ifelse(missing_type == "per_marker", " (average of all RILs)", paste0(" (", marker, ")"))))
    # save results
    if (!dir.exists(out_folder)) dir.create(out_folder)
    out_file <- paste0(out_folder, "/", marker, "_missing_", missing_type,".pdf")
    ggsave(filename = out_file, plot = plot_missing, device = "pdf")

  }
  
}

# get positions for markers on each chr
markers_pos <- data.frame()
for (chr in 1:10) {
  
  # load genotypic data from one chromosome at a time
  infile_name_chr <- gsub("chr[0-9]+", paste0("chr", chr), infile_name)
  geno_data <- fread(infile_name_chr, header = TRUE, data.table = FALSE)
  
  # get positions
  markers_pos_chr <- getMarkersPositions(geno_data = geno_data, list_of_SV_IDs = SVs)
  markers_pos <- rbind(markers_pos, markers_pos_chr)
}

# plot distribution of markers along each chomosome
plot_distribution <- ggplot(markers_pos, aes(x = pos, fill = marker)) + 
  geom_density(alpha = 0.5, position = "identity") +
  facet_wrap(~chrom, scales = "free_x") +
  scale_x_continuous(labels = function(x) x/1000000) +
  labs(x = "Position (Mb)",
       y = "Density",
       fill = "Marker type")

ggsave(filename = paste0(out_folder, "/SNPs_and_SVs_distribution_per_chr.pdf"),
       plot = plot_distribution, device = "pdf")