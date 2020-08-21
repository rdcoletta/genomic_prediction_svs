#### arguments for command line ----

args <- commandArgs(trailingOnly = TRUE)

# help
if (all(length(args) == 1 & args == "-h" | args == "--help")) {
  cat("
Description: this script filters Plink LD files to have only SNPs in LD with SVs or SNPs not in LD with SV,
             and also plots the distribution of R2 between SVs and SNPs.

Usage: Rscript distribution_ld_snps-svs.R [...]")
  quit()
}

# make sure the correct number of arguments are used
# you should provide 3 arguments
if (length(args) != 3) {
  stop("incorrect number of arguments provided.

Usage: Rscript distribution_ld_snps-svs.R [...]
       ")
}


# assign arguments to variables
plink.file <- args[1]
sv.file <- args[2]
out.dir.ld <- args[3]

if (!dir.exists(out.dir.ld)) dir.create(out.dir.ld, recursive = TRUE)

# plink.file <- "analysis/ld/ld_usda_rils_snp-sv_only.chr10.window-100kb.filter-0.25.ld"
# sv.file <- "data/usda_SVs_parents.sorted.hmp.txt"
# out.dir.ld <- "analysis/ld"



#### libraries ----

if(!require("data.table")) install.packages("data.table")
if(!require("ggplot2")) install.packages("ggplot2")
if(!require("doParallel")) install.packages("doParalell")


if (detectCores() > 10) {
  num.cores <- 10
} else {
  num.cores <- detectCores()
}

#### subsample ----


add.info.to.ld.df <- function(LD_results) {
  
  # add sv size
  LD_results$sv_size <- apply(LD_results, MARGIN = 1, FUN = function(info) {
    sv.index <- grep("^del|^dup|^ins|^inv|^tra", info[c("SNP_A", "SNP_B")], perl = TRUE)
    sv <- info[c("SNP_A", "SNP_B")][sv.index]
    sv <- unlist(strsplit(sv, split = ".", fixed = TRUE))
    sv.size <- as.numeric(sv[4]) - as.numeric(sv[3])
    return(sv.size)
  })
  
  # add R2 quarter
  LD_results$ld_quarter <- apply(LD_results, MARGIN = 1, FUN = function(info) {
    ld <- as.numeric(info["R2"])
    if (ld >= 0 & ld <= 0.25) return("0_to_0.25")
    if (ld > 0.25 & ld <= 0.5) return("0.25_to_0.5")
    if (ld > 0.5 & ld <= 0.75) return("0.5_to_0.75")
    if (ld > 0.75 & ld <= 1.1) return("0.75_to_1")
  })
  
  # add range of sv sizes
  LD_results$size_range <- apply(LD_results, MARGIN = 1, FUN = function(info) {
    if (all(!is.na(info))) {
      sv_size <- as.numeric(info["sv_size"])
      if (sv_size >= 0 & sv_size <= 10000) return("<10kb")
      if (sv_size > 10000 & sv_size <= 100000) return("10kb-100kb")
      if (sv_size > 100000 & sv_size <= 1000000) return("100kb-1Mb")
      if (sv_size > 1000000) return(">1Mb")
    } else {
      return("NA")
    }
  })
  
  return(LD_results)
}


# load sv info
SVs <- fread(sv.file, header = TRUE, data.table = FALSE)
SVs <- SVs[, c(1, 3, 4)]


for (chr in 1:10) {

  # subset SVs based on chr
  SVs.chr <- subset(SVs, chrom == chr)[1]
  SVs.chr <- as.character(SVs.chr[, 1])

  # load one chr at a time
  plink.file.chr <- gsub("chr[0-9]+", paste0("chr", chr), plink.file, perl = TRUE)
  # open table with LD among markers
  LD_results <- fread(plink.file.chr, header = TRUE, data.table = FALSE)

  # add column with distance between snp and closest SV boundary
  # LD_results$dist_to_sv2 <- abs(LD_results[, "BP_B"] - LD_results[, "BP_A"])
  LD_results$dist_to_sv <- apply(LD_results, MARGIN = 1, function(marker) {

    marker1 <- marker["SNP_A"]
    marker2 <- marker["SNP_B"]
    pos1 <- marker["BP_A"]
    pos2 <- marker["BP_B"]

    if (grepl("^del|^dup|^ins|^inv|^tra", marker1, perl = TRUE)) {
      sv.start <- as.numeric(unlist(strsplit(as.character(marker1), split = ".", fixed = TRUE))[3])
      sv.end <- as.numeric(unlist(strsplit(as.character(marker1), split = ".", fixed = TRUE))[4])
      snp.pos <- as.numeric(pos2)
    } else {
      sv.start <- as.numeric(unlist(strsplit(as.character(marker2), split = ".", fixed = TRUE))[3])
      sv.end <- as.numeric(unlist(strsplit(as.character(marker2), split = ".", fixed = TRUE))[4])
      snp.pos <- as.numeric(pos1)
    }

    snp.dist.start <- abs(sv.start - snp.pos)
    snp.dist.end <- abs(sv.end - snp.pos)

    # get closest distance
    if (snp.dist.start < snp.dist.end) {
      dist.to.sv <- snp.dist.start
    } else {
      dist.to.sv <- snp.dist.end
    }

    return(dist.to.sv)
  })


  cat("subsetting only SNPs with highest LD to SV\n")


  # create empty dataset
  LD_results_highest <- data.frame(matrix(nrow = 0, ncol = NCOL(LD_results)), stringsAsFactors = FALSE)
  colnames(LD_results_highest) <- colnames(LD_results)

  # create vector to make sure the same snp doesn't get picked twice
  snps.already.in.ld <- c()

  # get closest (highest LD) snps
  for (sv in SVs.chr) {

    # subset LD results to have only the SV being parsed
    snps_LD_with_sv <- LD_results[which(LD_results[, "SNP_A"] == sv | LD_results[, "SNP_B"] == sv ), ]

    if (length(snps.already.in.ld) > 0) {
      snps_LD_with_sv <- snps_LD_with_sv[which(!snps_LD_with_sv[, "SNP_A"] %in% snps.already.in.ld & !snps_LD_with_sv[, "SNP_B"] %in% snps.already.in.ld), ]
    }


    if (NROW(snps_LD_with_sv) > 0) {

      # select only SNP with highest LD with that SV
      snps_LD_with_sv <- snps_LD_with_sv[which(snps_LD_with_sv[, "R2"] == max(snps_LD_with_sv[, "R2"])), ]
      # if there are more than one SNP with the same R2, get the closest one to the SV
      if (NROW(snps_LD_with_sv) > 1) {
        snps_LD_with_sv <- snps_LD_with_sv[which(snps_LD_with_sv[, "dist_to_sv"] == min(snps_LD_with_sv[, "dist_to_sv"])), ]
      }
      # get closest snp
      snp.selected <- apply(snps_LD_with_sv, MARGIN = 1, function(row) {
        marker1 <- row["SNP_A"]
        marker2 <- row["SNP_B"]
        if (!grepl("^del|^dup|^ins|^inv|^tra", marker1, perl = TRUE)) {
          return(marker1)
        } else {
          return(marker2)
        }
      })
      snps.already.in.ld <- append(snps.already.in.ld, as.character(snp.selected))
      # add closest SNP in LD with SV into new df
      LD_results_highest <- rbind(LD_results_highest, snps_LD_with_sv)

    }
  }

  # check for duplicates
  LD_results_highest <- LD_results_highest[!duplicated(LD_results_highest[, c("SNP_A", "SNP_B")]), ]

  # add info
  LD_results_highest <- add.info.to.ld.df(LD_results_highest)

  # write filtered ld table
  outfile.highest <- gsub(".ld", ".highest-ld.ld", plink.file.chr, fixed = TRUE)
  outfile.highest <- rev(unlist(strsplit(outfile.highest, split = "/")))[1]
  fwrite(LD_results_highest, paste0(out.dir.ld, "/", outfile.highest), sep = "\t", quote = FALSE, row.names = FALSE, na = NA)


  cat("subsetting only SNPs not in LD to SV\n")


  # first get names of SNPs that are in LD with an SV (R2>0.8)
  snps.in.ld <- mclapply(1:NROW(LD_results), function(row, LD_results) {

    marker1 <- LD_results[row, "SNP_A"]
    marker2 <- LD_results[row, "SNP_B"]
    r2 <- as.numeric(LD_results[row, "R2"])

    if (!grepl("^del|^dup|^ins|^inv|^tra", marker1, perl = TRUE)) {
      snp <- marker1
    } else {
      snp <- marker2
    }

    if (r2 >= 0.8) return(snp)

  }, LD_results, mc.cores = num.cores)

  snps.in.ld.vector <- do.call(c, snps.in.ld)
  snps.in.ld.vector <- snps.in.ld.vector[!duplicated(snps.in.ld.vector)]

  # exclude such SNPs
  snps.to.keep <- which(!LD_results[, "SNP_A"] %in% snps.in.ld.vector & !LD_results[, "SNP_B"] %in% snps.in.ld.vector)
  LD_results_lowest <- LD_results[snps.to.keep, ]

  # remove duplicates
  LD_results_lowest <- LD_results_lowest[!duplicated(LD_results_lowest[, c("SNP_A", "SNP_B")]), ]

  # add info
  LD_results_lowest <- add.info.to.ld.df(LD_results_lowest)

  # write filtered ld table
  outfile.lowest <- gsub(".ld", ".not-in-ld.ld", plink.file.chr, fixed = TRUE)
  outfile.lowest <- rev(unlist(strsplit(outfile.lowest, split = "/")))[1]
  fwrite(LD_results_lowest, paste0(out.dir.ld, "/", outfile.lowest), sep = "\t", quote = FALSE, row.names = FALSE, na = NA)


}



#### plot distribution ----

plot_R2_dist <- function(ld_files_to_plot) {
  
  LD_results_all <- data.frame()
  for (chr in 1:10) {
    ld_file_chr <- ld_files_to_plot[grep(paste0("chr", chr, "."), ld_files_to_plot, fixed = TRUE)]
    ld_file_chr <- fread(ld_file_chr, header = TRUE, data.table = FALSE)
    LD_results_all <- rbind(LD_results_all, ld_file_chr)
  }
  
  # distribution of r2 of SNPs in LD with SVs
  plot_dist <- ggplot(LD_results_all, aes(x = R2)) +
    geom_histogram(fill = "#900721", binwidth = 0.01) +
    labs(title = paste0("LD between SVs and SNPs"),
         x = bquote("LD"~(r^2)),
         y = "Count") +
    coord_cartesian(xlim = c(0, 1)) +
    theme(title = element_text(size = 15),
          axis.title = element_text(size = 20),
          axis.text = element_text(size = 15))
  
  # summarize results
  LD_results_all$ld_quarter <- factor(LD_results_all$ld_quarter,
                                      levels = c("0_to_0.25", "0.25_to_0.5", "0.5_to_0.75", "0.75_to_1"))
  
  summary_ld <- aggregate(formula = cbind(sv_size, dist_to_sv) ~ ld_quarter, data = LD_results_all,
                          FUN = function(x) c(count = NROW(x), mean = mean(x), median = median(x)))
  summary_ld <- do.call(data.frame, summary_ld)
  summary_ld <- summary_ld[, c(1:4,6:7)]
  colnames(summary_ld) <- c("ld_quarter", "sv_count", "mean_sv_size", "median_sv_size", "mean_dist_to_sv", "median_dist_to_sv")
  
  # add counts per sv size range and ld quarter
  size_range_summary <- aggregate(LD_results_all$sv_size,
                                  by = list(LD_results_all$ld_quarter, LD_results_all$size_range),
                                  FUN = function(x) NROW(x))
  size_range_summary <- reshape(size_range_summary, v.names = "x", idvar = "Group.1", timevar = "Group.2", direction = "wide")
  colnames(size_range_summary) <- gsub("x.", "sv_", colnames(size_range_summary), fixed = TRUE)
  # match columns to summary df
  size_range_summary <- size_range_summary[order(match(size_range_summary$Group.1, summary_ld$ld_quarter)), ]
  
  # write summary table
  summary_ld <- cbind(summary_ld, size_range_summary[, -1])
  
  return(list(plot_dist, summary_ld))
  
}


# highest ld
files_highest_ld <- list.files(path = out.dir.ld, pattern = "highest-ld.ld", full.names = TRUE)
plot_highest_ld <- plot_R2_dist(files_highest_ld)
ggsave(plot_highest_ld[[1]], filename = paste0(out.dir.ld, "/dist-highest-LD_SNPs-SVs.png"), device = "png")
fwrite(plot_highest_ld[[2]], paste0(out.dir.ld, "/summary-highest-LD_SNPs-SVs.txt"), sep = "\t", quote = FALSE, row.names = FALSE, na = NA)


# not in ld
files_not_in_ld <- list.files(path = out.dir.ld, pattern = "not-in-ld.ld", full.names = TRUE)
plot_not_in_ld <- plot_R2_dist(files_not_in_ld)
ggsave(plot_not_in_ld[[1]], filename = paste0(out.dir.ld, "/dist-not-in-LD_SNPs-SVs.png"), device = "png")
fwrite(plot_not_in_ld[[2]], paste0(out.dir.ld, "/summary-not-in-LD_SNPs-SVs.txt"), sep = "\t", quote = FALSE, row.names = FALSE, na = NA)



