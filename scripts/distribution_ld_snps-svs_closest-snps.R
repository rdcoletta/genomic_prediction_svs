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
# you should provide 5 arguments
if (length(args) != 5) {
  stop("incorrect number of arguments provided.

Usage: Rscript distribution_ld_snps-svs_closest-snps.R [...]
       ")
}


# assign arguments to variables
plink.file <- args[1]
sv.file <- args[2]
snp.file <- args[3]
out.dir.ld <- args[4]
sv.type <- args[5]

if (!dir.exists(out.dir.ld)) dir.create(out.dir.ld, recursive = TRUE)

# plink.file <- "analysis/ld/ld_usda_rils_snp-sv_only.chr1.window-1kb.filter-0.25.ld"
# sv.file <- "data/usda_SVs_parents.sorted.hmp.txt"
# snp.file <- "analysis/ld/usda_rils_projected-SVs-SNPs.closest-snps-to-svs.filter-0.25.all.hmp.txt"
# out.dir.ld <- "analysis/ld/window-1kb_filter-0.25_closest-snps"
# sv.type <- "all"



#### libraries ----

if(!require("data.table")) install.packages("data.table")
if(!require("ggplot2")) install.packages("ggplot2")



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

if (sv.type == "all") SVs <- SVs
if (sv.type == "del") SVs <- SVs[grep("^del", SVs, perl = TRUE)]
if (sv.type == "ins") SVs <- SVs[grep("^ins", SVs, perl = TRUE)]
if (sv.type == "inv") SVs <- SVs[grep("^inv", SVs, perl = TRUE)]
if (sv.type == "dup") SVs <- SVs[grep("^dup", SVs, perl = TRUE)]

# load snps info
SNPs <- fread(snp.file, header = TRUE, data.table = FALSE)
SNPs <- SNPs[, c(1, 3, 4)]

for (chr in 1:10) {

  # subset SVs based on chr
  SVs.chr <- subset(SVs, chrom == chr)[1]
  SVs.chr <- as.character(SVs.chr[, 1])
  # subset SVs based on chr
  SNPs.chr <- subset(SNPs, chrom == chr)[1]
  SNPs.chr <- as.character(SNPs.chr[, 1])

  # load one chr at a time
  plink.file.chr <- gsub("chr[0-9]+", paste0("chr", chr), plink.file, perl = TRUE)
  # open table with LD among markers
  LD_results <- fread(plink.file.chr, header = TRUE, data.table = FALSE)
  
  # if (length(SVs.chr) != length(SNPs.chr)) stop(paste0("Number of SVs and closest SNPs in chr", chr, " doesn't match"))
  
  # filter ld df to have only SV and it's closest SNP
  LD_results_filtered <- data.frame(matrix(nrow = 0, ncol = NCOL(LD_results)), stringsAsFactors = FALSE)
  colnames(LD_results_filtered) <- colnames(LD_results)
  for (idx in 1:length(SVs.chr)) {
    
    row.to.keep <- which(LD_results$SNP_A == SVs.chr[idx] & LD_results$SNP_B == SNPs.chr[idx] |
                           LD_results$SNP_A == SNPs.chr[idx] & LD_results$SNP_B == SVs.chr[idx])
    LD_results_filtered <- rbind(LD_results_filtered, LD_results[row.to.keep, ])
    
  }
  
  # add column with distance between snp and closest SV boundary
  # LD_results_filtered$dist_to_sv2 <- abs(LD_results[, "BP_B"] - LD_results[, "BP_A"])
  LD_results_filtered$dist_to_sv <- apply(LD_results_filtered, MARGIN = 1, function(marker) {

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
  
  # add info
  LD_results_filtered <- add.info.to.ld.df(LD_results_filtered)

  # write filtered ld table
  outfile.highest <- gsub(".ld", ".highest-ld.ld", plink.file.chr, fixed = TRUE)
  outfile.highest <- rev(unlist(strsplit(outfile.highest, split = "/")))[1]
  fwrite(LD_results_filtered, paste0(out.dir.ld, "/", outfile.highest), sep = "\t", quote = FALSE, row.names = FALSE, na = NA)
  
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
         y = "Count",
         subtitle = paste0("(", NROW(LD_results_all), " closest SNPs)")) +
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
