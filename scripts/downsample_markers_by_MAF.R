library(data.table)
library(ggplot2)

usage <- function() {
  cat("
description: downsample SNPs in a way that SNPs MAF distribution are similar to SVs

usage: Rscript downsample_markers_by_MAF.R [summary_file] [svs_file] [outfile]

positional arguments:
  summary_file              path to TASSEL site summary file
  svs_file                  path to list of SV names file
  outfile                   name of output file

optional argument:
  --help                    show this helpful message
  --prop-downsample=VALUE   proportion of markers left after downsample (default: 0.25)
  --bin-size=VALUE          size of bins to separate the distribution and count number 
                            of markers (default: 0.01)
  --reps=VALUE              number of times the random sampling will be performed (default: 3)
  --seed=VALUE              value for set.seed (default: NULL; random number is selected)
  
"
  )
}

getArgValue <- function(arg) {
  
  # get value from command-line arguments
  arg <- unlist(strsplit(arg, "="))
  if (length(arg) == 1) return(TRUE)
  if (length(arg) > 1) return(arg[2])
  
}


#### command line options ----

# set default
prop_downsample <- "0.25"
bin_size <- "0.01"
reps <- "3"
seed <- NULL

args <- commandArgs(trailingOnly = TRUE)
if ("--help" %in% args) usage() & q(save = "no")

# assert to have the correct optional arguments
if (length(args) < 3) stop(usage(), "missing positional argument(s)")

if (length(args) > 3) {
  
  opt_args <- args[-1:-3]
  opt_args_allowed <- c("--prop-downsample", "--bin-size", "--reps", "--seed")
  opt_args_requested <- as.character(sapply(opt_args, function(x) unlist(strsplit(x, split = "="))[1]))
  if (any(!opt_args_requested %in% opt_args_allowed)) stop(usage(), "wrong optional argument(s)")
  
  # change default based on the argument provided
  for (argument in opt_args_allowed) {
    if (any(grepl(argument, opt_args_requested))) {
      arg_name <- gsub("-", "_", gsub("--", "", argument))
      arg_value <- getArgValue(opt_args[grep(argument, opt_args)])
      assign(arg_name, arg_value)
    }
  }
  
}

# make sure optional arguments are valid
if (suppressWarnings(!is.na(as.numeric(prop_downsample)))) {
  prop_downsample <- as.numeric(prop_downsample)
} else {
  stop("Optional argument '--prop-downsample' should be a number")
}

if (suppressWarnings(!is.na(as.numeric(bin_size)))) {
  bin_size <- as.numeric(bin_size)
} else {
  stop("Optional argument '--bin-size' should be a number")
}

if (suppressWarnings(!is.na(as.numeric(reps)))) {
  reps <- as.numeric(reps)
} else {
  stop("Optional argument '--reps' should be a number")
}

if (is.null(seed)) {
  seed <- ceiling(runif(1, 0, 1000000))
} else {
  if (suppressWarnings(!any(is.na(as.numeric(seed))))) {
    seed <- as.numeric(seed)
  } else {
    stop("Optional argument '--seed' should be a number")
  }
}

# get positional arguments
summary_file <- args[1]
svs_file <- args[2]
outfile <- args[3]



#### downsample markers ----

# load data
summary <- fread(summary_file, header = TRUE, data.table = FALSE)
SVs <- fread(svs_file, header = TRUE, data.table = FALSE)
SVs <- SVs[, 1]

# keep only MAF
summary <- summary[, c("Site Name", "Minor Allele Frequency")]
summary[which(summary[, "Site Name"] %in% SVs), "marker_type"] <- "SV"
summary[which(!summary[, "Site Name"] %in% SVs), "marker_type"] <- "SNP"

# get total number of markers 
n_markers_total <- NROW(summary)
n_svs_total <- NROW(summary[which(summary[, "Site Name"] %in% SVs), ])
n_snps_total <- NROW(summary[which(!summary[, "Site Name"] %in% SVs), ])
n_snps_downsample <- ceiling(n_snps_total * prop_downsample)

# # plot distribution before downsampling
# plot_before <- ggplot(summary, aes(x = `Minor Allele Frequency`)) +
#   # geom_histogram(fill = "#900721", binwidth = 0.01) +
#   geom_density(fill = "#900721") +
#   facet_grid(marker_type ~ ., scales = "free_y") +
#   labs(title = "MAF distribuition before downsampling",
#        subtitle = paste0("(", n_markers_total, " markers)"),
#        x = "MAF",
#        y = "Density") +
#   coord_cartesian(xlim = c(0, 0.5)) +
#   theme_bw()
# 
# ggsave(filename = gsub(".txt", ".before.pdf", outfile), plot = plot_before,
#        device = "pdf", width = 12, height = 7)

# analyze MAF distribution per bin
snps_downsampled <- vector(mode = "list", length = reps)
for (bin in seq(0, (0.5 -  bin_size), by = bin_size)) {
  
  # subset summary table
  summary_bin_svs <- subset(summary, marker_type == "SV"
                            & `Minor Allele Frequency` >= bin
                            & `Minor Allele Frequency` < bin + bin_size)
  summary_bin_snps <- subset(summary, marker_type == "SNP"
                             & `Minor Allele Frequency` >= bin
                             & `Minor Allele Frequency` < bin + bin_size)
  
  # get proportion of svs in bin
  prop_svs <- NROW(summary_bin_svs) / n_svs_total
  # get number of snps necessary to match proportion of SVs in that bin
  n_snps_downsample_bin <- ceiling(n_snps_downsample * prop_svs)
  
  # randomly downsample SNPs
  for (rep in 1:reps) {
    set.seed(seed * rep)
    if (length(summary_bin_snps$`Site Name`) >= n_snps_downsample_bin) {
      snps_downsampled[[rep]] <- append(snps_downsampled[[rep]],
                                        sample(summary_bin_snps$`Site Name`, size = n_snps_downsample_bin))
    } else {
      snps_downsampled[[rep]] <- append(snps_downsampled[[rep]], summary_bin_snps$`Site Name`)
    }
  }

}
rm(summary_bin_snps, summary_bin_svs)

# plot distribution after downsampling and save list of selected markers
for (rep in 1:reps) {
  
  # filter summary
  summary_downsample <- subset(summary, `Site Name` %in% snps_downsampled[[rep]])
  # add svs for plotting
  summary_downsample <- rbind(summary_downsample, subset(summary, marker_type == "SV"))
  
  # # plot
  # plot_after <- ggplot(summary_downsample, aes(x = `Minor Allele Frequency`)) +
  #   # geom_histogram(fill = "#900721", binwidth = 0.01) +
  #   geom_density(fill = "#900721") +
  #   facet_grid(marker_type ~ ., scales = "free_y") +
  #   labs(title = "MAF distribuition after downsampling",
  #        subtitle = paste0("(", NROW(summary_downsample), " markers - ", bin_size, " bins - rep", rep, ")"),
  #        x = "MAF",
  #        y = "Density") +
  #   coord_cartesian(xlim = c(0, 0.5)) +
  #   theme_bw()
  # 
  # ggsave(filename = gsub(".txt", paste0(".after.rep", rep, ".pdf"), outfile), plot = plot_after,
  #        device = "pdf", width = 12, height = 7)
  
  # save list of markers to sample
  sampled_markers_rep <- data.frame(summary_downsample$`Site Name`)
  fwrite(sampled_markers_rep, file = gsub(".txt", paste0(".rep", rep, ".txt"), outfile),
         quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
  
}



#### debug ----

# summary_file <- "analysis/ld_downsample/datasets/markers_chr10_SiteSummary.txt"
# svs_file <- "data/SVs_IDs_poly.txt"
# outfile <- "analysis/ld_downsample/datasets/markers_chr10_MAF-downsampled.txt"
# prop_downsample <- 0.25
# bin_size <- 0.01
# reps <- 3
# seed <- 123
