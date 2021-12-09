library(data.table)
library(ggplot2)

usage <- function() {
  cat("
description: plot SNPs and SVs MAF distribution

usage: Rscript distribution_maf_snps-svs.R [summary_file]

positional arguments:
  summary_file              path to TASSEL site summary file
  svs_file                  path to list of SV names file
  outfile                   name of output file

optional argument:
  --help                    show this helpful message
  
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

args <- commandArgs(trailingOnly = TRUE)
if ("--help" %in% args) usage() & q(save = "no")

# assert to have the correct optional arguments
if (length(args) != 3) stop(usage(), "missing positional argument(s)")

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

# plot distribution before downsampling
plot <- ggplot(summary, aes(x = `Minor Allele Frequency`)) +
  # geom_histogram(fill = "#900721", binwidth = 0.01) +
  geom_density(fill = "#900721") +
  facet_grid(marker_type ~ ., scales = "free_y") +
  labs(title = "MAF distribuition after downsampling and LD calculation",
       subtitle = paste0("(", n_markers_total, " markers)"),
       x = "MAF",
       y = "Density") +
  coord_cartesian(xlim = c(0, 0.5)) +
  theme_bw()

ggsave(filename = outfile, plot = plot, device = "pdf", width = 12, height = 7)


#### debug ----

# summary_file <- "analysis/ld_downsample/datasets/markers_chr10_SiteSummary.txt"
# svs_file <- "data/SVs_IDs_poly.txt"
# outfile <- "analysis/ld_downsample/datasets/markers_chr10_MAF-downsampled.txt"

