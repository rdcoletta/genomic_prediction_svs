library(data.table)
library(ggplot2)

usage <- function() {
  cat("
description: plot how many SNPs and SVs are missing per marker and per RIL, and also
             plot the distribution of these markers along the chromosome

usage: Rscript qc_snp-sv_hmp.R [infile_name] [sv_list] [outfile] [missing_threshold]

positional arguments:
  infile_name             hapmap file containing SNPs and SVs
  sv_list                 single-column file containing only SV IDs
  outfile                 output filename in with '.pdf' extension
  missing_threshold       amount of missing data to filter markers (number from 0 to 1)
  

optional argument:
  --help                  show this helpful message

"
  )
}



#### command line options ----

args <- commandArgs(trailingOnly = TRUE)
if ("--help" %in% args) usage() & q(save = "no")

# assert to have the correct optional arguments
if (length(args) != 4) stop(usage(), "missing positional argument(s)")

# get arguments
infile_name <- args[1]
sv_list <- args[2]
outfile <- args[3]
missing_threshold <- as.numeric(args[4])
# infile_name <- "analysis/qc/snp-sv_hmp/usda_rils_projected-SVs-SNPs_SiteSummary.txt"
# sv_list <- "data/SVs_IDs_poly.txt"
# outfile <- "analysis/qc/snp-sv_hmp/missing_markers_rils.pdf"
# missing_threshold <- 0.25



#### data qc ----

# load TASSEL summary
summary <- fread(infile_name, header = TRUE, data.table = FALSE)
# keep only interesting info
summary <- summary[, c("Site Name", "Chromosome", "Physical Position",
                       "Major Allele", "Minor Allele", "Minor Allele Frequency",
                       "Minor Allele Proportion", "Proportion Missing",
                       "Proportion Heterozygous")]

# load SV information
SVs <- fread(sv_list, header = FALSE, data.table = FALSE)
SVs <- SVs[, 1]

# add marker type info to summary table
summary[which(!summary[, "Site Name"] %in% SVs), "Marker Type"] <- "SNP"
summary[which(summary[, "Site Name"] %in% SVs), "Marker Type"] <- "SV"


# analyze snps and svs separately
for (marker in c("SNP", "SV")) {
  
  # separate data to have only SNPs or SVs
  summary_marker <- summary[summary[, "Marker Type"] == marker, ]
  
  prop_missing <- round(sum(summary_marker$`Proportion Missing` < missing_threshold) / NROW(summary_marker),
                             digits = 3)
  
  cat("Proportion of ", marker, "s < ", missing_threshold * 100, "% missing data: ", prop_missing * 100,
      "%\n", sep = "")
  
  # plot
  plot_missing <- ggplot(data = summary_marker, aes(x = `Proportion Missing`)) +
    geom_histogram(binwidth = 0.01) +
    coord_cartesian(xlim = c(0, 1)) +
    labs(title = paste0("Proportion of ", marker, "s < ", missing_threshold * 100,
                        "% missing data: ", prop_missing * 100, "%"),
         x = paste0("Percent missing data"),
         y = paste0("Number of ", marker, "s"))
  # save results
  outfile_marker <- gsub(".pdf", paste0(".", marker, "s.pdf"), outfile)
  ggsave(filename = outfile_marker, plot = plot_missing, device = "pdf")
  
}


# plot distribution of markers along each chomosome
plot_distribution <- ggplot(summary, aes(x = `Physical Position`, fill = `Marker Type`)) + 
  geom_density(alpha = 0.5, position = "identity") +
  facet_wrap(~`Chromosome`, scales = "free_x") +
  scale_x_continuous(labels = function(x) x/1000000) +
  labs(x = "Position (Mb)",
       y = "Density",
       fill = "Marker type")

ggsave(filename = gsub(".pdf", paste0(".", marker, "s.dist-chr.pdf"), outfile),
       plot = plot_distribution, device = "pdf")

# plot distribution of low missing markers along each chomosome
summary_low_missing <- subset(summary, `Proportion Missing` < missing_threshold)
plot_distribution_low <- ggplot(summary_low_missing, aes(x = `Physical Position`, fill = `Marker Type`)) + 
  geom_density(alpha = 0.5, position = "identity") +
  facet_wrap(~`Chromosome`, scales = "free_x") +
  scale_x_continuous(labels = function(x) x/1000000) +
  labs(x = "Position (Mb)",
       y = "Density",
       fill = "Marker type")

ggsave(filename = gsub(".pdf", paste0(".", marker, "s.dist-chr_low-missing.pdf"), outfile),
       plot = plot_distribution_low, device = "pdf")
