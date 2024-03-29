#### arguments for command line ----

args <- commandArgs(trailingOnly = TRUE)

# help
if (all(length(args) == 1 & args == "-h" | args == "--help")) {
  cat("
      Description: plot karyotypes of parents.
      Credits: code for plot was modified from Carles Hernandez-Ferrer's blog at
               https://carleshf.github.io/chromosome_karyotype_plot_with_ggplot2/


      Usage: Rscript plot_parents_karyotypes.R [chr_info] [centromere_info] [parents_hapmap] [output_folder]")
  quit()
}

# make sure the correct number of arguments are used
if (length(args) != 4) {
  stop("incorrect number of arguments provided.

       Usage: Rscript plot_parents_karyotypes.R [chr_info] [centromere_info] [parents_hapmap] [output_folder]
       ")
}

# assign arguments to variables
chr.info <- args[1]
cent.info <- args[2]
parents.file <- args[3]
out.folder <- args[4]

# chr.info <- "data/B73_RefGen_V4_chrm_info.txt"
# cent.info <- "data/centromeres_Schneider-2016-pnas_v4.bed"
# parents.file <- "data/usda_22kSNPs_7parents.sorted.diploid.hmp.txt"
# out.folder <- "analysis/qc/parents"


#### libraries ----

if(!require("data.table")) install.packages("data.table")
if(!require("ggplot2")) install.packages("ggplot2")



#### load chromosome and centromere positions ----

# chromosomes
chrms <- read.delim(chr.info, sep = "\t", header = TRUE)
chrms <- data.frame(chr = chrms$chr, start_pos = 0, end_pos = chrms$length)

# centromeres
centros <- read.delim(cent.info, sep = "\t", header = TRUE)
centros <- data.frame(chr = centros$chr, start_pos = centros$start_pos, end_pos = centros$end_pos)



#### plot karyotypes ----

# load data
geno.parents <- fread(parents.file, header = TRUE, data.table = FALSE)

# get names of parents
name.parents <- colnames(geno.parents)[12:NCOL(geno.parents)]

# create folder to store results
if (!dir.exists(out.folder)) dir.create(out.folder)


for (parent in name.parents) {

  # keep only necessary information for parent
  parent.info <- geno.parents[, c("rs#", "chrom", "pos", parent)]


  # transform genotypes in either homo, het, or missing
  parent.info[, parent] <- sapply(X = parent.info[, parent],
                                  USE.NAMES = FALSE,
                                  FUN = function(marker) {
                                    # get alleles
                                    allele <- unlist(strsplit(marker, split = ""))
                                    # check first if it's missing
                                    if (allele[1] == "N" | allele[2] == "N") {
                                      marker.type <- "missing"
                                    } else if (allele[1] == allele[2]) {
                                      # if it's not missing, check if alleles are the same...
                                      marker.type <- "homo"
                                    } else {
                                      # ...or not
                                      marker.type <- "het"
                                    }
                                    return(marker.type)
                                  })
  # rename columns for compatibility (chromosome column should be "chr" because
  # "centros" and "chrms" data have "chr" as a column as well)
  colnames(parent.info) <- c("markers", "chr", "pos", "marker_type")

  # subset df by marker type
  parent.info.homo <- subset(parent.info, marker_type == "homo")
  parent.info.het <- subset(parent.info, marker_type == "het")
  parent.info.missing <- subset(parent.info, marker_type == "missing")

  # count hets and missing (to add into plot)
  het.count <- NROW(parent.info.het)
  missing.count <- NROW(parent.info.missing)

  # plot karyotype
  karyo.plot <- ggplot() +
    geom_segment(data = chrms,
                 aes(x = 0, xend = 0, y = start_pos, yend = end_pos),
                 lineend = "round", color = "Gainsboro", size = 5) +
    geom_segment(data = centros,
                 aes(x = 0, xend = 0, y = start_pos, yend = end_pos),
                 lineend = "round", color = "DimGray", size = 5) +
    # geom_segment(data = parent.info.homo,
    #              aes(x = 0, xend = 0, y = pos, yend = pos + 1e6),  # increased size to be able to see the marker
    #              lineend = "butt", color = "#fdae61", size = 5, alpha = 0.3) +
    geom_segment(data = parent.info.het,
                 aes(x = 0, xend = 0, y = pos, yend = pos + 1e6),  # increased size to be able to see the marker
                 lineend = "butt", color = "darkblue", size = 5, alpha = 0.5) +
    geom_segment(data = parent.info.missing,
                 aes(x = 0, xend = 0, y = pos, yend = pos + 1e6),  # increased size to be able to see the marker
                 lineend = "butt", color = "firebrick", size = 5, alpha = 0.5) +
    scale_y_reverse(breaks = seq(0, 3.5e8, 0.50e8), labels = c(1, seq(50, 350, 50))) +
    scale_x_discrete(position = "top") +
    theme_minimal() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.caption = element_text(size = rel(1.1), color = "DimGray"),
          axis.text = element_text(size=rel(1.2)),
          axis.title = element_text(size=rel(1.2)),
          strip.text.x = element_text(size=rel(1.2))) +
    facet_grid(cols = vars(chr), switch = "y") +
    labs(caption = paste0("Parent ", parent, "\n\n",
                          "Het markers (blue): ", het.count, "\n",
                          "Missing markers (red): ", missing.count),
         x = "Chromosomes", y = "Genomic positions (Mb)")

  karyo.name <- paste0(out.folder, "/karyotype_", parent, ".pdf")
  ggsave(filename = karyo.name, plot = karyo.plot, device = "pdf")

}
