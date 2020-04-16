#### arguments for command line ----

args <- commandArgs(trailingOnly = TRUE)

# help
if (all(length(args) == 1 & args == "-h" | args == "--help")) {
  cat("
      Description: plot karyotypes for RILs of a biparental cross.
      Credits: code for plot was modified from Carles Hernandez-Ferrer's blog at
               https://carleshf.github.io/chromosome_karyotype_plot_with_ggplot2/


      Usage: Rscript plot_ril_karyotypes.R [chr_info] [centromere_info] [qc_folder]")
  quit()
}

# make sure the correct number of arguments are used
if (length(args) != 3) {
  stop("incorrect number of arguments provided.

       Usage: Rscript plot_ril_karyotypes.R [chr_info] [centromere_info] [qc_folder]
       ")
}

# assign arguments to variables
chr.info <- args[1]
cent.info <- args[2]
qc.folder <- args[3]

# chr.info <- "data/B73_RefGen_V4_chrm_info.txt"
# cent.info <- "data/centromeres_Schneider-2016-pnas_v4.bed"
# qc.folder <- "analysis/qc"


#### libraries ----

if(!require("data.table")) install.packages("data.table")
if(!require("ggplot2")) install.packages("ggplot2")




#### load chromosome and centromere positions ----

# chromosomes
chrms <- fread(chr.info, header = TRUE, data.table = FALSE)
chrms <- data.frame(chr = chrms$chr, start_pos = 0, end_pos = chrms$length)

# centromeres
centros <- read.delim(cent.info, sep = "\t", header = TRUE)
centros <- data.frame(chr = centros$chr, start_pos = centros$start_pos, end_pos = centros$end_pos)



#### plot karyotypes ----

# get names of RIL populations
RIL.pop.list <- list.dirs(path = qc.folder, full.names = FALSE, recursive = FALSE)
RIL.pop.list <- RIL.pop.list[grep("x", RIL.pop.list)]

for (RIL.pop in RIL.pop.list) {
  # load data with marker positions
  markers.filename <- paste0(qc.folder, "/", RIL.pop, "/recomb-freq_", RIL.pop, "_rils.txt")
  markers.infile <- fread(markers.filename, header = TRUE, data.table = FALSE)

  # load genotypic data for all RILs in the RIL population
  geno.data.filename <- paste0(qc.folder, "/", RIL.pop, "/geno-data_", RIL.pop, "_after_filtering.txt")
  geno.data.infile <- fread(geno.data.filename, header = TRUE, data.table = FALSE)

  # randomly select 5 RILs per population to plot karyotype
  set.seed(184)
  selected.RILs <- sample(colnames(geno.data.infile[-1]), size = 5, replace = FALSE)

  for (RIL in selected.RILs) {
    # merge information of RIL of interest with respective marker positions
    geno.data <- cbind(markers.infile[, c("marker", "chr", "pos")], geno.data.infile[, RIL])
    colnames(geno.data)[4] <- "geno"

    # select only AA genotype
    geno.data.AA <- subset(geno.data, geno == "AA")
    geno.data.BB <- subset(geno.data, geno == "BB")

    # get genotype frequency to add in the plot
    geno.freq.AA <- NROW(geno.data.AA) / (NROW(geno.data.AA) + NROW(geno.data.BB))
    geno.freq.BB <- NROW(geno.data.BB) / (NROW(geno.data.AA) + NROW(geno.data.BB))

    karyo.plot <- ggplot() +
      geom_segment(data = chrms,
                   aes(x = 0, xend = 0, y = start_pos, yend = end_pos),
                   lineend = "round", color = "Gainsboro", size = 5) +
      geom_segment(data = centros,
                   aes(x = 0, xend = 0, y = start_pos, yend = end_pos),
                   lineend = "round", color = "DimGray", size = 5) +
      geom_segment(data = geno.data.AA,
                   aes(x = 0, xend = 0, y = pos, yend = pos + 1e6),  # increased size to be able to see the marker
                   lineend = "butt", color = "darkblue", size = 5, alpha = 0.2) +
      geom_segment(data = geno.data.BB,
                   aes(x = 0, xend = 0, y = pos, yend = pos + 1e6),  # increased size to be able to see the marker
                   lineend = "butt", color = "firebrick", size = 5, alpha = 0.2) +
      scale_y_reverse(breaks = seq(0, 3.5e8, 0.50e8), labels = c(1, seq(50, 350, 50))) +
      scale_x_discrete(position = "top") +
      theme_minimal() +
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            plot.caption = element_text(size = rel(1.1), color = "DimGray"),
            axis.text = element_text(size=rel(2)),
            axis.title = element_text(size=rel(2)),
            strip.text.x = element_text(size=rel(2))) +
      facet_grid(cols = vars(chr), switch = "y") +
      labs(caption = paste0(RIL.pop, " - ", gsub("RIL_", "RIL ", RIL), "\n\n",
                            "AA freq (blue): ", round(geno.freq.AA, digits = 2), "\n",
                            "BB freq (red): ", round(geno.freq.BB, digits = 2)),
           x = "Chromosomes", y = "Genomic positions (Mb)")

    karyo.name <- paste0(qc.folder, "/", RIL.pop, "/karyotype_", RIL.pop, "_", gsub("RIL_", "ril-", RIL),".pdf")
    ggsave(filename = karyo.name, plot = karyo.plot, device = "pdf")
  }
}
