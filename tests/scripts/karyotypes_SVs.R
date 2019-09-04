#### intro ----

# plot distribution of SVs in chromosomes for each parent



#### libraries ----

library(data.table)
library(ggplot2)



#### plot karyotypes ----

# chromosomes
chrms <- read.delim("data/B73_RefGen_V4_chrm_info.txt", sep = "\t", header = TRUE)
chrms <- data.frame(chr = chrms$chr,
                    start_pos = 0,
                    end_pos = chrms$length)

# centromeres
centros <- read.delim("data/centromeres_Schneider-2016-pnas_v4.bed", sep = "\t", header = TRUE)
centros <- data.frame(chr = centros$chr,
                      start_pos = centros$start_pos,
                      end_pos = centros$end_pos)

# load data
geno.parents <- fread("tests/data/usda_SVs_7parents.sorted.hmp.txt",
                      header = TRUE, data.table = FALSE)

# get names of parents
name.parents <- colnames(geno.parents)[12:NCOL(geno.parents)]

# create folder to store results
if (!dir.exists("tests/analysis/SVs")) {
  dir.create("tests/analysis/SVs")
}

for (parent in name.parents) {
  
  # keep only necessary information for parent
  parent.info <- geno.parents[, c("rs", "chrom", "pos", parent)]
  
  # rename columns for compatibility (chromosome column should be "chr" because
  # "centros" and "chrms" data have "chr" as a column as well)
  colnames(parent.info) <- c("markers", "chr", "pos", "SV_type")
  
  # select only SVs
  parent.info.SVs <- subset(parent.info, SV_type == "TT")
  
  # count hets and missing (to add into plot)
  SV.count <- NROW(parent.info.SVs)
  
  # plot karyotype
  karyo.plot <- ggplot() +
    geom_segment(data = chrms,
                 aes(x = 0, xend = 0, y = start_pos, yend = end_pos),
                 lineend = "round", color = "Gainsboro", size = 5) +
    geom_segment(data = centros,
                 aes(x = 0, xend = 0, y = start_pos, yend = end_pos),
                 lineend = "round", color = "DimGray", size = 5) +
    geom_segment(data = parent.info.SVs,
                 aes(x = 0, xend = 0, y = pos, yend = pos + 1e6),  # increased size to be able to see the marker
                 lineend = "butt", color = "firebrick", size = 5, alpha = 0.1) +
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
    labs(caption = paste0(parent, " - ", SV.count, " SVs to B73"),
         x = "Chromosomes", y = "Genomic positions (Mb)")
  
  # print(karyo.plot)

  karyo.name <- paste0("tests/analysis/SVs/karyotype_SVs_", parent, ".png")
  ggsave(filename = karyo.name, plot = karyo.plot, device = "png")
  
}
