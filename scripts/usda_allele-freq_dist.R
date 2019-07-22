#### description ----

# extract the marker names from the recombination frequency tables (polymorphic, no seg distortion,
# etc.; `analysis/qc > cross folder > recomb-freq txt file`) and filter tassel's _SiteSummary
# tables. Then use this filtered table to make plots of distributions of allele frequencies.



#### libraries ----

library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)



#### create histograms for all crosses ----

# add folder names (i.e. cross names) into a vector
cross.list <- list.dirs("analysis/qc", full.names = FALSE, recursive = FALSE)

for (cross in cross.list) {
  
  # read in file with polymorphic markers without segregation distortion from rqtl analysis
  markers.filename <- paste0("analysis/qc/", cross, "/recomb-freq_", cross, "_rils.txt")
  markers.infile <- fread(markers.filename, header = TRUE, data.table = FALSE)
  markers <- markers.infile[, "marker"]
  
  # read in tassel summary table
  tassel.filename <- paste0("analysis/qc/", cross, "/", cross, "_SiteSummary.txt")
  tassel.infile <- fread(tassel.filename, header = TRUE, data.table = FALSE)
  
  # filter site summary table from tassel based on polymorphic markers from rqtl
  tassel.infile <- tassel.infile[which(tassel.infile[, "Site Name"] %in% markers), ]
  # keep only cokumna that may have useful info
  tassel.infile <- tassel.infile[,c("Site Name", "Chromosome", "Physical Position",
                                    "Number of Taxa", "Major Allele Frequency",
                                    "Minor Allele Frequency", "Proportion Heterozygous")]
  
  # plot distribution of allele frequencies
  allele.freq <- tassel.infile %>% 
    select(`Major Allele Frequency`, `Minor Allele Frequency`) %>%
    gather(key = "AF_type", value = "AF_value")
  
  ggplot(allele.freq, aes(x = AF_value, fill = "#053061")) +
    geom_histogram(binwidth = 0.07, show.legend = FALSE) +
    scale_x_continuous(name = "allele frequency", limits = c(0, 1), breaks = seq(0,1,0.25)) + 
    scale_fill_manual(values = "#053061") +
    labs(title = paste0(cross),
         subtitle = paste0("pop size = ", tassel.infile[1,"Number of Taxa"], " / markers = ",
                           NROW(tassel.infile)))
  
  # save plot
  figure_name <- paste0("analysis/qc/", cross, "/allele-freq_dist_", cross, ".png")
  ggsave(figure_name, device = "png")
}
