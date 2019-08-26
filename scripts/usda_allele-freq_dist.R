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

# create a list to store all markers for all populations
all.markers.list <- list()

# create a list to store the markers that have AF < 0.25 or > 0.75
extreme.markers.list <- list()

for (cross in cross.list) {
  
  # read in file with polymorphic markers without segregation distortion from rqtl analysis
  markers.filename <- paste0("analysis/qc/", cross, "/recomb-freq_", cross, "_rils.txt")
  markers.infile <- fread(markers.filename, header = TRUE, data.table = FALSE)
  markers <- markers.infile[, "marker"]
  pop.size <- markers.infile[1, "pop_size"]
  
  # read in tassel summary table
  tassel.filename <- paste0("analysis/qc/", cross, "/", cross, "_SiteSummary.txt")
  tassel.infile <- fread(tassel.filename, header = TRUE, data.table = FALSE)
  
  # filter site summary table from tassel based on polymorphic markers from rqtl
  tassel.infile <- tassel.infile[which(tassel.infile[, "Site Name"] %in% markers), ]
  # keep only cokumna that may have useful info
  tassel.infile <- tassel.infile[,c("Site Name", "Chromosome", "Physical Position",
                                    "Number of Taxa", "Major Allele Frequency",
                                    "Minor Allele Frequency", "Proportion Heterozygous")]
  
  # add all markers to list
  all.markers.list[[cross]] <- tassel.infile[, "Site Name"]
  
  # add markers with extreme allele frequencies to list
  AF.more.0.75 <- tassel.infile[which(tassel.infile[, "Major Allele Frequency"] > 0.75), "Site Name"]
  AF.less.0.25 <- tassel.infile[which(tassel.infile[, "Minor Allele Frequency"] < 0.25), "Site Name"]
  extreme.markers.list[[cross]] <- union(AF.more.0.75, AF.less.0.25)
  
  # plot distribution of allele frequencies
  allele.freq <- tassel.infile %>% 
    select(`Major Allele Frequency`, `Minor Allele Frequency`) %>%
    gather(key = "AF_type", value = "AF_value")
  
  ggplot(allele.freq, aes(x = AF_value, fill = "#053061")) +
    geom_histogram(binwidth = 0.07, show.legend = FALSE) +
    scale_x_continuous(name = "Allele Frequency", limits = c(0, 1), breaks = seq(0,1,0.25)) + 
    scale_fill_manual(values = "black") +
    labs(title = paste0(cross),
         subtitle = paste0("Population size = ", pop.size, " / markers = ", NROW(tassel.infile)),
         y = "Count") +
    theme(axis.text = element_text(size=rel(2)),
          axis.title = element_text(size=rel(2)))
  
  # save plot
  figure_name <- paste0("analysis/qc/", cross, "/allele-freq_dist_", cross, ".png")
  ggsave(figure_name, device = "png")
}



#### markers with AF < 0.25 or > 0.75 ----

extreme.markers <- unlist(extreme.markers.list, use.names = FALSE)

# total number of markers found
length(extreme.markers)
# 975

# total number of unique markers (i.e., without markers that show up in more than one population)
length(unique(extreme.markers))
# 910

# only 65 (975 - 910) markers have AF < 0.25 or > 0.75 in more than one population
