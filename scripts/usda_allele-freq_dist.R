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
    scale_fill_manual(values = "#053061") +
    labs(title = paste0(cross),
         subtitle = paste0("Population size = ", tassel.infile[1,"Number of Taxa"], " / markers = ",
                           NROW(tassel.infile)),
         y = "Count")
  
  # save plot
  figure_name <- paste0("analysis/qc/", cross, "/allele-freq_dist_", cross, ".png")
  ggsave(figure_name, device = "png")
}



#### markers with AF < 0.25 or > 0.75 ----

extreme.markers <- unlist(extreme.markers.list, use.names = FALSE)

# total number of markers found
length(extreme.markers)
# 1119

# total number of unique markers (i.e., without markers that show up in more than one population)
length(unique(extreme.markers))
# 1038

# only 81 markers have AF < 0.25 or > 0.75 in more than one population, and the majority of such
# markers are unique to 


# how many of the 1038 markers that have extreme AF in only one population are actually present in
# other population and has AF between 0.25 and 0.75?
for (i in 1:length(all.markers.list)) {
  print(
    round(sum(extreme.markers.unique %in% all.markers.list[[i]]) / length(extreme.markers.unique),
          digits = 2)
  )
}

# how many of the 81 markers with extreme AF in more than one population are present in each population?
for (i in 1:length(all.markers.list)) {
  print(
    round(sum(extreme.markers.duplicated %in% extreme.markers.list[[i]]) / length(extreme.markers.duplicated),
          digits = 2)
  )
}
