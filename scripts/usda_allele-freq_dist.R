#### arguments for command line ----

args <- commandArgs(trailingOnly = TRUE)

# help
if (all(length(args) == 1 & args == "-h" | args == "--help")) {
  cat("
      Description: extract the marker names from the recombination frequency tables (polymorphic, no seg distortion,
                   etc.) and filter tassel's site summary tables. Then use this filtered table to make plots of
                   distributions of allele frequencies.


      Usage: Rscript usda_allele-freq_dist.R [qc_folder]")
  quit()
}

# make sure the correct number of arguments are used
if (length(args) != 1) {
  stop("incorrect number of arguments provided.

       Usage: Rscript usda_allele-freq_dist.R [qc_folder]
       ")
}

# assign arguments to variables
qc.folder <- args[1]


# qc.folder <- "analysis/qc"



#### libraries ----

if(!require("data.table")) install.packages("data.table")
if(!require("dplyr")) install.packages("dplyr")
if(!require("tidyr")) install.packages("tidyr")
if(!require("ggplot2")) install.packages("ggplot2")



#### create histograms for all crosses ----

# add folder names (i.e. cross names) into a vector
cross.list <- list.dirs(qc.folder, full.names = FALSE, recursive = FALSE)

# create a list to store all markers for all populations
all.markers.list <- list()

# create a list to store the markers that have AF < 0.25 or > 0.75
extreme.markers.list <- list()

for (cross in cross.list) {

  # read in file with polymorphic markers without segregation distortion from rqtl analysis
  markers.filename <- paste0(qc.folder, "/", cross, "/recomb-freq_", cross, "_rils.txt")
  markers.infile <- fread(markers.filename, header = TRUE, data.table = FALSE)
  markers <- markers.infile[, "marker"]
  pop.size <- markers.infile[1, "pop_size"]

  # read in tassel summary table
  tassel.filename <- paste0(qc.folder, "/", cross, "/", cross, "_SiteSummary.txt")
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

  af.plot <- ggplot(allele.freq, aes(x = AF_value, fill = "#053061")) +
    geom_histogram(binwidth = 0.07, show.legend = FALSE) +
    scale_x_continuous(name = "Allele Frequency", limits = c(0, 1), breaks = seq(0,1,0.25)) +
    scale_fill_manual(values = "black") +
    labs(title = paste0(cross),
         subtitle = paste0("Population size = ", pop.size, " / markers = ", NROW(tassel.infile)),
         y = "Count") +
    theme(axis.text = element_text(size=rel(2)),
          axis.title = element_text(size=rel(2)))

  # save plot
  figure_name <- paste0(qc.folder, "/", cross, "/allele-freq_dist_", cross, ".pdf")
  ggsave(filename = figure_name, plot = af.plot, device = "pdf")
}


# # markers with AF < 0.25 or > 0.75
# extreme.markers <- unlist(extreme.markers.list, use.names = FALSE)
# # total number of markers found
# length(extreme.markers)
# # 975
# # total number of unique markers (i.e., without markers that show up in more than one population)
# length(unique(extreme.markers))
# # 910
# # only 65 (975 - 910) markers have AF < 0.25 or > 0.75 in more than one population
