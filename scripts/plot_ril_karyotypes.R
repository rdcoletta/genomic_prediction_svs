#### intro ----

# plot karyotypes for RILs of a biparental cross

# code for plot was modified from Carles Hernandez-Ferrer's blog:
# https://carleshf.github.io/chromosome_karyotype_plot_with_ggplot2/


#### libraries ----

library(ggplot2)
library(data.table)


#### get centromeres in v4 coordinates ----

# read in data
v2.to.v4 <- read.delim("data/centromeres_v2-to-v4.bed", sep = "\t", header = FALSE)
colnames(v2.to.v4) <- c("chr", "start", "end")
v3.to.v4 <- read.delim("data/centromeres_v3-to-v4.bed", sep = "\t", header = FALSE)
colnames(v3.to.v4) <- c("chr", "start", "end")

# data frame containing the rows for each chromosome
v3.to.v4.chromosomes <- data.frame(chr = 1:9,
                                   row_start = c(1, 198, 497, 697, 881, 1218, 1348, 1431, 1617),
                                   row_end = c(197, 496, 696, 880, 1217, 1347, 1430, 1616, 1814))
# i don't need to subset v3.to.v4 data by rows because there is only one chr in this df

# create empty data frame to store coordinates of centromeres in v4 assembly
centromeres.v4 <- data.frame(chr = as.numeric(),
                             start_pos = as.numeric(),
                             end_pos = as.numeric(),
                             centro_size = as.numeric())

# get min and max coordinates for each chromosome in v4
for (row in 1:NROW(v3.to.v4.chromosomes)) {
  
  # subset by row
  row_start <- v3.to.v4.chromosomes[row, "row_start"]
  row_end <- v3.to.v4.chromosomes[row, "row_end"]
  chr.info <- v2.to.v4[row_start:row_end,]
  
  # remove coordinates from other chromosomes
  curr.chr <- v3.to.v4.chromosomes[row, "chr"]
  chr.info <- subset(chr.info, subset = chr == curr.chr)
  
  # get min and max values
  v4.start <- min(chr.info$start)
  v4.end <- max(chr.info$end)
  v4.size <- v4.end - v4.start
  
  # write to centromeres df
  centromeres.v4 <- rbind(centromeres.v4, c(curr.chr, v4.start, v4.end, v4.size))
  
}

# get the coordinates for chr 10 now
v3.to.v4 <- subset(v3.to.v4, subset = chr == 10)
v4.start <- min(v3.to.v4$start)
v4.end <- max(v3.to.v4$end)
v4.size <- v4.end - v4.start
centromeres.v4 <- rbind(centromeres.v4, c(10, v4.start, v4.end, v4.size))

# fix column names
colnames(centromeres.v4) <- c("chr", "start_pos", "end_pos", "centro_size")

# write file
write.table(centromeres.v4, file = "data/centromeres_Schneider-2016-pnas_v4.bed",
            sep = "\t", quote = FALSE, row.names = FALSE)



#### load chromosome and centromere positions ----

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



#### load marker positions ----

# select a cross that has a RIL with high AA or BB genotype frequency
RIL.pop.list <- list.dirs(path = "analysis/qc/", full.names = FALSE, recursive = FALSE)
cross <- RIL.pop.list[1]
RIL <- "RIL_3"

# load data with marker positions
markers.filename <- paste0("analysis/qc/", cross, "/recomb-freq_", cross, "_rils.txt")
markers.infile <- fread(markers.filename, header = TRUE, data.table = FALSE)

# load genotypic data for all RILs in the cross
geno.data.filename <- paste0("analysis/qc/", cross, "/geno-data_", cross, "_after_filtering.txt")
geno.data.infile <- fread(geno.data.filename, header = TRUE, data.table = FALSE)

# merge information of RIL of interest with respective marker positions
geno.data <- cbind(markers.infile[, c("marker", "chr", "pos")], geno.data.infile[, RIL])
colnames(geno.data)[4] <- "geno"

# select only AA genotype
geno.data <- subset(geno.data, geno == "AA")



#### plot karyotypes ----

ggplot() +
  geom_segment(data = chrms,
               aes(x = 0, xend = 0, y = start_pos, yend = end_pos),
               lineend = "round", color = "Gainsboro", size = 5) +
  geom_segment(data = centros,
               aes(x = 0, xend = 0, y = start_pos, yend = end_pos),
               lineend = "round", color = "DimGray", size = 5) +
  geom_segment(data = geno.data,
               aes(x = 0, xend = 0, y = pos, yend = pos + 1e6),  # increased size to be able to see the marker
               lineend = "butt", color = "firebrick", size = 5, alpha = 0.3) +
  scale_y_reverse(breaks = seq(0, 3.5e8, 0.50e8), labels = c(1, seq(50, 350, 50))) +
  scale_x_discrete(position = "top") +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.caption = element_text(size = rel(1.1), color = "DimGray")) + 
  facet_grid(cols = vars(chr), switch = "y") +
  labs(caption = paste0(cross, " (", gsub("RIL_", "RIL ", RIL), ")"),
       x = "Chromosomes", y = "Genomic positions (Mb)")
