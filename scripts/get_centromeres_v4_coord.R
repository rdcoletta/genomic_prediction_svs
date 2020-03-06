#### arguments for command line ----

args <- commandArgs(trailingOnly = TRUE)

# help
if (all(length(args) == 1 & args == "-h" | args == "--help")) {
  cat("
      Description: convert centromere positions from v2/v3 genome assembly to v4.
      
      Usage: Rscript get_centromeres_v4_coord.R [bed_file_v2_to_v4] [bed_file_v3_to_v4] [output_name]")
  quit()
}

# make sure the correct number of arguments are used
if (length(args) != 3) {
  stop("incorrect number of arguments provided.
       
       Usage: Rscript get_centromeres_v4_coord.R [bed_file_v2_to_v4] [bed_file_v3_to_v4] [output_name]
       ")
}

# assign arguments to variables
cent.v2.v4 <- args[1]
cent.v3.v4 <- args[2]
out.file <- args[3]

# cent.v2.v4 <- "data/centromeres_v2-to-v4.bed"
# cent.v3.v4 <- "data/centromeres_v3-to-v4.bed"
# out.file <- "data/centromeres_Schneider-2016-pnas_v4.bed"

#### libraries ----

if(!require("data.table")) install.packages("data.table")


#### get centromeres in v4 coordinates ----

# read in data
v2.to.v4 <- read.delim(cent.v2.v4, sep = "\t", header = FALSE)
colnames(v2.to.v4) <- c("chr", "start", "end")
v3.to.v4 <- read.delim(cent.v3.v4, sep = "\t", header = FALSE)
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
write.table(centromeres.v4, file = out.file, sep = "\t", quote = FALSE, row.names = FALSE)

