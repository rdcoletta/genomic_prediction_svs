#### arguments for command line ----

args <- commandArgs(trailingOnly = TRUE)

# help
if (all(length(args) == 1 & args == "-h" | args == "--help")) {
  cat("
      Description: check if there are any differences between the v4 converted probes and the original probes from v2 assembly.

      Usage: Rscript check_refgen_SNPchip_v2tov4_probes.R [file_ref-markers_chr1] [v2-to-v4_probes]")
  quit()
}

# make sure the correct number of arguments are used
if (length(args) != 2) {
  stop("incorrect number of arguments provided.

       Usage: Rscript check_refgen_SNPchip_v2tov4_probes.R [file_ref-markers_chr1] [v2-to-v4_probes]
       ")
}

# assign arguments to variables
v2.file <- args[1]
v4.file <- args[2]

# v2.file <- "data/check_refgen_SNPchip/ref-markers_chr1.txt"
# v4.file <- "data/probes_v2-to-v4/SNP_positions_v2-to-v4_probes-100bp.txt"



#### check strand ----

# read data
v2.pos <- read.table(v2.file, header = TRUE, stringsAsFactors = FALSE)
v4.pos <- read.table(v4.file, header = TRUE, stringsAsFactors = FALSE)

# filter v4 data to have only chromosome 1 markers
v4.pos.chr1 <- v4.pos[which(v4.pos[, "chr_v4"] == 1), ]

# make sure that the same markers are in both datasets (because after correcting SNP positions in v4,
# some SNPs from chr1 went to other chromosome; so these were excluded from this analysis)
v2.pos <- v2.pos[v2.pos[, "pos"] %in% v4.pos.chr1[, "pos_v2"], ]
v4.pos.chr1 <- v4.pos.chr1[v4.pos.chr1[, "pos_v2"] %in% v2.pos[, "pos"], ]

# create a new data frame only with the respective SNPs in v2 and v4
df <- cbind(v2.pos["v2"], v4.pos.chr1["SNP"])
colnames(df) <- c("v2", "v4")

# count how many were different in v2 and v4
df_diff_alleles <- df[which(df$v2 != df$v4), ]

cat(NROW(df_diff_alleles), " marker alleles differ between v2 and v4\n", sep = "")
# 217

if (NROW(df_diff_alleles) > 0) {

  cat("Converting alleles that differ to complementary base...\n")
  # transform one of the columns into the complementary base
  df_diff_alleles$v2 <- sapply(df_diff_alleles$v2, FUN = function(base) {
    if (base == "A") {
      base <- "T"
    } else if (base == "T") {
      base <- "A"
    } else if (base == "C") {
      base <- "G"
    } else if (base == "G") {
      base <- "C"
    }
  })

  # now, check again how many markers are different
  cat(NROW(df[which(df_diff_alleles$v2 != df_diff_alleles$v4), ]), " marker alleles differ between v2 and v4\n", sep = "")
  # 0
}
