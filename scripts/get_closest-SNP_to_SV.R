# assign arguments to variables
args <- commandArgs(trailingOnly = TRUE)

hmp.file <- args[1]
outfile.name <- args[2]
missing.threshold <- as.numeric(args[3])
sv.type <- args[4]

# hmp.file <- "data/usda_SNPs-SVs_rils.not-in-SVs.projected.chr10.reseq-SNPs.hmp.txt"
# outfile.name <- "analysis/ld/closest_low-missing-data-SNPs_to_SVs.chr-10.txt"
# missing.threshold <- 0.25
# sv.type <- "all"

library(data.table)

# load data
hmp.chr <- fread(hmp.file, header = TRUE, data.table = FALSE)

# get list of svs
SVs <- hmp.chr[grep("^del|^ins|^dup|^inv|^tra", hmp.chr[, 1], perl = TRUE), 1]
# remove translocations
SVs <- SVs[grep("^tra", SVs, perl = TRUE, invert = TRUE)]

# keep only low missing data SNPs
SNPs.low.missing <- apply(hmp.chr[grep("^del|^ins|^dup|^inv|^tra", hmp.chr[, 1], perl = TRUE, invert = TRUE), ], MARGIN = 1, function(snp) {
  missing <- sum(snp[12:length(snp)] == "NN") / length(snp[12:length(snp)])
  if (missing < missing.threshold) return(as.character(snp[1]))
})
SNPs.low.missing <- do.call(c, SNPs.low.missing)

if (sv.type == "all") SVs <- SVs
if (sv.type == "del") SVs <- SVs[grep("^del", SVs, perl = TRUE)]
if (sv.type == "ins") SVs <- SVs[grep("^ins", SVs, perl = TRUE)]
if (sv.type == "inv") SVs <- SVs[grep("^inv", SVs, perl = TRUE)]
if (sv.type == "dup") SVs <- SVs[grep("^dup", SVs, perl = TRUE)]

# filter hmp to have only SVs and low missing data SNPs
hmp.chr.filtered <- hmp.chr[which(hmp.chr[, 1] %in% SVs | hmp.chr[, 1] %in% SNPs.low.missing), ]

# keep only 10 SNPs closer to SVs -- reduce total number of SNPs
match.idx <- which(hmp.chr.filtered[, 1] %in% SVs)
span <- seq(from = -10, to = 10)
extend.idx <- c(outer(match.idx, span, `+`))
extend.idx <- Filter(function(i) i > 0 & i < NROW(hmp.chr.filtered), extend.idx)
extend.idx <- sort(unique(extend.idx))
hmp.chr.filtered <- hmp.chr.filtered[extend.idx, ]


# create list to store snps to keep
snps.to.keep <- c()

for (sv in SVs) {

  # make sure to don't repeat snp
  hmp.chr.filtered <- hmp.chr.filtered[which(!hmp.chr.filtered[, 1] %in% snps.to.keep), ]

  # keep only that sv and all SNPs
  sv.row <- grep(sv, hmp.chr.filtered[, 1])
  snp.rows <- grep("^del|^ins|^dup|^inv|^tra", hmp.chr.filtered[, 1], perl = TRUE, invert = TRUE)
  rows.to.keep <- sort(c(sv.row, snp.rows))
  hmp.chr.filtered.sv.snps <- hmp.chr.filtered[rows.to.keep, ]

  # get sv info
  sv.row <- which(hmp.chr.filtered.sv.snps[, 1] == sv)
  sv.start <- as.numeric(unlist(strsplit(sv, split = ".", fixed = TRUE))[3])
  sv.end <- as.numeric(unlist(strsplit(sv, split = ".", fixed = TRUE))[4])

  # check the distance of up and downstream snps to sv
  snp.upstream <- hmp.chr.filtered.sv.snps[sv.row - 1, ]
  snp.downstream <- hmp.chr.filtered.sv.snps[sv.row + 1, ]

  if (NROW(snp.upstream) > 0) {
    snp.dist.up <- abs(sv.start - snp.upstream[1, 4])
  } else {
    # set to a very big number
    snp.dist.up <- 9999999999999
  }

  if (NROW(snp.downstream) > 0) {
    snp.dist.down <- abs(sv.end - snp.downstream[1, 4])
  } else {
    # set to a very big number
    snp.dist.down <- 9999999999999
  }

  # keep closest snp
  if (snp.dist.up < snp.dist.down) {
    snps.to.keep <- append(snps.to.keep, snp.upstream[1, 1])
  } else {
    snps.to.keep <- append(snps.to.keep, snp.downstream[1, 1])
  }

}

# write file
snps.to.keep <- data.frame(t(as.matrix(snps.to.keep)), stringsAsFactors = FALSE)
outfile.name <- gsub("txt", paste0(sv.type, ".txt"), outfile.name)
fwrite(snps.to.keep, outfile.name, sep = " ", row.names = FALSE, col.names = FALSE, quote = FALSE)
