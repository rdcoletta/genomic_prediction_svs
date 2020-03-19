#### arguments for command line ----

args <- commandArgs(trailingOnly = TRUE)

# help
if (all(length(args) == 1 & args == "-h" | args == "--help")) {
  cat("
      Description: plot karyotypes for RILs of a biparental cross using SVs.
      Credits: code for plot was modified from Carles Hernandez-Ferrer's blog at
      https://carleshf.github.io/chromosome_karyotype_plot_with_ggplot2/
      
      
      Usage: Rscript plot_ril_karyotypes_reseq-SNPs.R [...]")
  quit()
}

# make sure the correct number of arguments are used
if (length(args) != 7) {
  stop("incorrect number of arguments provided.
       
       Usage: Rscript plot_ril_karyotypes_reseq-SNPs.R [...]
       ")
}

# assign arguments to variables
chr.info <- args[1]
cent.info <- args[2]
cross <- args[3]
proj.folder <- args[4]
parents.folder <- args[5]
out.folder <- args[6]

if (args[7] == "--rils=random") {
  random.rils <- TRUE
} else if (grepl("--rils=", args[7])) {
  random.rils <- FALSE
  rils.list <- unlist(strsplit(args[7], split = "="))[2]
  rils.list <- unlist(strsplit(rils.list, split = ","))
} else {
  stop("Invalid list of rils. Make sure it's comma-separated")
}


# chr.info <- "data/B73_RefGen_V4_chrm_info.txt"
# cent.info <- "data/centromeres_Schneider-2016-pnas_v4.bed"
# cross <- "B73xLH82"
# proj.folder <- "analysis/projection_reseq-snps"
# parents.folder <- "data/reseq_snps"
# out.folder <- "analysis/qc/karyotypes"
# random.rils <- FALSE
# rils.list <- c("2011MOB_151")



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
cat("Plotting ", cross, "...\n", sep = "")

# load projected cross
data.filename <- list.files(path = proj.folder,
                            pattern = paste0(cross, ".poly.projected.hmp.txt"),
                            full.names = TRUE)
geno.data.cross <- fread(data.filename, header = TRUE, data.table = FALSE)

# load parental data
parents.filename <- list.files(path = paste0(parents.folder, "/", cross),
                               pattern = paste0("biomAP_parents_SNPs-reseq_SNP-chip.", cross, ".poly.hmp.txt"),
                               full.names = TRUE)
parents.data.cross <- fread(parents.filename, header = TRUE, data.table = FALSE)


# remove SNPs within SVs in parental file
parents.data.cross <- parents.data.cross[which(parents.data.cross[, 1] %in% geno.data.cross[, 1]), ]

# get parents names
parent1 <- unlist(strsplit(cross, "x"))[1]
parent2 <- unlist(strsplit(cross, "x"))[2]

# get parents column numbers in resequencing data
p1.col.reseq <- grep(parent1, colnames(parents.data.cross), fixed = TRUE)
p2.col.reseq <- grep(parent2, colnames(parents.data.cross), fixed = TRUE)


if (random.rils == TRUE) {
  # randomly select 3 RILs per population to plot karyotype
  set.seed(184)
  selected.RILs <- sample(colnames(geno.data.cross[12:NCOL(geno.data.cross)]), size = 3, replace = FALSE)
} else {
  selected.RILs <- rils.list
}

for (RIL in selected.RILs) {
  
  # get ril column number
  ril.col <- grep(RIL, colnames(geno.data.cross))
  
  cat("  RIL", RIL, "\n")
  
  # merge information of RIL of interest with respective marker positions
  geno.data <- cbind(geno.data.cross[, c(1,3,4)], geno.data.cross[, ril.col])
  colnames(geno.data) <- c("marker", "chr", "pos", "geno")
  
  # select only chip SNPs
  parents.SNPs <- parents.data.cross[grep("^snp", parents.data.cross[, 1], perl = TRUE, invert = TRUE), ]
  geno.data.SNPs <- geno.data[grep("^snp", geno.data$marker, perl = TRUE, invert = TRUE), ]
  # add parental info
  geno.data.SNPs <- cbind(geno.data.SNPs, parents.SNPs[, c(p1.col.reseq, p2.col.reseq)])
  colnames(geno.data.SNPs) <- c("marker", "chr", "pos", "geno", "parent1", "parent2")
  
  # select missing data and non-missing data
  geno.data.SNPs.not.missing <- subset(geno.data.SNPs, geno != "NN")
  # get proportion of non-missing data to add in the plot
  prop.SNPs.not.missing <- NROW(geno.data.SNPs.not.missing) / NROW(geno.data.SNPs)
  
  # check from which parent an allele came
  geno.data.SNPs.not.missing.p1 <- geno.data.SNPs.not.missing[which(geno.data.SNPs.not.missing[, "geno"] == geno.data.SNPs.not.missing[, "parent1"] &
                                                                      geno.data.SNPs.not.missing[, "geno"] != geno.data.SNPs.not.missing[, "parent2"]), ]
  geno.data.SNPs.not.missing.p2 <- geno.data.SNPs.not.missing[which(geno.data.SNPs.not.missing[, "geno"] != geno.data.SNPs.not.missing[, "parent1"] &
                                                                      geno.data.SNPs.not.missing[, "geno"] == geno.data.SNPs.not.missing[, "parent2"]), ]
  geno.data.SNPs.not.missing.rest <- geno.data.SNPs.not.missing[which(!geno.data.SNPs.not.missing[, "marker"] %in% geno.data.SNPs.not.missing.p1[, "marker"] &
                                                                        !geno.data.SNPs.not.missing[, "marker"] %in% geno.data.SNPs.not.missing.p2[, "marker"]), ]
  
  # find hets
  geno.data.SNPs.not.missing.het <- data.frame(matrix(nrow = 0, ncol = NCOL(geno.data.SNPs.not.missing.rest)))
  colnames(geno.data.SNPs.not.missing.het) <- colnames(geno.data.SNPs.not.missing.rest)
  for (snp in 1:NROW(geno.data.SNPs.not.missing.rest)) {
    # check if ril snp is a het or if it has a different allele from parents
    p1.alleles <- unlist(strsplit(as.character(geno.data.SNPs.not.missing.rest[snp, "parent1"]), split = ""))
    p2.alleles <- unlist(strsplit(as.character(geno.data.SNPs.not.missing.rest[snp, "parent2"]), split = ""))
    ril.alleles <- unlist(strsplit(as.character(geno.data.SNPs.not.missing.rest[snp, "geno"]), split = ""))
    # make sure to exclude hets in parental calls
    if (length(unique(p1.alleles)) == 1 & length(unique(p2.alleles)) == 1) {
      if (ril.alleles[1] != ril.alleles[2] & unique(p1.alleles) %in% ril.alleles & unique(p2.alleles) %in% ril.alleles) {
        geno.data.SNPs.not.missing.het <- rbind(geno.data.SNPs.not.missing.het, geno.data.SNPs.not.missing.rest[snp, ])
      }
    }
  }
  
  
  # plot karyotypes before SV projection (i.e. only parental SNPs)
  karyo.plot.before <- ggplot() +
    geom_segment(data = chrms,
                 aes(x = 0, xend = 0, y = start_pos, yend = end_pos),
                 lineend = "round", color = "Gainsboro", size = 5) +
    geom_segment(data = centros,
                 aes(x = 0, xend = 0, y = start_pos, yend = end_pos),
                 lineend = "round", color = "DimGray", size = 5) +
    geom_segment(data = geno.data.SNPs.not.missing.p1,
                 aes(x = 0, xend = 0, y = pos, yend = pos + 1e6),  # increased size to be able to see the marker
                 lineend = "butt", color = "firebrick", size = 5, alpha = 0.3) +
    geom_segment(data = geno.data.SNPs.not.missing.p2,
                 aes(x = 0, xend = 0, y = pos, yend = pos + 1e6),  # increased size to be able to see the marker
                 lineend = "butt", color = "#386cb0", size = 5, alpha = 0.3) +
    geom_segment(data = geno.data.SNPs.not.missing.het,
                 aes(x = 0, xend = 0, y = pos, yend = pos + 1e6),  # increased size to be able to see the marker
                 lineend = "butt", color = "#fec44f", size = 5, alpha = 0.5) +
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
    facet_grid(~chr, switch = "y") +
    labs(caption = paste0(cross, " - ", gsub("RIL_", "RIL ", RIL), "\n\n",
                          "Not missing: ", round(prop.SNPs.not.missing, digits = 3), "\n"),
         x = "Chromosomes", y = "Genomic positions (Mb)")
  
  
  # select all SVs
  parents.SVs <- parents.data.cross[grep("^snp", parents.data.cross[, 1], perl = TRUE), ]
  geno.data.all.SVs <- geno.data[grep("^snp", geno.data$marker, perl = TRUE), ]
  # add parental info
  geno.data.all.SVs <- cbind(geno.data.all.SVs, parents.SVs[, c(p1.col.reseq, p2.col.reseq)])
  colnames(geno.data.all.SVs) <- c("marker", "chr", "pos", "geno", "parent1", "parent2")
  
  # select only projected SVs
  geno.data.all.SVs.proj <- subset(geno.data.all.SVs, geno != "NN")
  geno.data.all.SVs.not.proj <- subset(geno.data.all.SVs, geno == "NN")
  
  # check from which parent an allele came
  geno.data.all.SVs.proj.p1 <- geno.data.all.SVs.proj[which(geno.data.all.SVs.proj[, "geno"] == geno.data.all.SVs.proj[, "parent1"] &
                                                              geno.data.all.SVs.proj[, "geno"] != geno.data.all.SVs.proj[, "parent2"] &
                                                              geno.data.all.SVs.proj[, "parent2"] != "NN"), ]
  geno.data.all.SVs.proj.p2 <- geno.data.all.SVs.proj[which(geno.data.all.SVs.proj[, "geno"] != geno.data.all.SVs.proj[, "parent1"] &
                                                              geno.data.all.SVs.proj[, "geno"] == geno.data.all.SVs.proj[, "parent2"] &
                                                              geno.data.all.SVs.proj[, "parent1"] != "NN"), ]
  geno.data.all.SVs.proj.rest <- geno.data.all.SVs.proj[which(!geno.data.all.SVs.proj[, "marker"] %in% geno.data.all.SVs.proj.p1[, "marker"] &
                                                                !geno.data.all.SVs.proj[, "marker"] %in% geno.data.all.SVs.proj.p2[, "marker"]), ]
  
  # find hets
  geno.data.all.SVs.proj.het <- apply(geno.data.all.SVs.proj.rest, MARGIN = 1, function(snp) {
    # check if ril snp is a het or if it has a different allele from parents
    p1.alleles <- unlist(strsplit(as.character(snp["parent1"]), split = ""))
    p2.alleles <- unlist(strsplit(as.character(snp["parent2"]), split = ""))
    ril.alleles <- unlist(strsplit(as.character(snp["geno"]), split = ""))
    # make sure to exclude hets in parental calls
    if (length(unique(p1.alleles)) == 1 & length(unique(p2.alleles)) == 1) {
      if (ril.alleles[1] != ril.alleles[2] & unique(p1.alleles) %in% ril.alleles & unique(p2.alleles) %in% ril.alleles) {
        return(snp)
        geno.data.all.SVs.proj.het <- rbind(geno.data.all.SVs.proj.het, geno.data.all.SVs.proj.rest[snp, ])
      }
    }
  })
  geno.data.all.SVs.proj.het <- data.frame(t(geno.data.all.SVs.proj.het), stringsAsFactors = FALSE)
  geno.data.all.SVs.proj.het$chr <- as.numeric(geno.data.all.SVs.proj.het$chr)
  geno.data.all.SVs.proj.het$pos <- as.numeric(geno.data.all.SVs.proj.het$pos)
  
  if (NROW(geno.data.all.SVs.proj.p1) == 0) {
    geno.data.all.SVs.proj.p1  <- rbind(geno.data.all.SVs.proj.p1, list(marker = NA, chr = NA, pos = NA, geno = NA, parent1 = NA, parent2 = NA))
  }
  if (NROW(geno.data.all.SVs.proj.p2) == 0) {
    geno.data.all.SVs.proj.p2  <- rbind(geno.data.all.SVs.proj.p2, list(marker = NA, chr = NA, pos = NA, geno = NA, parent1 = NA, parent2 = NA))
  }
  if (NROW(geno.data.all.SVs.proj.het) == 0) {
    geno.data.all.SVs.proj.het  <- rbind(geno.data.all.SVs.proj.het, list(marker = NA, chr = NA, pos = NA, geno = NA, parent1 = NA, parent2 = NA))
  }
  
  
  # plot karyotypes after SV projection 
  karyo.plot.after <- ggplot() +
    geom_segment(data = chrms,
                 aes(x = 0, xend = 0, y = start_pos, yend = end_pos),
                 lineend = "round", color = "Gainsboro", size = 5) +
    geom_segment(data = centros,
                 aes(x = 0, xend = 0, y = start_pos, yend = end_pos),
                 lineend = "round", color = "DimGray", size = 5) +
    geom_segment(data = geno.data.all.SVs.proj.p1,
                 aes(x = 0, xend = 0, y = pos, yend = pos + 1e6),  # increased size to be able to see the marker
                 lineend = "butt", color = "#008837", size = 5, alpha = 0.3) +
    geom_segment(data = geno.data.all.SVs.proj.p2,
                 aes(x = 0, xend = 0, y = pos, yend = pos + 1e6),  # increased size to be able to see the marker
                 lineend = "butt", color = "#7b3294", size = 5, alpha = 0.3) +
    geom_segment(data = geno.data.all.SVs.proj.het,
                 aes(x = 0, xend = 0, y = pos, yend = pos + 1e6),  # increased size to be able to see the marker
                 lineend = "butt", color = "#fec44f", size = 5, alpha = 0.5) +
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
    facet_grid(~chr, switch = "y") +
    labs(caption = paste0(cross, " - ", gsub("RIL_", "RIL ", RIL), "\n\n",
                          "Projected: ", round(NROW(geno.data.all.SVs.proj)/NROW(geno.data.all.SVs), digits = 3), "\n"),
         x = "Chromosomes", y = "Genomic positions (Mb)")
  
  # create directory for output if it doesn't exist
  if (!dir.exists(out.folder)) {
    dir.create(out.folder, recursive = TRUE)
  }
  
  # save plots
  karyo.name.before <- paste0(out.folder, "/", cross, "_", RIL,"_before-proj_reseq-snps.png")
  karyo.name.after <- paste0(out.folder, "/", cross, "_", RIL,"_after-proj_reseq-snps.png")
  ggsave(filename = karyo.name.before, plot = karyo.plot.before, device = "png")
  ggsave(filename = karyo.name.after, plot = karyo.plot.after, device = "png")
  
}

cat("Done!\n\n")






#### credits ----

# code for karyotype plot was adapted from Carles Hernandez-Ferrer's blog:
# https://carleshf.github.io/chromosome_karyotype_plot_with_ggplot2/
# assign arguments to variables
