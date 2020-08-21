# assign arguments to variables
args <- commandArgs(trailingOnly = TRUE)

ld.file <- args[1]
outfile.name <- args[2]

# ld.file <- "analysis/ld/ld_usda_rils_snp-sv_only.chr1.window-100kb.filter-0.25.ld"
# outfile.name <- "analysis/ld/LD-decay_SNPs-SVs.window-100kb.filter-0.25.png"

library(data.table)
library(ggplot2)



#### plot ld decay entire dataset ----

# load data
LD.results <- fread(ld.file, header = TRUE, data.table = FALSE)
for (chr in 2:10) {
  LD.results.chr <- fread(gsub("chr[0-9]+", paste0("chr", chr), ld.file, perl = TRUE), header = TRUE, data.table = FALSE)
  LD.results <- rbind(LD.results, LD.results.chr)
}

# add column with distance between sv and snp
LD.results$dist_to_sv <- abs(LD.results[, 5] - LD.results[, 2])

# plot LD decay for each chromosome (based on markers +/- 2000kb of a SV)
plot.LD.decay <- ggplot(LD.results, aes(x = abs(dist_to_sv), y = R2)) +
  facet_wrap(~CHR_A) +
  geom_point(size = 0.01, alpha = 0.01) +
  scale_x_continuous(labels = function(x) x/1000) +
  labs(title = "LD decay",
       x = "Distance between SNP and SV (kb)",
       y = bquote("LD"~(r^2)))

# save plot
ggsave(filename = outfile.name, plot = plot.LD.decay, device = "png")



#### plot ld decay for 10 random windows ----

# get SVs
SVs <- apply(LD.results, MARGIN = 1, function(marker) {
  sv <- marker[grep("^del|^ins|^dup|^inv|^tra", marker, perl = TRUE)]
  return(as.character(sv))
})
SVs <- unique(unlist(SVs))
# sample SVs
set.seed(7812)
SVs.samples <- sample(SVs, size = 10, replace = FALSE)

for (sv in SVs.samples) {
  
  # filter ld df
  LD.results.sample <- LD.results[which(LD.results$SNP_A == sv | LD.results$SNP_B == sv), ]
  
  # plot ld decay for sample
  plot.LD.decay.sample <- ggplot(LD.results.sample, aes(x = abs(dist_to_sv), y = R2)) +
    facet_wrap(~CHR_A) +
    geom_point(size = 0.7, alpha = 0.5) +
    scale_x_continuous(labels = function(x) x/1000) +
    labs(title = paste0("LD decay (", sv, ")"),
         x = "Distance between SNP and SV (kb)",
         y = bquote("LD"~(r^2)))
  
  # save plot
  sample.number <- grep(sv, SVs.samples)
  outfile.name.sample <- gsub(".png", paste0(".sv", sample.number, ".png"), outfile.name)
  ggsave(filename = outfile.name.sample, plot = plot.LD.decay.sample, device = "png")
  
}







