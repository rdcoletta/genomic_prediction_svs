# assign arguments to variables
args <- commandArgs(trailingOnly = TRUE)

dist.file <- args[1]
outfile.name <- args[2]

# dist.file <- "analysis/ld/closest_low-missing-data-SNPs_to_SVs.filter-0.25.all.distances.txt"
# outfile.name <- "analysis/ld/closest_low-missing-data-SNPs_to_SVs.filter-0.25.all.png"

library(data.table)
library(ggplot2)

# load data
dist.sv.snp <- fread(dist.file, header = TRUE, data.table = FALSE)

# plot distribution
dist.plot <- ggplot(dist.sv.snp, aes(x = distance)) + 
  geom_histogram(binwidth = 100) + 
  scale_x_continuous(limits = c(0, 10000), n.breaks = 10) +
  labs(x = "Distance of SV to its closest SNP",
       y = "Count")

ggsave(filename = outfile.name, plot = dist.plot, device = "png")

# count how many closest SNPs are within a certain distance of SV
for (distance in c(100, 500, 1000, 5000, 10000, 100000, 1000000)) {
  large.dist <- dist.sv.snp[, "distance"] > distance
  cat(sum(large.dist), " SNPs were more than ", distance/1000,"kb away from a SV\n", sep = "")
}