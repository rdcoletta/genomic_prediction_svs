#arguments for command line

args <- commandArgs(trailingOnly = TRUE)

# assign arguments to variables
eigenvec.file <- args[1]
outfile.pca.plot <- args[2]

# eigenvec.file <- "analysis/ld/pca_plink_usda_closest-snps_all.eigenvec"
# outfile.pca.plot <- "analysis/ld/pca_plink_usda_closest-snps_all.pdf" 


# libraries
library(data.table)
library(ggplot2)
library(RColorBrewer)


# load eigen vectors calculated from plink2
scores <- fread(eigenvec.file, header = TRUE, data.table = FALSE)
rownames(scores) <- scores[, "IID"]
scores <- scores[, 3:NCOL(scores)]

# visualize PCAs
scores.group <- cbind(scores, group = sapply(rownames(scores), function(x) unlist(strsplit(x, split = "-"))[1]))

# plot
getPalette <- colorRampPalette(brewer.pal(9, name = "Set1"))

pca.plot <- ggplot(data=scores.group, aes(x=PC1, y=PC2)) +
  geom_point(aes(color = group)) +
  theme(panel.background = element_blank()) +
  theme(legend.position="top", legend.title=element_blank(), legend.key=element_blank()) +
  theme(axis.line = element_line(colour = "black", size=.25)) +
  theme(axis.text=element_text(size=8)) +
  scale_color_manual(values = getPalette(13))


ggsave(filename = outfile.pca.plot, plot = pca.plot, device = "pdf")


