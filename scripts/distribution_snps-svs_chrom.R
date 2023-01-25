#### arguments for command line ----

args <- commandArgs(trailingOnly = TRUE)

# # help
# if (all(length(args) == 1 & args == "-h" | args == "--help")) {
#   cat("
# Description: this script merges SNPs and SVs hapmap files from usda parents to be used in
#              Tassel 5 when projecting SVs into RILs
#
# Usage: ")
#   quit()
# }

# make sure the correct number of arguments are used
# you should provide 2 arguments
if (length(args) != 2) {
  stop("incorrect number of arguments provided.
       
       Usage:
       ")
}

# assign arguments to variables
marker_info_file <- args[1]
plot_name <- args[2]


# file wiht format NAME CHR POS
# marker_info_file = "analysis/ld/window-100kb_filter-0.25/marker_info_highest-ld.txt"
# marker_info_file = "data/positions_projected_snps-svs.txt"
# plot_name <- "analysis/ld/window-100kb_filter-0.25/distribution_snps-svs_chrom.png"





#### libraries ----

if(!require("data.table")) install.packages("data.table")
if(!require("ggplot2")) install.packages("ggplot2")



#### plot ----

# # get the correct genotypic data
# marker_data <- data.frame(stringsAsFactors = FALSE)
# for (chr in 1:10) {
#   marker_data_chr <- gsub("chr[0-9]+", paste0("chr", chr), marker_info_file, perl = TRUE)
#   marker_data_chr <- fread(marker_data_chr, header = FALSE, data.table = FALSE)
#   marker_data <- rbind(marker_data, marker_data_chr)
# }
# marker_data_chr <- NULL

# get the correct genotypic data
# marker_data <- fread(marker_info_file, header = FALSE, data.table = FALSE)
marker_data <- fread(marker_info_file, header = TRUE, data.table = FALSE)


# add type of markers into main df
marker_type <- apply(marker_data, MARGIN = 1, function(marker) {
  
  if (grepl("^del|^dup|^ins|^inv|^tra", marker[1], perl = TRUE)) {
    sv_type <- unlist(strsplit(marker[1], split = ".", fixed = TRUE))[1]
    return(c("SV", toupper(sv_type)))
  } else {
    return(c("SNP", "SNP"))
  }
  
})
marker_data <- cbind(marker_data, t(marker_type), stringsAsFactors = FALSE)
colnames(marker_data) <- c("id", "chr", "pos", "marker", "type")

# plot by marker
dist_hist <- ggplot(marker_data, aes(x = pos, fill = marker)) + 
  geom_histogram(position = "identity", binwidth = 5000000, alpha = 0.4) +
  facet_wrap(~chr, scales = "free_x") +
  scale_x_continuous(labels = function(x) x/1000000) +
  scale_color_manual(values = c("firebrick", "black"), aesthetics = "fill") +
  labs(x = "Position (Mb)",
       y = "Count",
       fill = "Marker type")

ggsave(plot = dist_hist, filename = gsub("png", "marker.hist.png", plot_name), device = "png")

dist_densi <- ggplot(marker_data, aes(x = pos, fill = marker)) +
  geom_density(alpha = 0.5, position = "identity") +
  facet_wrap(~chr, scales = "free_x", nrow = 2) +
  scale_x_continuous(labels = function(x) x/1000000) +
  scale_y_continuous(labels = function(y) y*1000000) +
  scale_color_manual(values = c("firebrick", "black"), aesthetics = "fill") +
  labs(x = "Position (Mb)",
       y = expression(paste("Density (x", 10^{6}, ")")),
       fill = "Marker type") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggsave(plot = dist_densi, filename = gsub("png", "marker.densi.png", plot_name), device = "png")



# plot by type
dist_hist_type <- ggplot(marker_data, aes(x = pos, fill = type)) + 
  geom_histogram(position = "identity", binwidth = 5000000, alpha = 0.4) +
  facet_wrap(~chr, scales = "free_x") +
  scale_x_continuous(labels = function(x) x/1000000) +
  scale_color_manual(values = c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3"), aesthetics = "fill") +
  labs(x = "Position (Mb)",
       y = "Count",
       fill = "Marker type")

ggsave(plot = dist_hist_type, filename = gsub("png", "type.hist.png", plot_name), device = "png")


dist_densi_type <- ggplot(marker_data, aes(x = pos, fill = type)) + 
  geom_density(alpha = 0.5, position = "identity") +
  facet_wrap(~chr, scales = "free_x") +
  scale_x_continuous(labels = function(x) x/1000000) +
  scale_color_manual(values = c("#e41a1c", "#377eb8", "#4daf4a", "#984ea3"), aesthetics = "fill") +
  labs(x = "Position (Mb)",
       y = "Count",
       fill = "Marker type")

ggsave(plot = dist_densi_type, filename = gsub("png", "type.densi.png", plot_name), device = "png")
