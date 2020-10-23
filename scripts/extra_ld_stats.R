#### arguments for command line ----

args <- commandArgs(trailingOnly = TRUE)

# help
if (all(length(args) == 1 & args == "-h" | args == "--help")) {
  cat("
Description: 

Usage: ")
  quit()
}

# make sure the correct number of arguments are used
# you should provide 1 argument
if (length(args) != 1) {
  stop("incorrect number of arguments provided.

Usage: [...]
       ")
}


# assign arguments to variables
out.dir.ld <- args[1]

# out.dir.ld <- "analysis/ld/window-1000kb_filter-0.25"



#### libraries ----

if(!require("data.table")) install.packages("data.table")
if(!require("ggplot2")) install.packages("ggplot2")



###### other summaries ----

ld_files_to_plot <- list.files(path = out.dir.ld, pattern = "highest-ld.ld", full.names = TRUE)

LD_results_all <- data.frame()
for (chr in 1:10) {
  ld_file_chr <- ld_files_to_plot[grep(paste0("chr", chr, "."), ld_files_to_plot, fixed = TRUE)]
  ld_file_chr <- fread(ld_file_chr, header = TRUE, data.table = FALSE)
  LD_results_all <- rbind(LD_results_all, ld_file_chr)
}


LD_results_all$ld_quarter <- factor(LD_results_all$ld_quarter,
                                    levels = c("0_to_0.25", "0.25_to_0.5", "0.5_to_0.75", "0.75_to_1"))

LD_results_all$size_range <- factor(LD_results_all$size_range,
                                    levels = c("<10kb", "10kb-100kb", "100kb-1Mb", ">1Mb"))


# distribution of r2 of SNPs in LD with SVs by sv_type
for (type in c("del", "inv", "dup")) {
  
  sv_filter <- sort(unique(c(grep(type, LD_results_all[, "SNP_A"]), grep(type, LD_results_all[, "SNP_B"]))))
  LD_results_SVtype <- LD_results_all[sv_filter, ]
  
  plot_dist_SVtype <- ggplot(LD_results_SVtype, aes(x = R2)) +
    geom_histogram(fill = "#900721", binwidth = 0.01) +
    labs(title = paste0("LD between SVs and SNPs"),
         subtitle = paste0("(only ", type, ")"),
         x = bquote("LD"~(r^2)),
         y = "Count") +
    coord_cartesian(xlim = c(0, 1)) +
    theme(title = element_text(size = 15),
          axis.title = element_text(size = 20),
          axis.text = element_text(size = 15))
  
  ggsave(plot_dist_SVtype, filename = paste0(out.dir.ld, "/dist-highest-LD_SNPs-SVs_", type, "-only.png"), device = "png")
  
}

# distribution sv sizes
plot_SVsize <- ggplot(LD_results_all, aes(x = size_range)) + 
  geom_histogram(stat = "count", fill = "#900721") +
  labs(title = "Distribution SV sizes",
       y = "Count") +
  theme(title = element_text(size = 15),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 15))

ggsave(plot_SVsize, filename = paste0(out.dir.ld, "/dist_SV-sizes.png"), device = "png")


# plot frequency of SVs by LD quarter
LD_results_MAFbyLD <- apply(LD_results_all, MARGIN = 1, function(row) {
  
  if (grepl("^del|^dup|^ins|^inv|^tra", row["SNP_A"], perl = TRUE)) {
    maf <- row["MAF_A"]
  } else {
    maf <- row["MAF_B"]
  }
  
  return(c(maf, row["ld_quarter"]))
  
})
LD_results_MAFbyLD <- data.frame(t(LD_results_MAFbyLD), stringsAsFactors = FALSE)
colnames(LD_results_MAFbyLD) <- c("sv_maf", "ld_quarter")

LD_results_MAFbyLD$sv_maf <- as.numeric(LD_results_MAFbyLD$sv_maf)
LD_results_MAFbyLD$ld_quarter <- factor(LD_results_MAFbyLD$ld_quarter, 
                                        levels = c("0_to_0.25", "0.25_to_0.5", "0.5_to_0.75", "0.75_to_1"))

plot_SVmaf_LD <- ggplot(LD_results_MAFbyLD, aes(x = ld_quarter, y = sv_maf)) + 
  geom_boxplot() +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = bquote("LD"~(r^2)),
       y = "SV MAF") +
  theme(title = element_text(size = 15),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 15))

ggsave(plot_SVmaf_LD, filename = paste0(out.dir.ld, "/dist_SV-MAF_by_LD.png"), device = "png")

plot_SVmaf_LD_densi <- ggplot(LD_results_MAFbyLD, aes(x = sv_maf, color = ld_quarter)) + 
  geom_density() +
  scale_x_continuous(limits = c(0, 0.5)) +
  labs(x = bquote("LD"~(r^2)),
       y = "SV MAF") +
  theme(title = element_text(size = 15),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 15))

ggsave(plot_SVmaf_LD_densi, filename = paste0(out.dir.ld, "/dist_SV-MAF_by_LD_densi.png"), device = "png")


# distribution sv sizes and distance to sv by ld quarter
plot_SVsize_by_LD <- ggplot(LD_results_all[which(LD_results_all$sv_size < 100000), ], aes(x = ld_quarter, y = sv_size)) + 
  geom_boxplot()+
  scale_x_discrete(drop = FALSE) +
  labs(x = bquote("LD"~(r^2)),
       y = "SV size (<100kb)") +
  theme(title = element_text(size = 15),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 15))

ggsave(plot_SVsize_by_LD, filename = paste0(out.dir.ld, "/dist_SV-sizes_by_LD.png"), device = "png")

plot_distSV_by_LD <- ggplot(LD_results_all[which(LD_results_all$dist_to_sv < 1000000), ], aes(x = ld_quarter, y = dist_to_sv)) + 
  geom_boxplot() +
  scale_x_discrete(drop = FALSE) +
  labs(x = bquote("LD"~(r^2)),
       y = "Distance SNP to SV (<1Mb)") +
  theme(title = element_text(size = 15),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 15))

ggsave(plot_distSV_by_LD, filename = paste0(out.dir.ld, "/dist_dist-to-SV_by_LD.png"), device = "png")
