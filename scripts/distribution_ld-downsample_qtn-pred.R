library(data.table)
library(ggplot2)
library(dplyr)

usage <- function() {
  cat("
description: plot LD distribution between QTLs and predictors that have been downsampled according to their LD level.

usage: Rscript distribution_ld-downsample_qtn-pred.R [low_ld_file] [moderate_ld_file] [high_ld_file]

positional arguments:
  low_ld_file                path to low LD file
  moderate_ld_file           path to moderate LD file
  high_ld_file               path to high LD file
  qtn_list_file              path to file with list of QTNs
  plot_name                  name of plot

optional argument:
  --help                        show this helpful message

"
  )
}



#### command line options ----

args <- commandArgs(trailingOnly = TRUE)
if ("--help" %in% args) usage() & q(save = "no")

# assert to have the correct optional arguments
if (length(args) != 5) stop(usage(), "missing positional argument(s)")

# get positional arguments
low_ld_file <- args[1]
moderate_ld_file <- args[2]
high_ld_file <- args[3]
qtn_list_file <- args[4]
plot_name <- args[5]


#### plot ld distribution ----

# load LD data
ld_results <- data.frame()
for (ld in c("Low", "Moderate", "High")) {

  if (ld == "Low") ld_file <- low_ld_file
  if (ld == "Moderate") ld_file <- moderate_ld_file
  if (ld == "High") ld_file <- high_ld_file

  ld_file <- fread(ld_file, header = TRUE, data.table = FALSE)
  ld_results <- rbind(ld_results, cbind(ld_file, ld_category = ld))

}
rm(ld_file, ld)

# plot
ld_plot <- ggplot(ld_results, aes(x = R2)) +
  geom_histogram(fill = "#900721", binwidth = 0.01) +
  facet_grid(ld_category ~ .) +
  labs(title = paste0("LD between QTLs and predictors"),
       x = bquote("LD"~(r^2)),
       y = "Count") +
  coord_cartesian(xlim = c(0, 1)) +
  theme_bw()

ggsave(filename = paste0(plot_name, ".pdf"), plot = ld_plot, device = "pdf")


# load QTN list
qtn_list <- fread(qtn_list_file, header = FALSE, data.table = FALSE)
qtn_list <- as.character(qtn_list[, 1])

# create empty data frame
ld_results_highest <- data.frame(matrix(nrow = 0, ncol = NCOL(ld_results)), stringsAsFactors = FALSE)
colnames(ld_results_highest) <- colnames(ld_results)
# create vector to make sure the same snp doesn't get picked twice
snps_already_in_ld <- c()

# get markers in highest ld to each QTL
for (ld in unique(ld_results$ld_category)) {
  
  # subset by ld category
  ld_results_category <- subset(ld_results, ld_category == ld)
  
  # get closest (highest LD) snps
  for (qtn in qtn_list) {
    
    # subset LD results to have only the qtn being parsed
    snps_LD_with_sv <- ld_results_category[which(ld_results_category[, "SNP_A"] == qtn | ld_results_category[, "SNP_B"] == qtn ), ]
    
    if (length(snps_already_in_ld) > 0) {
      snps_LD_with_sv <- snps_LD_with_sv[which(!snps_LD_with_sv[, "SNP_A"] %in% snps_already_in_ld & !snps_LD_with_sv[, "SNP_B"] %in% snps_already_in_ld), ]
    }
    
    
    if (NROW(snps_LD_with_sv) > 0) {
      
      # select only SNP with highest LD with that qtn
      snps_LD_with_sv <- snps_LD_with_sv[which(snps_LD_with_sv[, "R2"] == max(snps_LD_with_sv[, "R2"])), ]
      # if there are more than one SNP with the same R2, get the first one
      # (in this case, it doesn't really matter which marker is in LD, only which is the maximum R2 value)
      if (NROW(snps_LD_with_sv) > 1) {
        snps_LD_with_sv <- snps_LD_with_sv[1, ]
      }
      # get marker name
      snp_selected <- apply(snps_LD_with_sv, MARGIN = 1, function(row) {
        marker1 <- row["SNP_A"]
        marker2 <- row["SNP_B"]
        if (!grepl(qtn, marker1, perl = TRUE)) {
          return(marker1)
        } else {
          return(marker2)
        }
      })
      snps_already_in_ld <- append(snps_already_in_ld, as.character(snp_selected))
      # add closest SNP in LD with qtn into new df
      ld_results_highest <- rbind(ld_results_highest, snps_LD_with_sv)
      
    }
  }
}
rm(ld_results_category, snps_already_in_ld, snps_LD_with_sv)

ld_plot_highest <- ggplot(ld_results_highest, aes(x = R2)) +
  geom_histogram(fill = "#900721", binwidth = 0.01) +
  facet_grid(ld_category ~ .) +
  labs(title = paste0("Predictors in highest LD to QTLs"),
       x = bquote("LD"~(r^2)),
       y = "Count") +
  coord_cartesian(xlim = c(0, 1)) +
  theme_bw()

ggsave(filename = paste0(plot_name, ".highest-ld.pdf"), plot = ld_plot_highest, device = "pdf")



#### plot MAF distribution ----

maf_qtns <- data.frame(maf = rep(NA, NROW(ld_results_highest)), ld_category = rep(NA, NROW(ld_results_highest)))
maf_qtns$maf <- c(ld_results_highest[which(ld_results_highest$SNP_A %in% qtn_list), "MAF_A"],
                  ld_results_highest[which(ld_results_highest$SNP_B %in% qtn_list), "MAF_B"])
maf_qtns$ld_category <- c(as.character(ld_results_highest[which(ld_results_highest$SNP_A %in% qtn_list), "ld_category"]),
                          as.character(ld_results_highest[which(ld_results_highest$SNP_B %in% qtn_list), "ld_category"]))
maf_qtns$ld_category <- factor(maf_qtns$ld_category, levels = c("Low", "Moderate", "High"))
plot_maf_qtns <-  ggplot(maf_qtns, aes(x = maf)) +
  # geom_histogram(fill = "#900721", binwidth = 0.01) +
  geom_density(fill = "#900721") +
  facet_grid(ld_category ~ .) +
  labs(title = paste0("Predictors in highest LD to QTLs"),
       x = "MAF QTLs",
       y = "Density") +
  coord_cartesian(xlim = c(0, 0.5)) +
  theme_bw()

ggsave(filename = paste0(plot_name, ".maf-qtls.pdf"), plot = plot_maf_qtns, device = "pdf")


maf_preds <- data.frame(maf = rep(NA, NROW(ld_results_highest)), ld_category = rep(NA, NROW(ld_results_highest)))
maf_preds$maf <- c(ld_results_highest[which(!ld_results_highest$SNP_A %in% qtn_list), "MAF_A"],
                  ld_results_highest[which(!ld_results_highest$SNP_B %in% qtn_list), "MAF_B"])
maf_preds$ld_category <- c(as.character(ld_results_highest[which(!ld_results_highest$SNP_A %in% qtn_list), "ld_category"]),
                          as.character(ld_results_highest[which(!ld_results_highest$SNP_B %in% qtn_list), "ld_category"]))
maf_preds$ld_category <- factor(maf_preds$ld_category, levels = c("Low", "Moderate", "High"))
plot_maf_preds <-  ggplot(maf_preds, aes(x = maf)) +
  # geom_histogram(fill = "#900721", binwidth = 0.01) +
  geom_density(fill = "#900721") +
  facet_grid(ld_category ~ .) +
  labs(title = paste0("Predictors in highest LD to QTLs"),
       x = "MAF predictors",
       y = "Density") +
  coord_cartesian(xlim = c(0, 0.5)) +
  theme_bw()

ggsave(filename = paste0(plot_name, ".maf-predictors.pdf"), plot = plot_maf_preds, device = "pdf")




#### debug ----

# low_ld_file <- "analysis/ld_downsample/pred_low/rep1/100-QTNs_from_SNP/pop1/qtn-pred.filtered.ld"
# moderate_ld_file <- "analysis/ld_downsample/pred_moderate/rep1/100-QTNs_from_SNP/pop1/qtn-pred.filtered.ld"
# high_ld_file <- "analysis/ld_downsample/pred_high/rep1/100-QTNs_from_SNP/pop1/qtn-pred.filtered.ld"
# qtn_list_file <- "analysis/ld_downsample/pred_low/rep1/100-QTNs_from_SNP/pop1/QTN_list.txt"
# plot_name <- "analysis/ld_downsample/ld_dist.rep1.pop1.100-QTNs_from_SNP"
