library(data.table)
library(ggplot2)
library(dplyr)
library(tidyr)
# install.packages("ggh4x")
# # for R.3.6 compatibility
# remotes::install_github("teunbrand/ggh4x@4e242295da442e0e412c7c3e713a383f603ceda0")  
library(ggh4x)

usage <- function() {
  cat("
description: plot prediction accuracy results after k-fold cross validation of simulated traits as a heatmap.

usage: Rscript plot_prediction_accuracy_heatmap.R [prediction_summary_file] [...]

positional arguments:
  prediction_summary_file       path to file with prediction summary

optional argument:
  --help                        show this helpful message

"
  )
}

getArgValue <- function(arg) {
  
  # get value from command-line arguments
  arg <- unlist(strsplit(arg, "="))
  if (length(arg) == 1) return(TRUE)
  if (length(arg) > 1) return(arg[2])
  
}



#### command line options ----

args <- commandArgs(trailingOnly = TRUE)
if ("--help" %in% args) usage() & q(save = "no")

# assert to have the correct optional arguments
if (length(args) != 1) stop(usage(), "wrong number of positional argument(s)")

# get positional arguments
prediction_summary_file <- args[1]



#### load prediction summary ----

# do causal variatnt 'both' separately
prediction_summary <- fread(prediction_summary_file, header = TRUE, data.table = FALSE)

# create output folder
outfolder <- unlist(strsplit(prediction_summary_file, "/"))
outfolder[length(outfolder)] <- paste0("plots_", gsub(".txt", "", outfolder[length(outfolder)]))
outfolder[length(outfolder)] <- gsub(".", "_", outfolder[length(outfolder)], fixed = TRUE)
outfolder <- paste0(outfolder, collapse = "/")
if (!dir.exists(outfolder)) dir.create(outfolder)



#### plot SNP only or SV only as causative variants ----

# set default
trait_var_source <- c("SNP", "SV")
trait_qtn_number <- c(10, 100)
trait_heritability <- c(0.3 ,0.7)

# create heatmap
results_plot <- prediction_summary %>% 
  filter(h2 %in% trait_heritability, var %in% trait_var_source, qtn %in% trait_qtn_number) %>%
  mutate(predictor = factor(predictor,
                            levels = c("all", "snp", "sv", "snp_ld", "snp_not_ld"),
                            labels = c("SNPs and SVs", "Only SNPs", "Only SVs",
                                       "SNPs in LD", "SNPs not in LD")),
         h2 = factor(h2, levels = c("0.3", "0.7"), labels = c("h2 = 0.3", "h2 = 0.7")),
         qtn = factor(qtn, levels = c("10", "100"), labels = c("10 QTLs", "100 QTLs"))) %>% 
  select(h2, var, qtn, predictor, cv, mean_accuracy_envs) %>%
  unite(h2:qtn, col = "sim_scenarios", sep = "/", remove = FALSE) %>% 
  unite(predictor:cv, col = "pred_scenarios", sep = "/", remove = FALSE) %>%
  ggplot(aes(x = pred_scenarios, y = sim_scenarios, fill = mean_accuracy_envs)) +
  facet_nested(h2 + qtn + var ~ predictor + cv, scales = "free", switch = "y",
               nest_line = element_line(linetype = 1, color = "gray80")) +
  geom_tile(aes(height = 1.5, width = 1.5)) +
  scale_fill_viridis_c(option = "B", limits = c(0, 1), name = "Accuracy") +
  theme_minimal() +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        strip.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        panel.grid = element_blank(),
        panel.spacing = unit(0.25, "lines"))

# save plot
plot_name <- paste0(outfolder, "/heatmap_pred-accuracy_vars-snp-sv.pdf")
ggsave(plot = results_plot, filename = plot_name, device = "pdf",
       width = 14, height = 8, units = "in")




#### plot both SNPs and SVs as causative variants ----

# set default
trait_var_source <- "both"
trait_qtn_number <- c(10, 100)
trait_heritability <- c(0.3 ,0.7)
snp_sv_ratio <- c(0.5, 0.8)
sv_effects <- c(0.1, 0.5)
sv_diff_dist <- FALSE  # not going to plot different GxE effects (no difference)

# create heatmap
results_plot <- prediction_summary %>% 
  filter(h2 %in% trait_heritability, var %in% trait_var_source, qtn %in% trait_qtn_number,
         ratio %in% snp_sv_ratio, sv_effect %in% sv_effects, diff_dist %in% sv_diff_dist) %>%
  mutate(predictor = factor(predictor,
                            levels = c("all", "snp", "sv", "snp_ld", "snp_not_ld"),
                            labels = c("SNPs and SVs", "Only SNPs", "Only SVs",
                                       "SNPs in LD", "SNPs not in LD")),
         h2 = factor(h2, levels = c("0.3", "0.7"), labels = c("h2 = 0.3", "h2 = 0.7")),
         qtn = factor(qtn, levels = c("10", "100"), labels = c("10 QTLs", "100 QTLs")),
         ratio = factor(ratio, levels = c("0.5", "0.8"), labels = c("50% SVs", "20% SVs")),
         # ratio = factor(ratio, levels = c("0.5", "0.8"), labels = c("50%\nSVs", "20%\nSVs")),
         sv_effect = factor(sv_effect, levels = c("0.1", "0.5"), labels = c("effect\nSV = SNP", "effect\nSV > SNP"))) %>% 
  select(h2, var, qtn, ratio, sv_effect, predictor, cv, mean_accuracy_envs) %>%
  unite(h2:sv_effect, col = "sim_scenarios", sep = "/", remove = FALSE) %>% 
  unite(predictor:cv, col = "pred_scenarios", sep = "/", remove = FALSE) %>%
  ggplot(aes(x = pred_scenarios, y = sim_scenarios, fill = mean_accuracy_envs)) +
  # facet_nested(h2 + qtn + ratio + sv_effect ~ predictor + cv, scales = "free", switch = "y",
  #              nest_line = element_line(linetype = 1, color = "gray80")) +
  facet_nested(h2 + qtn + sv_effect + ratio ~ predictor + cv, scales = "free", switch = "y",
               nest_line = element_line(linetype = 1, color = "gray80")) +
  geom_tile(aes(height = 1.5, width = 1.5)) +
  scale_fill_viridis_c(option = "B", limits = c(0, 1), name = "Accuracy") +
  theme_minimal() +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        strip.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        panel.grid = element_blank(),
        panel.spacing = unit(0.25, "lines"))

# save plot
plot_name <- paste0(outfolder, "/heatmap_pred-accuracy_vars-both.pdf")
ggsave(plot = results_plot, filename = plot_name, device = "pdf",
       width = 18, height = 24, units = "in")

#### debug ----

# prediction_summary_file <- "tests/test_prediction_results.pops15-4-18.summary.txt"
# prediction_summary_file <- "tests/test_prediction_results.15pops.summary.txt"
