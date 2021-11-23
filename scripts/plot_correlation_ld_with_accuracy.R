library(data.table)
library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)

usage <- function() {
  cat("
description: summarize and correlate LD between causative variants and prediction accuracy.

usage: Rscript plot_correlation_ld_with_accuracy.R [ld_pred_accuracy_file] [ld_summary_file] [folder_results] [...]

positional arguments:
  ld_pred_accuracy_file         path to file with summary of LD and prediction accuracy
  ld_summary_file               path to file with LD summary
  folder_results                path to folder to save plots

optional argument:
  --help                        show this helpful message

note: make sure that to provide the same number of `--trait-pops` and `--pred-iters`

"
  )
}




#### command line options ----

# set default
args <- commandArgs(trailingOnly = TRUE)
if ("--help" %in% args) usage() & q(save = "no")

# assert to have the correct optional arguments
if (length(args) != 3) stop(usage(), "missing positional argument(s)")

# get positional arguments
ld_pred_accuracy_file <- args[1]
ld_summary_file <- args[2]
folder_results <- args[3]



#### plot summary of LD between causative variants and predictors ----

# load files
results_accuracy_ld_pve <- fread(ld_pred_accuracy_file, header = TRUE, data.table = FALSE)
ld_summary_all_scenarios <- fread(ld_summary_file, header = TRUE, data.table = FALSE)
# results_accuracy_ld_pve <- results_accuracy_ld_pve %>%
#   select(!c(se, lowerCI, upperCI, causal_var, qtn_type, qtn_number, var_exp_se, chr, predictor_highest_ld))
results_accuracy_ld_pve <- results_accuracy_ld_pve %>%
  select(!c(se, lowerCI, upperCI, qtn_type, var_exp_se, chr))

# results_accuracy_ld_pve <- separate_rows(results_accuracy_ld_pve, r2, bp_distance, sep = ";")
# results_accuracy_ld_pve$r2 <- as.numeric(results_accuracy_ld_pve$r2)
# results_accuracy_ld_pve$bp_distance <- as.numeric(results_accuracy_ld_pve$bp_distance)
# results_accuracy_ld_pve$pop <- as.factor(results_accuracy_ld_pve$pop)



##### 1) overview ----

plot1 <- results_accuracy_ld_pve %>%
  filter(var %in% c("SNP", "SV")) %>%
  ggplot(aes(x = mean, y = r2)) +
  facet_grid(cv + predictor ~ qtn + h2 + var) +
  geom_point() +
  geom_smooth() +
  scale_x_discrete(drop = FALSE)

ggsave(plot = plot1, filename = paste0(folder_results, "/overview1.pdf"), device = "pdf")

plot2 <- results_accuracy_ld_pve %>%
  filter(var %in% c("SNP", "SV")) %>%
  ggplot(aes(x = r2, y = mean)) +
  facet_grid(cv + predictor ~ qtn + h2 + var) +
  geom_point(aes(color = var_exp_mean), size = 1) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1))

ggsave(plot = plot2, filename = paste0(folder_results, "/overview2.pdf"), device = "pdf")

plot3 <- results_accuracy_ld_pve %>%
  filter(var %in% c("SNP", "SV")) %>%
  ggplot(aes(x = r2, y = var_exp_mean)) +
  facet_grid(cv + predictor ~ qtn + h2 + var) +
  geom_point(aes(color = mean), size = 1) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1))
# plots above don't show much of a pattern -- it seems that qtn number, h2,
# cv, etc. are more important than just ld and/or qtn pve

ggsave(plot = plot3, filename = paste0(folder_results, "/overview3.pdf"), device = "pdf")



##### 2) closer look with less variables ----

plot4 <- results_accuracy_ld_pve %>%
  filter(predictor == "all", var %in% c("SNP", "SV"), qtn == 100) %>%
  ggplot(aes(x = r2, y = mean)) +
  facet_grid(cv ~ h2 + var) +
  geom_point(aes(color = var_exp_mean), size = 1) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1))

ggsave(plot = plot4, filename = paste0(folder_results, "/focused1.pdf"), device = "pdf")

plot5 <- results_accuracy_ld_pve %>%
  filter(predictor == "all", var %in% c("SNP", "SV"), qtn == 100) %>%
  ggplot(aes(x = r2, y = var_exp_mean)) +
  facet_grid(cv ~ h2 + var) +
  geom_point(aes(color = mean), size = 1) +
  coord_cartesian(xlim = c(0, 1))
# very hard to see any patterns with all these points over each other
# using transparency wouldn't help because the color is mapped to mean accuracy

ggsave(plot = plot5, filename = paste0(folder_results, "/focused2.pdf"), device = "pdf")



##### 3) transforming continuous variables to discrete----

plot6 <- results_accuracy_ld_pve %>%
  mutate(accuracy_bins = cut_width(mean, width = 0.1, boundary = 0, closed = "right")) %>%
  # filter(predictor == "all", var %in% c("SNP", "SV"), qtn == 100) %>%
  ggplot(aes(x = accuracy_bins, y = r2)) +
  geom_boxplot() +
  scale_x_discrete(drop = FALSE)

ggsave(plot = plot6, filename = paste0(folder_results, "/discrete1.pdf"), device = "pdf")

plot7 <- results_accuracy_ld_pve %>%
  mutate(accuracy_bins = cut_width(mean, width = 0.1, boundary = 0, closed = "right")) %>%
  filter(var %in% c("SNP", "SV")) %>%
  ggplot(aes(x = accuracy_bins, y = r2)) +
  facet_grid(qtn ~ h2) +
  geom_boxplot() +
  stat_summary(fun.data = function(x) c(y = 1.02, label = length(x)), geom = "text", size = 2) +
  scale_x_discrete(drop = FALSE)

ggsave(plot = plot7, filename = paste0(folder_results, "/discrete2.pdf"), device = "pdf")

plot8 <- results_accuracy_ld_pve %>%
  mutate(accuracy_bins = cut_width(mean, width = 0.1, boundary = 0, closed = "right")) %>%
  group_by(env, h2, qtn, var, ratio, sv_effect, diff_dist, pop, predictor, pred_iter, cv) %>%
  mutate(ld_NAs = sum(is.na(r2)),
         prop_high_ld = round(sum(r2 > 0.8, na.rm = TRUE) / (n() - ld_NAs), digits = 3)) %>%
  ungroup() %>%
  group_by(accuracy_bins) %>%
  mutate(prop_high_ld_mean = mean(prop_high_ld)) %>%
  ungroup() %>%
  filter(!is.na(accuracy_bins)) %>%
  ggplot(aes(x = accuracy_bins, y = r2)) +
  geom_boxplot(aes(fill = prop_high_ld_mean)) +
  stat_summary(fun.data = function(x) c(y = 1.02, label = length(x)), geom = "text", size = 2) +
  scale_x_discrete(drop = FALSE) +
  scale_fill_continuous(type = "viridis")

ggsave(plot = plot8, filename = paste0(folder_results, "/discrete3.pdf"), device = "pdf")

plot9 <- results_accuracy_ld_pve %>%
  mutate(accuracy_bins = cut_width(mean, width = 0.1, boundary = 0, closed = "right")) %>%
  group_by(env, h2, qtn, var, ratio, sv_effect, diff_dist, pop, predictor, pred_iter, cv) %>%
  mutate(ld_NAs = sum(is.na(r2)),
         prop_high_ld = round(sum(r2 > 0.8, na.rm = TRUE) / (n() - ld_NAs), digits = 3),
         prop_ld_NAs = ld_NAs / n()) %>%
  ungroup() %>%
  filter(!is.na(accuracy_bins)) %>%
  group_by(accuracy_bins) %>%
  mutate(prop_ld_NAs_mean = mean(prop_ld_NAs)) %>%
  ungroup() %>%
  ggplot(aes(x = accuracy_bins, y = prop_high_ld)) +
  geom_boxplot(aes(alpha = prop_ld_NAs_mean), fill = "firebrick") +
  coord_cartesian(ylim = c(0, 1)) +
  stat_summary(fun.data = function(x) c(y = 1.02, label = length(x)), geom = "text", size = 2) +
  scale_x_discrete(drop = FALSE) +
  scale_alpha_continuous(range = c(1, 0.5))

ggsave(plot = plot9, filename = paste0(folder_results, "/discrete4.pdf"), device = "pdf")

plot10 <- results_accuracy_ld_pve %>%
  mutate(accuracy_bins = cut_width(mean, width = 0.1, boundary = 0, closed = "right")) %>%
  filter(var %in% c("SNP", "SV")) %>%
  ggplot(aes(x = accuracy_bins, y = r2)) +
  facet_grid(cv + predictor ~ qtn + h2 + var) +
  geom_boxplot() +
  scale_x_discrete(drop = FALSE)

ggsave(plot = plot10, filename = paste0(folder_results, "/discrete5.pdf"), device = "pdf")

plot11 <- results_accuracy_ld_pve %>%
  group_by(h2, qtn, var, ratio, sv_effect, diff_dist, pop, predictor, pred_iter, cv) %>%
  summarize(accuracy_mean = mean(mean),
            qtn_maf_mean = mean(qtn_maf),
            qtn_effect_mean = mean(abs(qtn_effect)),
            var_exp_mean = mean(var_exp_mean),
            ld_NAs = sum(is.na(r2)),
            high_ld = round(sum(r2 > 0.8, na.rm = TRUE) / (n() - ld_NAs), digits = 3),
            bp_distance = mean(abs(bp_distance))) %>%
  ungroup() %>%
  arrange(h2, qtn, desc(var), ratio, sv_effect, diff_dist, pop, predictor, pred_iter, cv) %>%
  filter(var %in% c("SNP", "SV")) %>%
  ggplot(aes(x = accuracy_mean, y = high_ld)) +
  facet_grid(cv + predictor ~ qtn + h2 + var) +
  geom_point() +
  coord_cartesian(ylim = c(0, 1)) +
  scale_x_discrete(drop = FALSE)

ggsave(plot = plot11, filename = paste0(folder_results, "/discrete6.pdf"), device = "pdf")

plot12 <- results_accuracy_ld_pve %>%
  group_by(env, h2, qtn, var, ratio, sv_effect, diff_dist, pop, predictor, pred_iter, cv) %>%
  summarize(accuracy_mean = mean(mean),
            qtn_maf_mean = mean(qtn_maf),
            qtn_effect_mean = mean(abs(qtn_effect)),
            var_exp_mean = mean(var_exp_mean),
            ld_NAs = sum(is.na(r2)),
            high_ld = round(sum(r2 > 0.8, na.rm = TRUE) / (n() - ld_NAs), digits = 3),
            bp_distance = mean(abs(bp_distance))) %>%
  ungroup() %>%
  arrange(env, h2, qtn, desc(var), ratio, sv_effect, diff_dist, pop, predictor, pred_iter, cv) %>%
  mutate(accuracy_bins = cut_width(accuracy_mean, width = 0.1, boundary = 0, closed = "right")) %>%
  # filter(predictor == "all", var %in% c("SNP", "SV"), qtn == 100) %>%
  ggplot(aes(x = accuracy_bins, y = high_ld)) +
  geom_violin() +
  geom_point() +
  coord_cartesian(ylim = c(0, 1)) +
  scale_x_discrete(drop = FALSE)

ggsave(plot = plot12, filename = paste0(folder_results, "/discrete7.pdf"), device = "pdf")

plot13 <- results_accuracy_ld_pve %>%
  group_by(env, h2, qtn, var, ratio, sv_effect, diff_dist, pop, predictor, pred_iter, cv) %>%
  mutate(ld_NAs = sum(is.na(r2)),
         high_ld = round(sum(r2 > 0.8, na.rm = TRUE) / (n() - ld_NAs), digits = 3)) %>%
  ungroup() %>%
  arrange(env, h2, qtn, desc(var), ratio, sv_effect, diff_dist, pop, predictor, pred_iter, cv) %>%
  mutate(accuracy_bins = cut_width(mean, width = 0.1, boundary = 0, closed = "right")) %>%
  group_by(h2, qtn, accuracy_bins) %>%
  mutate(prop_high_ld_mean = mean(high_ld)) %>%
  ungroup() %>%
  filter(var %in% c("SNP", "SV")) %>%
  ggplot(aes(x = accuracy_bins, y = var_exp_mean)) +
  facet_grid(qtn ~ h2, scales = "free_y") +
  geom_violin(aes(fill = prop_high_ld_mean)) +
  scale_x_discrete(drop = FALSE) +
  scale_fill_continuous(type = "viridis")

ggsave(plot = plot13, filename = paste0(folder_results, "/discrete8.pdf"), device = "pdf")

plot14 <- results_accuracy_ld_pve %>%
  group_by(env, h2, qtn, var, ratio, sv_effect, diff_dist, pop, predictor, pred_iter, cv) %>%
  mutate(ld_NAs = sum(is.na(r2)),
         high_ld = round(sum(r2 > 0.8, na.rm = TRUE) / (n() - ld_NAs), digits = 3)) %>%
  ungroup() %>%
  arrange(env, h2, qtn, desc(var), ratio, sv_effect, diff_dist, pop, predictor, pred_iter, cv) %>%
  mutate(accuracy_bins = cut_width(mean, width = 0.1, boundary = 0, closed = "right")) %>%
  group_by(h2, qtn, accuracy_bins) %>%
  mutate(prop_high_ld_mean = mean(high_ld)) %>%
  ungroup() %>%
  filter(var %in% c("SNP", "SV")) %>%
  ggplot(aes(x = accuracy_bins, y = r2)) +
  facet_grid(qtn ~ h2, scales = "free_y") +
  geom_boxplot(aes(fill = prop_high_ld_mean)) +
  scale_x_discrete(drop = FALSE) +
  scale_fill_continuous(type = "viridis")
# legend: the proportion of qtns in that bin that has a r2 value higher than 0.8

ggsave(plot = plot14, filename = paste0(folder_results, "/discrete9.pdf"), device = "pdf")

plot15 <- results_accuracy_ld_pve %>%
  filter(var %in% c("SNP", "SV")) %>%
  group_by(env, h2, qtn, var, ratio, sv_effect, diff_dist, pop, predictor, pred_iter, cv) %>%
  mutate(ld_NAs = sum(is.na(r2)),
         high_ld = round(sum(r2 > 0.8, na.rm = TRUE) / (n() - ld_NAs), digits = 3)) %>%
  ungroup() %>%
  mutate(accuracy_bins = cut_width(mean, width = 0.1, boundary = 0, closed = "right")) %>%
  group_by(cv, predictor, qtn, h2, var, accuracy_bins) %>%
  mutate(prop_high_ld_mean = mean(high_ld)) %>%
  ungroup() %>%
  arrange(env, h2, qtn, desc(var), ratio, sv_effect, diff_dist, pop, predictor, pred_iter, cv) %>%
  filter(!is.na(accuracy_bins)) %>%
  ggplot(aes(x = accuracy_bins, y = r2)) +
  facet_grid(predictor + cv ~ var + qtn + h2, scales = "free_y") +
  geom_boxplot(aes(fill = prop_high_ld_mean)) +
  scale_x_discrete(drop = FALSE, labels = function(x) gsub(",", "\n", x)) +
  scale_fill_continuous(type = "viridis")

ggsave(plot = plot15, filename = paste0(folder_results, "/discrete10.pdf"), device = "pdf")



##### 4) weight metric (LD x PVE) -----

plot16 <- results_accuracy_ld_pve %>%
  mutate(weight_ld_pve = var_exp_mean * r2) %>%
  filter(var %in% c("SNP", "SV")) %>%
  ggplot(aes(x = mean, y = weight_ld_pve)) +
  facet_grid(qtn ~ h2, scales = "free_y") +
  geom_point() +
  geom_smooth()

ggsave(plot = plot16, filename = paste0(folder_results, "/ld-pve-weight1.pdf"), device = "pdf")

plot17 <- results_accuracy_ld_pve %>%
  mutate(weight_ld_pve = var_exp_mean * r2) %>%
  filter(var %in% c("SNP", "SV")) %>%
  ggplot(aes(x = mean, y = weight_ld_pve)) +
  facet_grid(cv + predictor ~ qtn + h2 + var) +
  geom_point() +
  geom_smooth()

ggsave(plot = plot17, filename = paste0(folder_results, "/ld-pve-weight2.pdf"), device = "pdf")



##### 5) plot LD distribution per scenario ----

# using only qtns in highest ld with predictors
plot18 <- results_accuracy_ld_pve %>%
  group_by(qtn, var, ratio, pop, predictor, pred_iter, causal_var, predictor_highest_ld) %>%
  summarize(r2 = mean(r2)) %>%
  ungroup() %>%
  ggplot(aes(x = r2)) +
  facet_grid(qtn + var ~ predictor, scales = "free_y") +
  geom_histogram() +
  labs(title = "Only predictor in highest LD to a QTN")

ggsave(plot = plot18, filename = paste0(folder_results, "/ld-dist1.pdf"), device = "pdf")

plot19 <- results_accuracy_ld_pve %>%
  group_by(qtn, var, ratio, pop, predictor, pred_iter, causal_var, predictor_highest_ld) %>%
  summarize(r2 = mean(r2)) %>%
  ungroup() %>%
  filter(var %in% c("SNP", "SV"), qtn == 100) %>%
  mutate(predictor = factor(predictor,
                            levels = c("all", "snp", "sv", "snp_ld", "snp_not_ld"),
                            labels = c("SNPs\nand SVs", "Only\nSNPs", "Only\nSVs",
                                       "SNPs\nin LD", "SNPs\nnot in LD"))) %>%
  ggplot(aes(x = r2)) +
  facet_grid(var ~ predictor) +
  geom_histogram() +
  theme_bw() +
  theme(title = element_text(size = 20),
        text = element_text(size = 20),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 25),
        strip.text.x = element_text(size = 20),
        strip.text.y = element_text(size = 20, angle = 0),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank()) +
  labs(x = bquote("LD"~(r^2)~"between QTNs and predictors"),
       y = "Count")

ggsave(plot = plot19, filename = paste0(folder_results, "/ld-dist2.pdf"), device = "pdf")

plot20 <- results_accuracy_ld_pve %>%
  group_by(qtn, var, ratio, pop, predictor, pred_iter, causal_var, predictor_highest_ld) %>%
  summarize(r2 = mean(r2)) %>%
  ungroup() %>%
  filter(var %in% c("SNP", "SV"), qtn == 100) %>%
  mutate(predictor = factor(predictor,
                            levels = c("all", "snp", "sv", "snp_ld", "snp_not_ld"),
                            labels = c("SNPs\nand SVs", "Only\nSNPs", "Only\nSVs",
                                       "SNPs\nin LD", "SNPs\nnot in LD"))) %>%
  ggplot(aes(x = r2)) +
  facet_grid(predictor ~ var) +
  geom_histogram() +
  theme_bw() +
  theme(title = element_text(size = 15),
        text = element_text(size = 15),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 20),
        strip.text.x = element_text(size = 15),
        strip.text.y = element_text(size = 15, angle = 0)) +
  labs(x = bquote("LD"~(r^2)~"between QTNs and predictors"),
       y = "Count")

ggsave(plot = plot20, filename = paste0(folder_results, "/ld-dist3.pdf"), device = "pdf")

plot21 <- results_accuracy_ld_pve %>%
  group_by(qtn, var, ratio, pop, predictor, pred_iter, causal_var, predictor_highest_ld) %>%
  summarize(r2 = mean(r2)) %>%
  ungroup() %>%
  filter(var %in% c("SNP"), qtn == 100, predictor %in% c("snp", "sv")) %>%
  mutate(predictor = factor(predictor,
                            levels = c("all", "snp", "sv", "snp_ld", "snp_not_ld"),
                            labels = c("all", "SNP markers", "SV markers",
                                       "snp_ld", "snp_not_ld"))) %>%
  ggplot(aes(x = r2)) +
  facet_grid(~ predictor) +
  geom_histogram() +
  coord_cartesian(ylim =  c(0, 150), xlim = c(0, 1)) +
  theme_bw() +
  theme(title = element_text(size = 20),
        text = element_text(size = 20),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 20),
        axis.title = element_text(size = 20),
        strip.text.x = element_text(size = 20),
        strip.text.y = element_text(size = 20, angle = 0),
        panel.spacing = unit(0.25, "inches"),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank()) +
  labs(x = bquote("LD"~(r^2)~"between QTLs and markers"),
       y = "Count",
       title = "SNPs as causative variants")
# 4hx10w inches

ggsave(plot = plot21, filename = paste0(folder_results, "/ld-dist4.pdf"), device = "pdf",
       width = 10, height = 4, units = "in")

plot22 <- results_accuracy_ld_pve %>%
  group_by(qtn, var, ratio, pop, predictor, pred_iter, causal_var, predictor_highest_ld) %>%
  summarize(r2 = mean(r2)) %>%
  ungroup() %>%
  filter(var %in% c("SV"), qtn == 100, predictor %in% c("snp", "sv")) %>%
  mutate(predictor = factor(predictor,
                            levels = c("all", "snp", "sv", "snp_ld", "snp_not_ld"),
                            labels = c("all", "SNP markers", "SV markers",
                                       "snp_ld", "snp_not_ld"))) %>%
  ggplot(aes(x = r2)) +
  facet_grid(~ predictor) +
  geom_histogram() +
  coord_cartesian(ylim =  c(0, 150), xlim = c(0, 1)) +
  theme_bw() +
  theme(title = element_text(size = 20),
        text = element_text(size = 20),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 20),
        axis.title = element_text(size = 20),
        strip.text.x = element_text(size = 20),
        strip.text.y = element_text(size = 20, angle = 0),
        panel.spacing = unit(0.25, "inches"),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank()) +
  labs(x = bquote("LD"~(r^2)~"between QTLs and markers"),
       y = "Count",
       title = "SVs as causative variants")
# 4hx10w inches

ggsave(plot = plot22, filename = paste0(folder_results, "/ld-dist5.pdf"), device = "pdf",
       width = 10, height = 4, units = "in")

plot23 <- results_accuracy_ld_pve %>%
  group_by(qtn, var, ratio, pop, predictor, pred_iter, causal_var, predictor_highest_ld) %>%
  summarize(r2 = mean(r2)) %>%
  ungroup() %>%
  group_by(qtn, var, ratio, pop, predictor, pred_iter) %>%
  summarize(prop_high_ld = mean(r2 > 0.8, na.rm = TRUE)) %>%
  ungroup() %>%
  ggplot(aes(x = var, y = prop_high_ld)) +
  facet_grid(qtn ~ predictor) +
  geom_boxplot() +
  coord_cartesian(ylim = c(0, 1)) +
  labs(title = "Proportion of predictor in highest LD to a QTN with r2 > 0.8")

ggsave(plot = plot23, filename = paste0(folder_results, "/ld-dist6.pdf"), device = "pdf")



##### 6) ld summary -----

# using all r2 values of qtns and predictors
plot24 <- ld_summary_all_scenarios %>%
  ggplot(aes(x = r2)) +
  facet_grid(qtn + var ~ pred_type, scales = "free_y") +
  geom_histogram() +
  labs(title = "LD between all predictors and all QTNs")

ggsave(plot = plot24, filename = paste0(folder_results, "/ld-summary1.pdf"), device = "pdf")

# average number of predictors in high ld with same QTN
plot25 <- ld_summary_all_scenarios %>%
  group_by(qtn, var, ratio, pop, pred_type, pred_iter, causal_var) %>%
  summarize(n_high_ld_per_qtn = sum(r2 > 0.8, na.rm = TRUE)) %>%
  ungroup() %>%
  ggplot(aes(x = var, y = n_high_ld_per_qtn)) +
  facet_grid(~ qtn ~ pred_type) +
  geom_boxplot()
# on average, a single SV causative variant tend to be tagged by more predictors than a SNP causative variant

ggsave(plot = plot25, filename = paste0(folder_results, "/ld-summary2.pdf"), device = "pdf")

# average proportion of QTNs being tagged
plot26 <- ld_summary_all_scenarios %>%
  group_by(qtn, var, pred_type) %>%
  summarize(prop_high_ld = mean(r2 > 0.8, na.rm = TRUE),
            se = sd(r2 > 0.8, na.rm = TRUE) / sqrt(n())) %>%
  ungroup() %>%
  ggplot(aes(x = var, y = prop_high_ld)) +
  facet_grid(qtn ~ pred_type) +
  geom_col() +
  geom_errorbar(aes(ymin = prop_high_ld - se, ymax = prop_high_ld + se),
                position = position_dodge(0.9), width = 0.2)
# when we use SNP predictors, ~3% of SNP causative variants are tagged but only ~1% of SV causative variants are tagged
# when we use SV predictors, the amount of causative variants tagged doesn't appear to depend on the causative variant type

ggsave(plot = plot26, filename = paste0(folder_results, "/ld-summary3.pdf"), device = "pdf")



#### debug ----

# ld_pred_accuracy_file <-"analysis/trait_sim/multi_env/pred-accuracy_qtn-effect_ld-pred-qtn.results.txt"
# ld_summary_file <- "analysis/trait_sim/multi_env/pred-accuracy_qtn-effect_ld-pred-qtn.ld-summary.txt"
# folder_results <- "analysis/trait_sim/multi_env/cor_ld-pred-accuracy"
