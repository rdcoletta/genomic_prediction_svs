library(data.table)
library(ggplot2)
library(grid)
library(dplyr)
library(doParallel)
library(tidyr)
library(ggh4x)
dir.create("figures")



#### figure 1 ----

# done in https://lucid.app/lucidchart/



#### figure 2 ----

##### a #####

ld_file <- "analysis/snp-chip_ld/usda_22kSNPs_rils.sorted.diploid.filtered.v4.window-1000kb.filter-0.25.ld.gz"
output_folder <- "figures"

# load data
ld_results <- fread(ld_file, header = TRUE, data.table = FALSE)
# if file is compressed, need to remove first line (which is copy of the header)
if (grepl(".gz", ld_file)) ld_results <- ld_results[-1,]

# make sure R2 is numeric
ld_results$R2 <- as.numeric(ld_results$R2)

# get distances between SNPs
ld_results$dist_markers <- abs(as.numeric(ld_results[, "BP_B"]) - as.numeric(ld_results[, "BP_A"]))

# set default window parameters -- max distance possible is 100mb
window_parameters <- data.frame(bp_start =    c(   1,  1001, 10001,  50001,  100001,  1000001,  10000000),
                                bp_stop =     c(1000, 10000, 50000, 100000, 1000000, 10000000, 100000000),
                                window_size = c( 100,  1000, 10000,  50000,  100000,  1000000,  10000000))

# adjust parameters based on maximum distance between markers
window_parameters <- window_parameters[which(window_parameters[, "bp_start"] < max(ld_results[, "dist_markers"])), ]

# create df to store results
df_plot <- data.frame(stringsAsFactors = FALSE)

for (row in 1:NROW(window_parameters)) {
  
  # set window values
  window_start <- window_parameters[row , "bp_start"]
  window_stop <- window_parameters[row, "bp_stop"]
  window_size <- window_parameters[row , "window_size"]
  
  # get ld for window
  ld_window <- ld_results[which(ld_results[, "dist_markers"] >= window_start & ld_results[, "dist_markers"] < window_stop), c("R2", "dist_markers")]
  
  # set df to store average of each window
  df_window <- data.frame(stringsAsFactors = FALSE)
  window_stop <- window_start + (window_size - 1)
  window_step <- window_size
  
  # get average R2 per window
  while (window_stop < max(ld_window[, "dist_markers"])) {
    
    # get all data points for each window
    window_avg <- ld_window[ld_window[, "dist_markers"] > window_start & ld_window[, "dist_markers"] <= window_stop, c("R2", "dist_markers")]
    window_avg <- data.frame(window_avg, n_markers = NROW(window_avg), bp_stop = window_stop, stringsAsFactors = FALSE)
    
    # append to main df
    df_window <- rbind(df_window, window_avg)
    
    # set up the start of next window
    window_start <- window_start + window_step
    window_stop <- window_start + (window_size - 1)
    
  }
  # get all data points for last window
  window_avg <- ld_window[ld_window[, "dist_markers"] > window_start & ld_window[, "dist_markers"] <= window_stop, c("R2", "dist_markers")]
  window_avg <- data.frame(window_avg, n_markers = NROW(window_avg), bp_stop = window_stop, stringsAsFactors = FALSE)
  # append to main df
  df_window <- rbind(df_window, window_avg)
  
  # append to main df
  df_plot <- rbind(df_plot, df_window)
  
}

# adjust scale for x axis
old_levels <- as.integer(levels(as.factor(df_plot[, "bp_stop"])))
new_levels <- c()
for (i in 1:length(old_levels)) {
  
  if (i == 1) {
    new_levels <- c(new_levels, paste0("0-", old_levels[i] / 1000))
  } else {
    new_levels <- c(new_levels, paste0(old_levels[i - 1] / 1000, "-", old_levels[i] / 1000))
  }
  
}

# plot ld decay
out_plot <- ggplot(df_plot, aes(x = as.factor(bp_stop), y = R2)) +
  geom_boxplot(width = 0.4, outlier.size = 0.01) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Physical distance (kb)",
       y = bquote("LD"~(r^2))) +
  scale_x_discrete(labels = new_levels) +
  stat_summary(fun.data = function(x) c(y = 1.03, label = round(length(x) / 1000, 1)), geom = "text", size = 1.5) +
  geom_boxplot(fill = "gray70", width = 0.8) +
  # labs(title = "LD decay") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        axis.text.x = element_text(size = 8, angle = 45, hjust = 1))
ggsave(filename = paste0(output_folder, "/fig_2a.png"), plot = out_plot,
       device = "png", dpi = 350, units = "in", width = 6, height = 4)


##### b #####

# get ld between snps and svs from multiple chr
ld_results <- data.frame()
for (chr in 1:10) {
  ld_file_chr <- paste0("analysis/ld/window-1kb_filter-0.25/ld_usda_rils_snp-sv_only.chr", chr, ".window-1kb.filter-0.25.highest-ld.ld")
  ld_file_chr <- fread(ld_file_chr, header = TRUE, data.table = FALSE)
  ld_results <- rbind(ld_results, ld_file_chr)
}
rm(ld_file_chr, chr)

# get histogram data
plot_data <- ggplot_build(ggplot(ld_results, aes(x = R2)) + geom_histogram(binwidth = 0.01))
plot_data <- data.frame(plot_data$data)[c("x", "ymin", "ymax")]
plot_data[NROW(plot_data), "ymin"] <- 2000
plot_data <- rbind(plot_data, data.frame(x = plot_data[NROW(plot_data), "x"], ymin = 0, ymax = 100))
plot_data[plot_data$ymin > 100, "y_axis_break"] <- "above"
plot_data[plot_data$ymin <= 100, "y_axis_break"] <- "below"

# plot distribution with y-axis split
dist_plot <- ggplot(plot_data, aes(x = x, y = ymin, xend = x, yend = ymax)) +
  geom_segment(size = 5, color = "gray50") +
  facet_grid(y_axis_break ~ ., scales = "free_y") +
  coord_cartesian(xlim = c(0, 1)) +
  labs(x = bquote("LD"~(r^2)),
       y = "Count") +
  theme_minimal() +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        strip.text.y = element_blank(),
        panel.grid = element_blank(),
        axis.ticks = element_line()) + 
  scale_x_continuous(breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1")) +
  annotation_custom(grob = linesGrob(gp = gpar(lwd = 2)), xmin = -Inf, xmax = -Inf, ymin = -Inf, ymax = +Inf) +
  annotation_custom(grob = linesGrob(gp = gpar(lwd = 2)), xmin = +Inf, xmax = +Inf, ymin = -Inf, ymax = +Inf) +
  annotation_custom(grob = linesGrob(gp = gpar(lwd = 2)), xmin = -Inf, xmax = +Inf, ymin = -5, ymax = -5) +
  annotation_custom(grob = linesGrob(gp = gpar(lwd = 2)), xmin = -Inf, xmax = +Inf, ymin = 4231, ymax = 4231)
ggsave(filename = paste0(output_folder, "/fig_2b.png"), plot = dist_plot,
       device = "png", dpi = 350, units = "in", width = 2, height = 4)

# inspired by:
#   https://stackoverflow.com/a/65242880
#   https://www.j4s8.de/post/2018-01-15-broken-axis-with-ggplot2/
# finished off on powerpoint:
#   - added white background and lines cutting y axis

# clean R environment
rm(list = ls())



#### figure 3 ----

##### a #####

prediction_summary_file <- "analysis/trait_sim/multi_env/prediction_results.summary.txt"
output_folder <- "figures"

# load data
prediction_summary <- fread(prediction_summary_file, header = TRUE, data.table = FALSE)

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
        strip.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.2, "in"),
        panel.grid = element_blank(),
        panel.spacing = unit(0.25, "lines"))
# save plot
ggsave(filename = paste0(output_folder, "/fig_3a.png"), plot = results_plot,
       device = "png", dpi = 350, units = "in", width = 8, height = 3.3)

# clean R environment
rm(list = ls())

##### b #####

# 1. As both heritability and number of QTLs increased, including both SNPs and SVs or only SVs as
#    predictors resulted in more similar accuracies as with SNPs alone

prediction_summary_file <- "analysis/trait_sim/multi_env/prediction_results.summary.txt"
output_folder <- "figures"

# load data
prediction_summary <- fread(prediction_summary_file, header = TRUE, data.table = FALSE)

# set default
trait_var_source <- c("SNP", "SV")
error_bars <- "CI_accuracy"

# filter summary
prediction_summary_plot <- prediction_summary %>% 
  filter(var %in% trait_var_source & predictor %in% c("all", "snp", "sv")) %>%
  group_by(h2, qtn, var, predictor) %>%
  summarize(mean_pred = mean(mean_accuracy_envs, na.rm = TRUE),
            se_pred = sd(mean_accuracy_envs, na.rm = TRUE) / n()) %>% 
  ungroup()

# reorder levels of marker type for better visualization
prediction_summary_plot$predictor <- factor(prediction_summary_plot$predictor,
                                            levels = c("snp", "sv", "all"),
                                            labels = c("SNPs only", "SVs only", "SNPs and SVs"))

# plot
results_plot <- prediction_summary_plot %>%
  mutate(h2_qtn = paste0("h2 = ",h2, "\n", qtn, " QTLs")) %>%
  ggplot(aes(x = h2_qtn, y = mean_pred)) +
  # facet_nested(~ var_label + var, nest_line = element_line(linetype = 1, color = "gray80")) +
  facet_nested(~ var, nest_line = element_line(linetype = 1, color = "gray80"),
               labeller = labeller(var = function(x) paste0(x, "s as causative variants"))) +
  geom_line(aes(color = predictor, group = predictor), size = 1) +
  geom_point(aes(color = predictor), size = 2) +
  coord_cartesian(ylim = c(0, 1)) +
  scale_x_discrete(drop = FALSE) +
  labs(x = "Trait architecture",
       y = "Prediction \naccuracy") +
  scale_color_manual(values = c("#A21E0DFF", "#F4A652FF", "#3F5B8AFF")) +
  guides(color = guide_legend("Predictors")) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 12),
        axis.title.x = element_text(vjust = -1.5),
        axis.title.y = element_text(vjust = 1.5),
        axis.text = element_text(size = 10),
        strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.position = "right",
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10)) +
  geom_errorbar(aes(ymin = mean_pred - se_pred, ymax = mean_pred + se_pred, color = predictor),
                width = 0.1)

# save plot
ggsave(filename = paste0(output_folder, "/fig_3b.png"),
       plot = results_plot, device = "png", dpi = 350, units = "in", width = 8, height = 2.2)

# fair_cols <- c("#38170B","#BF1B0B", "#FFC465", "#66ADE5", "#252A52")
# fair_ramp <- scales::colour_ramp(fair_cols)
# fair_ramp(seq(0, 1, 0.05))

# clean R environment
rm(list = ls())

##### c ----

# 2. when a trait was controlled by SVs exclusively, using SVs as predictors (either by themselves or
#    in conjunction with SNPs) resulted in higher accuracy than using only SNPs in all scenarios

prediction_summary_file <- "analysis/trait_sim/multi_env/prediction_results.summary.txt"
output_folder <- "figures"

# load data
prediction_summary <- fread(prediction_summary_file, header = TRUE, data.table = FALSE)

# set default
trait_var_source <- c("SNP", "SV")
error_bars <- "CI_accuracy"

# filter summary
prediction_summary_plot <- prediction_summary %>%
  filter(var %in% trait_var_source) %>%
  group_by(var, predictor) %>%
  summarize(mean_pred = mean(mean_accuracy_envs, na.rm = TRUE),
            se_pred = sd(mean_accuracy_envs, na.rm = TRUE) / n()) %>% 
  ungroup()

# reorder levels of marker type for better visualization
prediction_summary_plot$predictor <- factor(prediction_summary_plot$predictor,
                                            levels = c("all", "snp", "sv", "snp_ld", "snp_not_ld"),
                                            labels = c("SNPs\nand SVs", "SNPs\nonly", "SVs\nonly",
                                                       "SNPs\nin LD", "SNPs\nnot in LD"))

# # plot
# results_plot <- ggplot(data = prediction_summary_plot,
#        aes(x = predictor, y = mean_pred)) +
#   # facet_nested(~ var, nest_line = element_line(linetype = 1, color = "gray80"),
#   #              labeller = labeller(var = function(x) paste0(x, "s as causative variants"))) +
#   geom_line(aes(color = var, group = var), size = 1) +
#   geom_point(aes(color = var), size = 2) +
#   coord_cartesian(ylim = c(0, 1)) +
#   scale_x_discrete(drop = FALSE) +
#   labs(x = "Marker type used in prediciton",
#        y = "Prediction\naccuracy") +
#   scale_color_manual(values = c("grey70", "grey40")) +
#   guides(color = guide_legend("Causative variants")) +
#   theme_bw() +
#   theme(panel.grid = element_blank(),
#         axis.title = element_text(size = 12),
#         axis.title.x = element_text(vjust = -1.5),
#         axis.title.y = element_text(vjust = 1.5),
#         axis.text = element_text(size = 10),
#         strip.text.x = element_text(size = 12),
#         strip.text.y = element_text(size = 12),
#         strip.background = element_blank(),
#         strip.placement = "outside",
#         legend.title = element_text(size = 12),
#         legend.text = element_text(size = 10)) +
#   geom_errorbar(aes(ymin = mean_pred - se_pred, ymax = mean_pred + se_pred, color = var),
#                 width = 0.1)

# plot
results_plot <- ggplot(data = prediction_summary_plot,
       aes(x = predictor, y = mean_pred)) +
  facet_nested(~ var, nest_line = element_line(linetype = 1, color = "gray80"),
               labeller = labeller(var = function(x) paste0(x, "s as causative variants"))) +
  geom_bar(stat = "identity", position = "dodge", fill = "gray60") +
  # geom_line(aes(color = var, group = var), size = 1, show.legend = FALSE) +
  # geom_point(aes(color = var), size = 2, show.legend = FALSE) +
  coord_cartesian(ylim = c(0, 1)) +
  scale_x_discrete(drop = FALSE) +
  labs(x = "Marker type used in prediciton",
       y = "Prediction\naccuracy") +
  # scale_color_manual(values = c("grey70", "grey40")) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        # plot.margin = margin(5.5, 105.5, 5.5, 5.5),
        axis.title = element_text(size = 12),
        axis.title.x = element_text(vjust = -1.5),
        axis.title.y = element_text(vjust = 1.5),
        axis.text = element_text(size = 10),
        strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10)) +
  geom_errorbar(aes(ymin = mean_pred - se_pred, ymax = mean_pred + se_pred),
                width = 0.1)

# save plot
ggsave(filename = paste0(output_folder, "/fig_3c.png"),
       plot = results_plot, device = "png", dpi = 350, units = "in", width = 8, height = 2.2)

# clean R environment
rm(list = ls())

##### d ----

# 3. using SVs as predictors when SVs were the causative variants reduced the accuracy difference
#    between CV1 and CV2, especially at higher heritability and QTL numbers

prediction_summary_file <- "analysis/trait_sim/multi_env/prediction_results.summary.txt"
output_folder <- "figures"

# load data
prediction_summary <- fread(prediction_summary_file, header = TRUE, data.table = FALSE)

# set default
trait_var_source <- c("SNP", "SV")
error_bars <- "CI_accuracy"

# filter summary
prediction_summary_plot <- prediction_summary %>%
  filter(var %in% trait_var_source & predictor %in% c("all", "snp", "sv")) %>%
  group_by(h2, qtn, var, predictor) %>%
  summarize(diff_pred = diff(mean_accuracy_envs, na.rm = TRUE)) %>%
  ungroup()

# reorder levels of marker type for better visualization
prediction_summary_plot$predictor <- factor(prediction_summary_plot$predictor,
                                            levels = c("snp", "sv", "all"),
                                            labels = c("SNPs\nonly", "SVs\nonly", "SNPs\nand SVs"))

# plot
results_plot <- ggplot(data = prediction_summary_plot,
       aes(x = predictor, y = diff_pred)) +
  # facet_nested(~ var_label + var, nest_line = element_line(linetype = 1, color = "gray80")) +
  facet_nested(~ var, nest_line = element_line(linetype = 1, color = "gray80"),
               labeller = labeller(var = function(x) paste0(x, "s as causative variants"))) +
  geom_line(aes(linetype = as.factor(qtn), color = as.factor(h2), group = interaction(h2, qtn)), size = 1) +
  geom_point(aes(color = as.factor(h2)), size = 2) +
  coord_cartesian(ylim = c(0, 0.15)) +
  scale_x_discrete(drop = FALSE) +
  labs(x = "Marker type used in prediciton",
       y = "Prediction accuracy difference\n(CV2 - CV1)") +
  # scale_color_manual(values = c("#66ADE5FF", "#252A52FF")) +
  scale_color_manual(values = c("#DE8282FF", "#AD0000FF")) +
  scale_linetype_manual(values = c("dotted", "solid")) +
  guides(color = guide_legend("Heritability"),
         linetype = guide_legend("QTL number")) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 12),
        axis.title.x = element_text(vjust = -1.5),
        axis.title.y = element_text(size = 10, vjust = 1.5),
        axis.text = element_text(size = 10),
        strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.position = "right",
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10))

# save plot
ggsave(filename = paste0(output_folder, "/fig_3d.png"),
       plot = results_plot, device = "png", dpi = 350, units = "in", width = 8, height = 2.2)

# # filter summary
# prediction_summary_plot <- subset(prediction_summary,
#                                   var %in% trait_var_source
#                                   & predictor %in% c("all", "snp", "sv"))
# # prediction_summary_plot <- subset(prediction_summary,
# #                                   qtn == 100
# #                                   & h2 == 0.7
# #                                   & var %in% trait_var_source
# #                                   & predictor %in% c("all", "snp", "sv"))
# 
# # reorder levels of marker type for better visualization
# prediction_summary_plot$predictor <- factor(prediction_summary_plot$predictor,
#                                             levels = c("snp", "sv", "all"),
#                                             labels = c("SNPs\nonly", "SVs\nonly", "SNPs\nand SVs"))
# 
# # # add new label
# # prediction_summary_plot$var_label <- "Causative variants"
# 
# # plot
# results_plot <- ggplot(data = prediction_summary_plot,
#                        aes(x = predictor, y = mean_accuracy_envs)) +
#   # facet_nested(~ var_label + var, nest_line = element_line(linetype = 1, color = "gray80")) +
#   # geom_line(aes(color = cv, group = cv), size = 1) +
#   # geom_point(aes(color = cv), size = 2) +
#   facet_nested(~ var, nest_line = element_line(linetype = 1, color = "gray80"), 
#                labeller = labeller(var = function(x) paste0(x, "s as causative variants"))) +
#   geom_line(aes(linetype = as.factor(qtn), color = as.factor(h2), group = interaction(h2, qtn)), size = 1) +
#   geom_point(aes(color = as.factor(h2)), size = 2)
#   coord_cartesian(ylim = c(0, 1)) +
#   scale_x_discrete(drop = FALSE) +
#   labs(x = "Marker type used in prediciton",
#        y = "Prediction accuracy") +
#   scale_color_manual(values = c("#DE8282FF", "#AD0000FF")) +
#   guides(color = guide_legend("Cross-Validation")) +
#   theme_bw() +
#   theme(panel.grid = element_blank(),
#         axis.title = element_text(size = 12),
#         axis.title.x = element_text(vjust = -1.5),
#         axis.title.y = element_text(vjust = 1.5),
#         axis.text = element_text(size = 10),
#         strip.text.x = element_text(size = 12),
#         strip.text.y = element_text(size = 12),
#         strip.background = element_blank(),
#         strip.placement = "outside",
#         legend.position = "bottom",
#         legend.title = element_text(size = 12),
#         legend.text = element_text(size = 10)) +
#     geom_errorbar(aes(ymin = accuracy_lowerCI, ymax = accuracy_upperCI),
#                   position = position_dodge(0.9), width = 0.1)
# 
# # save plot
# ggsave(filename = paste0(output_folder, "/fig_3c.png"),
#        plot = results_plot, device = "png", dpi = 350, units = "in", width = 8, height = 4)

# clean R environment
rm(list = ls())



#### figure 4 ----

##### a #####

prediction_summary_file <- "analysis/trait_sim/multi_env/prediction_results.summary.txt"
output_folder <- "figures"

# load data
prediction_summary <- fread(prediction_summary_file, header = TRUE, data.table = FALSE)

# set default
trait_var_source <- "both"
snp_sv_ratio <- 0.5
sv_effects <- 0.1
sv_diff_dist <- FALSE
error_bars <- "CI_accuracy"

# filter summary
prediction_summary_plot <- prediction_summary %>% 
  filter(var %in% trait_var_source & predictor %in% c("all", "snp", "sv")
         & ratio %in% snp_sv_ratio & sv_effect %in% sv_effects) %>%
  group_by(h2, qtn, var, ratio, sv_effect, predictor) %>%
  summarize(mean_pred = mean(mean_accuracy_envs, na.rm = TRUE),
            se_pred = sd(mean_accuracy_envs, na.rm = TRUE) / n()) %>% 
  ungroup()

# reorder levels of marker type for better visualization
prediction_summary_plot$predictor <- factor(prediction_summary_plot$predictor,
                                            levels = c("snp", "sv", "all"),
                                            labels = c("SNPs \nonly", "SVs\nonly", "SNPs\nand SVs"))

# plot
results_plot <- prediction_summary_plot %>%
  mutate(h2_qtn = paste0("h2 = ",h2, "\n", qtn, " QTLs")) %>%
  ggplot(aes(x = h2_qtn, y = mean_pred)) +
  # facet_nested(~ ratio + sv_effect,
  #              nest_line = element_line(linetype = 1, color = "gray80"),
  #              labeller = labeller(ratio = function(x) paste0(as.numeric(x) * 100, "% SNPs - ", (1 - as.numeric(x)) * 100, "% SVs"),
  #                                  sv_effect = c("0.1" = "effect\nSVs = SNPs", "0.5" = "effect\nSVs > SNPs"))) +
  geom_line(aes(color = predictor, group = predictor), size = 1) +
  geom_point(aes(color = predictor), size = 2) +
  coord_cartesian(ylim = c(0, 1)) +
  scale_x_discrete(drop = FALSE) +
  labs(x = "Trait architecture",
       y = "Prediction accuracy") +
  scale_color_manual(values = c("#A21E0DFF", "#F4A652FF", "#3F5B8AFF")) +
  guides(color = guide_legend("Predictors")) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 12),
        axis.title.x = element_text(vjust = -1.5),
        axis.title.y = element_text(vjust = 1.5),
        axis.text = element_text(size = 10),
        strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12),
        strip.background = element_blank(),
        strip.placement = "outside",
        # legend.position = "right",
        legend.position = "bottom",
        # legend.key.size = unit(0.2, "in"),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 10, hjust = 0.5)) +
  geom_errorbar(aes(ymin = mean_pred - se_pred, ymax = mean_pred + se_pred, color = predictor),
                width = 0.1)

# save plot
ggsave(filename = paste0(output_folder, "/fig_4b.png"),
       plot = results_plot, device = "png", dpi = 350, units = "in", width = 4, height = 3)
# ggsave(filename = paste0(output_folder, "/fig_4b.png"),
#        plot = results_plot, device = "png", dpi = 350, units = "in", width = 8, height = 3)


# clean R environment
rm(list = ls())


##### b ----

prediction_summary_file <- "analysis/trait_sim/multi_env/prediction_results.summary.txt"
output_folder <- "figures"

# load data
prediction_summary <- fread(prediction_summary_file, header = TRUE, data.table = FALSE)

# set default
trait_var_source <- "both"
snp_sv_ratio <- 0.5
sv_effects <- c(0.1, 0.5)
sv_diff_dist <- FALSE
error_bars <- "CI_accuracy"

# filter summary
prediction_summary_plot <- prediction_summary %>% 
  filter(var %in% trait_var_source & ratio %in% snp_sv_ratio & sv_effect %in% sv_effects) %>%
  group_by(h2, qtn, ratio, sv_effect, predictor) %>%
  summarize(mean_pred = mean(mean_accuracy_envs, na.rm = TRUE),
            se_pred = sd(mean_accuracy_envs, na.rm = TRUE) / n()) %>% 
  ungroup()

# reorder levels of marker type for better visualization
prediction_summary_plot$predictor <- factor(prediction_summary_plot$predictor,
                                            levels = c("snp", "sv", "all", "snp_ld", "snp_not_ld"),
                                            labels = c("SNPs only", "SVs only", "SNPs and SVs", 
                                                       "SNPs in LD", "SNPs not in LD"))
prediction_summary_plot$sv_effect <- factor(prediction_summary_plot$sv_effect,
                                            levels = c("0.1", "0.5"),
                                            labels = c("SVs = SNPs", "SVs > SNPs"))

# plot
results_plot <- ggplot(data = prediction_summary_plot,
       aes(x = sv_effect, y = mean_pred)) +
  facet_nested(~ h2 + qtn, nest_line = element_line(linetype = 1, color = "gray80"),
               labeller = labeller(h2 = function(x) paste0("h2 = ", x),
                                   qtn = function(x) paste0(x, " QTLs"))) +
  geom_line(aes(color = predictor, group = predictor), size = 1) +
  geom_point(aes(color = predictor), size = 2) +
  coord_cartesian(ylim = c(0, 1)) +
  scale_x_discrete(drop = FALSE) +
  labs(x = "Effect sizes of causative variants",
       y = "Prediction accuracy") +
  scale_color_manual(values = c("#A21E0DFF", "#F4A652FF", "#3F5B8AFF", "gray30", "gray70")) +
  guides(color = guide_legend("Predictors")) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 12),
        axis.title.x = element_text(vjust = -1.5),
        axis.title.y = element_text(vjust = 1.5),
        axis.text = element_text(size = 10),
        axis.text.x = element_text(size = 10, angle = 30, hjust = 1),
        strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.position = "right",
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10))

# save plot
ggsave(filename = paste0(output_folder, "/fig_4c.png"),
       plot = results_plot, device = "png", dpi = 350, units = "in", width = 8, height = 3)

# clean R environment
rm(list = ls())


##### c -----

prediction_summary_file <- "analysis/trait_sim/multi_env/prediction_results.summary.txt"
output_folder <- "figures"

# load data
prediction_summary <- fread(prediction_summary_file, header = TRUE, data.table = FALSE)

# set default
trait_var_source <- "both"
snp_sv_ratio <- c(0.5, 0.8)
sv_effects <- c(0.1, 0.5)
sv_diff_dist <- FALSE
error_bars <- "CI_accuracy"

# filter summary
prediction_summary_plot <- prediction_summary %>% 
  filter(var %in% trait_var_source & predictor %in% c("all", "snp", "sv")
         & ratio %in% snp_sv_ratio & sv_effect %in% sv_effects) %>%
  group_by(ratio, predictor, cv) %>%
  summarize(mean_pred = mean(mean_accuracy_envs, na.rm = TRUE),
            se_pred = sd(mean_accuracy_envs, na.rm = TRUE) / n()) %>% 
  ungroup()

# reorder levels of marker type for better visualization
prediction_summary_plot$predictor <- factor(prediction_summary_plot$predictor,
                                            levels = c("snp", "sv", "all"),
                                            labels = c("SNPs\nonly", "SVs\nonly", "SNPs\nand SVs"))
prediction_summary_plot$ratio <- factor(prediction_summary_plot$ratio,
                                        levels = c("0.5", "0.8"),
                                        labels = c("0.5" = "50% SNPs\n50% SVs",
                                                   "0.8" = "80% SNPs\n20% SVs"))

# plot
results_plot <- ggplot(data = prediction_summary_plot,
       aes(x = predictor, y = mean_pred)) +
  facet_nested(~ cv, nest_line = element_line(linetype = 1, color = "gray80")) +
  geom_col(aes(fill = ratio), position = "dodge") +
  coord_cartesian(ylim = c(0, 1)) +
  scale_x_discrete(drop = FALSE) +
  labs(x = "Marker type used in prediciton",
       y = "Prediction accuracy") +
  scale_fill_manual(values = c("#96B1CDFF", "#4C75A7FF")) +
  guides(fill = guide_legend("Proportion of\ncausative variants")) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 12),
        axis.title.x = element_text(vjust = -1.5),
        axis.title.y = element_text(vjust = 1.5),
        axis.text = element_text(size = 10),
        strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12),
        strip.background = element_blank(),
        strip.placement = "outside",
        # legend.position = "right",
        legend.position = "bottom",
        legend.title = element_text(size = 10, hjust = 0.5),
        legend.text = element_text(size = 10, hjust = 0.5)) +
  geom_errorbar(aes(ymin = mean_pred - se_pred, ymax = mean_pred + se_pred, fill = ratio),
                position = position_dodge(0.9), width = 0.1)

# save plot
ggsave(filename = paste0(output_folder, "/fig_4d.png"),
       plot = results_plot, device = "png", dpi = 350, units = "in", width = 4, height = 3)
# ggsave(filename = paste0(output_folder, "/fig_4d.png"),
#        plot = results_plot, device = "png", dpi = 350, units = "in", width = 8, height = 3)

# # filter summary
# prediction_summary_plot <- subset(prediction_summary,
#                                   qtn == 100
#                                   & h2 == 0.7
#                                   & var %in% trait_var_source
#                                   & predictor %in% c("all", "snp", "sv")
#                                   & ratio %in% snp_sv_ratio
#                                   & sv_effect %in% sv_effects
#                                   & diff_dist == sv_diff_dist)
# 
# # reorder levels of marker type for better visualization
# prediction_summary_plot$predictor <- factor(prediction_summary_plot$predictor,
#                                             levels = c("snp", "sv", "all"),
#                                             labels = c("SNPs\nonly", "SVs\nonly", "SNPs\nand SVs"))
# 
# # # add new label
# # prediction_summary_plot$sv_effect_label <- "Effect sizes"
# # prediction_summary_plot$ratio_label <- "Trait architecture"
# # 
# # ggplot(data = prediction_summary_plot,
# #        aes(x = predictor, y = mean_accuracy_envs)) +
# #   facet_nested(~ sv_effect_label + sv_effect, nest_line = element_line(linetype = 1, color = "gray80")) +
# #   geom_line(aes(linetype = as.factor(ratio), color = cv, group = interaction(cv, ratio)), size = 1) +
# #   geom_point(aes(color = cv), size = 2) +
# #   coord_cartesian(ylim = c(0, 1)) +
# #   scale_x_discrete(drop = FALSE) +
# #   labs(x = "Marker type used in prediciton",
# #        y = "Prediction accuracy") +
# #   scale_color_manual(values = c("#DE8282FF", "#AD0000FF")) +
# #   guides(color = guide_legend("Cross-Validation")) +
# #   theme_bw() +
# #   theme(panel.grid = element_blank(),
# #         axis.title = element_text(size = 12),
# #         axis.title.x = element_text(vjust = -1.5),
# #         axis.title.y = element_text(vjust = 1.5),
# #         axis.text = element_text(size = 10),
# #         strip.text.x = element_text(size = 12),
# #         strip.text.y = element_text(size = 12),
# #         strip.background = element_blank(),
# #         strip.placement = "outside",
# #         legend.position = "bottom",
# #         legend.title = element_text(size = 12),
# #         legend.text = element_text(size = 10)) 
# 
# # plot
# results_plot <- ggplot(data = prediction_summary_plot,
#                        aes(x = predictor, y = mean_accuracy_envs)) +
#   facet_nested(~ ratio + sv_effect,
#                nest_line = element_line(linetype = 1, color = "gray80"),
#                labeller = labeller(ratio = function(x) paste0(as.numeric(x) * 100, "% SNPs - ", (1 - as.numeric(x)) * 100, "% SVs"),
#                                    sv_effect = c("0.1" = "effect\nSVs = SNPs", "0.5" = "effect\nSVs > SNPs"))) +
#   geom_line(aes(color = cv, group = cv), size = 1) +
#   geom_point(aes(color = cv), size = 2) +
#   coord_cartesian(ylim = c(0, 1)) +
#   scale_x_discrete(drop = FALSE) +
#   labs(x = "Marker type used in prediciton",
#        y = "Prediction accuracy") +
#   scale_color_manual(values = c("#DE8282FF", "#AD0000FF")) +
#   guides(color = guide_legend("Cross-Validation")) +
#   theme_bw() +
#   theme(panel.grid = element_blank(),
#         axis.title = element_text(size = 12),
#         axis.title.x = element_text(vjust = -1.5),
#         axis.title.y = element_text(vjust = 1.5),
#         axis.text = element_text(size = 10),
#         strip.text.x = element_text(size = 12),
#         strip.text.y = element_text(size = 12),
#         strip.background = element_blank(),
#         strip.placement = "outside",
#         legend.position = "bottom",
#         legend.title = element_text(size = 12),
#         legend.text = element_text(size = 10)) 
# 
# 
# # add error bars
# if (error_bars == "SE_accuracy") {
#   results_plot <- results_plot +
#     geom_errorbar(aes(ymin = mean_accuracy_envs - accuracy_se, ymax = mean_accuracy_envs + accuracy_se),
#                   position = position_dodge(0.9), width = 0.1)
# }
# 
# if (error_bars == "CI_accuracy") {
#   results_plot <- results_plot +
#     geom_errorbar(aes(ymin = accuracy_lowerCI, ymax = accuracy_upperCI),
#                   position = position_dodge(0.9), width = 0.1)
# }
# 
# if (error_bars == "SE_error") {
#   results_plot <- results_plot +
#     geom_errorbar(aes(ymin = mean_accuracy_envs - mean_SE, ymax = mean_accuracy_envs + mean_SE),
#                   position = position_dodge(0.9), width = 0.1)
# }
# 
# if (error_bars == "CI_error") {
#   results_plot <- results_plot +
#     geom_errorbar(aes(ymin = mean_lowerCI, ymax = mean_upperCI),
#                   position = position_dodge(0.9), width = 0.1)
# }
# 
# # save plot
# ggsave(filename = paste0(output_folder, "/fig_4b.png"),
#        plot = results_plot, device = "png", dpi = 350, units = "in", width = 8, height = 4)

# clean R environment
rm(list = ls())



#### figure 5 ----

# files
ld_pred_accuracy_file <-"analysis/trait_sim/multi_env/pred-accuracy_qtn-effect_ld-pred-qtn.results.txt"
ld_summary_file <- "analysis/trait_sim/multi_env/pred-accuracy_qtn-effect_ld-pred-qtn.ld-summary.txt"
output_folder <- "figures"

# load files
results_accuracy_ld_pve <- fread(ld_pred_accuracy_file, header = TRUE, data.table = FALSE)
results_accuracy_ld_pve <- results_accuracy_ld_pve %>%
  select(!c(se, lowerCI, upperCI, qtn_type, var_exp_se, chr))
ld_summary_all_scenarios <- fread(ld_summary_file, header = TRUE, data.table = FALSE)

##### a #####

ld_accuracy <- results_accuracy_ld_pve %>%
  mutate(mean = if_else(mean < 0, true = 0, false = mean),
         accuracy_bins = cut_width(mean, width = 0.1, boundary = 0, closed = "right")) %>% 
  group_by(env, h2, qtn, var, ratio, sv_effect, diff_dist, pop, predictor, pred_iter, cv) %>%
  mutate(ld_NAs = sum(is.na(r2)),
         prop_high_ld = round(sum(r2 > 0.9, na.rm = TRUE) / (n() - ld_NAs), digits = 3)) %>%
  ungroup() %>% 
  group_by(accuracy_bins) %>%
  mutate(prop_high_ld_mean = mean(prop_high_ld)) %>%
  ungroup() %>% 
  filter(!is.na(accuracy_bins)) %>%
  mutate(accuracy_bins = factor(accuracy_bins,
                                levels = c("[0,0.1]", "(0.1,0.2]", "(0.2,0.3]", "(0.3,0.4]", "(0.4,0.5]",
                                           "(0.5,0.6]", "(0.6,0.7]", "(0.7,0.8]", "(0.8,0.9]", "(0.9,1]"),
                                labels = c("0-0.1", "0.1-0.2", "0.2-0.3", "0.3-0.4", "0.4-0.5",
                                           "0.5-0.6", "0.6-0.7", "0.7-0.8", "0.8-0.9", "0.9-1"))) %>% 
  ggplot(aes(x = accuracy_bins, y = r2)) +
  geom_boxplot(aes(fill = prop_high_ld_mean), outlier.size = 1) +
  stat_summary(fun.data = function(x) c(y = 1.03, label = length(x)), geom = "text", size = 2) +
  scale_x_discrete(drop = FALSE) +
  scale_fill_viridis_c(option = "B", begin = 0, limits = c(0, 1),
                       name = "Proportion of markers\nin high LD") +
  labs(x = "Prediction accuracy",
       y = bquote("LD"~(r^2))) +
  theme_bw() +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        axis.text.x = element_text(angle = 30, hjust = 1),
        legend.title = element_text(size = 10, hjust = 0.5, vjust = 0.75),
        legend.text = element_text(size = 8),
        legend.position = "bottom",
        panel.grid = element_blank())
ggsave(filename = paste0(output_folder, "/fig_5a.png"),
       plot = ld_accuracy, device = "png", dpi = 350, units = "in", width = 4, height = 4)

##### b #####

ld_dist <- results_accuracy_ld_pve %>%
  # group_by(qtn, var, ratio, pop, predictor, pred_iter, causal_var, predictor_highest_ld) %>%
  group_by(var, ratio, pop, predictor, pred_iter, causal_var, predictor_highest_ld) %>%
  summarize(r2 = mean(r2)) %>%
  ungroup() %>% 
  filter(predictor %in% c("snp", "sv") & var %in% c("SNP", "SV")) %>% 
  mutate(predictor = factor(predictor,
                            levels = c("snp", "sv"),
                            labels = c("SNPs only", "SVs only")),
         # qtn = factor(qtn, levels = c("10", "100"), labels = c("10 QTLs", "100 QTLs")),
         var = factor(var, levels = c("SNP", "SV"), labels = c("SNPs", "SVs")),
         predictor_label = "Predictors",
         causal_var_label = "Causative variants") %>% 
  # ggplot(aes(x = r2, fill = interaction(var, qtn, predictor))) +
  ggplot(aes(x = r2, fill = interaction(var, predictor))) +
  # facet_nested(causal_var_label + qtn + var ~ predictor_label + predictor, scales = "free_y",
  #              nest_line = element_line(linetype = 1, color = "gray80")) +
  facet_nested(causal_var_label + var ~ predictor_label + predictor, scales = "free_y",
               nest_line = element_line(linetype = 1, color = "gray80")) +
  geom_histogram(show.legend = FALSE) +
  coord_cartesian(xlim = c(0, 1)) +
  scale_x_continuous(labels = c("0", "0.25", "0.5", "0.75", "1")) +
  # scale_fill_manual(values = c("gray60", "gray60", "gray60", "firebrick",
  #                              "gray60", "gray60", "gray60", "gray60")) +
  scale_fill_manual(values = c("gray60", "firebrick",
                               "gray60", "gray60")) +
  labs(x = bquote("LD"~(r^2)),
       y = "Count") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 12),
        axis.title.x = element_text(vjust = -1.5),
        axis.title.y = element_text(vjust = 1.5),
        axis.text = element_text(size = 10),
        strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.position = "right",
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10))
ggsave(filename = paste0(output_folder, "/fig_5b.png"),
       plot = ld_dist, device = "png", dpi = 350, units = "in", width = 4, height = 4)

# clean R environment
rm(list = ls())



#### figure 6 ----

prediction_summary_file <- "analysis/ld_downsample/prediction_results.reps1-10.pops1-3.summary.txt"
output_folder <- "figures"

# load data
prediction_summary <- fread(prediction_summary_file, header = TRUE, data.table = FALSE)

# reorder levels of marker type for better visualization
prediction_summary$predictor <- factor(prediction_summary$predictor,
                                       levels = c("low", "moderate", "high"),
                                       labels = c("Low", "Moderate", "High"))

# set default
trait_var_source <- c("SNP", "SV", "both")
trait_heritability <- c(0.3 ,0.7)
error_bars <- "CI_accuracy"

# filter summary
prediction_summary_plot <- subset(prediction_summary, qtn == 100 & h2 %in% trait_heritability & var %in% trait_var_source)

# adjust labels
prediction_summary_plot$h2 <- as.character(sapply(prediction_summary_plot$h2, function(x) {
  label <- bquote(h^2 ~ "=" ~ .(x))
  return(label)
}))
prediction_summary_plot$h2_label <- "Heritability"
prediction_summary_plot$var_label <- "Causative variants"

# plot
results_plot <- ggplot(data = prediction_summary_plot,
                       aes(x = predictor, y = mean_accuracy_envs, fill = cv)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(y = 0.05, label = paste0("n=", pops)), position = position_dodge(0.9), size = 2) +
  facet_nested(h2_label + h2 ~ var_label + var,
               nest_line = element_line(linetype = 1, color = "gray80"),
               labeller = labeller(h2 = label_parsed,
                                   var = c("both" = "SNPs and SVs",
                                           "SNP" = "Only SNPs",
                                           "SV" = "Only SVs"))) +
  coord_cartesian(ylim = c(0, 1)) +
  scale_x_discrete(drop = FALSE) +
  labs(x = "Marker type used in prediciton",
       y = "Prediction accuracy") +
  scale_fill_manual(values = c("#DE8282FF", "#AD0000FF")) +
  guides(fill = guide_legend("Cross-Validation")) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 12),
        axis.title.x = element_text(vjust = -1.5),
        axis.title.y = element_text(vjust = 1.5),
        axis.text = element_text(size = 10),
        strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.position = "bottom",
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10)) 

# add error bars
if (error_bars == "SE_accuracy") {
  results_plot <- results_plot +
    geom_errorbar(aes(ymin = mean_accuracy_envs - accuracy_se, ymax = mean_accuracy_envs + accuracy_se),
                  position = position_dodge(0.9), width = 0.2)
}

if (error_bars == "CI_accuracy") {
  results_plot <- results_plot +
    geom_errorbar(aes(ymin = accuracy_lowerCI, ymax = accuracy_upperCI),
                  position = position_dodge(0.9), width = 0.2)
}

if (error_bars == "SE_error") {
  results_plot <- results_plot +
    geom_errorbar(aes(ymin = mean_accuracy_envs - mean_SE, ymax = mean_accuracy_envs + mean_SE),
                  position = position_dodge(0.9), width = 0.2)
}

if (error_bars == "CI_error") {
  results_plot <- results_plot +
    geom_errorbar(aes(ymin = mean_lowerCI, ymax = mean_upperCI),
                  position = position_dodge(0.9), width = 0.2)
}

# save plot
ggsave(filename = paste0(output_folder, "/fig_6.png"),
       plot = results_plot, device = "png", dpi = 350, units = "in", width = 8, height = 4)

# clean R environment
rm(list = ls())



#### supp fig 1 ----

# done in https://lucid.app/lucidchart/



#### supp fig 2 ----

output_folder <- "supp_materials"

# projection rate and accuracy
proj_summary_sv <- fread("analysis/projection_svs-snps/projection_svs_summary.txt",
                         header = TRUE, data.table = FALSE)
proj_summary_snp <- fread("analysis/projection_svs-snps/projection_reseq-snps_summary.txt",
                          header = TRUE, data.table = FALSE)
proj_summary <- rbind(data.frame(proj_summary_snp[, c(1, 3, 4)], marker = "SNP"),
                      data.frame(proj_summary_sv[, c(1, 5, 7)], marker = "SV"))
rm(proj_summary_sv, proj_summary_snp)

proj_rate_plot <- ggplot(proj_summary, aes(x = family, y = avg_percent_projected, fill = marker)) +
  geom_col(position = "dodge", width = 0.5) +
  labs(x = "Family", y = "Projection rate (%)", fill = "Marker type") +
  scale_y_continuous(labels = function(x) x*100, limits = c(0, 1))+
  scale_x_discrete(labels = function(cross) gsub("x", "\nx\n", cross)) +
  scale_fill_manual(values = c("gray70", "gray30")) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        axis.text.x = element_text(size = 8))
ggsave(filename = paste0(output_folder, "/supp_fig_2a.png"),
       plot = proj_rate_plot, device = "png", dpi = 350, units = "in", width = 8, height = 3)

proj_accu_plot <- ggplot(proj_summary, aes(x = family, y = proj_accuracy, fill = marker)) +
  geom_col(position = "dodge", width = 0.5) +
  labs(x = "Family", y = "Projection accuracy (%)", fill = "Marker type") +
  scale_y_continuous(labels = function(x) x*100, limits = c(0, 1))+
  scale_x_discrete(labels = function(cross) gsub("x", "\nx\n", cross)) +
  scale_fill_manual(values = c("gray70", "gray30")) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        axis.text.x = element_text(size = 8))
ggsave(filename = paste0(output_folder, "/supp_fig_2b.png"),
       plot = proj_accu_plot, device = "png", dpi = 350, units = "in", width = 8, height = 3)

# clean R environment
rm(list = ls())


#### supp fig 3 (msi) ----

# ran the code below on HPC due to size of file
marker_info_file <- "data/positions_projected_snps-svs.txt"
output_folder <- "supp_materials"

# load data
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
dist_densi <- ggplot(marker_data, aes(x = pos, fill = marker)) +
  geom_density(alpha = 0.5, position = "identity") +
  facet_wrap(~chr, nrow = 10, labeller = labeller(chr = function(x) paste0("chr", x))) +
  scale_x_continuous(labels = function(x) x/1000000) +
  scale_y_continuous(labels = function(y) y*1000000) +
  scale_color_manual(values = c("firebrick", "black"), aesthetics = "fill") +
  labs(x = "Position (Mb)",
       y = expression(paste("Density (x", 10^{6}, ")")),
       fill = "Marker type") +
  theme_bw() +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 10),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        legend.position = "bottom",
        panel.grid = element_blank())

ggsave(filename = "supp_fig_3.pdf",
       plot = dist_densi, device = "pdf", units = "in", width = 4, height = 10)




#### supp fig 4 ----

prediction_summary_file <- "analysis/trait_sim/multi_env/prediction_results.summary.txt"
output_folder <- "supp_materials"

# load data
prediction_summary <- fread(prediction_summary_file, header = TRUE, data.table = FALSE)

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
        strip.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        panel.grid = element_blank(),
        panel.spacing = unit(0.25, "lines"))
# save plot
ggsave(filename = paste0(output_folder, "/supp_fig_4.png"), plot = results_plot,
       device = "png", dpi = 350, units = "in", width = 8, height = 11)

# clean R environment
rm(list = ls())


#### supp fig 5 ----

ld_pred_accuracy_file <-"analysis/trait_sim/multi_env/pred-accuracy_qtn-effect_ld-pred-qtn.results.txt"
ld_summary_file <- "analysis/trait_sim/multi_env/pred-accuracy_qtn-effect_ld-pred-qtn.ld-summary.txt"
output_folder <- "supp_materials"

# load files
results_accuracy_ld_pve <- fread(ld_pred_accuracy_file, header = TRUE, data.table = FALSE)
results_accuracy_ld_pve <- results_accuracy_ld_pve %>%
  select(!c(se, lowerCI, upperCI, qtn_type, var_exp_se, chr))
ld_summary_all_scenarios <- fread(ld_summary_file, header = TRUE, data.table = FALSE)

ld_dist <- results_accuracy_ld_pve %>%
  group_by(qtn, var, ratio, pop, predictor, pred_iter, causal_var, predictor_highest_ld) %>%
  summarize(r2 = mean(r2)) %>%
  ungroup() %>% 
  mutate(predictor = factor(predictor,
                            levels = c("all", "snp", "sv", "snp_ld", "snp_not_ld"),
                            labels = c("SNPs and SVs", "Only SNPs", "Only SVs",
                                       "SNPs in LD", "SNPs not in LD")),
         qtn = factor(qtn, levels = c("10", "100"), labels = c("10 QTLs", "100 QTLs")),
         var = factor(var, levels = c("both", "SNP", "SV"), labels = c("SNPs and SVs", "SNPs", "SVs")),
         predictor_label = "Predictors",
         causal_var_label = "Causative variants") %>% 
  ggplot(aes(x = r2, fill = interaction(var, qtn, predictor))) +
  # facet_grid(qtn + var ~ predictor, scales = "free_y") +
  facet_nested(causal_var_label + qtn + var ~ predictor_label + predictor, scales = "free_y",
               nest_line = element_line(linetype = 1, color = "gray80")) +
  geom_histogram(show.legend = FALSE) +
  coord_cartesian(xlim = c(0, 1)) +
  scale_x_continuous(breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1")) +
  # scale_fill_manual(values = c("gray60", "gray60", "gray60", "gray60", "gray60", "gray60",
  #                              "gray60", "gray60", "gray60", "gray60", "gray60", "firebrick",
  #                              "gray60", "gray60", "gray60", "gray60", "gray60", "gray60",
  #                              "gray60", "gray60", "gray60", "gray60", "gray60", "gray60",
  #                              "gray60", "gray60", "gray60", "gray60", "gray60", "gray60")) +
  scale_fill_manual(values = rep("gray60", 30)) +
  labs(x = bquote("LD"~(r^2)), y = "Count") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.title = element_text(size = 12),
        axis.title.x = element_text(vjust = -1.5),
        axis.title.y = element_text(vjust = 1.5),
        axis.text = element_text(size = 10),
        strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size = 12),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.position = "right",
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10))
ggsave(filename = paste0(output_folder, "/supp_fig_5.png"),
       plot = ld_dist, device = "png", dpi = 350, units = "in", width = 8, height = 8)

# clean R environment
rm(list = ls())
