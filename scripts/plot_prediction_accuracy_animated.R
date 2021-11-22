library(data.table)
library(ggplot2)
library(dplyr)
library(tidyr)
# install.packages("ggh4x")
# # for R.3.6 compatibility
# remotes::install_github("teunbrand/ggh4x@4e242295da442e0e412c7c3e713a383f603ceda0")  
library(ggh4x)
library(gganimate)


# prediction_summary_file <- "tests/test_prediction_results.pops15-4-18.summary.txt"
prediction_summary_file <- "tests/test_prediction_results.15pops.summary.txt"


#### load prediction summary ----

# do causal variatnt 'both' separately
prediction_summary <- fread(prediction_summary_file, header = TRUE, data.table = FALSE)

# # reorder levels of marker type for better visualization
# prediction_summary$predictor <- factor(prediction_summary$predictor,
#                                        levels = c("all", "snp", "sv", "snp_ld", "snp_not_ld"),
#                                        labels = c("SNPs\nand SVs", "Only\nSNPs", "Only\nSVs",
#                                                   "SNPs\nin LD", "SNPs\nnot in LD"))
# 
# # create output folder
# outfolder <- unlist(strsplit(prediction_summary_file, "/"))
# outfolder[length(outfolder)] <- paste0("plots_", gsub(".txt", "", outfolder[length(outfolder)]))
# outfolder[length(outfolder)] <- gsub(".", "_", outfolder[length(outfolder)], fixed = TRUE)
# outfolder <- paste0(outfolder, collapse = "/")
# if (!dir.exists(outfolder)) dir.create(outfolder)



#### plot SNP only or SV only as causative variants ----

# set default
trait_var_source <- c("SNP", "SV")
trait_qtn_number <- c(10, 100)
trait_heritability <- c(0.3 ,0.7)

prediction_summary %>%
  filter(h2 %in% trait_heritability, var %in% trait_var_source, qtn %in% trait_qtn_number) %>%
  select(h2, var, qtn, predictor, cv, mean_accuracy_envs) %>%
  unite(h2:qtn, col = "sim_scenarios", sep = "/") %>%
  unite(predictor:cv, col = "pred_scenarios", sep = "/") %>%
  ggplot(aes(x = pred_scenarios, y = sim_scenarios, fill = mean_accuracy_envs)) +
  geom_tile() +
  scale_fill_viridis_c(option = "B", limits = c(0, 1))

prediction_summary %>% 
  filter(h2 %in% trait_heritability, var %in% trait_var_source, qtn %in% trait_qtn_number) %>%
  select(h2, var, qtn, predictor, cv, mean_accuracy_envs) %>%
  unite(h2:qtn, col = "sim_scenarios", sep = "/", remove = FALSE) %>% 
  unite(predictor:cv, col = "pred_scenarios", sep = "/", remove = FALSE) %>%
  ggplot(aes(x = pred_scenarios, y = sim_scenarios, fill = mean_accuracy_envs)) +
  facet_nested(h2 + qtn + var ~ predictor + cv, scales = "free", switch = "y",
               nest_line = element_line(linetype = 1, color = "gray80")) +
  geom_tile(aes(height = 1.5, width = 1.5)) +
  scale_fill_viridis_c(option = "B", limits = c(0, 1)) +
  theme_minimal() +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        panel.spacing = unit(0.25, "lines"))




#### plot both SNPs and SVs as causative variants ----

# set default
trait_var_source <- "both"
trait_qtn_number <- c(10, 100)
trait_heritability <- c(0.3 ,0.7)
snp_sv_ratio <- c(0.5, 0.8)
sv_effects <- c(0.1, 0.5)
sv_diff_dist <- c(FALSE, TRUE)

prediction_summary %>%
  filter(h2 %in% trait_heritability, var %in% trait_var_source, qtn %in% trait_qtn_number,
         ratio %in% snp_sv_ratio, sv_effect %in% sv_effects, diff_dist %in% sv_diff_dist) %>%
  select(h2, var, qtn, ratio, sv_effect, diff_dist, predictor, cv, mean_accuracy_envs) %>%
  unite(h2:diff_dist, col = "sim_scenarios", sep = "/") %>%
  unite(predictor:cv, col = "pred_scenarios", sep = "/") %>%
  ggplot(aes(x = pred_scenarios, y = sim_scenarios, fill = mean_accuracy_envs)) +
  geom_tile() +
  scale_fill_viridis_c(option = "B", limits = c(0, 1))

prediction_summary %>% 
  filter(h2 %in% trait_heritability, var %in% trait_var_source, qtn %in% trait_qtn_number,
         ratio %in% snp_sv_ratio, sv_effect %in% sv_effects, diff_dist %in% sv_diff_dist) %>%
  select(h2, var, qtn, ratio, sv_effect, diff_dist, predictor, cv, mean_accuracy_envs) %>%
  unite(h2:diff_dist, col = "sim_scenarios", sep = "/", remove = FALSE) %>% 
  unite(predictor:cv, col = "pred_scenarios", sep = "/", remove = FALSE) %>%
  ggplot(aes(x = pred_scenarios, y = sim_scenarios, fill = mean_accuracy_envs)) +
  facet_nested(h2 + qtn + var + ratio + sv_effect + diff_dist ~ predictor + cv, scales = "free", switch = "y",
               nest_line = element_line(linetype = 1, color = "gray80")) +
  geom_tile(aes(height = 1.5, width = 1.5)) +
  scale_fill_viridis_c(option = "B", limits = c(0, 1)) +
  theme_minimal() +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        panel.spacing = unit(0.25, "lines"))



########

# # install.packages("gganimate")
# library(gganimate)
# 
# p <- ggplot(iris, aes(x = Petal.Width, y = Petal.Length)) + 
#   geom_point()
# 
# anim <- p + 
#   transition_states(Species,
#                     transition_length = 2,
#                     state_length = 1)
# 
# anim_save(animation = anim, filename = "tests/animation_test.gif", renderer = gifski_renderer(loop = TRUE))
# anim_save(animation = anim, filename = "tests/animation_test.no-loop.gif", renderer = gifski_renderer(loop = FALSE))
# 
# anim <- ggplot(iris, aes(x = Petal.Width, y = Petal.Length)) + 
#   geom_point(aes(colour = Species), size = 2) + 
#   transition_states(Species,
#                     transition_length = 2,
#                     state_length = 1)
# 
# anim + 
#   enter_fade() + 
#   exit_shrink()




# library(gapminder)
# # Make a ggplot, but add frame=year: one image per year
# ggplot(gapminder, aes(gdpPercap, lifeExp, size = pop, colour = country)) +
#   geom_point(alpha = 0.7, show.legend = FALSE) +
#   scale_colour_manual(values = country_colors) +
#   scale_size(range = c(2, 12)) +
#   scale_x_log10() +
#   facet_wrap(~continent) +
#   # Here comes the gganimate specific bits
#   labs(title = 'Year: {frame_time}', x = 'GDP per capita', y = 'life expectancy') +
#   transition_time(year) +
#   ease_aes('linear')




#### animated heatmap ----

trait_var_source <- c("SNP", "SV")
trait_qtn_number <- c(10, 100)
trait_heritability <- c(0.3 ,0.7)

heatmap_plot <- prediction_summary %>% 
  filter(h2 %in% trait_heritability, var %in% trait_var_source, qtn %in% trait_qtn_number) %>%
  mutate(predictor = factor(predictor,
                            levels = c("all", "snp", "sv", "snp_ld", "snp_not_ld"),
                            labels = c("SNPs and SVs", "Only SNPs", "Only SVs",
                                       "SNPs in LD", "SNPs not in LD"))) %>%
  select(h2, var, qtn, predictor, cv, mean_accuracy_envs) %>%
  unite(c(h2, qtn), col = "h2_qtn", sep = "_", remove = FALSE) %>%
  unite(predictor:cv, col = "pred_scenarios", sep = "/", remove = FALSE) %>%
  filter(h2_qtn %in% c("0.3_10", "0.7_100")) %>% 
  ggplot(aes(x = pred_scenarios, y = var, fill = mean_accuracy_envs)) +
  facet_nested(var ~ predictor + cv, scales = "free", switch = "y",
               nest_line = element_line(linetype = 1, color = "gray80")) +
  geom_tile(aes(height = 1.5, width = 1.5)) +
  scale_fill_viridis_c(option = "B", limits = c(0, 1), name = "Accuracy") +
  theme_minimal() +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        panel.grid = element_blank(),
        # panel.spacing = unit(0.25, "lines"),
        panel.spacing = unit(1.5, "lines"),
        # strip.text = element_text(size = 15),
        # legend.title = element_text(size = 12),
        # legend.text = element_text(size = 10),
        strip.text = element_text(size = 50),
        legend.title = element_text(size = 40),
        legend.text = element_text(size = 35),
        legend.key.size = unit(6, 'lines'), 
        legend.key.height = unit(6, 'lines'),
        legend.key.width = unit(6, 'lines')) 

anim_heatmap_plot <- heatmap_plot +
  transition_states(states = h2_qtn, transition_length = 4, state_length = 1, wrap = FALSE)
  # transition_states(states = h2_qtn, transition_length = 1, state_length = 1, wrap = FALSE)

anim_heatmap <- animate(anim_heatmap_plot, nframes = 100, fps = 10, width = 4000, height = 1000, renderer = gifski_renderer(loop = FALSE))
anim_save("animated_heatmap.15pops.gif", animation = anim_heatmap)

# anim_heatmap <- animate(anim_heatmap_plot, nframes = 150, fps = 25, width = 4000, height = 1000, renderer = ffmpeg_renderer(format = "mp4"))
# anim_save("animated_heatmap.mp4", animation = anim_heatmap)





#### animated bar plot 1 ----

# set default
trait_var_source <- c("SNP", "SV")
trait_qtn_number <- 100
trait_heritability <- 0.7

barplot <- prediction_summary %>% 
  filter(h2 %in% trait_heritability, var %in% trait_var_source, qtn %in% trait_qtn_number) %>%
  mutate(predictor = factor(predictor,
                            levels = c("all", "snp", "sv", "snp_ld", "snp_not_ld"),
                            labels = c("SNPs\nand SVs", "Only\nSNPs", "Only\nSVs",
                                       "SNPs\nin LD", "SNPs\nnot in LD"))) %>%
  ggplot(aes(x = predictor, y = mean_accuracy_envs, fill = cv)) +
  geom_bar(stat = "identity", position = "dodge") +
  coord_cartesian(ylim = c(0, 1)) +
  scale_x_discrete(drop = FALSE) +
  labs(x = "Marker type used in prediciton",
       y = "Prediction accuracy") +
  scale_fill_manual(values = c("#DE8282FF", "#AD0000FF")) +
  guides(fill = guide_legend("Cross-Validation")) +
  theme_bw() +
  theme(title = element_text(size = 30),
        text = element_text(size = 30),
        axis.text.y = element_text(size = 30),
        strip.text.x = element_text(size = 30),
        legend.position = "top",
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 17)) +
  geom_errorbar(aes(ymin = accuracy_lowerCI, ymax = accuracy_upperCI),
                position = position_dodge(0.9), width = 0.2) 

anim_barplot <- barplot +
  transition_states(states = var, transition_length = 1, state_length = 1, wrap = FALSE)

anim_bars <- animate(anim_barplot, nframes = 300, fps = 50, width = 1000, height = 1000, renderer = gifski_renderer(loop = FALSE))
anim_save("animated_bars.15pops.gif", animation = anim_bars)



#### animated bar plot 2 ----

trait_var_source <- "both"
trait_qtn_number <- 100
trait_heritability <- 0.7
snp_sv_ratio <- 0.5
sv_effects <- c(0.1, 0.5)
sv_diff_dist <- FALSE

barplot_both <- prediction_summary %>% 
  filter(h2 %in% trait_heritability, var %in% trait_var_source, qtn %in% trait_qtn_number,
         ratio %in% snp_sv_ratio, sv_effect %in% sv_effects, diff_dist %in% sv_diff_dist) %>% 
  mutate(predictor = factor(predictor,
                            levels = c("all", "snp", "sv", "snp_ld", "snp_not_ld"),
                            labels = c("SNPs\nand SVs", "Only\nSNPs", "Only\nSVs",
                                       "SNPs\nin LD", "SNPs\nnot in LD"))) %>%
  ggplot(aes(x = predictor, y = mean_accuracy_envs, fill = cv)) +
  geom_bar(stat = "identity", position = "dodge") +
  coord_cartesian(ylim = c(0, 1)) +
  scale_x_discrete(drop = FALSE) +
  labs(x = "Marker type used in prediciton",
       y = "Prediction accuracy") +
  scale_fill_manual(values = c("#DE8282FF", "#AD0000FF")) +
  guides(fill = guide_legend("Cross-Validation")) +
  theme_bw() +
  theme(title = element_text(size = 30),
        text = element_text(size = 30),
        axis.text.y = element_text(size = 30),
        strip.text.x = element_text(size = 30),
        legend.position = "top",
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 17)) +
  geom_errorbar(aes(ymin = accuracy_lowerCI, ymax = accuracy_upperCI),
                position = position_dodge(0.9), width = 0.2) 

anim_barplot_both <- barplot_both +
  transition_states(states = sv_effect, transition_length = 1, state_length = 1, wrap = FALSE)

anim_bars_both <- animate(anim_barplot_both, nframes = 300, fps = 50, width = 1000, height = 1000, renderer = gifski_renderer(loop = FALSE))
anim_save("animated_bars_both.15pops.gif", animation = anim_bars_both)




trait_var_source <- "both"
trait_qtn_number <- 100
trait_heritability <- 0.7
snp_sv_ratio <- c(0.5, 0.8)
sv_effects <- 0.5
sv_diff_dist <- FALSE

barplot_ratio <- prediction_summary %>% 
  filter(h2 %in% trait_heritability, var %in% trait_var_source, qtn %in% trait_qtn_number,
         ratio %in% snp_sv_ratio, sv_effect %in% sv_effects, diff_dist %in% sv_diff_dist) %>% 
  mutate(predictor = factor(predictor,
                            levels = c("all", "snp", "sv", "snp_ld", "snp_not_ld"),
                            labels = c("SNPs\nand SVs", "Only\nSNPs", "Only\nSVs",
                                       "SNPs\nin LD", "SNPs\nnot in LD"))) %>%
  ggplot(aes(x = predictor, y = mean_accuracy_envs, fill = cv)) +
  geom_bar(stat = "identity", position = "dodge") +
  coord_cartesian(ylim = c(0, 1)) +
  scale_x_discrete(drop = FALSE) +
  labs(x = "Marker type used in prediciton",
       y = "Prediction accuracy") +
  scale_fill_manual(values = c("#DE8282FF", "#AD0000FF")) +
  guides(fill = guide_legend("Cross-Validation")) +
  theme_bw() +
  theme(title = element_text(size = 30),
        text = element_text(size = 30),
        axis.text.y = element_text(size = 30),
        strip.text.x = element_text(size = 30),
        legend.position = "top",
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 17)) +
  geom_errorbar(aes(ymin = accuracy_lowerCI, ymax = accuracy_upperCI),
                position = position_dodge(0.9), width = 0.2) 

anim_barplot_ratio <- barplot_ratio +
  transition_states(states = ratio, transition_length = 1, state_length = 1, wrap = FALSE)

anim_bars_ratio <- animate(anim_barplot_ratio, nframes = 300, fps = 50, width = 1000, height = 1000, renderer = gifski_renderer(loop = FALSE))
anim_save("animated_bars_ratio.15pops.gif", animation = anim_bars_ratio)

