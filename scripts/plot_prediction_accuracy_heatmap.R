library(data.table)
library(ggplot2)
library(dplyr)
library(tidyr)
# install.packages("ggh4x")
# # for R.3.6 compatibility
# remotes::install_github("teunbrand/ggh4x@4e242295da442e0e412c7c3e713a383f603ceda0")  
library(ggh4x)

prediction_summary_file <- "tests/test_prediction_results.pops15-4-18.summary.txt"


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
