library(data.table)
library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)

getMarkersHighestLD <- function(ld_summary) {
  
  # create empty df to store results
  highest_ld_results <- data.frame()
  
  # get predictor in highest LD with a causal variant
  for (var in unique(ld_summary$causal_var)) {
    
    # keep only respective causal variant
    ld_summary_causal <- subset(ld_summary, causal_var == var)
    
    # if all r2 are NA (i.e. no predictor in ld)...
    if (sum(is.na(ld_summary_causal$r2)) == NROW(ld_summary_causal)) {
      # ...get closest marker position
      pred_highest_ld <- ld_summary_causal[which.min(abs(ld_summary_causal$bp_distance)), ]
    } else {
      # if not, get the predictor with highest ld
      pred_highest_ld <- ld_summary_causal[which(ld_summary_causal$r2 == max(ld_summary_causal$r2, na.rm = TRUE)), ]
      # if more than one predictor with same highest ld...
      if (NROW(pred_highest_ld) > 1) {
        # ...get the closest one
        pred_highest_ld <- pred_highest_ld[which.min(abs(pred_highest_ld$bp_distance)), ]
      }
    }
    
    # append to results df
    highest_ld_results <- rbind(highest_ld_results, pred_highest_ld)
    
  }
  # change second column name
  colnames(highest_ld_results)[2] <- "predictor_highest_ld"
  
  return(highest_ld_results)
  
}




pred_accuracy_file <- "analysis/trait_sim/multi_env/full_prediction_results.iter1.pops1-5.txt"
folder_results <- "analysis/trait_sim/multi_env/no_gxe/additive_model/equal_effects"
trait_qtn_number <- c("10", "100")
trait_var_source <- c("SNP", "SV")
trait_pops <- 1:5
prediction_marker_type <- c("all", "sv", "snp", "snp_ld", "snp_not_ld")
predictor_datasets <- 1

# pred_accuracy_file <- "tests/test_full_prediction_results.txt"
# folder_results <- "tests"
# trait_qtn_number <- 100
# trait_var_source <-"SNP"
# trait_pops <- 1
# prediction_marker_type <- "all"
# predictor_datasets <- 1




# load prediction accuracy results
pred_accuracy <- fread(pred_accuracy_file, header = TRUE, data.table = FALSE)

# create empty df to store results
pred_accuracy_ld_summary <- data.frame()

for (qtns in trait_qtn_number) {
  for (causal_var_type in trait_var_source) {
    for (population in trait_pops) {
      for (predictors in prediction_marker_type) {
        for (n_dataset in predictor_datasets) {
          
          # get ld summary file
          ld_summary_file <- paste0(folder_results, "/", qtns, "-QTNs_from_",
                                    causal_var_type, "/ld_causative-vars_predictors/pop",
                                    population, "/ld_summary.", predictors,
                                    "_markers.pred-iter", n_dataset, ".causal-pop",
                                    population, ".txt")
          # load ld file
          ld_summary <- fread(ld_summary_file, header = TRUE, data.table = FALSE)
          
          # get predictors in highest ld to causative variants
          highest_ld_results <- getMarkersHighestLD(ld_summary)
          
          
          # # write results
          # outfile_highest_ld <- paste0(folder_results, "/", qtns, "-QTNs_from_",
          #                              causal_var_type, "/ld_causative-vars_predictors/pop",
          #                              population, "/pred_highest_ld_causal-vars.txt")
          # outfile_highest_ld <- "analysis/trait_sim/multi_env/pred_highest_ld_causal-vars.example.txt"
          # fwrite(highest_ld_results, outfile_highest_ld, quote = FALSE, sep = "\t", na = NA, row.names = FALSE)
          # highest_ld_results <- fread("tests/pred_highest_ld_causal-vars.example.txt", header = TRUE, data.table = FALSE)

          
          
          ############# SHOULD I COUNT HOW MANY IN HIGH LD?
         
          # count how many predictors in high ld (>0.8) or moderate ld (0.5-0.8) and
          # causative variants, and get the mean/median distance between them
          ld_high <- highest_ld_results$r2 >= 0.8
          ld_moderate <- highest_ld_results$r2 >= 0.5 & highest_ld_results$r2 < 0.8

          highest_ld_results_summary <- data.frame(high_ld_n = sum(ld_high, na.rm = TRUE),
                                                   high_ld_perc = round(sum(ld_high, na.rm = TRUE) / NROW(highest_ld_results), 2),
                                                   high_ld_bp_mean = mean(abs(highest_ld_results[which(ld_high), "bp_distance"])),
                                                   high_ld_bp_median = median(abs(highest_ld_results[which(ld_high), "bp_distance"])),
                                                   moderate_ld_n = sum(ld_moderate, na.rm = TRUE),
                                                   moderate_ld_perc = round(sum(ld_moderate, na.rm = TRUE) / NROW(highest_ld_results), 2),
                                                   moderate_ld_bp_mean = mean(abs(highest_ld_results[which(ld_moderate), "bp_distance"])),
                                                   moderate_ld_bp_median = median(abs(highest_ld_results[which(ld_moderate), "bp_distance"])))
          
          ############# OR SHOULD I USE THE ACTUAL R2 VALUES?
          
          # highest_ld_results_filtered <- highest_ld_results[!is.na(highest_ld_results$r2), ]
          # highest_ld_results_summary <- data.frame(r2 = paste0(highest_ld_results_filtered$r2, collapse = ";"),
          #                                          bp_distance = paste0(highest_ld_results_filtered$bp_distance, collapse = ";"))
          
          
          
          
          # keep only trait combination of interest
          pred_accuracy_subset <- subset(pred_accuracy, qtn == qtns & var == causal_var_type
                                         & pop == population & predictor == predictors)
          
          # append to prediction accuracy results
          pred_accuracy_ld_summary <- rbind(pred_accuracy_ld_summary,
                                            cbind(pred_accuracy_subset, highest_ld_results_summary))
          
          
        }
      }
    }
  }
}

# # write final summary
# outfile <- paste0(folder_results, "/pred_accuracy_ld_summary.txt")
# outfile <- paste0("analysis/trait_sim/multi_env/pred_accuracy_ld_summary.iter1.pops1-5.txt")
# fwrite(pred_accuracy_ld_summary, outfile, quote = FALSE, sep = "\t", na = NA, row.names = FALSE)

#### for testing...
# pred_accuracy_ld_summary <- fread("tests/pred_accuracy_ld_summary.iter1.pops1-5.txt", header = TRUE, data.table = FALSE)
####

# pred_accuracy_ld_summary <- separate_rows(pred_accuracy_ld_summary, r2, bp_distance, sep = ";")
# pred_accuracy_ld_summary$r2 <- as.numeric(pred_accuracy_ld_summary$r2)
# pred_accuracy_ld_summary$bp_distance <- as.numeric(pred_accuracy_ld_summary$bp_distance)
# pred_accuracy_ld_summary$pop <- as.factor(pred_accuracy_ld_summary$pop)



pred_accuracy_ld_summary %>% 
  filter(predictor == "all") %>%
  ggplot(aes(x = high_ld_perc, y = mean_accuracy_envs)) + 
  facet_grid(cv + var ~ qtn + h2) +
  geom_point(aes(color = predictor)) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1))


pred_accuracy_ld_summary %>% 
  ggplot(aes(x = high_ld_perc, y = mean_accuracy_envs)) + 
  facet_grid(cv + var ~ qtn + h2) +
  geom_point(aes(color = predictor)) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1))

pred_accuracy_ld_summary %>% 
  ggplot(aes(x = high_ld_perc, y = mean_accuracy_envs, color = as.factor(h2), shape = as.factor(qtn))) +
  geom_point() +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  # geom_smooth(method = "lm")
  geom_smooth()

pred_accuracy_ld_summary %>% 
  ggplot(aes(x = high_ld_perc, y = mean_accuracy_envs)) +
  facet_grid(qtn ~ h2) +
  geom_point() +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
  # geom_smooth(method = "lm")
  geom_smooth()


# figure: perhaps a heatmap with all possible combinations (x axis is trait combinations, y axis is predictor type)

#### correlation accounting for causative variant
cor_pred_accuracy_ld <- pred_accuracy_ld_summary %>%
  group_by(gxe, h2, qtn, var, ratio, sv_effect, diff_dist, predictor, cv) %>%
  summarize(pops = n(),
            cor_accuracy_high_ld = cor(mean_accuracy_envs, high_ld_perc,
                                       method = "pearson", use = "complete.obs")) %>% 
  ungroup()

heatmap_data <- cor_pred_accuracy_ld  %>% 
  select(-c(pops, ratio, sv_effect, diff_dist)) %>%
  mutate_at(vars(-cor_accuracy_high_ld), as.character) %>% 
  mutate(gxe = paste0(gxe, "_GxE"),
         h2 = paste0("h2_", h2),
         qtn = paste0(qtn, "_QTNs"),
         var = paste0(var, "_as_causal")) 

ggplot(heatmap_data, aes(x = var, y = predictor, fill = cor_accuracy_high_ld)) + 
  facet_grid(cv ~ gxe + qtn + h2) +
  geom_tile() +
  scale_fill_gradientn(colours = brewer.pal(7, "RdBu"), limits = c(-1, 1),
                       na.value = "grey50", aesthetics = "fill")

ggplot(heatmap_data, aes(x = var, y = predictor, fill = cor_accuracy_high_ld)) + 
  facet_grid(cv ~ qtn + h2) +
  geom_tile() +
  scale_fill_gradientn(colours = brewer.pal(7, "RdBu"), limits = c(-1, 1),
                       na.value = "grey50", aesthetics = "fill")

#### correlation not accounting for causative variant
cor_pred_accuracy_ld <- pred_accuracy_ld_summary %>%
  group_by(gxe, h2, qtn, ratio, sv_effect, diff_dist, predictor, cv) %>%
  summarize(cor_accuracy_high_ld = cor(mean_accuracy_envs, high_ld_perc,
                                       method = "pearson", use = "complete.obs")) %>% 
  ungroup()

heatmap_data <- cor_pred_accuracy_ld  %>% 
  select(-c(ratio, sv_effect, diff_dist)) %>%
  mutate_at(vars(-cor_accuracy_high_ld), as.character) %>% 
  mutate(gxe = paste0(gxe, "_GxE"),
         h2 = paste0("h2_", h2),
         qtn = paste0(qtn, "_QTNs")) 

ggplot(heatmap_data, aes(x = cv, y = predictor, fill = cor_accuracy_high_ld)) + 
  facet_grid(gxe ~ qtn + h2) +
  geom_tile() +
  scale_fill_gradientn(colours = brewer.pal(7, "RdBu"), limits = c(-1, 1),
                       na.value = "grey50", aesthetics = "fill")

ggplot(heatmap_data, aes(x = cv, y = predictor, fill = cor_accuracy_high_ld)) + 
  facet_grid(qtn ~ h2) +
  geom_tile() +
  scale_fill_gradientn(colours = brewer.pal(7, "RdBu"), limits = c(-1, 1),
                       na.value = "grey50", aesthetics = "fill")


#### correlation not accounting for causative variant nor predictor type
cor_pred_accuracy_ld <- pred_accuracy_ld_summary %>%
  group_by(gxe, h2, qtn, ratio, sv_effect, diff_dist, cv) %>%
  summarize(cor_accuracy_high_ld = cor(mean_accuracy_envs, high_ld_perc,
                                       method = "pearson", use = "complete.obs")) %>% 
  ungroup()

heatmap_data <- cor_pred_accuracy_ld  %>% 
  select(-c(ratio, sv_effect, diff_dist)) %>%
  mutate_at(vars(-cor_accuracy_high_ld), as.character) %>% 
  mutate(gxe = paste0(gxe, "_GxE"),
         h2 = paste0("h2_", h2),
         qtn = paste0(qtn, "_QTNs")) 

ggplot(heatmap_data, aes(x = qtn, y = h2, fill = cor_accuracy_high_ld)) +
  facet_grid(cv ~ gxe) +
  geom_tile() +
  scale_fill_gradientn(colours = brewer.pal(7, "RdBu"), limits = c(-1, 1),
                       na.value = "grey50", aesthetics = "fill")

ggplot(heatmap_data, aes(x = qtn, y = h2, fill = cor_accuracy_high_ld)) +
  facet_grid(~ cv) +
  geom_tile() +
  scale_fill_gradientn(colours = brewer.pal(7, "RdBu"), limits = c(-1, 1),
                       na.value = "grey50", aesthetics = "fill")
