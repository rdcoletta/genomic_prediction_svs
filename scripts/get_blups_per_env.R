library(data.table)
library(tidyr)
library(dplyr)
library(ggplot2)
library(asreml)
library(asremlPlus)
library(gtools)

usage <- function() {
  cat("
description: get BLUPs for each environment in simplePHENOTYPES.

usage: Rscript get_blups_per_env.R [sim_trait_filename] [outfolder]

positional arguments:
  sim_trait_filename      simulated traits from simplePHENOTYPES
  outfolder               output folder name

optional argument:
  --help                  show this helpful message

credits:
  part of this script was modified from Fernandes et al. 2018 (TAG, doi:10.1007/s00122-017-3033-y),
  available at https://github.com/samuelbfernandes/Trait-assisted-GS
"
  )
}



#### command line options ----

args <- commandArgs(trailingOnly = TRUE)
if ("--help" %in% args) usage() & q(save = "no")

# assert to have the correct optional arguments
if (length(args) != 2) stop(usage(), "missing positional argument(s)")

# get positional arguments
sim_trait_filename <- args[1]
outfolder <- args[2]
# sim_trait_filename <- "analysis/test_prediction/multi_env/with_gxe/100qtns_SVs_equal_0.5h2_pop1/Simulated_Data_3_Reps_Herit_0.5_0.5_0.5_0.5_0.5_0.5_0.5_0.5_0.5_0.5_0.5_0.5_0.5_0.5_0.5_0.5_0.5_0.5_0.5_0.5.txt"
# outfolder <- "analysis/test_prediction/multi_env/with_gxe/100qtns_SVs_equal_0.5h2_pop1"


#### get BLUPs for each environment ----

# load data
pheno_pop <- fread(sim_trait_filename, header = TRUE, data.table = FALSE)

# adjust column names
colnames(pheno_pop)[1] <- "genotype"
colnames(pheno_pop)[NCOL(pheno_pop)] <- "rep"
colnames(pheno_pop)[2:(NCOL(pheno_pop) - 1)] <- paste0("env", 1:length(2:(NCOL(pheno_pop) - 1)))
# transform df into long format
pheno_pop <- pivot_longer(pheno_pop, !c(genotype, rep), names_to = "environment", values_to = "trait_value") %>%
  mutate(rep = paste0("rep", rep)) %>% 
  relocate(rep, .after = environment) %>%
  arrange(environment, genotype) %>% 
  mutate(genotype = factor(genotype),
         environment = factor(environment, levels = mixedsort(unique(environment))),
         rep = factor(rep, levels = mixedsort(unique(rep))))

# # histogram of raw data
# ggplot(pheno_pop, aes(x = trait_value)) +
#   geom_histogram(binwidth = 1) +
#   labs(x = "Simulated trait values",
#        y = "Count")
# 
# # stripchart of raw data across environments
# stripchart(pheno_pop$trait_value ~ pheno_pop$environment,
#            xlab = "Environment",
#            ylab = "Simulated trait",
#            vertical = TRUE,
#            method = "jitter",
#            bg = "cadetblue",
#            pch = 21)

# create empty data frame
results_1st_stage <- data.frame(environment = character(),
                      genotype = character(),
                      predicted_values = numeric())

for (env in levels(pheno_pop$environment)) {
  
  # filter by environment
  pheno_pop_env <- droplevels(subset(pheno_pop, environment == env))
  
  # run mixed modes
  sim_trait_env <- asreml(fixed = trait_value ~ genotype,
                          random = ~ rep,
                          maxit = 300,
                          data = pheno_pop_env)
  sim_trait_env <- update.asreml(sim_trait_env)
  
  # get AIC and BIC
  print(cbind(environment = env, infoCriteria.asreml(sim_trait_env)))
  
  # predicted values
  predicted_values_env <- predict(sim_trait_env, classify = "genotype", sed = TRUE, maxit = 1000)
  predicted_values_env <- data.frame(predicted_values_env$pvals[, c(1, 2)])
  
  # append results to main data frame
  results_1st_stage <- rbind(results_1st_stage, data.frame(environment = env, predicted_values_env))
  
}

# readjust columns
colnames(results_1st_stage) <- c("environment", "genotype", "predicted_trait_value")
results_1st_stage <- results_1st_stage[, c("genotype", "environment", "predicted_trait_value")]

# save adjusted means for 2nd stage analysis
fwrite(results_1st_stage, file = paste0(outfolder, "/blups_1st_stage.txt"),
       quote = FALSE, sep = "\t", na = NA, row.names = FALSE)


# more info about ASReml predict function:
# https://asreml.kb.vsni.co.uk/knowledge-base/predict-asreml/
