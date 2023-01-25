library(data.table)
library(tidyr)
library(dplyr)
library(asreml)
library(gtools)
# library(ggplot2)

usage <- function() {
  cat("
description: get BLUPs for each environment in simplePHENOTYPES.

usage: Rscript get_blups_per_env.R [sim_trait_filename] [outfolder]

positional arguments:
  sim_trait_filename      simulated traits from simplePHENOTYPES
  outfolder               output folder name

optional argument:
  --help                  show this helpful message
  --envs-weight           add weights for each environment based on the inverse
                          of the variance of the predictions

credits:
  part of this script was modified from Fernandes et al. 2018
  (TAG, doi:10.1007/s00122-017-3033-y), available at
  https://github.com/samuelbfernandes/Trait-assisted-GS
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

# set default
envs_weight <- FALSE

args <- commandArgs(trailingOnly = TRUE)
if ("--help" %in% args) usage() & q(save = "no")

# assert to have the correct optional arguments
if (length(args) < 2) stop(usage(), "missing positional argument(s)")

if (length(args) > 2) {

  opt_args <- args[-1:-2]
  opt_args_allowed <- c("--envs-weight")
  opt_args_requested <- as.character(sapply(opt_args, function(x) unlist(strsplit(x, split = "="))[1]))
  if (any(!opt_args_requested %in% opt_args_allowed)) stop(usage(), "wrong optional argument(s)")

  # change default based on the argument provided
  for (argument in opt_args_allowed) {
    if (any(grepl(argument, opt_args_requested))) {
      arg_name <- gsub("-", "_", gsub("--", "", argument))
      arg_value <- getArgValue(opt_args[grep(argument, opt_args)])
      assign(arg_name, arg_value)
    }
  }

}

# get positional arguments
sim_trait_filename <- args[1]
outfolder <- args[2]



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
results_1st_stage <- data.frame()

# get blups for each environment
for (env in levels(pheno_pop$environment)) {

  # filter data by environment
  pheno_pop_env <- droplevels(subset(pheno_pop, environment == env))

  # run mixed modes
  sim_trait_env <- asreml(fixed = trait_value ~ genotype,
                          random = ~ rep,
                          maxit = 100,
                          data = pheno_pop_env,
                          trace = FALSE)
  # update model until converge and % change in parameters is lower than 0.1
  while (!sim_trait_env$converge | any(summary(sim_trait_env)$varcomp$`%ch` > 0.1, na.rm = TRUE)) {
    cat("updating model...\n")
    sim_trait_env <- update.asreml(sim_trait_env, trace = FALSE)
    cat("...done\n")
  }

  # get predicted values
  predicted_values_env <- predict(sim_trait_env, classify = "genotype", sed = TRUE, maxit = 100)
  if (!envs_weight) {
    predicted_values_env <- data.frame(predicted_values_env$pvals[, c(1, 2)])
  } else {
    predicted_values_env <- data.frame(predicted_values_env$pvals[, c(1, 2, 3)])
  }

  # append results to main data frame
  results_1st_stage <- rbind(results_1st_stage, data.frame(environment = env, predicted_values_env))

}

# readjust columns
colnames(results_1st_stage)[3] <- "predicted_trait_value"
if (!envs_weight) {
  results_1st_stage <- results_1st_stage[, c("genotype", "environment", "predicted_trait_value")]
} else {
  # calculate weights based on the inverse of the variance of the predictions
  results_1st_stage[, 4] <- 1 / (results_1st_stage[, 4]^2)
  colnames(results_1st_stage)[4] <- "weight"
  results_1st_stage <- results_1st_stage[, c("genotype", "environment", "predicted_trait_value", "weight")]
}

# save adjusted means for 2nd stage analysis
if (!envs_weight) {
  fwrite(results_1st_stage, file = paste0(outfolder, "/blups_1st_stage.txt"),
         quote = FALSE, sep = "\t", na = NA, row.names = FALSE)
} else {
  fwrite(results_1st_stage, file = paste0(outfolder, "/blups_1st_stage_weighted.txt"),
         quote = FALSE, sep = "\t", na = NA, row.names = FALSE)
}



#### debug ----

# sim_trait_filename <- "analysis/test_prediction/multi_env/with_gxe/100qtns_SVs_equal_0.5h2_pop1/Simulated_Data_3_Reps_Herit_0.5_0.5_0.5_0.5_0.5_0.5_0.5_0.5_0.5_0.5_0.5_0.5_0.5_0.5_0.5_0.5_0.5_0.5_0.5_0.5.txt"
# outfolder <- "analysis/test_prediction/multi_env/with_gxe/100qtns_SVs_equal_0.5h2_pop1"
# envs_weight <- TRUE
