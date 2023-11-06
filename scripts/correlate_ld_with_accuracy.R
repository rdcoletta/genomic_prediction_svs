library(data.table)
library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)

usage <- function() {
  cat("
description: summarize and correlate LD between causative variants and prediction accuracy.

usage: Rscript correlate_ld_with_accuracy.R [folder_results] [pred_accuracy_file] [outfile] [...]

positional arguments:
  folder_results                path to folder with results of cross validation
  pred_accuracy_file            path to file with summary of prediction accuracy
  outfile                       name of output file

optional argument:
  --help                        show this helpful message
  --trait-pops=[LIST]           comma-separated list of pop_ns of QTNs used for simulating traits
  --pred-iters=[LIST]           comma-separated list of predictors iterations used in genomic prediction

note: make sure that to provide the same number of `--trait-pops` and `--pred-iters`

"
  )
}

getArgValue <- function(arg) {

  # get value from command-line arguments
  arg <- unlist(strsplit(arg, "="))
  if (length(arg) == 1) return(TRUE)
  if (length(arg) > 1) return(arg[2])

}

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



#### command line options ----

# set default
qtn_opts <- c(10, 100)
h2_opts <- c(0.3, 0.7)
predictor_opts <- c("all", "sv", "snp", "snp_ld", "snp_not_ld")
ratio_opts <- c(0.5, 0.8)
sv_effect_opts <- c(0.1, 0.5)
trait_pops <- "1,2,3"
pred_iters <- "1,2,3"

args <- commandArgs(trailingOnly = TRUE)
if ("--help" %in% args) usage() & q(save = "no")

# assert to have the correct optional arguments
if (length(args) < 3) stop(usage(), "missing positional argument(s)")

if (length(args) > 3) {

  opt_args <- args[-1:-3]
  opt_args_allowed <- c("--trait-pops", "--pred-iters")
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

# adjust format from optional arguments
trait_pops <- unlist(strsplit(trait_pops, split = ","))
pred_iters <- unlist(strsplit(pred_iters, split = ","))

if (length(trait_pops) != length(pred_iters)) stop("Number of QTN pops and marker iterations doesn't match")

# get positional arguments
folder_results <- args[1]
pred_accuracy_file <- args[2]
outfile <- args[3]



#### summarize LD between causative variants and predictors ----

# load prediction accuracy results
pred_accuracy <- fread(pred_accuracy_file, header = TRUE, data.table = FALSE)

# create empty df to store results
ld_summary_all_scenarios <- data.frame(stringsAsFactors = FALSE)
results_accuracy_ld_pve <- data.frame(stringsAsFactors = FALSE)

# get results for traits controlled by SNPs only or SVs only
trait_ratio <- NA
trait_sv_effect <- NA
trait_diff_dist <- NA
for (trait_var in c("SNP", "SV")) {
  print(trait_var)
  for (trait_qtn in qtn_opts) {
    for (pred_type in predictor_opts) {
      for (i in 1:length(trait_pops)) {

        # get qtn pops and pred iters
        trait_pop <- trait_pops[i]
        pred_n <- pred_iters[i]

        # load ld summary file
        ld_summary_file <- paste0(folder_results, "/", trait_qtn, "-QTNs_from_",
                                  trait_var, "/ld_causative-vars_predictors/pop",
                                  trait_pop, "/ld_summary.", pred_type,
                                  "_markers.pred-iter", pred_n, ".causal-pop",
                                  trait_pop, ".txt")
        ld_summary <- fread(ld_summary_file, header = TRUE, data.table = FALSE)

        # combine all ld summaries into one data frame
        ld_summary_all_scenarios <- rbind(ld_summary_all_scenarios,
                                          cbind(data.frame(qtn = trait_qtn, var = trait_var,
                                                           ratio = trait_ratio, pop = trait_pop,
                                                           pred_type = pred_type, pred_iter = pred_n),
                                                ld_summary))

        for (trait_h2 in h2_opts) {

          # get accuracy for each particular simulated and predicted scenario
          pred_accuracy_scenario <- pred_accuracy %>%
            filter(h2 == trait_h2, qtn == trait_qtn, var == trait_var, is.na(ratio),
                   is.na(sv_effect), is.na(diff_dist), pop == trait_pop,
                   predictor == pred_type, pred_iter == pred_n) %>%
            pivot_longer(!(h2:mean_upperCI), names_to = "stat", values_to = "vals") %>%
            select(!(mean_accuracy_envs:mean_upperCI)) %>%
            separate(stat, into = c("env", "stat"), sep = "_") %>%
            pivot_wider(names_from = "stat", values_from = "vals") %>%
            as.data.frame()

          # load pve per qtn summary file
          pve_qtn_summary_file <- paste0(folder_results, "/", trait_qtn, "-QTNs_from_",
                                         trait_var, "/", trait_h2, "-heritability/pop",
                                         trait_pop, "/summary_var_explained_per_qtn.txt")
          pve_qtn_summary <- fread(pve_qtn_summary_file, header = TRUE, data.table = FALSE)

          # get predictors in highest ld to causative variants
          highest_ld_results <- getMarkersHighestLD(ld_summary)

          # change column name to facilitate merging
          colnames(pve_qtn_summary)[2] <- "causal_var"
          # merge ld results with qtn effects
          qtn_pred_ld_pve <- merge(x = pve_qtn_summary, y = highest_ld_results, by = "causal_var", all = TRUE)
          # sort columns
          qtn_pred_ld_pve <- qtn_pred_ld_pve[order(qtn_pred_ld_pve$qtn_number, qtn_pred_ld_pve$env), ]
          # merge pred accuracy with ld results and qtn effects
          accuracy_ld_pve_scenario <- merge(pred_accuracy_scenario, qtn_pred_ld_pve, by = "env", all = TRUE)
          # append to final df
          results_accuracy_ld_pve <- rbind(results_accuracy_ld_pve, accuracy_ld_pve_scenario)

        }
      }
    }
  }
}

# get results for traits controlled by both SNPs and SVs
trait_var <- "both"
print(trait_var)
for (trait_qtn in qtn_opts) {
  for (pred_type in predictor_opts) {
    for (trait_ratio in ratio_opts) {
      for (i in 1:length(trait_pops)) {

        # get qtn pops and pred iters
        trait_pop <- trait_pops[i]
        pred_n <- pred_iters[i]

        # load ld filename
        ld_summary_file <- paste0(folder_results, "/", trait_qtn, "-QTNs_from_",
                                  trait_var, "/SNP-SV-ratio_", trait_ratio,
                                  "/ld_causative-vars_predictors/pop", trait_pop,
                                  "/ld_summary.", pred_type, "_markers.pred-iter", pred_n,
                                  ".causal-pop", trait_pop, ".txt")
        ld_summary <- fread(ld_summary_file, header = TRUE, data.table = FALSE)

        # combine all ld summaries into one data frame
        ld_summary_all_scenarios <- rbind(ld_summary_all_scenarios,
                                          cbind(data.frame(qtn = trait_qtn, var = trait_var,
                                                           ratio = trait_ratio, pop = trait_pop,
                                                           pred_type = pred_type, pred_iter = pred_n),
                                                ld_summary))

        for (trait_h2 in h2_opts) {
          for (trait_sv_effect in sv_effect_opts) {
            for (trait_diff_dist in c(FALSE, TRUE)) {

              # get accuracy for each particular simulated and predicted scenario
              pred_accuracy_scenario <- pred_accuracy %>%
                filter(h2 == trait_h2, qtn == trait_qtn, var == trait_var, ratio == trait_ratio,
                       sv_effect == trait_sv_effect, diff_dist == trait_diff_dist,
                       pop == trait_pop, predictor == pred_type, pred_iter == pred_n) %>%
                pivot_longer(!(h2:mean_upperCI), names_to = "stat", values_to = "vals") %>%
                select(!(mean_accuracy_envs:mean_upperCI)) %>%
                separate(stat, into = c("env", "stat"), sep = "_") %>%
                pivot_wider(names_from = "stat", values_from = "vals") %>%
                as.data.frame()

              # load pve per qtn summary file
              pve_qtn_summary_file <- paste0(folder_results, "/", trait_qtn, "-QTNs_from_",
                                             trait_var, "/SNP-SV-ratio_", trait_ratio,
                                             "/effects_SNP-0.1_SV-",  trait_sv_effect, "/",
                                             trait_h2, "-heritability/pop", trait_pop,
                                             "/summary_var_explained_per_qtn.txt")
              pve_qtn_summary <- fread(pve_qtn_summary_file, header = TRUE, data.table = FALSE)

              # get predictors in highest ld to causative variants
              highest_ld_results <- getMarkersHighestLD(ld_summary)

              # change column name to facilitate merging
              colnames(pve_qtn_summary)[2] <- "causal_var"
              # merge ld results with qtn effects
              qtn_pred_ld_pve <- merge(x = pve_qtn_summary, y = highest_ld_results, by = "causal_var", all = TRUE)
              # sort columns
              qtn_pred_ld_pve <- qtn_pred_ld_pve[order(qtn_pred_ld_pve$qtn_number, qtn_pred_ld_pve$env), ]
              # merge pred accuracy with ld results and qtn effects
              accuracy_ld_pve_scenario <- merge(pred_accuracy_scenario, qtn_pred_ld_pve, by = "env", all = TRUE)
              # append to final df
              results_accuracy_ld_pve <- rbind(results_accuracy_ld_pve, accuracy_ld_pve_scenario)

            }
          }
        }
      }
    }
  }
}

# sort columns
ld_summary_all_scenarios <- ld_summary_all_scenarios %>%
  arrange(qtn, var, ratio, pop, pred_type, pred_iter, chr)

results_accuracy_ld_pve <- results_accuracy_ld_pve %>%
  arrange(env, h2, qtn, desc(var), ratio, sv_effect, diff_dist, pop, predictor, pred_iter, cv, qtn_number)

# write final results
outfile_ld_summary <- gsub(".txt", ".ld-summary.txt", outfile)
fwrite(ld_summary_all_scenarios, outfile_ld_summary, quote = FALSE, sep = "\t", na = NA, row.names = FALSE)

outfile_cor <- gsub(".txt", ".results.txt", outfile)
fwrite(results_accuracy_ld_pve, outfile_cor, quote = FALSE, sep = "\t", na = NA, row.names = FALSE)



#### debug ----

# folder_results <- "analysis/trait_sim/multi_env/with_gxe/additive_model/equal_effects"
# pred_accuracy_file <- "analysis/trait_sim/multi_env/prediction_results.full.txt"
# outfile <- "analysis/trait_sim/multi_env/pred-accuracy_qtn-effect_ld-pred-qtn.txt"
# trait_pops <- c(15, 4, 18) #c(15, 4, 18, 16, 14, 20, 2, 8, 17, 11, 6, 1, 13, 7, 5, 10, 3, 9, 12, 19)
# pred_iters <- c(5, 19, 10) #c(5, 19, 10, 7, 16, 4, 17, 9, 1, 3, 18, 11, 15, 6, 2, 12, 20, 13, 8, 14)
