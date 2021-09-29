library(data.table)
library(dplyr)
library(tidyr)

usage <- function() {
  cat("
description: summarize prediction accuracy results after k-fold cross validation of simulated traits.

usage: Rscript summarize_prediction_accuracies.R [folder_base] [outfile_name] [...]

positional arguments:
  folder_base                   path to folder with results of cross validation
  outfile_name                  output filename

optional argument:
  --help                        show this helpful message
  --trait-pops=[LIST]           comma-separated list of populations of QTNs used for simulating traits
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

getAccuracy <- function(folder_accuracy, h2, qtn, var, ratio, sv_effect,
                        diff_dist, pop, predictor, pred_iter) {
  
  if (is.null(ratio)) ratio <- NA
  if (is.null(sv_effect)) sv_effect <- NA
  if (is.null(diff_dist)) diff_dist <- NA
  
  # create empty df to store results of CVs
  df_results <- data.frame(stringsAsFactors = FALSE)
  
  for (cv in c("CV1", "CV2")) {
    
    # try to access file with results
    pred_accuracy <- try(fread(paste0(folder_accuracy, "/prediction_accuracy.", cv, ".txt"),
                               header = TRUE, data.table = FALSE))
    
    if (class(pred_accuracy) != "try-error") {
      
      # get avg accuracy across all environments
      pred_accuracy_means <- data.frame(t(rowMeans(pred_accuracy[, -1])))
      colnames(pred_accuracy_means) <- c("mean_accuracy_envs", "mean_se", "mean_lowerCI", "mean_upperCI")
      
      # keep accuracy about each environment as well -- just need to reformat
      pred_accuracy_envs <- pred_accuracy %>%
        pivot_longer(-stat, names_to = "env", values_to = "vals") %>%
        arrange(env) %>% 
        unite(env:stat, col = "env_info", sep = "_") %>%
        t() %>% as.data.frame(stringsAsFactors = FALSE)
      colnames(pred_accuracy_envs) <- apply(pred_accuracy_envs[1, ], MARGIN = 1, function(x) as.character(gsub("_CI", "CI", x)))
      pred_accuracy_envs <- data.frame(apply(pred_accuracy_envs[-1, ], MARGIN = c(1, 2), function(x) as.numeric(x)))
      rownames(pred_accuracy_envs) <- NULL
      
      # add metadata
      df_results <- rbind(df_results,
                          data.frame(h2 = h2, qtn = qtn, var = var,
                                     ratio = ratio, sv_effect = sv_effect,
                                     diff_dist = diff_dist, pop = pop,
                                     predictor = predictor, pred_iter = pred_iter,
                                     cv = cv, pred_accuracy_means, pred_accuracy_envs))
      
    }
  }
  
  # only return values if both CV1 and CV2 results exist
  if (NROW(df_results) == 2) return(df_results)
  
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
if (length(args) < 2) stop(usage(), "missing positional argument(s)")

if (length(args) > 2) {

  opt_args <- args[-1:-2]
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
folder_base <- args[1]
outfile_name <- args[2]



#### summarize prediction results ----

# create empty df to store results
prediction_results <- data.frame(stringsAsFactors = FALSE)

# get results for traits controlled by SNPs only or SVs only
ratio <- NULL
sv_effect <- NULL
diff_dist <- NULL
for (var in c("SNP", "SV")) {
  for (qtn in qtn_opts) {
    for (h2 in h2_opts) {
      for (predictor in predictor_opts) {
        for (i in 1:length(trait_pops)) {
          
          pop <- trait_pops[i]
          pred_iter <- pred_iters[i]
          
          # get folder with prediction accuracy results
          folder_accuracy <- paste0(folder_base, "/with_gxe/additive_model/equal_effects/",
                                    qtn, "-QTNs_from_", var, "/", h2, "-heritability/pop", pop,
                                    "/prediction_iter", pred_iter, "/", predictor, "_markers")
          
          # get accuracy results
          prediction_results <- rbind(prediction_results,
                                      getAccuracy(folder_accuracy, h2 = h2, qtn = qtn,
                                                  var = var, ratio = ratio, sv_effect = sv_effect,
                                                  diff_dist = diff_dist, pop = pop,
                                                  predictor = predictor, pred_iter = pred_iter))
          
        }
      }
    }
  }
}

# get results for traits controlled by both SNPs and SVs
var <- "both"
for (qtn in qtn_opts) {
  for (h2 in h2_opts) {
    for (predictor in predictor_opts) {
      for (ratio in ratio_opts) {
        for (sv_effect in sv_effect_opts) {
          for (diff_dist in c(FALSE, TRUE)) {
            for (i in 1:length(trait_pops)) {
              
              pop <- trait_pops[i]
              pred_iter <- pred_iters[i]
              
              # get folder with prediction accuracy results
              folder_accuracy <- paste0(folder_base, "/with_gxe/additive_model/equal_effects/", qtn,
                                        "-QTNs_from_", var, "/SNP-SV-ratio_", ratio, "/effects_SNP-0.1_SV-", 
                                        sv_effect, "/", h2, "-heritability/pop", pop,
                                        "/prediction_iter", pred_iter, "/", predictor, "_markers")
              
              # get accuracy results
              prediction_results <- rbind(prediction_results,
                                          getAccuracy(folder_accuracy, h2 = h2, qtn = qtn,
                                                      var = var, ratio = ratio, sv_effect = sv_effect,
                                                      diff_dist = diff_dist, pop = pop, predictor = predictor,
                                                      pred_iter = pred_iter))
              
            }
          }
        }
      }
    }
  }
}

# summarize results
prediction_summary <- prediction_results %>%
  group_by(h2, qtn, var, ratio, sv_effect, diff_dist, predictor, cv) %>%
  summarize(reps = n(), across(mean_accuracy_envs:env5_upperCI, ~ mean(.x, na.rm = TRUE)))

# round numbers
prediction_results[, 11:NCOL(prediction_results)] <- apply(prediction_results[, 11:NCOL(prediction_results)],
                                                           MARGIN = 2, function(x) round(x, digits = 4))
prediction_summary[, 10:NCOL(prediction_summary)] <- apply(prediction_summary[, 10:NCOL(prediction_summary)],
                                                           MARGIN = 2, function(x) round(x, digits = 4))

# write results
fwrite(prediction_results, file = gsub(".txt", ".full.txt", outfile_name),
       quote = FALSE, sep = "\t", row.names = FALSE, na = NA)
fwrite(prediction_summary, file = gsub(".txt", ".summary.txt", outfile_name),
       quote = FALSE, sep = "\t", row.names = FALSE, na = NA)




#### debug ----

# folder_base <- "analysis/trait_sim/multi_env"
# outfile_name <- "analysis/trait_sim/multi_env/test_prediction_results.pops15-4-18.txt"
# trait_pops <- c(15, 4, 18) #, 16, 14, 20, 2, 8, 17, 11, 6, 1, 13, 7, 5, 10, 3, 9, 12, 19)
# pred_iters <- c(5, 19, 10) #, 7, 16, 4, 17, 9, 1, 3, 18, 11, 15, 6, 2, 12, 20, 13, 8, 14)