library(data.table)
library(ggplot2)
library(dplyr)

usage <- function() {
  cat("
description: plot prediction accuracy results after k-fold cross validation of simulated traits.

usage: Rscript select_random_markers.R [results_folder] [outfile_name] [...]

positional arguments:
  results_folder                path to folder with results of cross validation
  outfile_name                  plot name

optional argument:
  --help                        show this helpful message
  --trait-qtn=[LIST]            comma-separated list of QTNs used for simulating traits
  --trait-var-source=[LIST]     comma-separated list of causative variants used for simulating traits
  --trait-h2=[LIST]             comma-separated list of heritabilities used for simulating traits
  --pred-marker-n=[LIST]        comma-separated list of number of markers used in genomic prediction
  --pred-marker-type=[LIST]     comma-separated list of type of markers used in genomic prediction
  --error-bars=[VALUE]          plot error bars as standard error (SE; default) or confidence intervals (CI)

"
  )
}

getArgValue <- function(arg) {

  # get value from command-line arguments
  arg <- unlist(strsplit(arg, "="))
  if (length(arg) == 1) return(TRUE)
  if (length(arg) > 1) return(arg[2])

}

getPredictionAccuracies <- function(trait_qtn_number = c("3", "25", "75"),
                                    trait_var_source = c("SNP", "SV", "both"),
                                    trait_heritability = c("0.2", "0.5" ,"0.9"),
                                    prediction_marker_number = c("max-number", "1000", "50"),
                                    prediction_marker_type = c("all", "sv", "snp", "snp_ld", "snp_not_ld"),
                                    dir_with_results = "analysis/trait_sim/additive_model/geom_series_effects") {

  # create empty df to store results
  df_results <- data.frame(stringsAsFactors = FALSE)

  # get first folder names within the folder provided
  initial_dirs <- list.dirs(dir_with_results, recursive = FALSE)

  for (qtn in trait_qtn_number) {

    for (var in trait_var_source) {

      dir_qtn_var <- initial_dirs[grepl(qtn, initial_dirs) & grepl(var, initial_dirs)]
      dir_qtn_var <- list.dirs(dir_qtn_var, recursive = FALSE)

      for (marker_n in prediction_marker_number) {

        dir_pred_n <- dir_qtn_var[grepl(marker_n, dir_qtn_var)]
        dir_pred_n <- list.dirs(dir_pred_n, recursive = FALSE)

        for (marker_type in prediction_marker_type) {

          dir_pred_type <- dir_pred_n[grepl(marker_type, dir_pred_n)]
          dir_pred_type <- list.dirs(dir_pred_type, recursive = FALSE)

          for (h2 in trait_heritability) {

            dir_h2 <- dir_pred_type[grepl(h2, dir_pred_type)]
            dir_h2 <- list.dirs(dir_h2, recursive = TRUE)

            # keep only dirs with reps
            dir_reps <- dir_h2[grepl("rep", dir_h2)]

            # get file names of all reps
            reps <- list.files(dir_reps, pattern = "CV_Results", full.names = TRUE)
            reps <- reps[grep(".txt", reps)]

            for (file in reps) {

              # get rep number
              rep <- unlist(strsplit(file, split = "/"))
              rep <- rep[grep("rep", rep)]
              rep <- as.numeric(gsub("rep", "", rep))

              # extract CV result from rep
              cv_result <- fread(file, header = TRUE, data.table = FALSE)

              # add results to table along with metadata
              df_results <- rbind(df_results, data.frame(trait_qtn = qtn,
                                                         trait_var = var,
                                                         trait_h2 = h2,
                                                         pred_marker_n = marker_n,
                                                         pred_marker_type = marker_type,
                                                         rep = rep,
                                                         cv_result))

            }
          }
        }
      }
    }
  }

  # order final df
  df_results <- df_results[order(df_results$trait_qtn,
                                 df_results$trait_var,
                                 df_results$trait_h2,
                                 df_results$pred_marker_n,
                                 df_results$pred_marker_type,
                                 df_results$rep), ]

  return(df_results)

}


#### command line options ----

# set default
trait_qtn <- "3,25,75"
trait_var_source <- "SNP,SV,both"
trait_h2 <- "0.2,0.5,0.9"
pred_marker_n <- "max-number,1000,50"
pred_marker_type <- "all_markers,sv_markers,snp_markers,snp_ld_markers,snp_not_ld_markers"
error_bars <- "SE"

args <- commandArgs(trailingOnly = TRUE)
if ("--help" %in% args) usage() & q(save = "no")

# assert to have the correct optional arguments
if (length(args) < 2) stop(usage(), "missing positional argument(s)")

if (length(args) > 2) {

  opt_args <- args[-1:-2]
  opt_args_allowed <- c("--trait-qtn", "--trait-var-source", "--trait-h2", "--pred-marker-n",
                        "--pred-marker-type", "--error-bars")
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

if (!error_bars %in% c("SE", "CI")) {
  stop("Optional argument '--error-bars' should be either 'SE', or 'CI'")
}


# transform arguments into lists
trait_qtn <- unlist(strsplit(trait_qtn, split = ","))
trait_var_source <- unlist(strsplit(trait_var_source, split = ","))
trait_h2 <- unlist(strsplit(trait_h2, split = ","))
pred_marker_n <- unlist(strsplit(pred_marker_n, split = ","))
pred_marker_type <- unlist(strsplit(pred_marker_type, split = ","))

# get positional arguments
results_folder <- args[1]
outfile_name <- args[2]
# results_folder <- "analysis/trait_sim/additive_model/geom_series_effects"
# outfile_name <- "analysis/trait_sim/additive_model/geom_series_effects/prediction_accuracy_sim-traits_single-env.txt"
# results_folder <- "analysis/trait_sim_mult-env/additive_model/geom_series_effects"
# outfile_name <- "analysis/trait_sim_mult-env/additive_model/geom_series_effects/prediction_accuracy_sim-traits_mult-env.txt"


#### summarize prediction results ----

# get prediction accuracies according to the genetic architecture
prediction_results <- getPredictionAccuracies(trait_qtn_number = trait_qtn,
                                              trait_var_source = trait_var_source,
                                              trait_heritability = trait_h2,
                                              prediction_marker_number = pred_marker_n,
                                              prediction_marker_type = pred_marker_type,
                                              dir_with_results = results_folder)

# get mean prediction accuracy across folds and add as another column
# same for SE and CI bounds (or CI bounds should be calculated based on all reps?)

prediction_summary <- prediction_results %>%
  group_by(trait_qtn, trait_var, trait_h2, pred_marker_n, pred_marker_type) %>%
  summarise(reps = n(),
            mean_cor = mean(r_avg),
            mean_SE = mean(r_SE),
            mean_lowerCI = mean(r_lowerCI),
            mean_upperCI = mean(r_upperCI))

# write results
fwrite(prediction_results, file = outfile_name, quote = FALSE, sep = "\t", row.names = FALSE, na = NA)
fwrite(prediction_summary, file = gsub(".txt", ".summary.txt", outfile_name), quote = FALSE, sep = "\t", row.names = FALSE, na = NA)



#### plot prediction results ----

# reorder levels of marker type for better visualization
prediction_summary$pred_marker_type <- factor(prediction_summary$pred_marker_type,
                                              levels = c("snp_markers", "snp_not_ld_markers", "snp_ld_markers",
                                                         "sv_markers","all_markers"),
                                              labels = c("Only SNPs", "SNPs not in LD", "SNPs in LD",
                                                         "Only SVs","SNPs and SVs"))


# create plots for each scenario
for (qtn in trait_qtn) {
  for (h2 in trait_h2) {
    for (n_markers in pred_marker_n) {
      
      # filter summary
      prediction_summary_plot <- prediction_summary[prediction_summary$trait_qtn == qtn &
                                                      prediction_summary$trait_h2 == h2 & 
                                                      prediction_summary$pred_marker_n == n_markers, ]
    
      # plot
      results_plot <- ggplot(data = prediction_summary_plot, aes(x = pred_marker_type, y = mean_cor, fill = pred_marker_type)) +
        geom_bar(stat = "identity", show.legend = FALSE) +
        annotate("text", label = "(causative variant)", x = 3, y = 1.0, vjust = -0.5) +
        facet_grid(cols = vars(trait_var)) +
        coord_cartesian(ylim = c(0, 1)) +
        scale_x_discrete(drop = FALSE) +
        labs(title = paste0(qtn, " QTNs, ", h2, " heritability"),
             subtitle = paste0("(", gsub("-", " ", n_markers), " markers in prediction)"),
             x = "Marker type used in prediciton",
             y = "Prediction accuracy") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              axis.title = element_text(size = 15),
              title = element_text(size = 15),
              axis.text.y = element_text(size = 15),
              text = element_text(size = 10), strip.text.x = element_text(size = 15)) +
        scale_colour_viridis_d(option = "B", begin = 0.3, end = 0.8, aesthetics = "fill", drop = FALSE)
        # scale_fill_manual(values = colorRampPalette(c("#900721", "#FFD75F"))(5))
        # scale_fill_brewer(type = "div", palette = "RdYlBu")
      
      # add error bars
      if (error_bars == "SE") {
        results_plot <- results_plot +
          geom_errorbar(aes(ymin = mean_cor - mean_SE, ymax = mean_cor + mean_SE), width = 0.2)
      }
      
      if (error_bars == "CI") {
        results_plot <- results_plot +
          geom_errorbar(aes(ymin = mean_lowerCI, ymax = mean_upperCI), width = 0.2)
      }
      
      # save plot
      plot_name <- gsub(".txt", paste0(".qtns_", qtn, ".h2_", h2, ".nmarkers_", n_markers, ".pdf"), outfile_name)
      ggsave(results_plot, filename = plot_name, device = "pdf")
      
    }
  }
}
