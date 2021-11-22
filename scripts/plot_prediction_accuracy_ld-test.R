library(data.table)
library(ggplot2)
library(dplyr)
library(RColorBrewer)

usage <- function() {
  cat("
description: plot prediction accuracy results after k-fold cross validation of simulated traits with different LD levels.

usage: Rscript plot_prediction_accuracy_ld-test.R [prediction_summary_file] [...]

positional arguments:
  prediction_summary_file       path to file with prediction summary

optional argument:
  --help                        show this helpful message
  --error-bars=[VALUE]          error bars can take one of the following values:
                                'CI_accuracy' (confidence intervals of the mean prediction accuracy; default),
                                'CI_error' (confidence intervals of the mean prediction error),
                                'SE_accuracy' (standard error of the mean prediction accuracy), or
                                'SE_error' (standard error of the mean prediction error)

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
error_bars <- "CI_accuracy"

args <- commandArgs(trailingOnly = TRUE)
if ("--help" %in% args) usage() & q(save = "no")

# assert to have the correct optional arguments
if (length(args) < 1) stop(usage(), "missing positional argument(s)")

if (length(args) > 1) {

  opt_args <- args[-1:-2]
  opt_args_allowed <- c("--error-bars")
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
prediction_summary_file <- args[1]



#### load prediction summary ----

# do causal variatnt 'both' separately
prediction_summary <- fread(prediction_summary_file, header = TRUE, data.table = FALSE)

# reorder levels of marker type for better visualization
prediction_summary$predictor <- factor(prediction_summary$predictor,
                                       levels = c("low", "moderate", "high"),
                                       labels = c("Low", "Moderate", "High"))

# create output folder
outfolder <- unlist(strsplit(prediction_summary_file, "/"))
outfolder[length(outfolder)] <- paste0("plots_", gsub(".txt", "", outfolder[length(outfolder)]))
outfolder[length(outfolder)] <- gsub(".", "_", outfolder[length(outfolder)], fixed = TRUE)
outfolder <- paste0(outfolder, collapse = "/")
if (!dir.exists(outfolder)) dir.create(outfolder)



#### plot SNP only or SV only as causative variants ----

# set default
trait_var_source <- "SNP"
trait_qtn_number <- 100
# trait_var_source <- c("SNP", "SV")
# trait_qtn_number <- c(10, 100)
trait_heritability <- c(0.3 ,0.7)

# create plots for each scenario
for (qtn_number in trait_qtn_number) {
  for (heritability in trait_heritability) {
    
    # filter summary
    prediction_summary_plot <- subset(prediction_summary, qtn == qtn_number
                                      & h2 == heritability
                                      & var %in% trait_var_source)
    
    # plot
    results_plot <- ggplot(data = prediction_summary_plot,
                           aes(x = predictor, y = mean_accuracy_envs, fill = cv)) +
      geom_bar(stat = "identity", position = "dodge") +
      # annotate("text", label = "(causative variant)", x = 3, y = 1.0, vjust = -0.5, size = 7) +
      facet_grid(cols = vars(var)) +
      coord_cartesian(ylim = c(0, 1)) +
      scale_x_discrete(drop = FALSE) +
      labs(title = paste0(qtn_number, " QTLs, ", heritability, " heritability"),
           x = "LD level of marker to QTL",
           y = "Prediction accuracy") +
      scale_fill_manual(values = c("#DE8282FF", "#AD0000FF")) +
      guides(fill = guide_legend("Cross-Validation")) +
      theme_bw() +
      theme(title = element_text(size = 20),
            text = element_text(size = 20),
            axis.title = element_text(size = 20),
            axis.title.x = element_text(vjust = -1.5),
            axis.title.y = element_text(vjust = 1.5),
            axis.text.y = element_text(size = 20),
            strip.text.x = element_text(size = 20),
            legend.position = "bottom",
            legend.title = element_text(size = 15),
            legend.text = element_text(size = 12)) 
    
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
    
    # print(results_plot)
    # Sys.sleep(3)
    
    # save plot
    plot_name <- paste0(outfolder, "/pred_accuracy.qtns_", qtn_number, ".h2_", heritability, ".pdf")
    ggsave(results_plot, filename = plot_name, device = "pdf", width = 15, height = 10)
    
  }
}

# View(subset(prediction_summary, gxe == "with"
#             & qtn  %in% c(10, 100)
#             & h2 %in% c(0.3, 0.7)
#             & var == "SV"
#             & cv == "CV1"),
#      "summary_cv1")
# 
# View(subset(prediction_summary, gxe == "with"
#             & qtn  %in% c(10, 100)
#             & h2 %in% c(0.3, 0.7)
#             & var == "SV"
#             & cv == "CV2"),
#      "summary_cv2")
# 
# View(subset(prediction_summary, gxe == "with"
#             & qtn  %in% c(10, 100)
#             & h2 %in% c(0.3, 0.7)
#             & var == "SV"),
#      "summary_cvs")
# 
# View(subset(prediction_summary, gxe == "with"
#             & qtn  %in% c(10, 100)
#             & h2 %in% c(0.3, 0.7)
#             & var == "SV"),
#      "summary_cvs")



#### plot both SNPs and SVs as causative variants ----

# set default
trait_var_source <- "both"
trait_qtn_number <- c(10, 100)
trait_heritability <- c(0.3 ,0.7)
snp_sv_ratio <- c(0.5, 0.8)
marker_effect <- 0.1
sv_effects <- c(0.1, 0.5)


sv_diff_dist <- FALSE
for (qtn_number in trait_qtn_number) {
  for (heritability in trait_heritability) {
    for (var_ratio in snp_sv_ratio) {
      
      # filter summary
      prediction_summary_plot <- subset(prediction_summary, qtn == qtn_number
                                        & h2 == heritability
                                        & ratio == var_ratio
                                        & sv_effect %in% sv_effects
                                        & diff_dist == sv_diff_dist)
      
      # plot
      results_plot <- ggplot(data = prediction_summary_plot,
                             aes(x = predictor, y = mean_accuracy_envs, fill = cv)) +
        geom_bar(stat = "identity", position = "dodge") +
        annotate("text", label = "(effect size)", x = 3, y = 1.0, vjust = -0.5, size = 7) +
        facet_grid(cols = vars(sv_effect),
                   labeller = labeller(sv_effect = c("0.1" = "SV = SNP", "0.5" = "SV = 5x SNP"))) +
        coord_cartesian(ylim = c(0, 1)) +
        scale_x_discrete(drop = FALSE) +
        labs(title = paste0(qtn_number, " QTNs, ", heritability, " heritability, ", var_ratio, " SNP/SV ratio"),
             # subtitle = "(SNPs and SVs as causative variants)",
             x = "Marker type used in prediciton",
             y = "Prediction accuracy") +
        scale_fill_manual(values = c("#DE8282FF", "#AD0000FF")) +
        guides(fill = guide_legend("Cross-Validation")) +
        theme_bw() +
        theme(title = element_text(size = 20),
              text = element_text(size = 20),
              axis.title = element_text(size = 20),
              axis.title.x = element_text(vjust = -1.5),
              axis.title.y = element_text(vjust = 1.5),
              axis.text.y = element_text(size = 20),
              strip.text.x = element_text(size = 20),
              legend.position = "bottom",
              legend.title = element_text(size = 15),
              legend.text = element_text(size = 12)) 
      
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
      
      # print(results_plot)
      # Sys.sleep(3)
      
      # save plot
      plot_name <- paste0(outfolder, "/pred_accuracy.qtns_", qtn_number, ".h2_", heritability,
                          ".both-causal.", var_ratio, "snp-", (1 - var_ratio), "sv.pdf")
      ggsave(results_plot, filename = plot_name, device = "pdf", width = 15, height = 10)
      
      
    }
  }
}

# View(subset(prediction_summary, gxe == "with"
#             & qtn %in% c(100)
#             & h2 %in% c(0.3, 0.7)
#             & ratio %in% c(0.5, 0.8)
#             & sv_effect %in% c(0.1)
#             & diff_dist == FALSE
#             & var == "both"),
#      "summary_cvs_qtns100_ratio0.5")
# 
# View(subset(prediction_summary, gxe == "with"
#             & qtn %in% c(100)
#             & h2 %in% c(0.7)
#             & ratio %in% c(0.5, 0.8)
#             & sv_effect %in% c(0.1, 0.5)
#             & diff_dist == FALSE
#             & var == "both"),
#      "summary_cvs_qtns100_h2-0.7_ratio0.5")
# 
# View(subset(prediction_summary, gxe == "with"
#             & qtn %in% c(100)
#             & h2 %in% c(0.3, 0.7)
#             & ratio %in% c(0.8)
#             & sv_effect %in% c(0.1, 0.5)
#             & diff_dist == FALSE
#             & cv == "CV1"
#             & var == "both"),
#      "summary_cv1_qtns100_ratio0.8")
# 
# View(subset(prediction_summary, gxe == "with"
#             & qtn %in% c(10, 100)
#             & h2 %in% c(0.3, 0.7)
#             & ratio %in% c(0.5, 0.8)
#             & sv_effect %in% c(0.1, 0.5)
#             & diff_dist == FALSE
#             & var == "both"
#             & cv == "CV2"),
#      "summary_cv2")
# 
# View(subset(prediction_summary, gxe == "with"
#             & qtn %in% c(10, 100)
#             & h2 %in% c(0.3, 0.7)
#             & ratio %in% c(0.5, 0.8)
#             & sv_effect %in% c(0.1, 0.5)
#             & diff_dist == FALSE
#             & var == "both"
#             & cv == "CV1"),
#      "summary_cv1")
# 
# View(subset(prediction_summary, gxe == "with"
#             & qtn %in% c(10, 100)
#             & h2 %in% c(0.3, 0.7)
#             & ratio %in% c(0.5, 0.8)
#             & sv_effect %in% c(0.1, 0.5)
#             & diff_dist == FALSE
#             & var == "both"),
#      "summary_cvs")



sv_effects <- 0.5
sv_diff_dist <- c(FALSE, TRUE)
for (qtn_number in trait_qtn_number) {
  for (heritability in trait_heritability) {
    for (var_ratio in snp_sv_ratio) {
      
      # filter summary
      prediction_summary_plot <- subset(prediction_summary, qtn == qtn_number
                                        & h2 == heritability
                                        & ratio == var_ratio
                                        & sv_effect %in% sv_effects
                                        & diff_dist %in% sv_diff_dist)
      
      # plot
      results_plot <- ggplot(data = prediction_summary_plot,
                             aes(x = predictor, y = mean_accuracy_envs, fill = cv)) +
        geom_bar(stat = "identity", position = "dodge") +
        annotate("text", label = "GxE effects", x = 3, y = 1.0, vjust = -0.5, size = 7) +
        facet_grid(cols = vars(diff_dist),
                   labeller = labeller(diff_dist = c("FALSE" = "SV = SNP", "TRUE" = "SV = 2xSNP"))) +
        coord_cartesian(ylim = c(0, 1)) +
        scale_x_discrete(drop = FALSE) +
        labs(title = paste0(qtn_number, " QTNs, ", heritability, " heritability, ", var_ratio, " SNP/SV ratio, SV effect = 5x SNP effect"),
             # subtitle = "(SNPs and SVs as causative variants)",
             x = "Marker type used in prediciton",
             y = "Prediction accuracy") +
        scale_fill_manual(values = c("#DE8282FF", "#AD0000FF")) +
        guides(fill = guide_legend("Cross-Validation")) +
        theme_bw() +
        theme(title = element_text(size = 20),
              text = element_text(size = 20),
              axis.title = element_text(size = 20),
              axis.title.x = element_text(vjust = -1.5),
              axis.title.y = element_text(vjust = 1.5),
              axis.text.y = element_text(size = 20),
              strip.text.x = element_text(size = 20),
              legend.position = "bottom",
              legend.title = element_text(size = 15),
              legend.text = element_text(size = 12)) 
      
      
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
      
      # print(results_plot)
      # Sys.sleep(3)
      
      # save plot
      plot_name <- paste0(outfolder, "/pred_accuracy.qtns_", qtn_number, ".h2_", heritability,
                          ".both-causal.diff-gxe-effects.", var_ratio, "snp-", (1 - var_ratio), "sv.pdf")
      ggsave(results_plot, filename = plot_name, device = "pdf", width = 15, height = 10)
      
    }
  }
}



#### choosing colors ----

# library(paletteer)
# library(tinter)
# paletteer_c(`"viridis::inferno"`, n = 20)
# for (col in c("#360961FF", "#8C2369FF", "#D84D3EFF", "#F1731DFF", "#FCA108FF")) {
#   print(prismatic::color(rev(tinter(col, direction = "tints"))[c(1,3)]))
# }



#### debug ----

# prediction_summary_file <- "analysis/ld_downsample/test_prediction_results.reps1-10.pops1-3.summary.txt"
