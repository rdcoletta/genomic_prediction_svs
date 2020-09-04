library(data.table)
library(ggplot2)

usage <- function() {
  cat("
description: plot LD decay from PLINK files

usage: Rscript plot_ld_decay.R [ld_filename] [output_prefix] [...]

positional arguments:
  ld_filename             input ld file from PLINK
  output_prefix           name of ld decay plot

optional argument:
  --help                  show this helpful message
  --average-ld            plot average ld from 100bp windows (default) instead of scatter plot with all points
  --window-size=VALUE     control window size in bp (default: 100)
  --window-step=VALUE     control window step in bp (default: 100)
  --decay-random-sv       plot ld decay from three randomly chosen SVs

"
  )
}


getArgValue <- function(arg) {
  
  # get value from command-line arguments
  arg <- unlist(strsplit(arg, "="))
  if (length(arg) == 1) return(TRUE)
  if (length(arg) > 1) return(as.integer(arg[2]))
  
}

plot_ld_decay <- function(data, x, y, point_size, point_alpha) {
  
  # make ld decay plot
  plot <- ggplot(data, aes_string(x = data[, x], y = data[, y])) +
    geom_point(size = point_size, alpha = point_alpha) +
    coord_cartesian(ylim = c(0, 1)) +
    scale_x_continuous(labels = function(x) x/1000) +
    labs(title = "LD decay",
         x = "Physical distance (kb)",
         y = bquote("LD"~(r^2)))
  
  return(plot)
  
}

getWindowAvg <- function(data, col_r2, col_dist, win_start, win_stop) {
  
  # get r2 value average within a window
  window_avg <- data[data[, col_dist] > win_start & data[, col_dist] <= win_stop, c(col_r2, col_dist)]
  window_avg <- data.frame(t(colMeans(window_avg)), stringsAsFactors = FALSE)
  
  return(window_avg)
  
}

filterByRandomSV <- function(data, sv_name = NULL) {
  
  # get SVs
  svs <- apply(data, MARGIN = 1, function(marker) {
    return(as.character(marker[grep("^del|^ins|^dup|^inv|^tra", marker, perl = TRUE)])) 
  })
  svs <- unique(unlist(svs))
  
  if (length(svs) > 0) {
    
    # filter data based on random SVs
    sv <- ifelse(is.null(sv_name), sample(svs, size = 1, replace = FALSE), sv_name)
    data_sample <- data[which(data$SNP_A == sv | data$SNP_B == sv), ]
    
    return(data_sample)
    
  } else {
    
    stop("no SVs in the dataset")
    
  }
}



#### command line options ----

# set default
avg_ld <- FALSE
window_size <- 100
window_step <- 100
decay_random_sv <- FALSE

# get arguments
args <- commandArgs(trailingOnly = TRUE)

# assert to have the correct optional arguments
if ("--help" %in% args) usage() & q(save = "no")

if (length(args) < 2) stop(usage(), "missing positional argument(s)")

if (length(args) > 2) {
  
  opt_args <- args[-1:-2]
  opt_args_allowed <- c("--average-ld", "--window-size", "--window-step", "--decay-random-sv")
  opt_args_requested <- as.character(sapply(opt_args, function(x) unlist(strsplit(x, split = "="))[1]))
  if (any(!opt_args_requested %in% opt_args_allowed)) stop(usage(), "wrong optional argument(s)")
  
  # change default based on the argument provided
  if (any(grepl("--average-ld", opt_args_requested))) avg_ld <- getArgValue(opt_args[grep("--average-ld", opt_args)])
  if (any(grepl("--window-size", opt_args_requested))) window_size <- getArgValue(opt_args[grep("--window-size", opt_args)])
  if (any(grepl("--window-step", opt_args_requested))) window_step <- getArgValue(opt_args[grep("--window-step", opt_args)])
  if (any(grepl("--decay-random-sv", opt_args_requested))) decay_random_sv <- getArgValue(opt_args[grep("--decay-random-sv", opt_args)])
  
}



#### plot ld decay ----

ld_file <- args[1]
out_prefix <- args[2]
# ld_file <- "analysis/ld/ld_usda_rils_snp-sv_only.chr1.window-100kb.filter-0.25.ld"
# ld_file <- "analysis/snp-chip_ld/usda_22kSNPs_rils.sorted.diploid.filtered.v4.window-1000kb.filter-0.25.ld.gz"
# ld_file <- "analysis/div-panel_ld/SNP55K_maize282_AGP3_20190419_no-scaff.window-1000kb.filter-0.25.ld.gz"
# ld_file <- "analysis/ld/ld_usda_rils_snp-sv_only.chr1.window-1kb.filter-0.25.ld"
# ld_file <- "/scratch.global/della028/hirsch_lab/genomic_prediction/ld/ld_usda_rils_snp-snp_only.chr1.window-1kb.filter-0.25.ld.gz"
# out_prefix <- "analysis/ld/LD-decay_SNPs-SVs.chr1.window-100kb.filter-0.25"

# load data
ld_results <- fread(ld_file, header = TRUE, data.table = FALSE)
# if file is compressed, need to remove first line (which is copy of the header)
if (grepl(".gz", ld_file)) ld_results <- ld_results[-1,]

# add column with distance between markers
ld_results$dist_markers <- abs(as.numeric(ld_results[, 5]) - as.numeric(ld_results[, 2]))
# make sure R2 is numeric
ld_results$R2 <- as.numeric(ld_results$R2)


if (avg_ld) {
  
  cat("Running rolling window")
  
  # set starting window
  window_start <- 1
  window_stop <- window_start + (window_size - 1)
  # set df to store average of each window
  df_window_avgs <- data.frame(stringsAsFactors = FALSE)
  
  while (window_stop <= max(ld_results$dist_markers)) {
    
    # get r2 average for each window
    window_avg <- getWindowAvg(ld_results, "R2", "dist_markers", window_start, window_stop)
    df_window_avgs <- rbind(df_window_avgs, window_avg)
    # set up the start of next window
    window_start <- window_start + window_step
    window_stop <- window_start + (window_size - 1)
    
  }
  
  # get r2 average for last window
  window_avg <- getWindowAvg(ld_results, "R2", "dist_markers", window_start, window_stop)
  df_window_avgs <- rbind(df_window_avgs, window_avg)
  
  # plot average r2
  out_plot <- plot_ld_decay(data = df_window_avgs, x = "dist_markers", y = "R2", point_size = 0.2, point_alpha = 0.2) +
    geom_smooth()
  
} else {
  
  # plot ld decay without averaging r2 values
  out_plot <- plot_ld_decay(data = ld_results, x = "dist_markers", y = "R2", point_size = 0.01, point_alpha = 0.01)
  
  if (decay_random_sv) {
    
    cat("Selecting random SVs")
    
    for (i in 1:3) {
      
      ld_results_sample <- filterByRandomSV(ld_results, sv_name = NULL)
      
      if (NROW(ld_results_sample) > 0) {
        
        # plot ld decay of random SVs
        out_plot <- plot_ld_decay(data = ld_results_sample, x = "dist_markers", y = "R2", point_size = 0.7, point_alpha = 0.5) +
          labs(subtitle = sv)
        # save plot
        ggsave(filename = paste0(out_prefix, ".", sv, ".png"), plot = out_plot, device = "png")
        
      }
    }
  }
}

# save plot
if (decay_random_sv == FALSE) ggsave(filename = paste0(out_prefix, ".png"), plot = out_plot, device = "png")
