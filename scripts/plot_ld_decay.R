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
  --unequal-windows       plot distribution of ld values across different window sizes (i.e. smaller windows for
                          closer markers and larger windows for distant markers); this is the default behavior
  --average-ld            plot average ld from 100bp windows (default)
  --window-size=VALUE     control window size in bp (default: 100)
  --window-step=VALUE     control window step in bp (default: 100)
  --raw-ld                plot scatter plot with all ld values by physical distance
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

removeLargeSVs <- function(plink_data, max_size) {

  plink_data$sv_size <- apply(plink_data, MARGIN = 1, function(marker) {

    # check which marker ("A" or "B") is a SV with name formatted as 'type.chr.start.end'
    if (grepl("^del|^dup|^ins|^inv|^tra", marker["SNP_A"], perl = TRUE)) {
      sv_start <- as.numeric(unlist(strsplit(as.character(marker["SNP_A"]), split = ".", fixed = TRUE))[3])
      sv_end <- as.numeric(unlist(strsplit(as.character(marker["SNP_A"]), split = ".", fixed = TRUE))[4])
    } else {
      sv_start <- as.numeric(unlist(strsplit(as.character(marker["SNP_B"]), split = ".", fixed = TRUE))[3])
      sv_end <- as.numeric(unlist(strsplit(as.character(marker["SNP_B"]), split = ".", fixed = TRUE))[4])
    }

    # return SV size SVs
    sv_size <- abs(sv_end - sv_start)
    return(sv_size)

  })

  # filter out SVs larger than max_size
  plink_data_filtered <- plink_data[which(plink_data[, "sv_size"] < max_size), ]

  return(plink_data_filtered)

}

distSNPtoSV <- function(plink_data) {

  distance <- apply(plink_data, MARGIN = 1, function(marker) {

    # extract marker info
    marker1 <- marker["SNP_A"]
    marker2 <- marker["SNP_B"]
    pos1 <- marker["BP_A"]
    pos2 <- marker["BP_B"]

    # check which marker ("A" or "B") is a SV with name formatted as 'type.chr.start.end'
    if (grepl("^del|^dup|^ins|^inv|^tra", marker1, perl = TRUE)) {
      sv_start <- as.numeric(unlist(strsplit(as.character(marker1), split = ".", fixed = TRUE))[3])
      sv_end <- as.numeric(unlist(strsplit(as.character(marker1), split = ".", fixed = TRUE))[4])
      snp_pos <- as.numeric(pos2)
    } else {
      sv_start <- as.numeric(unlist(strsplit(as.character(marker2), split = ".", fixed = TRUE))[3])
      sv_end <- as.numeric(unlist(strsplit(as.character(marker2), split = ".", fixed = TRUE))[4])
      snp_pos <- as.numeric(pos1)
    }

    # calculate distance from snp to sv
    snp_dist_start <- abs(sv_start - snp_pos)
    snp_dist_end <- abs(sv_end - snp_pos)

    # return distance to closest SV boundary
    if (snp_dist_start < snp_dist_end) {
      dist_to_sv <- snp_dist_start
    } else {
      dist_to_sv <- snp_dist_end
    }
    return(dist_to_sv)
  })

  return(distance)

}

plotLDdecay <- function(data, x, y, point_size, point_alpha) {

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

plotLDdecayBoxplot <- function(data, x, y) {

  # make ld decay plot
  plot <- ggplot(data, aes_string(x = as.factor(data[, x]), y = data[, y])) +
    geom_boxplot() +
    coord_cartesian(ylim = c(0, 1)) +
    labs(title = "LD decay",
         x = "Physical distance (kb)",
         y = bquote("LD"~(r^2)))

  return(plot)

}

avgR2byWindow <- function(data, col_r2, col_dist, window_size, window_step, window_start, all_points = FALSE) {

  window_stop <- window_start + (window_size - 1)

  # set df to store average of each window
  df_window_avgs <- data.frame(stringsAsFactors = FALSE)

  while (window_stop <= max(data[, col_dist])) {

    window_avg <- data[data[, col_dist] > window_start & data[, col_dist] <= window_stop, c(col_r2, col_dist)]

    if (all_points) {
      # get all data points for each window
      window_avg <- data.frame(window_avg, n_markers = NROW(window_avg), bp_stop = window_stop, stringsAsFactors = FALSE)
    } else {
      # get r2 average for each window
      window_avg <- data.frame(t(colMeans(window_avg)), n_markers = NROW(window_avg), bp_stop = window_stop, stringsAsFactors = FALSE)
    }

    df_window_avgs <- rbind(df_window_avgs, window_avg)

    # set up the start of next window
    window_start <- window_start + window_step
    window_stop <- window_start + (window_size - 1)

  }

  return(df_window_avgs)

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
unequal_windows <- TRUE
avg_ld <- FALSE
window_size <- 100
window_step <- 100
raw_ld <- FALSE
decay_random_sv <- FALSE

# get arguments
args <- commandArgs(trailingOnly = TRUE)

# assert to have the correct optional arguments
if ("--help" %in% args) usage() & q(save = "no")

if (length(args) < 2) stop(usage(), "missing positional argument(s)")

if (length(args) > 2) {

  opt_args <- args[-1:-2]
  opt_args_allowed <- c("--unequal-windows", "--average-ld", "--window-size", "--window-step", "--raw-ld", "--decay-random-sv")
  opt_args_requested <- as.character(sapply(opt_args, function(x) unlist(strsplit(x, split = "="))[1]))
  if (any(!opt_args_requested %in% opt_args_allowed)) stop(usage(), "wrong optional argument(s)")

  # change default based on the argument provided
  if (any(grepl("--unequal-windows", opt_args_requested))) unequal_windows <- getArgValue(opt_args[grep("--unequal-windows", opt_args)])
  if (any(grepl("--average-ld", opt_args_requested))) avg_ld <- getArgValue(opt_args[grep("--average-ld", opt_args)])
  if (any(grepl("--window-size", opt_args_requested))) window_size <- getArgValue(opt_args[grep("--window-size", opt_args)])
  if (any(grepl("--window-step", opt_args_requested))) window_step <- getArgValue(opt_args[grep("--window-step", opt_args)])
  if (any(grepl("--raw-ld", opt_args_requested))) raw_ld <- getArgValue(opt_args[grep("--raw-ld", opt_args)])
  if (any(grepl("--decay-random-sv", opt_args_requested))) decay_random_sv <- getArgValue(opt_args[grep("--decay-random-sv", opt_args)])

}

if (unequal_windows) avg_ld <- FALSE; window_size <- NULL; window_step <- NULL


#### plot ld decay ----

ld_file <- args[1]
out_prefix <- args[2]
# ld_file <- "analysis/ld/ld_usda_rils_snp-sv_only.chr1.window-1000kb.filter-0.25.ld"
# ld_file <- "analysis/ld/ld_usda_rils_snp-sv_only.chr2.window-100kb.filter-0.25.ld"
# ld_file <- "analysis/snp-chip_ld/usda_22kSNPs_rils.sorted.diploid.filtered.v4.window-1000kb.filter-0.25.ld.gz"
# ld_file <- "analysis/div-panel_ld/SNP55K_maize282_AGP3_20190419_no-scaff.window-1000kb.filter-0.25.ld.gz"
# ld_file <- "analysis/ld/ld_usda_rils_snp-sv_only.chr1.window-1kb.filter-0.25.ld"
# ld_file <- "/scratch.global/della028/hirsch_lab/genomic_prediction/ld/ld_usda_rils_snp-snp_only.chr1.window-1kb.filter-0.25.ld.gz"
# out_prefix <- "analysis/ld/LD-decay_SNPs-SVs.chr1.window-100kb.filter-0.25"
# out_prefix <- "test_ld_decay"

# load data
ld_results <- fread(ld_file, header = TRUE, data.table = FALSE)
# if file is compressed, need to remove first line (which is copy of the header)
if (grepl(".gz", ld_file)) ld_results <- ld_results[-1,]

# make sure R2 is numeric
ld_results$R2 <- as.numeric(ld_results$R2)

# check for SVs in dataset
if (any(grepl("^del|^ins|^dup|^inv|^tra", ld_results[, "SNP_A"], perl = TRUE)) | any(grepl("^del|^ins|^dup|^inv|^tra", ld_results[, "SNP_B"], perl = TRUE))) {

  # remove very large SVs (>1Mb) as they can have misleading results if a SNP them
  ld_results <- removeLargeSVs(ld_results, max_size = 1000000)
  # get distance from SNP to closest SV boundary
  ld_results$dist_markers <- distSNPtoSV(ld_results)

} else {

  # get distances between SNPs
  ld_results$dist_markers <- abs(as.numeric(ld_results[, "BP_B"]) - as.numeric(ld_results[, "BP_A"]))

}
# ld_results$dist_markers <- abs(as.numeric(ld_results[, "BP_B"]) - as.numeric(ld_results[, "BP_A"]))

if (unequal_windows) {

  # set default window parameters -- max distance possible is 100mb
  window_parameters <- data.frame(bp_start =    c(   1,  1001, 10001,  50001,  100001,  1000001,  10000000),
                                  bp_stop =     c(1000, 10000, 50000, 100000, 1000000, 10000000, 100000000),
                                  window_size = c( 100,  1000, 10000,  50000,  100000,  1000000,  10000000))

  # adjust parameters based on maximum distance between markers
  window_parameters <- window_parameters[which(window_parameters[, "bp_start"] < max(ld_results[, "dist_markers"])), ]

  df_plot <- data.frame(stringsAsFactors = FALSE)
  for (row in 1:NROW(window_parameters)) {

    window_start <- window_parameters[row , "bp_start"]
    window_stop <- window_parameters[row, "bp_stop"]
    window_size <- window_parameters[row , "window_size"]

    # get ld for window
    ld_window <- ld_results[which(ld_results[, "dist_markers"] > window_start & ld_results[, "dist_markers"] <= window_stop), c("R2", "dist_markers")]

    df_window <- avgR2byWindow(ld_window, "R2", "dist_markers",
                               window_size = window_size,
                               window_step = window_size,
                               window_start = window_start,
                               all_points = TRUE)

    df_plot <- rbind(df_plot, df_window)

  }

  # out_plot <- plotLDdecayBoxplot(data = df_plot, x = "bp_stop", y = "R2") +
  #   scale_x_discrete(labels = function(x) as.numeric(x)/1000)

  # adjust scale for x axis
  old_levels <- as.integer(levels(as.factor(df_plot[, "bp_stop"])))
  new_levels <- c()
  for (i in 1:length(old_levels)) {

    if (i == 1) {
      new_levels <- c(new_levels, paste0("0-", old_levels[1] / 1000))
    } else {
      new_levels <- c(new_levels, paste0(old_levels[i - 1] / 1000, "-", old_levels[i] / 1000))
    }

  }
  
  n_per_bin <- function(x) {
    return(c(y = 1.02, label = length(x)))
  }
  
  out_plot <- plotLDdecayBoxplot(data = df_plot, x = "bp_stop", y = "R2") +
    scale_x_discrete(labels = new_levels) +
    stat_summary(fun.data = n_per_bin, geom = "text", size = 2)

} else if (avg_ld) {

  cat("Running rolling window")

  # get r2 average for each window
  df_window_avgs <- avgR2byWindow(ld_results, "R2", "dist_markers", window_size, window_step, window_start = 1)

  # plot average r2
  out_plot <- plotLDdecay(data = df_window_avgs, x = "dist_markers", y = "R2", point_size = 0.2, point_alpha = 0.2) +
    geom_smooth()

} else if (raw_ld) {

  # plot ld decay without averaging r2 values
  out_plot <- plotLDdecay(data = ld_results, x = "dist_markers", y = "R2", point_size = 0.01, point_alpha = 0.01)

  if (decay_random_sv) {

    cat("Selecting random SVs")

    for (i in 1:3) {

      ld_results_sample <- filterByRandomSV(ld_results, sv_name = NULL)

      if (NROW(ld_results_sample) > 0) {

        # plot ld decay of random SVs
        out_plot <- plotLDdecay(data = ld_results_sample, x = "dist_markers", y = "R2", point_size = 0.7, point_alpha = 0.5) +
          labs(subtitle = sv)
        # save plot
        ggsave(filename = paste0(out_prefix, ".", sv, ".png"), plot = out_plot, device = "png")

      }
    }
  }
}

# save plot  --- MIGHT NEED TO MOVE THIS UP IN THE AVG_LD AND RAW_LD LOOP
if (decay_random_sv == FALSE) ggsave(filename = paste0(out_prefix, ".png"), plot = out_plot, device = "png", width = 20)
