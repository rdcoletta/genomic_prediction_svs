library(data.table)
library(ggplot2)
library(dplyr)
library(tidyr)
library(gtools)

usage <- function() {
  cat("
description: plot percent variance explained per QTN for simulated traits.

usage: Rscript plot_sim_pve_qtns.R [dir_sim_traits] [sv_list] [...]

positional arguments:
  dir_sim_traits          path to folder containing simulated traits from simplePHENOTYPES
  sv_list                 single-column file containing only SV IDs

optional argument:
  --help                  show this helpful message

"
  )
}



#### command line options ----

args <- commandArgs(trailingOnly = TRUE)
if ("--help" %in% args) usage() & q(save = "no")

# assert to have the correct optional arguments
if (length(args) != 2) stop(usage(), "missing positional argument(s)")

# get positional arguments
sim_folder <- args[1]
sv_list <- args[2]



#### variance explained by qtn ----

# load info about QTN names
file_qtn_info <- list.files(sim_folder, pattern = "Selected_QTNs", full.names = TRUE)
qtn_info <- fread(file_qtn_info, header = TRUE, data.table = FALSE)
# reorder markers based on position -- to match order displayed on PVE files
qtn_info <- qtn_info[order(qtn_info$chr, qtn_info$pos), ]
# get QTN names only
qtn_names <- qtn_info[, 2]

# load SV information
SVs <- fread(sv_list, header = FALSE, data.table = FALSE)
SVs <- SVs[, 1]
# identify which QTNs are SVs
qtn_types <- case_when(qtn_names %in% SVs ~ "SV",
                       !qtn_names %in% SVs ~ "SNP")
# get MAF for each QTN
qtn_mafs <- qtn_info[, "maf"]

# find files with percent variation explained
files_var_qtn <- list.files(sim_folder, pattern = "PVE", full.names = TRUE)
# create empty df to store results
pve_qtns <- data.frame()

# load file with initial effect sizes of each QTN to calculate correlation with PVE later
file_qtn_effect <- list.files(sim_folder, pattern = "Log_Sim", full.names = TRUE)
qtn_effect <- scan(file_qtn_effect, what = "character", sep = "\n")
n_log_lines_effects <- grep("$Trait_2", qtn_effect, fixed = TRUE)[1] - grep("$Trait_1", qtn_effect, fixed = TRUE)[1] - 1

for (env in 1:length(files_var_qtn)) {

  # get inital qtn effects for this environment
  env_idx_log <- grep(paste0("$Trait_", env), qtn_effect, fixed = TRUE)[1]
  qtn_effects_env <- qtn_effect[(env_idx_log + 1):(env_idx_log + n_log_lines_effects)]
  qtn_effects_env <- sapply(qtn_effects_env, function(x) {
    x <- unlist(strsplit(x, split = "] +", perl = TRUE))[2]
    x <- as.numeric(unlist(strsplit(x, " +", perl = TRUE)))
    return(x)
  })
  if (is.matrix(qtn_effects_env)) qtn_effects_env <- as.list(qtn_effects_env)
  names(qtn_effects_env) <- NULL
  qtn_effects_env <- do.call("c", qtn_effects_env)

  # find file from that environemnt
  var_qtn_env <- files_var_qtn[grep(paste0("Trait_", env, ".txt"), files_var_qtn)]
  var_qtn_env <- fread(var_qtn_env, header = FALSE, data.table = FALSE)
  # add qtn names and type
  var_qtn_env <- rbind(c("names", qtn_names), c("type", qtn_types), c("maf", qtn_mafs), c("effect", qtn_effects_env), var_qtn_env)
  # format data frame
  var_qtn_env <- data.frame(env, t(var_qtn_env[, -1]), stringsAsFactors = FALSE, row.names = NULL)
  colnames(var_qtn_env)[1:6] <- c("env", "qtn_name", "qtn_type", "qtn_maf", "qtn_effect", "qtn_number")
  colnames(var_qtn_env)[7:NCOL(var_qtn_env)] <- paste0("rep", 1:length(7:NCOL(var_qtn_env)))
  var_qtn_env <- pivot_longer(var_qtn_env, !c(env, qtn_name, qtn_type, qtn_maf, qtn_effect, qtn_number), names_to = "rep", values_to = "var") %>%
    mutate(env = paste0("env", env)) %>%
    select(env, rep, qtn_name, qtn_type, qtn_maf, qtn_effect, qtn_number, var)
  var_qtn_env$qtn_number <- as.numeric(gsub("QTN_", "", var_qtn_env$qtn_number))
  var_qtn_env$qtn_effect <- as.numeric(var_qtn_env$qtn_effect)
  var_qtn_env$var <- as.numeric(var_qtn_env$var)
  # append to main df
  pve_qtns <- rbind(pve_qtns, var_qtn_env)

}

# calculate correlation of initial qtn effects and final variance explained
pearson <- round(cor(pve_qtns$qtn_effect, pve_qtns$var, method = "pearson", use = "complete.obs"), digits = 2)
spearman <- round(cor(pve_qtns$qtn_effect, pve_qtns$var, method = "spearman", use = "complete.obs"), digits = 2)

# reorder envs to plot
pve_qtns$env <- factor(pve_qtns$env, levels = mixedsort(unique(pve_qtns$env)))

# plot percent variance explained by qtns in each environment
plot_pve_qtns <- ggplot(pve_qtns, aes(x = as.factor(qtn_number), y = var, color = qtn_type)) +
  geom_violin() +
  facet_wrap(~env, scales = "free_y") +
  labs(caption = paste0("Correlation between QTN effects and variance explained: ",
                        pearson, " (Pearson) / ",
                        spearman, " (Spearman)"),
       x = "QTN",
       y = "variance explained") +
  theme(axis.text.x = element_text(size = 5))

ggsave(filename = paste0(sim_folder, "/var_explained_by_each_qtn.pdf"),
       plot = plot_pve_qtns, device = "pdf", width = 15)

# plot percent variance explained from all qtns in each environment
# (separated by qtn type)
total_pve <- pve_qtns %>%
  group_by(env, rep, qtn_type) %>%
  summarize(var_all_qtns = sum(var))

plot_total_pve <- ggplot(total_pve, aes(x = env, y = var_all_qtns, color = qtn_type)) +
  geom_boxplot() +
  labs(x = "environment",
       y = "variance explained (all QTNs)")

ggsave(filename = paste0(sim_folder, "/var_explained_all_qtns.pdf"),
       plot = plot_total_pve, device = "pdf", width = 12)

# # plot percent variance explained by MAF
# plot_pve_maf <- ggplot(pve_qtns, aes(x = as.numeric(qtn_maf), y = as.numeric(var), color = rep)) +
#   geom_point() +
#   facet_wrap(~env, nrow = length(unique(pve_qtns$env)), scales = "free_y") +
#   labs(x = "MAF", y = "PVE")
# ggsave(filename = paste0(sim_folder, "/var_explained_by_qtn_maf_v2.pdf"),
#        plot = plot_pve_maf, device = "pdf", width = 10)

# get average variance across reps
pve_qtns_avg <- pve_qtns %>% 
  group_by(env, qtn_name, qtn_type, qtn_maf, qtn_effect, qtn_number) %>% 
  summarize(var_exp_mean = mean(var), var_exp_se = sd(var)/sqrt(n())) %>% 
  ungroup() %>% 
  arrange(env, qtn_number)

# write summary of PVE per QTNs (averaged across reps)
fwrite(x = pve_qtns_avg, file = paste0(sim_folder, "/summary_var_explained_per_qtn.txt"),
       quote = FALSE, sep = "\t", na = NA, row.names = FALSE)




#### debug ----

# sim_folder <- "analysis/test_prediction/multi_env/with_gxe/100qtns_SVs_equal_0.5h2_pop1"
# sv_list <- "data/SVs_IDs_poly.txt"
