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
  dir_sim_traits              path to folder containing simulated traits from simplePHENOTYPES
  sv_list                     single-column file containing only SV IDs

optional argument:
  --help                      show this helpful message
  --qtn-effect=VALUE          type of QTN effect ('add', 'dom', 'both') from simulated traits 
  --info-file-format=VALUE    format of file containing information of additive/dominance QTNs ('wide', 'long') 

"
  )
}

getArgValue <- function(arg) {
  
  # get value from command-line arguments
  arg <- unlist(strsplit(arg, "="))
  if (length(arg) == 1) return(TRUE)
  if (length(arg) > 1) return(arg[2])
  
}

loadQTNs <- function(sim_folder, pattern, format, qtn_effect) {
  
  # load file
  file_qtn_info <- list.files(sim_folder, pattern = pattern, full.names = TRUE, ignore.case = TRUE)
  qtn_info <- fread(file_qtn_info, header = TRUE, data.table = FALSE)
  
  # format columns
  if (format == "long") {
    colnames(qtn_info)[3] <- "env"
    qtn_info$env <- gsub("trait_", "env", qtn_info$env)
    qtn_info$type <- qtn_effect
    colnames(qtn_info)[4] <- "effect"
  }
  if (format == "wide") {
    qtn_info <- qtn_info %>% 
      pivot_longer(cols = contains("_eff_t"), names_to = "env", values_to = "effect") %>% 
      mutate(type = qtn_effect, env = gsub(paste0(qtn_effect, "_eff_t"), "env", env)) %>% 
      relocate(type, env, effect, .before = snp) %>%
      as.data.frame()
  } 
  
  return(qtn_info)
  
}



#### command line options ----

# set default
qtn_effect <- "add"
info_file_format <- "long"

args <- commandArgs(trailingOnly = TRUE)
if ("--help" %in% args) usage() & q(save = "no")

# assert to have the correct optional arguments
if (length(args) < 2) stop(usage(), "missing positional argument(s)")

if (length(args) > 2) {
  
  opt_args <- args[-1:-2]
  opt_args_allowed <- c("--qtn-effect", "--info-file-format")
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

# assert optional arguments
if (!qtn_effect %in% c("add", "dom", "both")) {
  stop("Optional argument '--qtn-effect' should be 'add', 'dom' or 'both'")
}

if (!info_file_format %in% c("long", "wide")) {
  stop("Optional argument '--info-file-format' should be 'long' or 'wide'")
}


# get positional arguments
sim_folder <- args[1]
sv_list <- args[2]



#### load info about QTNs ----

if (qtn_effect == "add") {
  qtn_info <- loadQTNs(sim_folder = sim_folder, pattern = "Additive_QTNs", format = info_file_format, qtn_effect = qtn_effect)
}

if (qtn_effect == "dom") {
  qtn_info <- loadQTNs(sim_folder = sim_folder, pattern = "Dominance_QTNs", format = info_file_format, qtn_effect = qtn_effect)
}

if (qtn_effect == "both") {
  # load add and dom qtns separately
  qtn_info_add <- loadQTNs(sim_folder = sim_folder, pattern = "Additive_QTNs", format = info_file_format, qtn_effect = "add")
  qtn_info_dom <- loadQTNs(sim_folder = sim_folder, pattern = "Dominance_QTNs", format = info_file_format, qtn_effect = "dom")
  # merge add and dom files
  qtn_info <- rbind(qtn_info_add, qtn_info_dom)
  rm(qtn_info_add, qtn_info_dom)
}

# load SV information
SVs <- fread(sv_list, header = FALSE, data.table = FALSE)
SVs <- SVs[, 1]

# identify which QTNs are SVs
qtn_info$var <- case_when(qtn_info$snp %in% SVs ~ "SV", !qtn_info$snp %in% SVs ~ "SNP")
qtn_info <- relocate(var, .before = snp, .data = qtn_info)



#### variance explained by qtn ----

# find files with percent variation explained
files_var_qtn <- list.files(sim_folder, pattern = "PVE", full.names = TRUE)
# create empty df to store results
pve_qtns <- data.frame()
# count number of environments
if (qtn_effect == "both") {
  n_envs <- length(files_var_qtn) / 2
} else {
  n_envs <- length(files_var_qtn)
}

for (environment in 1:n_envs) {
  
  # get inital qtn effects for this environment
  qtn_info_env <- subset(qtn_info, env == paste0("env", environment))

  # find file from that environemnt
  if (qtn_effect == "both") {
    
    # load add qtns first
    var_qtn_env_add <- files_var_qtn[grep(paste0("ADD_QTNs_Trait_", environment, ".txt"), files_var_qtn)]
    var_qtn_env_add <- fread(var_qtn_env_add, header = FALSE, data.table = FALSE)
    n_add_qtns <- NCOL(var_qtn_env_add) - 1
    # then dom qtns
    var_qtn_env_dom <- files_var_qtn[grep(paste0("DOM_QTNs_Trait_", environment, ".txt"), files_var_qtn)]
    var_qtn_env_dom <- fread(var_qtn_env_dom, header = FALSE, data.table = FALSE)
    n_dom_qtns <- NCOL(var_qtn_env_dom) - 1
    # merge tables
    var_qtn_env <- cbind(var_qtn_env_add, var_qtn_env_dom[, -1])
    var_qtn_env[1, -1] <- paste0("QTN_", 1:(n_add_qtns + n_dom_qtns))
    rm(var_qtn_env_add, n_add_qtns, var_qtn_env_dom, n_dom_qtns)
    
  } else {
    
    var_qtn_env <- files_var_qtn[grep(paste0("Trait_", environment, ".txt"), files_var_qtn)]
    var_qtn_env <- fread(var_qtn_env, header = FALSE, data.table = FALSE)
    
  }
  
  # add qtn names and type
  var_qtn_env <- rbind(c("names", qtn_info_env$snp), c("type", qtn_info_env$type), c("var", qtn_info_env$var),
                       c("maf", qtn_info_env$maf), c("effect", qtn_info_env$effect), var_qtn_env)
  # format data frame
  var_qtn_env <- data.frame(env = environment, t(var_qtn_env[, -1]), stringsAsFactors = FALSE, row.names = NULL)
  colnames(var_qtn_env)[1:7] <- c("env", "qtn_name", "qtn_type", "var", "qtn_maf", "qtn_effect", "qtn_number")
  colnames(var_qtn_env)[8:NCOL(var_qtn_env)] <- paste0("rep", 1:length(8:NCOL(var_qtn_env)))
  var_qtn_env <- var_qtn_env %>% 
    pivot_longer(contains("rep"), names_to = "rep", values_to = "pve") %>%
    mutate(env = paste0("env", env),
           qtn_number = as.numeric(gsub("QTN_", "", qtn_number)),
           qtn_effect = as.numeric(qtn_effect),
           qtn_maf = as.numeric(qtn_maf),
           pve = as.numeric(pve)) %>% 
    relocate(rep, .after = env)
  # append to main df
  pve_qtns <- rbind(pve_qtns, var_qtn_env)
  
}

# calculate correlation of initial qtn effects and final variance explained
pearson <- round(cor(pve_qtns$qtn_effect, pve_qtns$pve, method = "pearson", use = "complete.obs"), digits = 2)
spearman <- round(cor(pve_qtns$qtn_effect, pve_qtns$pve, method = "spearman", use = "complete.obs"), digits = 2)

# reorder envs to plot
pve_qtns$env <- factor(pve_qtns$env, levels = mixedsort(unique(pve_qtns$env)))
pve_qtns$var <- factor(pve_qtns$var, levels = c("SV", "SNP"))
pve_qtns$qtn_type <- factor(pve_qtns$qtn_type, levels = c("add", "dom"))


# plot percent variance explained by qtns in each environment
plot_pve_qtns <- ggplot(pve_qtns, aes(x = as.factor(qtn_number), y = pve, color = interaction(var, qtn_type))) +
  geom_violin() +
  geom_text(aes(y = 0, label = qtn_type), vjust = 1.5, size = 2) +
  facet_wrap(~env, scales = "free_y") +
  labs(x = "QTN", y = "variance explained",
       caption = paste0("Correlation between QTN effects and variance explained: ",
                        pearson, " (Pearson) / ", spearman, " (Spearman)")) +
  theme(axis.text.x = element_text(size = 5))

ggsave(filename = paste0(sim_folder, "/var_explained_by_each_qtn.pdf"),
       plot = plot_pve_qtns, device = "pdf", width = 15)

# plot percent variance explained from all qtns in each environment
# (separated by qtn type)
total_pve <- pve_qtns %>%
  group_by(env, rep, qtn_type, var) %>%
  summarize(var_all_qtns = sum(pve))

plot_total_pve <- ggplot(total_pve, aes(x = env, y = var_all_qtns, color = interaction(var, qtn_type))) +
  geom_boxplot() +
  labs(x = "environment", y = "variance explained (all QTNs)")

ggsave(filename = paste0(sim_folder, "/var_explained_all_qtns.pdf"),
       plot = plot_total_pve, device = "pdf", width = 12)

# # plot percent variance explained by MAF
# plot_pve_maf <- ggplot(pve_qtns, aes(x = as.numeric(qtn_maf), y = as.numeric(pve), color = rep)) +
#   geom_point() +
#   facet_wrap(~env, nrow = length(unique(pve_qtns$env)), scales = "free_y") +
#   labs(x = "MAF", y = "PVE")
# ggsave(filename = paste0(sim_folder, "/var_explained_by_qtn_maf_v2.pdf"),
#        plot = plot_pve_maf, device = "pdf", width = 10)

# get average variance across reps
pve_qtns_avg <- pve_qtns %>% 
  group_by(env, qtn_name, qtn_type, var, qtn_maf, qtn_effect, qtn_number) %>% 
  summarize(var_exp_mean = mean(pve), var_exp_se = sd(pve)/sqrt(n())) %>% 
  ungroup() %>% 
  arrange(env, qtn_number)

# write summary of PVE per QTNs (averaged across reps)
fwrite(x = pve_qtns_avg, file = paste0(sim_folder, "/summary_var_explained_per_qtn.txt"),
       quote = FALSE, sep = "\t", na = NA, row.names = FALSE)




#### debug ----

# {
#   sim_folder <- "tests/trait_sim_hybrids/10-add-QTNs_0-dom-QTNs_from_both/0.7-heritability/pop1"
#   qtn_effect <- "add"
#   info_file_format <- "long"
# }
# {
#   sim_folder <- "tests/trait_sim_hybrids/5-add-QTNs_5-dom-QTNs_from_both/0.7-heritability/pop1"
#   qtn_effect <- "both"
#   info_file_format <- "long"
# }
# {
#   sim_folder <- "tests/trait_sim_hybrids/0-add-QTNs_10-dom-QTNs_from_both/0.7-heritability/pop1"
#   qtn_effect <- "dom"
#   info_file_format <- "long"
# }

# {
#   sim_folder <- "tests/trait_sim_hybrids/10-add-QTNs_0-dom-QTNs_from_SNP/0.7-heritability/pop1"
#   qtn_effect <- "add"
#   info_file_format <- "wide"
# }
# {
#   sim_folder <- "tests/trait_sim_hybrids/5-add-QTNs_5-dom-QTNs_from_SNP/0.7-heritability/pop1"
#   qtn_effect <- "both"
#   info_file_format <- "wide"
# }
# {
#   sim_folder <- "tests/trait_sim_hybrids/0-add-QTNs_10-dom-QTNs_from_SNP/0.7-heritability/pop1"
#   qtn_effect <- "dom"
#   info_file_format <- "wide"
# }
# sv_list <- "data/SVs_IDs_poly.txt"
