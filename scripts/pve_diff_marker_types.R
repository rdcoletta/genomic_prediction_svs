library(data.table)
library(asreml)
library(cvTools)
library(plyr)
library(rrBLUP)
library(MASS)
library(reshape)
library(tidyr)
library(gtools)
if(!require("GAPIT3")) {
  source("http://zzlab.net/GAPIT/GAPIT.library.R")
  source("http://zzlab.net/GAPIT/gapit_functions.txt")
}
suppressWarnings(suppressMessages(library(doParallel)))


usage <- function() {
  cat("
description: simulate trait based on user-defined genetic architecture.

usage: Rscript genomic_prediction_from_blups.R [markers_file] [blups_file] [output_folder] [...]

positional arguments:
  markers_file            hapmap file with genotypic data
  blups_file              file with BLUPs in long format
  output_folder           output folder name

optional argument:
  --help                  show this helpful message
  --cv-type=VALUE         choose 'CV2 or 'CV1' cross-validation scheme (default: CV2)
  --n-folds=VALUE         number of folds to be used in cross validation
  --cv-iter=VALUE         number of cross validataion iterations
  --total-envs=VALUE      choose how many environments will be used (default: all)
  --seed=VALUE            value for set.seed (default: NULL; random number is selected)
  --envs-weight           incorporate weights (inverse of the variance of the predictions
                          from first-stage analysis) into the multi-environment model

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
cv_type <- "CV2"
n_folds <- "5"
cv_iter <- "3"
total_envs <- NULL
seed <- NULL
envs_weight <- FALSE

args <- commandArgs(trailingOnly = TRUE)
if ("--help" %in% args) usage() & q(save = "no")

# assert to have the correct optional arguments
if (length(args) < 3) stop(usage(), "missing positional argument(s)")

if (length(args) > 3) {
  
  opt_args <- args[-1:-3]
  opt_args_allowed <- c("--cv-type", "--n-folds", "--cv-iter", "--total-envs", "--seed", "--envs-weight")
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

# make sure optional arguments are valid
if (!cv_type %in% c("CV2", "CV1")) {
  stop("Optional argument '--cv-type' should be either 'CV2' or 'CV1'")
}

if (suppressWarnings(!is.na(as.numeric(n_folds)))) {
  n_folds <- as.numeric(n_folds)
} else {
  stop("Optional argument '--n-folds' should be a number")
}

if (suppressWarnings(!is.na(as.numeric(cv_iter)))) {
  cv_iter <- as.numeric(cv_iter)
} else {
  stop("Optional argument '--cv-iter' should be a number")
}

if (!is.null(total_envs)) {
  if (suppressWarnings(!is.na(as.numeric(total_envs)))) {
    total_envs <- as.numeric(total_envs)
  } else {
    stop("Optional argument '--total-envs' should be a number")
  }
}

if (is.null(seed)) {
  seed <- ceiling(runif(1, 0, 1000000))
} else {
  if (suppressWarnings(!any(is.na(as.numeric(seed))))) {
    seed <- as.numeric(seed)
  } else {
    stop("Optional argument '--seed' should be a number")
  }
}

# get positional arguments
markers_file <- args[1]
blups_file <- args[2]
output_folder <- args[3]



#### load data ----

# "header = FALSE" for compatibility with GAPIT
hapmap <- fread(markers_file, header = FALSE, data.table = FALSE)
# numericalize hapmap
hapmap <- GAPIT.HapMap(G = hapmap, SNP.effect = "Add", SNP.impute = "Major")
#   hm$GT: vector with line names
#   hm$GD: matrix where each row is a line and each column is a marker (the numbers of the cells are the numeric genotypes)
#   hm$GI: data frame with marker information (name, chromosome and position)

# length(hapmap$GI$SNP) == NCOL(markers)
# length(hapmap$GT) == NROW(markers)
markers <- data.frame(hapmap$GD - 1)
colnames(markers) <- hapmap$GI$SNP
rownames(markers) <- hapmap$GT
rm(hapmap)

# load and format phenotypic data
means <- fread(blups_file, header = TRUE, data.table = FALSE)

if (envs_weight) {
  # create a separate data frame with environmental weights
  weights <- aggregate(weight ~ environment, FUN = mean, data = means)
  weights <- weights[mixedorder(weights$environment), ]
  # remove weight column for compatibility
  means <- means[, colnames(means) != "weight"]
}

# transform data to wide format
means <- pivot_wider(means, names_from = environment, values_from = predicted_trait_value)
means <- data.frame(means)

# filter total number of environments
if (is.null(total_envs)) total_envs <- NCOL(means) - 1
means <- means[, c(1:(total_envs + 1))]
if (envs_weight) weights <- weights[1:total_envs, ]


# using only phenotypes with snp information
means <- means[means$genotype %in% rownames(markers), ]
means <- means[match(rownames(markers), means$genotype), ]
# all(rownames(markers) == means$genotype)

# changing names for numbers, but save ids for later as well
geno_names <- means$genotype
geno_ids <- 1:length(geno_names)
means$genotype <- geno_ids
means <- transform(means, genotype = factor(genotype))



#### create relationship matrix ----

# separating svs from snps
markers_svs <- markers[, grep("^del|^dup|^ins|^inv|^tra", colnames(markers), perl = TRUE)]
markers_snps <- markers[, grep("^del|^dup|^ins|^inv|^tra", colnames(markers), perl = TRUE, invert = TRUE)]

# create relationship matrix for each marker type
kmatrix_svs <- A.mat(markers_svs, return.imputed = FALSE)
kmatrix_snps <- A.mat(markers_snps, return.imputed = FALSE)
rm(markers, markers_svs, markers_snps)

# inverting relationship matrix
A_svs <- ginv(kmatrix_svs)
A_snps <- ginv(kmatrix_snps)
colnames(A_svs) <- geno_ids
rownames(A_svs) <- geno_ids
colnames(A_snps) <- geno_ids
rownames(A_snps) <- geno_ids
rm(kmatrix_svs, kmatrix_snps)

# changing inverted A matrix format to be used in asreml
A_svs[lower.tri(A_svs)] <- NA
A_snps[lower.tri(A_snps)] <- NA
A_svs <- na.omit(reshape2::melt(A_svs))
A_snps <- na.omit(reshape2::melt(A_snps))
rownames(A_svs) <- NULL
rownames(A_snps) <- NULL
ginv_svs <- data.frame(A_svs[,2], A_svs[,1], A_svs[,3])
ginv_snps <- data.frame(A_snps[,2], A_snps[,1], A_snps[,3])
colnames(ginv_svs)<- c("Row", "Column", "GINV")
colnames(ginv_snps)<- c("Row", "Column", "GINV")
attr(ginv_svs,"rowNames") <- geno_ids
attr(ginv_snps,"rowNames") <- geno_ids
rm(A_svs, A_snps)



#### estimate variance ----

# format inverted A matrix for asreml compatibility
attr(ginv_svs, "INVERSE") <- TRUE
attr(ginv_snps, "INVERSE") <- TRUE
attr(ginv_svs, "rowNames") <- as.character(attr(ginv_svs, "rowNames"))
attr(ginv_snps, "rowNames") <- as.character(attr(ginv_snps, "rowNames"))

# create empty dfs
G_init <- data.frame(stringsAsFactors = FALSE)
R_init <- data.frame(stringsAsFactors = FALSE)
for (env in 1:total_envs) {
  
  # subset env data
  means_env <- means[, c("genotype", paste0("env", env))]
  # run model
  model_env <- asreml(fixed = means_env[, 2] ~ 1,
                      random = ~ vm(genotype, source = ginv_snps) + vm(genotype, source = ginv_svs),
                      residual = ~ id(units),
                      workspace = "128mb",
                      maxit = 100,
                      trace = TRUE,
                      data = means)
  # get variance estimates for env
  G_init_env <- data.frame(Component = c(paste0("vm(genotype, source = ginv_snps):trait!trait_env", env), 
                                         paste0("vm(genotype, source = ginv_svs):trait!trait_env", env)),
                           Value = c(as.numeric(summary(model_env)$varcomp["vm(genotype, source = ginv_snps)", "component"]),
                                     as.numeric(summary(model_env)$varcomp["vm(genotype, source = ginv_svs)", "component"])),
                           Constraint = c(as.character(summary(model_env)$varcomp["vm(genotype, source = ginv_snps)", "bound"]),
                                          as.character(summary(model_env)$varcomp["vm(genotype, source = ginv_svs)", "bound"])))
  R_init_env <- data.frame(Component = paste0("units:trait!trait_env", env),
                           Value = as.numeric(summary(model_env)$varcomp["units!R", "component"]),
                           Constraint = as.character(summary(model_env)$varcomp["units!R", "bound"]))
  # add var estimates to main df
  G_init <- rbind(G_init, G_init_env)
  R_init <- rbind(R_init, R_init_env)
  
}
# add one R param to table
R_init <- rbind(data.frame(Component = "units:trait!R", Value = 1, Constraint = "F"), R_init)

# run multivariate model with snps only 
model_snps <- asreml(fixed = means[, 2:ncol(means)] ~ trait,
                     random = ~ vm(genotype, source = ginv_snps):diag(trait),
                     residual = ~ id(units):diag(trait),
                     workspace = "500mb",
                     maxit = 100,
                     trace = TRUE,
                     data = means,
                     G.param = G_init,
                     R.param = R_init)
# predict phenotypes
pred_snps <- predict(model_snps, classify = "trait:genotype", pworkspace = "128mb")
pred_snps <- data.frame(pred_snps$pvals)
pred_snps <- pred_snps[, c("genotype", "trait", "predicted.value")]
pred_snps <- pivot_wider(pred_snps, names_from = trait, values_from = predicted.value)
pred_snps <- pred_snps[match(means$genotype, pred_snps$genotype), ]
# calculate predictive ability
pred_ability_snps <- mean(diag(cor(means[, 2:ncol(means)], pred_snps[, 2:ncol(pred_snps)])))

# run multivariate model with svs only
model_svs <- asreml(fixed = means[, 2:ncol(means)] ~ trait,
                    random = ~ vm(genotype, source = ginv_svs):diag(trait),
                    residual = ~ id(units):diag(trait),
                    workspace = "500mb",
                    maxit = 100,
                    trace = TRUE,
                    data = means,
                    G.param = G_init,
                    R.param = R_init)
# predict phenotypes
pred_svs <- predict(model_svs, classify = "trait:genotype", pworkspace = "128mb")
pred_svs <- data.frame(pred_svs$pvals)
pred_svs <- pred_svs[, c("genotype", "trait", "predicted.value")]
pred_svs <- pivot_wider(pred_svs, names_from = trait, values_from = predicted.value)
pred_svs <- pred_svs[match(means$genotype, pred_svs$genotype), ]
# calculate predictive ability
pred_ability_svs <- mean(diag(cor(means[, 2:ncol(means)], pred_svs[, 2:ncol(pred_svs)])))

# run multivariate model with both snps and svs
model_all <- asreml(fixed = means[, 2:ncol(means)] ~ trait,
                    random = ~ vm(genotype, source = ginv_snps):diag(trait) + vm(genotype, source = ginv_svs):diag(trait),
                    residual = ~ id(units):diag(trait),
                    workspace = "500mb",
                    maxit = 100,
                    trace = TRUE,
                    data = means,
                    G.param = G_init,
                    R.param = R_init)
# predict phenotypes
pred_all <- predict(model_all, classify = "trait:genotype", pworkspace = "128mb")
pred_all <- data.frame(pred_all$pvals)
pred_all <- pred_all[, c("genotype", "trait", "predicted.value")]
pred_all <- pivot_wider(pred_all, names_from = trait, values_from = predicted.value)
pred_all <- pred_all[match(means$genotype, pred_all$genotype), ]
# calculate predictive ability
pred_ability_all <- mean(diag(cor(means[, 2:ncol(means)], pred_all[, 2:ncol(pred_all)])))

# calculate pve (h2?)
mean_var_snps <- mean(model_all$vparameters[grep("vm(genotype, source = ginv_snps)", names(model_all$vparameters), fixed = TRUE)])
mean_var_svs <- mean(model_all$vparameters[grep("vm(genotype, source = ginv_svs)", names(model_all$vparameters), fixed = TRUE)])
mean_var_res <- mean(model_all$vparameters[grep("units:trait!trait_env", names(model_all$vparameters), fixed = TRUE)])
mean_PVE_snps <- mean_var_snps / (mean_var_snps + mean_var_svs + mean_var_res)
mean_PVE_svs <- mean_var_svs / (mean_var_snps + mean_var_svs + mean_var_res)
h2_all <- (mean_var_snps + mean_var_svs) / (mean_var_snps + mean_var_svs + mean_var_res)

# summarize all results
summary <- data.frame(parameter = c("loglik", "h2", "pred_ability"),
                      snps_only = c(model_snps$loglik, mean_PVE_snps, pred_ability_snps),
                      svs_only = c(model_svs$loglik, mean_PVE_svs, pred_ability_svs),
                      snps_and_svs = c(model_all$loglik, h2_all, pred_ability_all),
                      stringsAsFactors = FALSE)
summary <- pivot_longer(summary, -parameter, names_to = "marker", values_to = "value")

# make sure folder exists
if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)

# write results
fwrite(summary, file = paste0(output_folder, "/summary_PVE_marker_types.txt"),
       quote = FALSE, sep = "\t", na = NA, row.names = FALSE)



#### debug ----

# markers_file <- "analysis/trait_sim/datasets/iter1/usda_rils.all_markers.adjusted-n-markers.hmp.txt"
# blups_file <- "analysis/trait_sim/multi_env/with_gxe/additive_model/equal_effects/100-QTNs_from_SV/0.7-heritability/pop1/blups_1st_stage.txt"
# output_folder <- "analysis/trait_sim/multi_env/with_gxe/additive_model/equal_effects/100-QTNs_from_SV/0.7-heritability/pop1/prediction_iter1/all_markers"
# cv_type <- "CV1"
# n_folds <- 5
# cv_iter <- 3
# total_envs <- 5
# seed <- 2021
# envs_weight <- FALSE