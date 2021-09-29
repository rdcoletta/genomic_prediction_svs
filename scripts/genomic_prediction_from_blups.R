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
means <- means[match(means$genotype, rownames(markers)), ]
# all(rownames(markers) == means$genotype)

# changing names for numbers, but save ids for later as well
geno_names <- means$genotype
geno_ids <- 1:length(geno_names)
means$genotype <- geno_ids
means <- transform(means, genotype = factor(genotype))



#### create relationship matrix ----

kmatrix <- A.mat(markers, return.imputed = FALSE)
rm(markers)

# inverting relationship matrix
A <- ginv(kmatrix)
colnames(A) <- geno_ids
rownames(A) <- geno_ids
rm(kmatrix)

# changing inverted A matrix format to be used in asreml
A[lower.tri(A)] <- NA
A <- na.omit(reshape2::melt(A))
rownames(A) <- NULL
ginv <- data.frame(A[,2], A[,1], A[,3])
colnames(ginv)<- c("Row", "Column", "GINV")
attr(ginv,"rowNames") <- geno_ids
rm(A)



#### create groups for 5-fold CV ----

# i'll do a few extra CV iterations than what user supplied in the command line
# to account for possible problems during model convergence by asreml

k_folds <- list()
for(i in 1:(cv_iter * 3)){

  # set seed for reproducibility
  set.seed(seed + i)

  # divide genotype into folds according to CV scheme
  if (cv_type == "CV1") {
    folds <- cvFolds(nlevels(means$genotype), type = "random", K = n_folds, R = 1)
    sample <- cbind(folds$which, folds$subsets)
    cv <- split(sample[, 2], f = sample[, 1])
  }

  if (cv_type == "CV2") {
    folds <- cvFolds(nlevels(means$genotype), type = "random", K = n_folds, R = total_envs)
    sample <- cbind(folds$which, folds$subsets)
    # different environments in the same fold should have different genotypes selected
    cv <- list()
    for (env in 1:total_envs) {
      cv_env <- split(sample[, 1 + env], f = sample[, 1])
      names(cv_env) <- paste0("fold", names(cv_env), "_env", env)
      cv <- append(cv, cv_env)
    }
    rm(cv_env, env)
  }

  # add folds to iteration
  k_folds[[i]] <- cv

}
rm(i, folds, cv, sample)



#### genomic prediction ----

# format inverted A matrix for asreml compatibility
attr(ginv, "INVERSE") <- TRUE
attr(ginv, "rowNames") <- as.character(attr(ginv, "rowNames"))

# start list to store results
accuracy <- list()
gebv <- geno_ids

# run k-fold cross-validation if CV1 or CV2
for(j in 1:cv_iter) {

  cat(paste0("\n---- ", n_folds, "-fold ", cv_type," iteration ", j, " ---\n"))

  extra_cv_iter <- 0
  models_converged <- FALSE
  while (models_converged == FALSE) {

    registerDoParallel(cores = n_folds)
    results_folds <- try(foreach(i = 1:n_folds, .combine = c) %dopar% {

      cat("  fold ", i, ":\n", sep = "")

      # get pheno data
      means_test <- means
      # mask genotypes in test population according to CV scheme
      if (cv_type == "CV1") {
        means_test[which(means[, "genotype"] %in% k_folds[[j+extra_cv_iter]][[i]]), 2:NCOL(means)] <- NA
      }
      if (cv_type == "CV2") {
        for (env in 1:total_envs) {
          means_test[which(means[, "genotype"] %in% k_folds[[j+extra_cv_iter]][[paste0("fold", i, "_env", env)]]), 1 + env] <- NA
        }
      }

      # run mixed model
      cat("    fitting model...\n")

      if (!envs_weight) {
        # if not using weights, data can be used in the wide format
        envs_data <- as.matrix(means_test[, 2:NCOL(means)])
        # get initial values from a simpler model
        model_init <- asreml(fixed = envs_data ~ trait,
                             random = ~ diag(trait):vm(genotype, source = ginv),
                             residual = ~ id(units):diag(trait),
                             workspace = "128mb",
                             maxit = 100,
                             trace = FALSE,
                             data = means_test)
        G_init <- data.frame(Component = names(model_init$G.param$`trait:vm(genotype, source = ginv)`$trait$initial),
                             Value = as.numeric(model_init$G.param$`trait:vm(genotype, source = ginv)`$trait$initial),
                             Constraint = model_init$G.param$`trait:vm(genotype, source = ginv)`$trait$con)
        R_init <- data.frame(Component = names(model_init$R.param$`units:trait`$trait$initial),
                             Value = as.numeric(model_init$R.param$`units:trait`$trait$initial),
                             Constraint = model_init$R.param$`units:trait`$trait$con)
        # run more complex model with initial values
        model <- asreml(fixed = envs_data ~ trait,
                        random = ~ diag(trait):vm(genotype, source = ginv),
                        residual = ~ id(units):us(trait),
                        workspace = "128mb",
                        maxit = 100,
                        trace = FALSE,
                        data = means_test,
                        G.param = G_init,
                        R.param = R_init)
      } else {
        # if using weights, data need to be in long format
        means_test <- pivot_longer(means_test, -genotype, names_to = "environment", values_to = "sim_trait")
        means_test <- means_test[order(means_test$genotype, means_test$environment), ]
        means_test$environment <- as.factor(means_test$environment)
        # add weights back to data frame with simulated phenotypes
        means_test$weight <- weights[match(means_test$environment, weights$environment), "weight"]
        # get initial values from a simpler model
        model_init <- asreml(fixed = sim_trait ~ environment,
                             random = ~ diag(environment):vm(genotype, source = ginv),
                             residual = ~ id(units):diag(environment),
                             weights = weight,
                             asmv = environment,
                             family = asr_gaussian(dispersion = 1),
                             workspace = "128mb",
                             maxit = 100,
                             trace = FALSE,
                             data = means_test)
        G_init <- data.frame(Component = names(model_init$G.param$`environment:vm(genotype, source = ginv)`$environment$initial),
                             Value = as.numeric(model_init$G.param$`environment:vm(genotype, source = ginv)`$environment$initial),
                             Constraint = model_init$G.param$`environment:vm(genotype, source = ginv)`$environment$con)
        R_init <- data.frame(Component = names(model_init$R.param$`units:environment`$environment$initial),
                             Value = as.numeric(model_init$R.param$`units:environment`$environment$initial),
                             Constraint = model_init$R.param$`units:environment`$environment$con)
        # run more complex model with initial values
        model <- asreml(fixed = sim_trait ~ environment,
                        random = ~ diag(environment):vm(genotype, source = ginv),
                        residual = ~ id(units):us(environment),
                        weights = weight,
                        asmv = environment,
                        family = asr_gaussian(dispersion = 1),
                        workspace = "128mb",
                        maxit = 100,
                        trace = FALSE,
                        data = means_test,
                        G.param = G_init,
                        R.param = R_init)
      }

      # update model until converge -- do it max 5 times
      try <- 1
      while (!model$converge & try <= 5) {
        cat("    updating model until converge...\n")
        model <- update.asreml(model)
        try <- try + 1
      }
      if (!model$converge & try == 6) stop("Model didn't converge")

      # predict phenotypes
      cat("    running predictions...\n")
      if (!envs_weight) {
        pred <- predict(model, classify = "trait:genotype", pworkspace = "128mb", trace = FALSE)
        pred <- data.frame(pred$pvals)
        pred <- pred[, c("genotype", "trait", "predicted.value")]
        pred <- pivot_wider(pred, names_from = trait, values_from = predicted.value)
      } else {
        pred <- predict(model, classify = "environment:genotype", pworkspace = "128mb", trace = FALSE)
        pred <- data.frame(pred$pvals)
        pred <- pred[, c("genotype", "environment", "predicted.value")]
        pred <- pivot_wider(pred, names_from = environment, values_from = predicted.value)
      }

      # return predicted phenotypes for that fold according to CV scheme
      if (cv_type == "CV1") {
        if (all(k_folds[[j+extra_cv_iter]][[i]] == pred[k_folds[[j+extra_cv_iter]][[i]], "genotype"])) {
          # cv1 table will be in wide format
          pred_fold <- pred[k_folds[[j+extra_cv_iter]][[i]], ]
        } else {
          stop("genotypes don't match")
        }
      }
      if (cv_type == "CV2") {
        pred_fold <- data.frame()
        for (env in 1:total_envs) {
          fold_env <- names(k_folds[[j+extra_cv_iter]])[grep(paste0("fold", i, "_env", env), names(k_folds[[j+extra_cv_iter]]))]
          fold_values <- pred[k_folds[[j+extra_cv_iter]][[fold_env]], c("genotype", paste0("env", env))]
          colnames(fold_values) <- c("genotype", "pred_values")
          # cv2 table will be in long format to keep track of genotypes masked in each env
          pred_fold <- rbind(pred_fold, cbind(fold = paste0("fold", i), environment = paste0("env", env), fold_values))
        }
      }

      cat("    done!\n\n")

      return(list(pred_fold))

    })
    stopImplicitCluster()

    # if any of the folds had singularity/convergence issues, try again with
    # an another set of folds (i.e. another CV iteration)
    if (class(results_folds) == "try-error") {

      models_converged <- FALSE
      extra_cv_iter <- extra_cv_iter + cv_iter
      if (extra_cv_iter >= cv_iter * 3) stop("ASREML couldn't converge models after 3 extra CV iterations")

    } else {

      # if there was no issue, proceed to calculate accuracy across folds
      models_converged <- TRUE

    }
  }

  # combine predicted phenotypes from all folds
  pred_pheno <- do.call(rbind, results_folds)
  if (cv_type == "CV1") {
    # sort by genotype
    pred_pheno <- pred_pheno[order(pred_pheno$genotype), ]
  }
  if (cv_type == "CV2") {
    pred_pheno_merged <- data.frame()
    # combine one environment at a time
    for (env in 1:total_envs) {
      pred_pheno_env <- subset(pred_pheno, environment == paste0("env", env))
      # sort by genotype
      pred_pheno_env <- pred_pheno_env[order(pred_pheno_env$genotype), c("genotype", "pred_values")]
      colnames(pred_pheno_env) <- c("genotype", paste0("env", env))
      if (env == 1) {
        # keep genotype column
        pred_pheno_merged <- pred_pheno_env
      } else {
        # skip genotype column
        pred_pheno_merged <- cbind(pred_pheno_merged, pred_pheno_env[, 2, drop = FALSE])
      }
    }
    # overwrite original table for compatibility with the rest of the script
    pred_pheno <- pred_pheno_merged
    rm(pred_pheno_merged, pred_pheno_env)
  }

  # calculate prediction accuracy
  accuracy[[j]] <- list()
  for (env in colnames(means)[2:NCOL(means)]) {
    accuracy[[j]] <- append(accuracy[[j]], cor(pred_pheno[, env], means[, env]))
  }
  accuracy[[j]] <- unlist(accuracy[[j]])

  # also get the prediction values in each iteration
  colnames(pred_pheno)[-1] <- paste0(colnames(pred_pheno)[-1], "_iter", j)
  if (all(geno_ids == pred_pheno[, "genotype"])) {
    gebv <- cbind(gebv, pred_pheno[, 2:NCOL(pred_pheno)])
  } else {
    stop("genotypes don't match")
  }

}

# transform into a data frame
accuracy <- do.call(rbind, accuracy)

# calculate summary stats after many iterations of CV
summary <- data.frame(apply(accuracy, MARGIN = 2, function(x) {
  # mean
  mean <- mean(x)
  # std error
  se <- sd(x)/sqrt(length(x))
  # 95% confidence interval
  t_value <- qt(p = 0.05/2, df = NROW(accuracy) - 1, lower.tail = FALSE)
  lower_CI <- mean - (t_value * se)
  upper_CI <- mean + (t_value * se)
  return(c(mean = mean, se = se, lower_CI = lower_CI, upper_CI = upper_CI))
}))
colnames(summary) <- colnames(means)[2:NCOL(means)]




#### write results ----

# make sure folder exists
if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)

# prediction accuracy
summary <- cbind(stat = rownames(summary), summary)
rownames(summary) <- NULL
fwrite(summary, file = paste0(output_folder, "/prediction_accuracy.", cv_type, ".txt"),
       quote = FALSE, sep = "\t", na = NA, row.names = FALSE)

# GEBVs
colnames(gebv)[1] <- "genotype"
if (any(geno_ids != gebv$genotype)) gebv[intersect(geno_ids, gebv$genotype), ]
gebv$genotype <- geno_names
fwrite(gebv, file = paste0(output_folder, "/GEBVs.", cv_type, ".txt"),
       quote = FALSE, sep = "\t", na = NA, row.names = FALSE)

# get mean predicted performance of each genotype at each environment
gebv_top <- pivot_longer(data = gebv, cols = -genotype, names_to = c("env", "iter"), names_sep = "_", values_to = "pred_pheno")
gebv_top <- aggregate(pred_pheno ~ genotype + env, FUN = mean, data = gebv_top)

gebv_ranks <- gebv$genotype
for (i in 1:total_envs) {

  # sort by best to worst performance
  gebv_top_env <- subset(gebv_top, env == paste0("env", i))
  gebv_top_env <- gebv_top_env[order(gebv_top_env$pred_pheno, decreasing = TRUE), c("genotype", "pred_pheno")]
  rownames(gebv_top_env) <- NULL
  fwrite(gebv_top_env, file = paste0(output_folder, "/GEBVs_best-env", i, ".", cv_type, ".txt"),
         quote = FALSE, sep = "\t", na = NA, row.names = FALSE)

  # also save ranks of genotypes by environments
  gebv_ranks <- cbind(gebv_ranks, match(gebv$genotype, gebv_top_env$genotype))
  colnames(gebv_ranks)[NCOL(gebv_ranks)] <- paste0("env", i)

}
colnames(gebv_ranks)[1] <- "genotype"
fwrite(gebv_ranks, file = paste0(output_folder, "/GEBVs_ranks.", cv_type, ".txt"),
       quote = FALSE, sep = "\t", na = NA, row.names = FALSE)



#### debug ----

markers_file <- "analysis/test_prediction/linux/usda_rils.all_markers.adjusted-n-markers.iter1.hmp.txt"
# blups_file <- "analysis/test_prediction/linux/blups_1st_stage.h2-0.5.txt"
# blups_file <- "analysis/test_prediction/linux/blups_1st_stage.h2-0.9.txt"
# blups_file <- "analysis/test_prediction/linux/blups_1st_stage.qtns200.h2-0.3.txt"
# blups_file <- "analysis/test_prediction/linux/blups_1st_stage.no-gxe.qtns-10.h2-0.7.txt"
blups_file <- "analysis/test_prediction/linux/blups_1st_stage_weighted.no-gxe.qtns-10.h2-0.7.txt"
output_folder <- "analysis/test_prediction/linux/prediction_all_markers"
# markers_file <- "analysis/trait_sim/datasets/usda_rils.all_markers.adjusted-n-markers.hmp.txt"
# blups_file <- "analysis/test_prediction/multi_env/with_gxe/100qtns_SVs_equal_0.5h2_pop1/blups_1st_stage.txt"
# output_folder <- "analysis/test_prediction/multi_env/with_gxe/100qtns_SVs_equal_0.5h2_pop1/prediction_all_markers"
# blups_file <- "analysis/test_prediction/multi_env/no_gxe/100qtns_SVs_equal_0.5h2_pop1/blups_1st_stage.txt"
# output_folder <- "analysis/test_prediction/multi_env/no_gxe/100qtns_SVs_equal_0.5h2_pop1/prediction_all_markers"
# blups_file <- "analysis/test_prediction/multi_env/with_gxe/100qtns_SVs_equal_0.5h2_pop1/blups_1st_stage_weighted.txt"
# output_folder <- "analysis/test_prediction/multi_env/with_gxe/100qtns_SVs_equal_0.5h2_pop1/prediction_all_markers_weighted"
n_folds <- 5
cv_iter <- 3
total_envs <- 5
seed <- 2021
cv_type <- "CV1"
envs_weight <- FALSE
envs_weight <- TRUE

# # CV1
#                 env1        env2      env3        env4        env5
# mean     0.759196107 0.760904541 0.7600242 0.768950437 0.746356322
# se       0.005113508 0.003153316 0.0062245 0.003008618 0.004207666
# lower_CI 0.744998733 0.752149531 0.7427423 0.760597175 0.734673968
# upper_CI 0.773393481 0.769659550 0.7773062 0.777303699 0.758038676

# # CV1 weighted
#                 env1        env2        env3        env4        env5
# mean     0.757290588 0.760251005 0.759442276 0.767393203 0.746022087
# se       0.005099289 0.003198881 0.006219393 0.002903152 0.004249904
# lower_CI 0.743132692 0.751369486 0.742174472 0.759332760 0.734222460
# upper_CI 0.771448484 0.769132524 0.776710079 0.775453645 0.757821713


# # CV2
#                 env1       env2        env3        env4        env5
# mean     0.783199219 0.76027952 0.766580095 0.787550657 0.770312076
# se       0.005195691 0.00771124 0.005664869 0.003605817 0.001882068
# lower_CI 0.768773668 0.73886969 0.750851898 0.777539305 0.765086618
# upper_CI 0.797624769 0.78168936 0.782308293 0.797562008 0.775537534

# # diag() structure ~ 2min per cv iteration (my Mac)
# # fa(1) structure ~ 10min per cv iteration (my Mac)

# # CV1 -- FA structure
#                 env1        env2         env3        env4        env5
# mean     0.785480816 0.777680454 0.7859339233 0.802942808 0.789777850
# se       0.003570838 0.003508044 0.0002324807 0.003130235 0.002592054
# lower_CI 0.770116739 0.762586561 0.7849336393 0.789474492 0.778625141
# upper_CI 0.800844894 0.792774347 0.7869342072 0.816411123 0.800930559

# accuracy pheno iter 1 h2 0.9 (cv 1 diag): -0.06422884 -0.06739743 -0.07900174 -0.07144830 -0.06569532
# accuracy pheno iter 1 h2 0.9 (cv 2 diag):
# accuracy pheno iter 1 h2 0.9 (cv 1 fa1): singularity...
