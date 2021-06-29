library(data.table)
library(asreml)
library(cvTools)
library(plyr)
library(dplyr)
library(gtools)
library(rrBLUP)
library(MASS)
library(reshape)
library(tidyr)
suppressWarnings(suppressMessages(library(GAPIT3)))
suppressWarnings(suppressMessages(library(doParallel)))

usage <- function() {
  cat("
description: simulate trait based on user-defined genetic architecture.

usage: Rscript genomic_prediction_from_blups.R [markers_file] [sim_trait_filename] [output_folder] [...]

positional arguments:
  markers_file            hapmap file with genotypic data
  sim_trait_filename      simulated traits from simplePHENOTYPES
  output_folder           output folder name

optional argument:
  --help                  show this helpful message
  --n-folds=VALUE         number of folds to be used in cross validation
  --cv-iter=VALUE         number of cross validataion iterations
  --total-envs=VALUE      choose how many environments will be used (default: all)

credits:
  part of this script was modified from Fernandes et al. 2018 (TAG, doi:10.1007/s00122-017-3033-y),
  available at https://github.com/samuelbfernandes/Trait-assisted-GS, and from Dias et al. 2018
  (Heredity, doi:10.1038/s41437-018-0053-6; script kindly provided by Kaio Olímpio Das Graças Dias)
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
n_folds <- "5"
cv_iter <- "3"
total_envs <- NULL

args <- commandArgs(trailingOnly = TRUE)
if ("--help" %in% args) usage() & q(save = "no")

# assert to have the correct optional arguments
if (length(args) < 3) stop(usage(), "missing positional argument(s)")

if (length(args) > 3) {
  
  opt_args <- args[-1:-3]
  opt_args_allowed <- c("--n-folds", "--cv-iter", "--total-envs")
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
# get positional arguments
markers_file <- args[1]
sim_trait_filename <- args[2]
output_folder <- args[3]

# markers_file <- "analysis/trait_sim/datasets/usda_rils.all_markers.adjusted-n-markers.hmp.txt"
# sim_trait_filename <- "analysis/test_prediction/multi_env/with_gxe/100qtns_SVs_equal_0.5h2_pop1/Simulated_Data_3_Reps_Herit_0.5_0.5_0.5_0.5_0.5_0.5_0.5_0.5_0.5_0.5_0.5_0.5_0.5_0.5_0.5_0.5_0.5_0.5_0.5_0.5.txt"
# output_folder <- "analysis/test_prediction/multi_env/with_gxe/100qtns_SVs_equal_0.5h2_pop1/prediction_all_markers"
# n_folds <- 5
# cv_iter <- 5
# total_envs <- 5



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
pheno_pop <- fread(sim_trait_filename, header = TRUE, data.table = FALSE)
# adjust column names
colnames(pheno_pop)[1] <- "genotype"
colnames(pheno_pop)[NCOL(pheno_pop)] <- "rep"
colnames(pheno_pop)[2:(NCOL(pheno_pop) - 1)] <- paste0("env", 1:length(2:(NCOL(pheno_pop) - 1)))
# filter total number of environments
pheno_pop <- pheno_pop[, c("genotype", paste0("env", 1:total_envs), "rep")]
# transform df into long format
pheno_pop <- pivot_longer(pheno_pop, !c(genotype, rep), names_to = "environment", values_to = "trait_value") %>%
  mutate(rep = paste0("rep", rep)) %>% 
  relocate(rep, .after = environment) %>%
  arrange(environment, genotype) %>% 
  mutate(genotype = factor(genotype),
         environment = factor(environment, levels = mixedsort(unique(environment))),
         rep = factor(rep, levels = mixedsort(unique(rep))))
# make sure it's data.frame, not tibble
pheno_pop <- data.frame(pheno_pop)
# using only phenotypes with snp information
pheno_pop <- pheno_pop[pheno_pop$genotype %in% rownames(markers), ]
pheno_pop$genotype <- droplevels(pheno_pop$genotype)

# changing names for numbers, but save ids for later as well
geno_names <- pheno_pop$genotype
pheno_pop$genotype <- match(pheno_pop$genotype, levels(geno_names))
pheno_pop <- transform(pheno_pop, genotype = factor(genotype))



#### create relationship matrix ----

kmatrix <- A.mat(markers, return.imputed = FALSE)
kmatrix <- kmatrix[match(rownames(kmatrix), levels(geno_names)), match(colnames(kmatrix), levels(geno_names))]
rm(markers)
# kmatrix<-as.matrix(read.table("kmatrix.txt", h=T)) ##### importing previously calculated kinship matrix.

# inverting relationship matrix
A <- ginv(kmatrix)
colnames(A) <- match(colnames(kmatrix), levels(geno_names))
rownames(A) <- match(rownames(kmatrix), levels(geno_names))
rm(kmatrix)

# changing inverted A matrix format to be used in asreml
A[lower.tri(A)] <- NA
A <- na.omit(reshape2::melt(A))
rownames(A) <- NULL
#ginv <- read.csv("ginv.csv")
ginv <- data.frame(A[,2], A[,1], A[,3])
colnames(ginv)<- c("Row", "Column", "GINV")
attr(ginv,"rowNames") <- 1:nlevels(geno_names)
rm(A)
#write.csv(ginv, "ginv.csv", row.names = FALSE)

# # pheno_pop over location
# pheno_pop2 <- ddply(pheno_pop,.(geno), summarize, Y=mean(Y,na.rm=T), M=mean(M,na.rm=T),  h1=mean(h1,na.rm=T), h2=mean(h2,na.rm=T), h3=mean(h3,na.rm=T), h4=mean(h4,na.rm=T)  )



#### create groups for 5-fold CV ----

sort <- list()
for(a in 1:cv_iter){
  for(j in 1:n_folds){
    folds <- cvFolds(nlevels(pheno_pop$genotype), type ="random", K = 5, R = 1)
    Sample <- cbind(folds$which, folds$subsets)
    cv <- split(Sample[,2], f=Sample[,1])
  }
  sort[[a]] <- cv  
}
rm(a, folds, j, cv, Sample)




##### genomic prediction ----

# format inverted A matrix for asreml compatibility
attr(ginv, "INVERSE") <- TRUE
attr(ginv, "rowNames") <- as.character(attr(ginv, "rowNames"))

# start list to store results
accuracy <- list()
# gebv <- geno_ids
gebv <- 1:nlevels(geno_names)

for(j in 1:cv_iter){
  
  cat(paste0("\n---- k-fold CV iteration ", j, " ---\n"))
  
  rm <- list()
  registerDoParallel(cores = n_folds)
  results_folds <- foreach(i = 1:n_folds, .combine = c) %dopar% {
    
    # get pheno data
    test <- pheno_pop
    # mask genotypes in test population (CV2 method)
    test[which(pheno_pop[, "genotype"] %in% sort[[j]][[i]]), "trait_value"] <- NA
    
    # run mixed model
    # cols <- as.matrix(test[, 2:NCOL(pheno_pop)])
    resp  <- test[, "trait_value"]
    model <- asreml(fixed = resp ~ environment,   #### try: resp ~ 1 (no fixed effect)
                    random = ~ fa(environment, 2):vm(genotype, source = ginv) + environment:rep,
                    residual = ~ dsum(~ id(units)|environment),
                    workspace = "128mb",
                    maxit = 100,
                    data = test)
    
    ##### ERROR: The estimation was aborted; too many exceptions.
    
    
    if(!convergence(model, tol = 0.01)){   #### use a while loop until get convergence and parameters less than 0.1 change
      model <- update.asreml(model)}   
    if(!convergence(model, tol = 0.01)){
      model <- update.asreml(model)}   
    if(!model$converge){print ('ERROR')}
    if(!model$converge){print ('ERROR')}
    
    
    # predict phenotypes
    x <- as.matrix(coef(model, pattern = ~ genotype))  #### try: x2 <- predict(model, classify = "environment:genotype")
    rownames(x) <- sapply(rownames(x), function(name) {
      trait <- unlist(strsplit(name, split = "_"))[2]
      trait <- unlist(strsplit(trait, split = ":"))[1]
      name <- unlist(strsplit(name, split = "_"))[3]
      return(paste0(trait, "_", name))
    })
    
    # format table
    x <- data.frame(treatment = as.character(rownames(x)), effect = x)
    rownames(x) <- NULL
    x <- separate(x, col = treatment, into = c("env", "geno"), sep = "_")
    x <- pivot_wider(x, names_from = env, values_from = effect)
    x <- data.frame(x)
    
    # rm[[i]] <- cbind(sort[[j]][[i]], x[, 1][sort[[j]][[i]]])
    # rm[[i]] <- cbind(sort[[j]][[i]], x[sort[[j]][[i]], ])
    
    # return predicted phenotypes for that fold
    rm <- cbind(sort[[j]][[i]], x[sort[[j]][[i]], ])
    return(list(rm))
    
  }
  stopImplicitCluster()
  
  # combine predicted phenotypes from all folds
  ram <- rbind(results_folds[[1]], results_folds[[2]], results_folds[[3]], results_folds[[4]], results_folds[[5]])
  ram <- ram[order(ram[, 1]), ]
  
  # calculate prediction accuracy
  accuracy[[j]] <- list()
  for (env in colnames(pheno_pop)[2:NCOL(pheno_pop)]) {
    accuracy[[j]] <- append(accuracy[[j]], cor(ram[, env], pheno_pop[, env]))
  }
  accuracy[[j]] <- unlist(accuracy[[j]])
  
  # also get the prediction values in each iteration
  colnames(ram)[3:NCOL(ram)] <- paste0(colnames(ram)[3:NCOL(ram)], "_iter", j)
  gebv <- cbind(gebv, ram[, 3:NCOL(ram)])
  
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
colnames(summary) <- colnames(pheno_pop)[2:NCOL(pheno_pop)]



##### write results ----

# make sure direcotry exists
if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)

# prediction accuracy
fwrite(summary, file = paste0(output_folder, "/prediction_accuracy.txt"),
       quote = FALSE, sep = "\t", na = NA, row.names = FALSE)

# GEBVs
fwrite(gebv, file = paste0(output_folder, "/GEBVs.txt"),
       quote = FALSE, sep = "\t", na = NA, row.names = FALSE)
