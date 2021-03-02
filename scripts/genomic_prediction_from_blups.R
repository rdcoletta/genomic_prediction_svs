####
#### code from Fernandes et al (2018, TAG)
####

library(data.table)
library(asreml)
library(cvTools)
library(plyr)
library(rrBLUP)
library(MASS)
library(reshape)
library(tidyr)
suppressWarnings(suppressMessages(library(GAPIT3)))
suppressWarnings(suppressMessages(library(doParallel)))

usage <- function() {
  cat("
description: simulate trait based on user-defined genetic architecture.

usage: Rscript genomic_prediction.R [markers_file] [...]

positional arguments:
  markers_file            hapmap file with genotypic data
  ...
  ...
  ...

optional argument:
  --help                  show this helpful message

"
  )
}



#### command line options ----

args <- commandArgs(trailingOnly = TRUE)
if ("--help" %in% args) usage() & q(save = "no")

# assert to have the correct optional arguments
if (length(args) != 6) stop(usage(), "missing positional argument(s)")

# get positional arguments
markers_file <- args[1]
blups_file <- args[2]
output_folder <- args[3]
n_folds <- as.numeric(args[4])
cv_iter <- as.numeric(args[5])
total_envs <- as.numeric(args[6])

# markers_file <- "analysis/trait_sim/datasets/usda_rils.all_markers.adjusted-n-markers.hmp.txt"
# blups_file <- "analysis/test_prediction/multi_env/with_gxe/100qtns_SVs_equal_0.5h2_pop1/blups_1st_stage.txt"
# output_folder <- "analysis/test_prediction/multi_env/with_gxe/100qtns_SVs_equal_0.5h2_pop1/prediction_all_markers"
# n_folds <- 5
# cv_iter <- 5
# total_envs <- 5


if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)


### TRY WITH DIAG STRUCTURE FIRST, THEN UNSTRUCTURED



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
means <- pivot_wider(means, names_from = environment, values_from = predicted_trait_value)
means <- data.frame(means)

# filter total number of environments
means <- means[, c(1:(total_envs + 1))]

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
# kmatrix<-as.matrix(read.table("kmatrix.txt", h=T)) ##### importing previously calculated kinship matrix.

# inverting relationship matrix
A <- ginv(kmatrix)
colnames(A) <- geno_ids
rownames(A) <- geno_ids
rm(kmatrix)

# changing inverted A matrix format to be used in asreml
A[lower.tri(A)] <- NA
A <- na.omit(reshape2::melt(A))
rownames(A) <- NULL
#ginv <- read.csv("ginv.csv")
ginv <- data.frame(A[,2], A[,1], A[,3])
colnames(ginv)<- c("Row", "Column", "GINV")
attr(ginv,"rowNames") <- geno_ids
rm(A)
#write.csv(ginv, "ginv.csv", row.names = FALSE)

# # means over location
# means2 <- ddply(means,.(geno), summarize, Y=mean(Y,na.rm=T), M=mean(M,na.rm=T),  h1=mean(h1,na.rm=T), h2=mean(h2,na.rm=T), h3=mean(h3,na.rm=T), h4=mean(h4,na.rm=T)  )



#### create groups for 5-fold CV ----

sort <- list()
for(a in 1:cv_iter){
  for(j in 1:n_folds){
    folds <- cvFolds(nlevels(means$genotype), type ="random", K = 5, R = 1)
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
gebv <- geno_ids

for(j in 1:cv_iter){
  
  cat(paste0("\n---- k-fold CV iteration ", j, " ---\n"))
  
  rm <- list()
  registerDoParallel(cores = n_folds)
  results_folds <- foreach(i = 1:n_folds, .combine = c) %dopar% {
    
    # get pheno data
    test <- means
    # mask genotypes in test population (CV2 method)
    test[which(means[, "genotype"] %in% sort[[j]][[i]]), 2:NCOL(means)] <- NA
    
    # run mixed model
    cols <- as.matrix(test[, 2:NCOL(means)])
    model <- asreml(fixed = cols ~ trait,
                    random = ~ diag(trait):vm(genotype, source = ginv),
                    residual = ~ units:diag(trait),
                    workspace = "128mb",
                    maxit = 100,
                    data = test)
    
    # predict phenotypes
    x <- as.matrix(coef(model, pattern = ~ genotype))
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
  for (env in colnames(means)[2:NCOL(means)]) {
    accuracy[[j]] <- append(accuracy[[j]], cor(ram[, env], means[, env]))
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
colnames(summary) <- colnames(means)[2:NCOL(means)]



##### write results ----

# prediction accuracy
fwrite(summary, file = paste0(output_folder, "/prediction_accuracy.txt"),
       quote = FALSE, sep = "\t", na = NA, row.names = FALSE)

# GEBVs
fwrite(gebv, file = paste0(output_folder, "/GEBVs.txt"),
       quote = FALSE, sep = "\t", na = NA, row.names = FALSE)
