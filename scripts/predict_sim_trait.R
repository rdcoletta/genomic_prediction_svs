library(data.table)
library(rrBLUP)
suppressWarnings(suppressMessages(library(GAPIT3)))
suppressWarnings(suppressMessages(library(doParallel)))
suppressWarnings(suppressMessages(library(foreach)))

usage <- function() {
  cat("
description: predict simplePHENOTYPE simulated traits with GBLUP.

usage: Rscript predict_sim_trait.R [infile_name] [dir_sim_traits] [outfolder] [...]

positional arguments:
  infile_name             hapmap file containing SNPs and/or SVs
  dir_sim_traits          path to folder containing simulated traits from simplePHENOTYPES
  outfolder               name of folder to save results

optional argument:
  --help                  show this helpful message
  --multiple-envs         add this option if want to predict simulated traits across multiple environments
  --number-folds=VALUE    number of folds for k-fold cross validation (default: 5)
  --SNP-impute=VALUE      impute missing data based on the 'Major' (default), 'Minor' or 'Middle' allele frequencies
  --seed=VALUE            value for set.seed (default: NULL; random number is selected)
  
  
credits: the main function used in this script was modified from Alex Lipka's script 'k.fold.CV.Function.to.Read.In.v.1.5.R'
    
"
  )
}

getArgValue <- function(arg) {
  
  # get value from command-line arguments
  arg <- unlist(strsplit(arg, "="))
  if (length(arg) == 1) return(TRUE)
  if (length(arg) > 1) return(arg[2])
  
}

GBLUPkfoldCV <- function(pheno = NULL,
                         geno = NULL,
                         traitname = "Name_of_Trait",
                         multiple_envs = NULL,
                         path.for.results = getwd(),
                         number.of.folds = 5,
                         seed.number = 999,
                         SNP.effect = "Add",
                         SNP.impute = "Major") {

  #### numericalize genotypic data ----
  
  # make sure "header = FALSE" when reading in the hapmap for compatibility with GAPIT
  hm <- GAPIT.HapMap(G = geno, SNP.effect = SNP.effect, SNP.impute = SNP.impute)
  # hm$GT: vector with line names
  # hm$GD: matrix where each row is a line and each column is a marker (the numbers of the cells are the numeric genotypes)
  # hm$GI: data frame with marker information (name, chromosome and position)
  
  
  #### remove markers with MAF < 0.05 ----
  
  # get total number of lines
  ns <- nrow(hm$GD)
  # get the sum of the allele scores for each SNP (i.e. get the sum for each column)
  ss <- apply(hm$GD, MARGIN = 2, sum)
  # combine two situations: one where the allele coded as "2" is major; one where "0" is coded as major.
  # need to multiply by 0.5 to scale values between 0 and 1 (since GAPIT use 0, 1 and 2 notation for marker values)
  maf.matrix <- rbind((0.5 * ss/ns), (1 - (0.5 * ss/ns)))
  # copy the minor allele frequencies for all SNPs
  # for each column, get the minimum value
  maf <- apply(maf.matrix, MARGIN = 2, min)
  # find out which SNPs have MAF < 0.05
  snps.below.0.05.maf <- which(maf < 0.05)
  # remove these SNPs from hm$GD (if there are SNPs with MAF < 0.05)
  if (length(snps.below.0.05.maf) > 0) {
    hm.GD.without.snps.below.0.05.maf <- hm$GD[, -snps.below.0.05.maf]
  } else {
    hm.GD.without.snps.below.0.05.maf <- hm$GD
  }
  
  
  #### QC dataset with GAPIT -----
  
  # first, set phenotypic covariates to 1
  CV <- pheno[, 1:2]
  CV[, 2] <- 1
  colnames(CV) <- c("taxa", "overall")
  # then qc
  GK <- cbind(hm$GT, hm.GD.without.snps.below.0.05.maf)
  qc <- GAPIT.QC(Y = pheno, GT = hm$GT, CV = CV, GK = GK)
  
  
  #### calculate kinship matrix in rrBLUP ----
  
  # # debug
  # all(qc$Y[, 1] == qc$CV[, 1])
  # all(qc$Y[, 1] == qc$GK[, 1])
  # all(rownames(qc$Y) == rownames(qc$CV))
  # all(rownames(qc$Y) == rownames(qc$GK))
  
  # correct row names
  rownames(qc$GK) <- rownames(qc$Y)
  # make sure all dfs have the same sorted order of rows
  qc$Y <- qc$Y[order(as.integer(rownames(qc$Y))), ]
  qc$GK <- qc$GK[order(as.integer(rownames(qc$GK))), ]
  qc$CV <- qc$CV[order(as.integer(rownames(qc$CV))), ]
  
  # get only phenotypic values (after QC) by excluding first column
  y <- as.matrix(qc$Y[-1])
  # get only genotypic values (after QC) by excluding first column and transform it to a matrix
  G <- as.numeric(qc$GK[, -1])
  G <- matrix(G, nrow(y), ncol(qc$GK[, -1]))
  # transform genotypic values from 0,1,2 scale to -1,0,1 scale for compatibility with rrBLUP
  G <- G - 1
  # get covariate values
  cv <- (as.matrix(qc$CV[, -1]))
  # get line names
  taxa.names <- qc$CV[, 1]
  
  # calculate realized relationship matrix
  # (shrink = TRUE will give a better estimation of relationship when using low-density markers)
  A1 <- A.mat(G, shrink = TRUE)
  # inbreeding coeff = mean(diag(A1)) - 1
  
  
  #### run k-fold cross validation ----
  
  # randomly sort the number of lines, and subdivide them into subgroups (i.e. number.of.folds)
  sample.size <- NROW(y)
  sequence.sample <- rep(1:sample.size)
  set.seed(seed.number)
  random.sample <- sample(1:sample.size, replace = FALSE)
  increment <- ceiling(length(random.sample) / number.of.folds)
  # so each subgroup will have "increment" number of lines
  
  # have a "for" loop, start it at 0, and end it at k (number of folds - 1).
  # this is done so that the for loop will work correctly.
  k <- number.of.folds - 1
  
  # register number of cores for parallelizing as the number of folds
  registerDoParallel(number.of.folds)
  
  # predict phenotype for each fold
  prediction.results <- foreach(i=0:k, .combine = c) %dopar% {
    
    # get group of lines that will have their phenotype predicted (i.e. 1 subgroup)
    pred <- random.sample[((increment*i) + 1):min(((increment*i) + increment), sample.size)]
    # make sure it has the same order of rownames of trait file
    order.pred <- match(rownames(y), pred)
    order.pred <- order.pred[!is.na(order.pred)]
    pred <- pred[order.pred]
    
    # get group of lines that will be used to train the model (i.e. all the other subgroups)
    train <- random.sample[-(((increment*i) + 1):min(((increment*i) + increment), sample.size))] 
    # make sure it has the same order of rownames of trait file
    order.train <- match(rownames(y), train)
    order.train <- order.train[!is.na(order.train)]
    train <- train[order.train]
    
    if (multiple_envs == FALSE) {
      
      #### prediction for single environment ----
      
      # mask phenotypes from subgroup that will have their phenotype predicted
      yNA <- y
      yNA[pred] <- NA
      # create df with masked phenotypes, genotype id (or line id) and covariate
      data1 <- data.frame(y = yNA, gid = 1:length(y), cv = cv)
      
      # get all covariate names
      the.cv.names <- NULL
      for(j in 1:ncol(cv)) the.cv.names <- c(the.cv.names, paste("CV_", j, sep = ""))
      # fix data1 column names
      colnames(data1) <- c("y", "gid", the.cv.names)
      
      # name rows of kinship matrix with line numbers
      rownames(A1) <- 1:nrow(A1)
      # predict breeding values of masked subgroup
      ans1 <- kin.blup(data1, K=A1, geno="gid", pheno="y", covariate = the.cv.names)
      
      # get correlation between predicted BV (ans1$g[pred]) and observed (y[pred)
      r.gy <- cor(ans1$g[pred], y[pred, ])
      vector.of.observed.values <- y[pred, ]
      vector.of.predicted.values <- ans1$g[pred]
      vector.of.taxa.names <- as.character(taxa.names[pred])
      
    } else {
      
      #### prediction for multiple environments ----
      
      # create new data frame with training population for a particular fold
      y_train <- c()
      env <- c()
      gid <- c()
      for (col in 1:NCOL(y)) {
        y_train <- append(y_train, y[rownames(y) %in% train, col])
        env <- append(env, rep(col, NROW(y[rownames(y) %in% train, col])))
        # gid <- append(gid, 1:NROW(y[rownames(y) %in% train, col]))
        gid <- append(gid, rownames(y)[rownames(y) %in% train])
      }
      data1 <- data.frame(y = y_train, env = env, gid = gid)
      
      # name rows of kinship matrix with line numbers
      rownames(A1) <- 1:nrow(A1)
      # predict breeding values with environment as fixed effect
      ans1 <- kin.blup(data1, K = A1, geno = "gid", pheno = "y", fixed = "env")
      
      # get correlation between predicted BV and observed
      r.gy <- cor(ans1$g[pred], rowMeans(y[rownames(y) %in% pred, ]))
      # using mean observed values across environments (should I do one environment at a time???)
      vector.of.observed.values <- rowMeans(y[rownames(y) %in% pred, ])
      vector.of.predicted.values <- ans1$g[pred]
      vector.of.taxa.names <- as.character(taxa.names[pred])
    }
    
    results <- list(ans1 = ans1,
                    r.gy = r.gy,
                    vector.of.observed.values = vector.of.observed.values,
                    vector.of.predicted.values = vector.of.predicted.values,
                    vector.of.taxa.names = vector.of.taxa.names)
    
    return(list(results))
    
  }
  
  # clean parallel clusters
  stopImplicitCluster()
  
  # make data available for downstream analysis
  r.gy <- c()
  vector.of.observed.values <- c()
  vector.of.predicted.values <- c()
  vector.of.taxa.names <- c()
  for (i in 1:number.of.folds) {
    r.gy <- append(r.gy, prediction.results[[i]]$r.gy)
    vector.of.observed.values <- append(vector.of.observed.values, prediction.results[[i]]$vector.of.observed.values)
    vector.of.predicted.values <- append(vector.of.predicted.values, prediction.results[[i]]$vector.of.predicted.values)
    vector.of.taxa.names <- append(vector.of.taxa.names, prediction.results[[i]]$vector.of.taxa.names)
  }
  
  
  #### fit a linear regression model ----
  
  # Observed.Phenotype_i = beta_0 + beta_1*Predicrted.Phenotype_i + epsilon_i
  SLR.model <- lm(vector.of.observed.values ~ vector.of.predicted.values)
  
  # plot the regression
  pdf(paste0(path.for.results, number.of.folds, "-fold_CV_Results_", traitname, "Obs.vs.Pred.pdf"))
  plot(vector.of.observed.values ~ vector.of.predicted.values, col = "blue", xlab = "Predicted values", ylab = "Observed Values")
  abline(SLR.model)
  legend("topleft", paste("Intercept = ", round(coef(SLR.model)["(Intercept)"],2), ", Slope = ", round(coef(SLR.model)["vector.of.predicted.values"],2), sep = ""))
  dev.off()
  
  #### write output ----
  
  # once the loop is over, output the values of the correlation coefficients,
  # as well as their means, standard deviations , 
  r.gy <- c(r.gy, mean(r.gy), sd(r.gy))
  r.gy.output <- t(as.matrix(r.gy))
  
  # create a vector of column names for the k-fold cross validation output
  colname.r.gy <- NULL
  for(i in 1:number.of.folds) colname.r.gy <- c(colname.r.gy, paste("r_CV",i,sep = ""))
  colnames(r.gy.output) <- c(colname.r.gy, "r_avg", "r_std")
  
  # also add standard errors
  r.gy.output <- cbind(r.gy.output, r_SE = as.numeric(r.gy.output[1, "r_std"])/sqrt(number.of.folds))
  # and 95% confidence interval bounds
  t_value <- qt(p = 0.05/2, df = number.of.folds - 1, lower.tail = FALSE)
  r.gy.output <- cbind(r.gy.output,
                       r_lowerCI = as.numeric(r.gy.output[1, "r_avg"]) - (t_value * as.numeric(r.gy.output[1, "r_SE"])),
                       r_upperCI = as.numeric(r.gy.output[1, "r_avg"]) + (t_value * as.numeric(r.gy.output[1, "r_SE"])))
  
  # write correlation table
  write.table(r.gy.output, paste0(path.for.results, number.of.folds, "-fold_CV_Results_", traitname, ".txt"),
              quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
  
  # write out the observed and predicted values
  the.observed.predicted.and.taxa.names <- cbind(as.character(vector.of.taxa.names),
                                                 vector.of.observed.values,
                                                 vector.of.predicted.values)
  
  colnames(the.observed.predicted.and.taxa.names) <- c("SampleID",
                                                       paste("Observed_", traitname, sep = ""),
                                                       paste("Predicted_", traitname, sep = ""))
  
  write.table(the.observed.predicted.and.taxa.names,  paste0(path.for.results, "Obs.and.Predicted.Trait.values", traitname, ".txt"),
              quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
  
}


#### command line options ----

# set default
multiple_envs <- FALSE
number_folds <- 5
SNP_impute <- "Major"
seed <- NULL


args <- commandArgs(trailingOnly = TRUE)
if ("--help" %in% args) usage() & q(save = "no")

# assert to have the correct optional arguments
if (length(args) < 3) stop(usage(), "missing positional argument(s)")

if (length(args) > 3) {
  
  opt_args <- args[-1:-3]
  opt_args_allowed <- c("--multiple-envs", "--number-folds", "--SNP-impute", "--seed")
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
if (suppressWarnings(!is.na(as.integer(number_folds)))) {
  number_folds <- as.integer(number_folds) 
} else {
  stop("Optional argument '--number-folds' should be an integer")
}

if (!SNP_impute %in% c("Major", "Minor", "Middle")) {
  stop("Optional argument '--SNP-impute' should be 'Major', 'Minor' or 'Middle'")
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
infile_name <- args[1]
sim_trait_filename <- args[2]
outfolder <- args[3]
# infile_name <- "analysis/trait_sim/datasets/usda_rils.sv_markers.adjusted-n-markers.hmp.txt"
# # sim_trait_filename <- "analysis/trait_sim/additive_model/geom_series_effects/3-QTNs_from_SNP/0.2-heritability/rep1/Simulated_Data_1_Reps_Herit_0.2.txt"
# sim_trait_filename <- "analysis/trait_sim_mult-env/additive_model/geom_series_effects/3-QTNs_from_SV/0.9-heritability/rep1/Simulated_Data_1_Reps_Herit_0.9_0.9_0.9_0.9_0.9_0.9_0.9_0.9_0.9_0.9_0.9_0.9_0.9_0.9_0.9_0.9_0.9_0.9_0.9_0.9.txt"
# # outfolder <- "analysis/trait_sim/additive_model/geom_series_effects/3-QTNs_from_SNP/gs_with_max-number_markers/all_markers/0.2-heritability/rep1"
# outfolder <- "analysis/trait_sim_mult-env/additive_model/geom_series_effects/3-QTNs_from_SV/gs_with_max-number_markers/sv_markers/0.9-heritability/rep1"
# # multiple_envs <- FALSE
# multiple_envs <- TRUE
# SNP_impute <- "Major"
# seed <- 999
# number_folds <- 5


#### simulate trait ----

# read data with "header = FALSE" for compatibility with GAPIT.Hapmap() function
geno_data <- fread(infile_name, header = FALSE, data.table = FALSE)

# read in file with simulated trait replicated 50 times
sim_trait <- fread(sim_trait_filename, header = TRUE, data.table = FALSE)
# remove column with rep number
sim_trait <- sim_trait[, colnames(sim_trait) != "Rep"]


# retrieve name of trait being parsed (heritability and replicate number)
trait_name <- colnames(sim_trait)[2]


# check if directory already exists, create a new one if it doesn't
if (!dir.exists(outfolder)) dir.create(outfolder, recursive = TRUE)
# make sure output folder ends with "/", otherwise simplePHENOTYPE may mess up the directory name
if (!grepl("/$", outfolder, perl = TRUE)) outfolder <- paste0(outfolder, "/")

# generate kinship matrix and run k-fold cross validation with GBLUP
GBLUPkfoldCV(pheno = sim_trait,
             geno = geno_data,
             traitname = trait_name,
             multiple_envs = multiple_envs,
             path.for.results = outfolder,
             number.of.folds = number_folds,
             seed.number = seed,
             SNP.effect = "Add",
             SNP.impute = SNP_impute)
