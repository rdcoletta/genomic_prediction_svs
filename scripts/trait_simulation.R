library(data.table)
library(randcorr)
suppressWarnings(suppressMessages(library(simplePHENOTYPES)))
suppressWarnings(suppressMessages(library(GAPIT3)))

usage <- function() {
  cat("
description: simulate trait based on user-defined genetic architecture.

usage: Rscript trait_simulation.R [infile_name] [sv_list] [out_folder] [...]

positional arguments:
  infile_name             hapmap file containing SNPs and SVs
  sv_list                 single-column file containing only SV IDs
  out_folder              path to folder to save plots

optional argument:
  --help                  show this helpful message
  --causal-variant=VALUE  marker type controlling trait variation ('SNP', 'SV' or 'both')
  --rep=VALUE             number of reps (default: 50)
  --ntraits=VALUE         number of traits to simulate (default: 1)
  --h2=VALUE              heritability (default: 0.2; comma-separated list of values also allowed)
  --model=VALUE           additive, dominant or epistatic genetic model ('A', 'D' or 'E'), or any
                          combination of these models ('AE', 'DE' or 'ADE')
  --add-QTN-num=VALUE     number of additive loci controlling the trait (default: 3; comma-separated
                          list of values also allowed)
  --add-effect=VALUE      size of additive effect. If number is provided (default: 0.5), a geometric
                          series of effects starting at 0.5 (default) is applied to the loci. If
                          'norm' or 'unif' is provided, then additive effects will be
                          sampled from a normal or uniform distribution, respectively. Finally,
                          if 'equal' is provided, additive effects across all loci will be 0.5.
  --architecture=VALUE    genetic architecture to be simulated ('pleiotropic' for traits being
                          controlled by the same QTNs, 'partially' for traits being controlled by
                          pleiotropic and trait-specific QTNs, 'LD' for traits being exclusively
                          controlled by different QTNs)
  --cor-matrix=FILE       file containing a correlation matrix among environments to be simulated.
                          If '--ntraits' > 1 and a file is not provided, a random correlation matrix
                          will be generated, with number of rows = number of traits.
  --seed=VALUE            value for set.seed (default: NULL; random number is selected)
         
"
  )
}

getArgValue <- function(arg) {
  
  # get value from command-line arguments
  arg <- unlist(strsplit(arg, "="))
  if (length(arg) == 1) return(TRUE)
  if (length(arg) > 1) return(arg[2])
  
}

simulate_trait <- function(geno_data = NULL,
                           list_of_SV_IDs = NULL,
                           source_trait_variation = c("SNPs", "SVs", "both"),
                           # experiments
                           rep = 50,
                           # architecture
                           ntraits = 1,
                           h2 = 0.2,
                           model = c("A", "D", "E", "AE", "DE", "ADE"),
                           add_QTN_num = 3,
                           add_effect = NULL,
                           architecture = "pleiotropic",
                           cor_matrix = NULL,
                           seed = 2020,
                           # ouput
                           out_folder = getwd()) {
  
  # make sure list of SVs is provided
  if (is.null(list_of_SV_IDs)) stop("No list of SV IDs provided!")
  
  # separate data into SNPs and SVs
  geno_data_SVs <- geno_data[which(geno_data[, 1] %in% list_of_SV_IDs), ]
  geno_data_SNPs <- geno_data[which(!geno_data[, 1] %in% list_of_SV_IDs), ]
  
  # filter data based on sorce of which marker type is the causative variant 
  if (source_trait_variation == "SV") geno2trait_sim <- geno_data_SVs
  if (source_trait_variation == "SNP") geno2trait_sim <- geno_data_SNPs
  if (source_trait_variation == "both") {
    
    # sample dataset based on the marker type that has the lowest amount of markers
    n_markers_to_sample <- min(NROW(geno_data_SVs), NROW(geno_data_SNPs))
    set.seed(seed)
    SVs_to_keep <- sort(sample(1:NROW(geno_data_SVs), size = n_markers_to_sample, replace = FALSE))
    set.seed(seed)
    SNPs_to_keep <- sort(sample(1:NROW(geno_data_SNPs), size = n_markers_to_sample, replace = FALSE))
    
    # filter datasets for each marker type separately
    geno_data_SVs <- geno_data_SVs[SVs_to_keep, ]
    geno_data_SNPs <- geno_data_SNPs[SNPs_to_keep, ]
    
    # combine and order dataset where number of SVs = number of SNPs
    geno2trait_sim <- rbind(geno_data_SVs, geno_data_SNPs)
    geno2trait_sim <- geno2trait_sim[order(geno2trait_sim$chr, geno2trait_sim$pos), ]
    
  }
  
  # select QTNs
  set.seed(seed)
  marker_qtns <- sort(sample(1:NROW(geno2trait_sim), size = add_QTN_num, replace = FALSE))
  geno2trait_sim <- geno2trait_sim[marker_qtns, ]
  
  # check which type of effect will be used
  if (is.numeric(add_effect)) {
    
    # if add_effect is numeric, use geometric series
    sim_method <- "geometric"
    add_effect_value <- add_effect
    
  } else {
    
    # otherwise...
    if (add_effect == "norm") {
      
      # generate a normal distribution to draw additive effects
      sim_method <- "custom"
      dist_mean <- 0.5
      set.seed(seed)
      add_effect_value <- rnorm(n = add_QTN_num, mean = dist_mean, sd = dist_mean/3)
      
      # make sure the effects are 0 and 1
      add_effect_value <- sapply(add_effect_value, function(i) {
        if (i > 1) i = 1
        if (i < 0) i = 0
        return(i)
      })
      
      # transform to list for compatibility to simplePHENOTYPES
      add_effect_value <- list(add_effect_value)
      
    }
    
    if (add_effect == "unif") {
      
      # generate an uniform distribution to draw additive effects
      sim_method <- "custom"
      set.seed(seed)
      add_effect_value <- runif(n = add_QTN_num, min = 0, max = 1)
      # transform to list for compatibility to simplePHENOTYPES
      add_effect_value <- list(add_effect_value)
    }
    
    if (add_effect == "equal") {
      # generate an uniform distribution and sample one additive effect to be the same across all loci
      sim_method <- "custom"
      # set.seed(seed)
      # add_effect_value <- sample(runif(n = add_QTN_num, min = 0, max = 1), size = 1)
      add_effect_value <- rep(0.5, times = add_QTN_num)
      # transform to list for compatibility to simplePHENOTYPES
      add_effect_value <- list(add_effect_value)
    }
    
  }
  
  
  # adjust some parameters if multiple traits
  if (ntraits > 1) {
    
    # make sure architecture option is provided if number of traits is bigger than 1
    if(is.null(architecture)) stop("Provide genetic architecture (pleiotropic, partially, or LD) when more than 1 trait is simulated")
    # set the same additive effects for all traits
    add_effect_value <- rep(add_effect_value, ntraits)
    
    # set the same heritability for all traits
    h2 <- rep(h2, ntraits)
    
  }
  
  #### if causative variant = "both", resample SNPs again???
  
  cat("--- Simulating traits with ", add_QTN_num, " QTNs and ", h2[1], " heritability (seed number: ", seed, ") ---\n", sep = "")
  
  # simulate trait for single environment
  create_phenotypes(geno_obj = geno2trait_sim,
                    rep = 1,
                    # architecture
                    ntraits = ntraits,
                    h2 = h2,
                    model = model,
                    add_QTN_num = add_QTN_num,
                    add_effect = add_effect_value,
                    architecture = architecture,
                    sim_method = sim_method,
                    vary_QTN = FALSE,
                    cor = cor_matrix,
                    seed = seed,
                    # ouput
                    export_gt = TRUE,
                    home_dir = getwd(),
                    output_dir = out_folder,
                    out_geno = NULL,
                    quiet = TRUE,
                    # numericalization
                    SNP_effect = "Add",
                    SNP_impute = "Major")
  
}




#### command line options ----

# set default
causal_variant <- "SNP"
rep <- "50"
ntraits <- "1"
h2 <- "0.2"
model <- "A"
add_QTN_num <- "3"
add_effect <- "0.5"
architecture <- "pleiotropic"
cor_matrix <- NULL
seed <- NULL

args <- commandArgs(trailingOnly = TRUE)
if ("--help" %in% args) usage() & q(save = "no")

# assert to have the correct optional arguments
if (length(args) < 3) stop(usage(), "missing positional argument(s)")

if (length(args) > 3) {
  
  opt_args <- args[-1:-3]
  opt_args_allowed <- c("--causal-variant", "--rep", "--ntraits", "--h2", "--model", "--add-QTN-num", 
                        "--add-effect", "--architecture", "--seed",  "--cor-matrix")
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
if (!causal_variant %in% c("SNP", "SV", "both")) {
  stop("Optional argument '--causal-variant' should be 'SNP', 'SV' or 'both'")
}

if (!model %in% c("A", "D", "E", "AE", "DE", "ADE")) {
  stop("Optional argument '--model' should be 'A', 'D', 'E', 'AE', 'DE' or 'ADE'")
}

if (!architecture %in% c("pleiotropic", "partially", "LD")) {
  stop("Optional argument '--architecture' should be 'pleiotropic', 'partially' or 'LD'")
}

if (suppressWarnings(!is.na(as.integer(rep)))) {
  rep <- as.integer(rep) 
} else {
  stop("Optional argument '--rep' should be an integer")
}

if (suppressWarnings(!is.na(as.integer(ntraits)))) {
  ntraits <- as.integer(ntraits) 
} else {
  stop("Optional argument '--ntraits' should be an integer")
}

if (suppressWarnings(!is.na(as.numeric(add_effect)))) {
  add_effect <- as.numeric(add_effect) 
} else {
  if (!add_effect %in% c("norm", "unif", "equal")) {
    stop("Optional argument '--add-effect' should be a either a number or 'norm', 'unif', or 'equal'")
  }
}

h2 <- unlist(strsplit(h2, split = ","))
if (suppressWarnings(!any(is.na(as.numeric(h2))))) {
  h2 <- as.numeric(h2) 
} else {
  stop("Optional argument '--h2' should be a number or a comma-separated list of numbers")
}

add_QTN_num <- unlist(strsplit(add_QTN_num, split = ","))
if (suppressWarnings(!any(is.na(as.integer(add_QTN_num))))) {
  add_QTN_num <- as.integer(add_QTN_num) 
} else {
  stop("Optional argument '--add-QTN-num' should be an integer or a comma-separated list of integers")
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

if (ntraits > 1) {
  if (is.null(cor_matrix)) {
    set.seed(seed)
    cor_matrix <- apply(randcorr(ntraits), MARGIN = c(1,2), function(x) as.numeric(as.character(x)))
  } else {
    cor_matrix <- fread(cor_matrix, header = FALSE, data.table = FALSE, stringsAsFactors = FALSE)
    cor_matrix <- as.matrix(cor_matrix)
  }
}

# print conditions selected
cat("\nNumber of traits to be simulated: ", ntraits, "\n", sep = "")
cat("Number of QTNs: ", paste0(add_QTN_num, collapse = ", "), "\n", sep = "")
cat("Causal variant: ", causal_variant, "\n", sep = "")
cat("Effect size: ", add_effect, "\n", sep = "")
cat("Genetic model: ", model, "\n", sep = "")
cat("Architecture: ", architecture, "\n", sep = "")
cat("Heritability: ", paste0(h2, collapse = ", "), "\n", sep = "")
cat("Number of reps: ", rep, "\n", sep = "")
cat("Seed number: ", seed, "\n\n", sep = "")

# get positional arguments
infile_name <- args[1]
sv_list <- args[2]
out_folder <- args[3]
# infile_name <- "data/test_usda_rils_projected-SVs-only.poly.hmp.txt"
# sv_list <- "data/test_SVs_IDs.txt"
# out_folder <- "analysis/trait_sim/additive_model/geom_series_effects/3-QTNs_from_SNP/0.2-heritability"
# out_folder <- "analysis/trait_sim_mult-env/additive_model/geom_series_effects/3-QTNs_from_SNP/0.2-heritability"
# causal_variant <- "SNP"
# # causal_variant <- "SV"
# rep <- 3
# # ntraits <- 1
# ntraits <- 5
# # h2 <- 0.2
# h2 <- 0.5
# model <- "A"
# add_QTN_num <- 3
# # add_QTN_num <- 25
# # add_effect <- 0.5
# add_effect <- "equal"
# architecture <- "pleiotropic"
# seed <- 2020
# set.seed(seed); cor_matrix <- apply(randcorr(ntraits), MARGIN = c(1,2), function(x) as.numeric(as.character(x)))
# # out_folder <- paste0("analysis/trait_sim_mult-env/additive_model/geom_series_effects/", add_QTN_num,
# #                      "-QTNs_from_", causal_variant, "/", h2, "-heritability")
# out_folder <- paste0("analysis/trait_sim_mult-env/additive_model/equal_effects/", add_QTN_num,
#                      "-QTNs_from_", causal_variant, "/", h2, "-heritability")


#### simulate trait ----

# load genotypic data
geno_data <- fread(infile_name, header = TRUE, data.table = FALSE)

# load SV information
SVs <- fread(sv_list, header = FALSE, data.table = FALSE)
SVs <- SVs[, 1]

# check if directory already exists, create a new one if it doesn't
if (!dir.exists(out_folder)) dir.create(out_folder, recursive = TRUE)

for (rep_number in 1:rep) {
  
  # create folder for rep
  out_folder_rep <- paste0(out_folder, "/rep", rep_number)
  # get rep specific seed number
  seed_rep <- seed + rep_number
  
  cat ("rep #", rep_number, "\n", sep = "")
  
  # simulate traits with different QTNs per experiment
  simulate_trait(geno_data = geno_data,
                 list_of_SV_IDs = SVs,
                 source_trait_variation = causal_variant,
                 # experiments
                 rep = 1,
                 # architecture
                 ntraits = ntraits,
                 h2 = h2,
                 model = model,
                 add_QTN_num = add_QTN_num,
                 add_effect = add_effect,
                 architecture = architecture,
                 cor_matrix = cor_matrix,
                 seed = seed_rep,
                 # ouput
                 out_folder = out_folder_rep)
  
}

