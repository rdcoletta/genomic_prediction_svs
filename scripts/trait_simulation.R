library(data.table)
library(randcorr)
suppressWarnings(suppressMessages(library(simplePHENOTYPES)))
suppressWarnings(suppressMessages(library(GAPIT3)))

usage <- function() {
  cat("
description: simulate trait based on user-defined genetic architecture.

usage: Rscript trait_simulation.R [infile_name] [sv_list] [out_folder] [...]

positional arguments:
  infile_name               hapmap file containing SNPs and SVs
  sv_list                   single-column file containing only SV IDs
  out_folder                path to folder to save plots

optional argument:
  --help                    show this helpful message
  --causal-variant=VALUE    marker type controlling trait variation ('SNP', 'SV' or 'both')
  --SNP-SV-ratio=VALUE      if both SNPs and SVs can be the causal variants, this option controls
                            the SNP/SV ratio to be sampled (default: 0.5)
  --pops=VALUE              number of populations to be simulated with different QTNs (default: 10)
  --reps=VALUE              number of reps (default: 3)
  --envs=VALUE              number of environments to simulate (default: 1)
  --h2=VALUE                heritability (default: 0.5; comma-separated list of values also allowed)
  --model=VALUE             additive, dominant or epistatic genetic model ('A', 'D' or 'E'), or any
                            combination of these models ('AE', 'DE' or 'ADE')
  --add-QTN-num=VALUE       number of additive loci controlling the trait (default: 10)
  --add-effect-type=VALUE   type of additive effect among QTNs. If 'equal' is provided (default),
                            additive effects across all loci will be '--marker-effect'. If 'geometric'
                            is provided, a geometric series of effects starting at '--marker-effect'
                            is applied to the loci
  --marker-effect=VALUE     size of marker effect (default: 0.1)
  --SV-effect=VALUE         if both SNPs and SVs can be the causal variants, this option allows
                            providing a different effect size for SVs (default: NULL; same as
                            '--marker-effect')
  --architecture=VALUE      genetic architecture to be simulated ('pleiotropic' for traits being
                            controlled by the same QTNs, 'partially' for traits being controlled by
                            pleiotropic and trait-specific QTNs, 'LD' for traits being exclusively
                            controlled by different QTNs)
  --gen-cor-matrix=VALUE    if option is not provided, no genetic correlation matrix will be used.
                            If '/path/to/file', then the file should contain a genetic correlation
                            matrix among environments to be simulated (no column or row names). The
                            lower the genetic correlation among environments, higher the chance of
                            having GxE interactions. If '1', a perfect positive correlation matrix
                            (i.e. all environments have correlation = 1) will be generated with
                            'number of rows' = 'number of envs'
  --res-cor-matrix=VALUE    if option is not provided, a random residual correlation matrix will be
                            generated with 'number of rows' = 'number of envs'.  If '/path/to/file',
                            then the file should contain a residual correlation matrix among
                            environments to be simulated (no column or row names)
  --gxe                     add this option to simulate GxE by adding random effects (drawn from a
                            normal distribution) to the marker effects in each environment. This
                            option also overrides any genetic correlation matrix provided
  --diff-dist-gxe           add this option to draw random GxE effects from different normal
                            distributions depending on marker type (SNP or SV). The normal distribution
                            used to sample effects will have mean 0 for both SNPs and SVs, but the
                            standard deviation for SVs will be multiplied by the value of '--SV-effect'
                            option.
  --diff-env-mean           add this option to force simulated traits to have different means in
                            each environment
  --QTN-variance            add this option to write files with percent variance explained per QTN
  --seed=VALUE              value for set.seed (default: NULL; random number is selected)

"
  )
}

getArgValue <- function(arg) {

  # get value from command-line arguments
  arg <- unlist(strsplit(arg, "="))
  if (length(arg) == 1) return(TRUE)
  if (length(arg) > 1) return(arg[2])

}

createGeomSeries <- function(qtn_number, largest_effect) {

  list_effects <- rep(largest_effect, times = qtn_number)
  order_effects <- 1:qtn_number
  geom_series <- list(list_effects ^ order_effects)

  return(geom_series)

}

simulateTraits <- function(geno_data = NULL,
                           list_of_SV_IDs = NULL,
                           source_trait_variation = c("SNPs", "SVs", "both"),
                           SNP_SV_ratio = 0.5,
                           # experiments
                           pops = 10,
                           reps = 3,
                           # architecture
                           envs = 1,
                           h2 = 0.5,
                           model = c("A", "D", "E", "AE", "DE", "ADE"),
                           add_QTN_num = 3,
                           add_effect_type = c("geometric", "equal"),
                           marker_effect = 0.1,
                           SV_effect = NULL,
                           architecture = "pleiotropic",
                           diff_env_mean = FALSE,
                           gxe = FALSE,
                           diff_dist_gxe = FALSE,
                           gen_cor_matrix = NULL,
                           res_cor_matrix = NULL,
                           seed = 2020,
                           # ouput
                           QTN_variance = QTN_variance,
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

    # get how many SNPs and SVs will be selected
    n_SNPs_to_sample <- ceiling(add_QTN_num * SNP_SV_ratio)
    n_SVs_to_sample <- ceiling(add_QTN_num * (1 - SNP_SV_ratio))

    # sample dataset
    set.seed(seed)
    SVs_to_keep <- sort(sample(1:NROW(geno_data_SVs), size = n_SVs_to_sample, replace = FALSE))
    set.seed(seed)
    SNPs_to_keep <- sort(sample(1:NROW(geno_data_SNPs), size = n_SNPs_to_sample, replace = FALSE))

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
  if (add_effect_type == "geometric") {

    # create a geometric series where first QTN effect is 'marker_effect"
    add_effect_value <- createGeomSeries(qtn_number = add_QTN_num, largest_effect = marker_effect)

  }

  if (add_effect_type == "equal") {

    # make sure all QTNs have the same 'marker_effect' value
    add_effect_value <- rep(marker_effect, times = add_QTN_num)
    # transform to list for compatibility to simplePHENOTYPES
    add_effect_value <- list(add_effect_value)

    # adjust effect size if SVs have different effects than SNPs
    if (!is.null(SV_effect)) {

      # find which qtns are SVs
      # sv_qtns <- grep("^del|^dup|^ins|^inv|^tra", geno2trait_sim[, 1], perl = TRUE)
      sv_qtns <- which(geno2trait_sim[, 1] %in% list_of_SV_IDs)
      # modify effect size of svs
      if (length(sv_qtns) > 0) add_effect_value[[1]][sv_qtns] <- SV_effect

    }
  }

  # #### FOR TESTING
  # # make QTNs to have 3 different effects
  # qtns_1 <- 1:ceiling(add_QTN_num/3)
  # qtns_2 <- (length(qtns_1) + 1):(length(qtns_1) + floor(add_QTN_num/3))
  # qtns_3 <- (length(qtns_1) + length(qtns_2) + 1):(add_QTN_num)
  # add_effect_value[[1]][qtns_2] <- 0.4
  # add_effect_value[[1]][qtns_3] <- 0.7

  # adjust some parameters if multiple traits
  if (envs > 1) {

    # make sure architecture option is provided if number of traits is bigger than 1
    if(is.null(architecture)) stop("Provide genetic architecture (pleiotropic, partially, or LD) when more than 1 trait is simulated")
    # set the same additive effects for all traits
    add_effect_value <- rep(add_effect_value, envs)

    if (diff_env_mean) {

      # scale effect sizes among environments so they have different means
      for (env in 1:envs) {
        # use the absolute mean correlation among environments as the scaling factor
        mean_cor_env <- abs(mean(res_cor_matrix[env, -env]))
        add_effect_value[[env]] <- round(add_effect_value[[env]] * mean_cor_env, digits = 4)
      }

    }

    # set the same heritability for all traits
    h2 <- rep(h2, envs)

    # if GxE, change effect sizes for each environment
    if (gxe) {

      # select QTNs to have constant effects across envs -- 10% QTNs without GxE
      n_constant_qtns <- ceiling(0.1 * add_QTN_num)
      constant_seed <- seed + seed
      set.seed(constant_seed)
      constant_qtns <- sort(sample(1:add_QTN_num, size = n_constant_qtns, replace = FALSE))

      # select QTNs that can change signs across environments -- 10% QTNs can change signs
      n_diff_sign_qtns <- ceiling(0.1 * add_QTN_num)
      diff_sign_seed <- seed * envs
      set.seed(diff_sign_seed)
      diff_sign_qtns <- sort(sample((1:add_QTN_num)[-constant_qtns], size = n_diff_sign_qtns, replace = FALSE))

      for (env in 1:length(add_effect_value)) {

        # select QTNs to have no effect in this env -- 10% QTNs with no effect at all
        n_zero_qtns <- ceiling(0.1 * add_QTN_num)
        zero_seed <- (seed + seed) * env
        set.seed(zero_seed)
        zero_qtns <- sort(sample((1:add_QTN_num)[-constant_qtns], size = n_zero_qtns, replace = FALSE))

        # draw specific effects from a normal distribution
        gxe_seed <- seed + env
        set.seed(gxe_seed)
        gxe_effect <- rnorm(n = add_QTN_num, mean = 0, sd = 0.3)
        # if requested, draw GxE effects for SVs from a different distribution
        if (diff_dist_gxe & !is.null(SV_effect)) {
          if (SV_effect != marker_effect & length(sv_qtns) > 0) {
            set.seed(gxe_seed)
            gxe_effect[sv_qtns] <- rnorm(n = length(sv_qtns), mean = 0, sd = (SV_effect / marker_effect) * 0.3)
          }
        }
        gxe_effect <- round(gxe_effect, digits = 4)
        # remove gxe effects for qtns with zero or constant effects
        gxe_effect <- gxe_effect[-c(zero_qtns, constant_qtns)]
        # add effects for each qtn
        add_effect_value[[env]][zero_qtns] <- 0
        add_effect_value[[env]][-c(zero_qtns, constant_qtns)] <- add_effect_value[[env]][-c(zero_qtns, constant_qtns)] + gxe_effect
        add_effect_value[[env]] <- sapply(add_effect_value[[env]], function(i) {
          if (i > 1) i = 1
          if (i < -1) i = -1
          return(i)
        })

        #### DO THIS IF QTN CAN CHANGE SIGN ACCORDING TO PREVIOUS ENVS
        # # allow only some QTNs to change the sign of their effects in different envs
        # if (env > 1) {
        #
        #   # first, get all QTNs that changed sign
        #   all_diff_sign_qtns <- sort(c(which(add_effect_value[[env]] > 0 & add_effect_value[[env - 1]] < 0),
        #                                which(add_effect_value[[env]] < 0 & add_effect_value[[env - 1]] > 0)))
        #   # select QTNs to don't change sign from previous env
        #   n_diff_sign_qtns <- ceiling(0.05 * add_QTN_num)
        #   diff_sign_seed <- (seed + env) * env
        #   set.seed(diff_sign_seed)
        #   keep_sign_qtns <- sort(sample(all_diff_sign_qtns, size = (length(all_diff_sign_qtns) - n_diff_sign_qtns), replace = FALSE))
        #   # add back original sign
        #   add_effect_value[[env]][keep_sign_qtns] <- add_effect_value[[env]][keep_sign_qtns] * (-1)
        #
        # }
        ####

      }

      # make sure only selected QTNs are allowed to change signs
      for (env in 1:length(add_effect_value)) {

        add_effect_value[[env]][-diff_sign_qtns] <- sapply(add_effect_value[[env]][-diff_sign_qtns], function(qtn, effect) {
          # change signs according to initial 'marker_effect' provided
          if (effect > 0) if (qtn < 0) qtn = qtn * (-1)
          if (effect < 0) if (qtn > 0) qtn = qtn * (-1)
          return(qtn)
        }, effect = marker_effect)

      }

      # make sure to disable the use of any genetic correlation among environments
      gen_cor_matrix <- NULL

    }

  }

  cat("--- Simulating traits with ", add_QTN_num, " QTNs and ", h2[1], " heritability (seed number: ", seed, ") ---\n", sep = "")

  if (envs == 1) {
    cat("QTN effects: ", paste0(add_effect_value[[1]], collapse = " "), "\n", sep = "")
  } else {
    cat("QTN effects: env1 - ", paste0(add_effect_value[[1]], collapse = " "), "\n", sep = "")
    for (i in 2:length(add_effect_value)) {
      cat("             env", i, " - ", paste0(add_effect_value[[i]], collapse = " "), "\n", sep = "")
    }
  }

  # simulate trait for single environment
  create_phenotypes(geno_obj = geno2trait_sim,
                    rep = reps,
                    # architecture
                    ntraits = envs,
                    h2 = h2,
                    model = model,
                    add_QTN_num = add_QTN_num,
                    add_effect = add_effect_value,
                    architecture = architecture,
                    sim_method = "custom",
                    vary_QTN = FALSE,
                    cor = gen_cor_matrix,
                    cor_res = res_cor_matrix,
                    seed = seed,
                    # ouput
                    export_gt = TRUE,
                    home_dir = getwd(),
                    output_dir = out_folder,
                    out_geno = NULL,
                    QTN_variance = QTN_variance,
                    quiet = TRUE,
                    # numericalization
                    SNP_effect = "Add",
                    SNP_impute = "Major")

}




#### command line options ----

# set default
causal_variant <- "SNP"
SNP_SV_ratio <- NULL
pops <- "10"
reps <- "3"
envs <- "1"
h2 <- "0.5"
model <- "A"
add_QTN_num <- "10"
add_effect_type <- "equal"
marker_effect <- "0.1"
SV_effect <- NULL
architecture <- "pleiotropic"
gxe <- FALSE
diff_dist_gxe <- FALSE
gen_cor_matrix <- NULL
res_cor_matrix <- NULL
diff_env_mean <- FALSE
QTN_variance <- FALSE
seed <- NULL

args <- commandArgs(trailingOnly = TRUE)
if ("--help" %in% args) usage() & q(save = "no")

# assert to have the correct optional arguments
if (length(args) < 3) stop(usage(), "missing positional argument(s)")

if (length(args) > 3) {

  opt_args <- args[-1:-3]
  opt_args_allowed <- c("--causal-variant", "--SNP-SV-ratio", "--pops", "--reps", "--envs", "--h2",
                        "--model", "--add-QTN-num", "--add-effect-type", "--marker-effect",
                        "--SV-effect", "--architecture", "--gen-cor-matrix", "--res-cor-matrix",
                        "--gxe", "--diff-dist-gxe", "--diff-env-mean", "--QTN-variance", "--seed")
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

if (!is.null(SNP_SV_ratio)) {
  if (suppressWarnings(!is.na(as.numeric(SNP_SV_ratio)))) {
    SNP_SV_ratio <- as.numeric(SNP_SV_ratio)
  } else {
    stop("Optional argument '--SNP-SV-ratio' should be a number")
  }
} else {
  # if 'SNP_SV_ratio' is NULL, set default to 0.5 if 'causal_variant' is both
  if (causal_variant == "both") SNP_SV_ratio <- 0.5
}

if (suppressWarnings(!is.na(as.integer(pops)))) {
  pops <- as.integer(pops)
} else {
  stop("Optional argument '--pops' should be an integer")
}

if (suppressWarnings(!is.na(as.integer(reps)))) {
  reps <- as.integer(reps)
} else {
  stop("Optional argument '--reps' should be an integer")
}

if (suppressWarnings(!is.na(as.integer(envs)))) {
  envs <- as.integer(envs)
} else {
  stop("Optional argument '--envs' should be an integer")
}

h2 <- unlist(strsplit(h2, split = ","))
if (suppressWarnings(!any(is.na(as.numeric(h2))))) {
  h2 <- as.numeric(h2)
} else {
  stop("Optional argument '--h2' should be a number or a comma-separated list of numbers")
}

if (!model %in% c("A", "D", "E", "AE", "DE", "ADE")) {
  stop("Optional argument '--model' should be 'A', 'D', 'E', 'AE', 'DE' or 'ADE'")
}

if (suppressWarnings(!any(is.na(as.integer(add_QTN_num))))) {
  add_QTN_num <- as.integer(add_QTN_num)
} else {
  stop("Optional argument '--add-QTN-num' should be an integer")
}

if (!add_effect_type %in% c("equal", "geometric")) {
  stop("Optional argument '--add-effect-type' should be a either 'equal' or 'geometric'")
}

if (suppressWarnings(!is.na(as.numeric(marker_effect)))) {
  marker_effect <- as.numeric(marker_effect)
} else {
  stop("Optional argument '--marker-effect' should be an integer")
}

if (!is.null(SV_effect)) {
  if (suppressWarnings(!is.na(as.numeric(SV_effect)))) {
    SV_effect <- as.numeric(SV_effect)
  } else {
    stop("Optional argument '--SV-effect' should be a number")
  }
}

if (!architecture %in% c("pleiotropic", "partially", "LD")) {
  stop("Optional argument '--architecture' should be 'pleiotropic', 'partially' or 'LD'")
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

if (envs > 1) {
  # genetic correlation
  if (!is.null(gen_cor_matrix)) {
    if (gen_cor_matrix == "1") {
      gen_cor_matrix <- matrix(data = rep(1, times = envs * envs), nrow = envs, ncol = envs)
    } else {
      gen_cor_matrix <- fread(gen_cor_matrix, header = FALSE, data.table = FALSE, stringsAsFactors = FALSE)
      # remove header, if it exists
      if (nrow(gen_cor_matrix) > ncol(gen_cor_matrix)) gen_cor_matrix <- apply(gen_cor_matrix[-1, ], c(1,2), as.numeric)
      gen_cor_matrix <- as.matrix(gen_cor_matrix)
    }
  }
  # residual correlation
  if (is.null(res_cor_matrix)) {
    set.seed(seed)
    res_cor_matrix <- apply(randcorr(envs), MARGIN = c(1,2), function(x) as.numeric(as.character(x)))
  } else {
    res_cor_matrix <- fread(res_cor_matrix, header = FALSE, data.table = FALSE, stringsAsFactors = FALSE)
    # remove header, if it exists
    if (nrow(res_cor_matrix) > ncol(res_cor_matrix)) res_cor_matrix <- apply(res_cor_matrix[-1, ], c(1,2), as.numeric)
    res_cor_matrix <- as.matrix(res_cor_matrix)
  }
}

if (gxe & envs == 1) {
  gen_cor_matrix <- NULL
  stop("Optional argument '--gxe' only valid if '--envs' > 1")
}

# print conditions selected
cat("\nNumber of environments to be simulated: ", envs, "\n", sep = "")
cat("Number of QTNs per environment: ", add_QTN_num, "\n", sep = "")
cat("Causal variant: ", ifelse(causal_variant == "both", paste0("SNPs and SVs (SNP/SV ratio: ", SNP_SV_ratio, ")"), causal_variant), "\n", sep = "")
if (causal_variant != "both") {
  cat("QTNs effect size: ", marker_effect, ifelse(add_effect_type == "geometric", " (geometric series)\n", "\n"), sep = "")
} else {
  cat("QTNs effect size: ", marker_effect, " (SNPs) and ", SV_effect, " (SVs)\n", sep = "")
}
cat("Genetic model: ", model, "\n", sep = "")
cat("Architecture: ", architecture, "\n", sep = "")
cat("Heritability: ", paste0(h2, collapse = ", "), "\n", sep = "")
cat("Number of populations: ", pops, "\n", sep = "")
cat("Number of reps: ", reps, "\n", sep = "")
cat("Seed number: ", seed, "\n", sep = "")
if (envs > 1) {
  # cat("Environments have different means: ", diff_env_mean, "\n", sep = "")
  if (!is.null(gen_cor_matrix)) {
    cat("Mean genetic correlation among environments: ", round(mean(gen_cor_matrix[upper.tri(gen_cor_matrix)]), digits = 2), "\n", sep = "")
  }
  cat("Mean residual correlation among environments: ", round(mean(res_cor_matrix[upper.tri(res_cor_matrix)]), digits = 2), "\n", sep = "")
}
cat("\n")

# get positional arguments
infile_name <- args[1]
sv_list <- args[2]
out_folder <- args[3]



#### simulate trait ----

# load genotypic data
geno_data <- fread(infile_name, header = TRUE, data.table = FALSE)

# load SV information
SVs <- fread(sv_list, header = FALSE, data.table = FALSE)
SVs <- SVs[, 1]

# check if directory already exists, create a new one if it doesn't
if (!dir.exists(out_folder)) dir.create(out_folder, recursive = TRUE)

for (pop_number in 1:pops) {

  # create folder for rep
  out_folder_pop <- paste0(out_folder, "/pop", pop_number)
  # get rep specific seed number
  seed_pop <- seed + pop_number

  cat ("population #", pop_number, "\n", sep = "")

  # simulate traits with different QTNs per experiment
  simulateTraits(geno_data = geno_data,
                 list_of_SV_IDs = SVs,
                 source_trait_variation = causal_variant,
                 SNP_SV_ratio = SNP_SV_ratio,
                 # experiments
                 pops = pops,
                 reps = reps,
                 # architecture
                 envs = envs,
                 h2 = h2,
                 model = model,
                 add_QTN_num = add_QTN_num,
                 add_effect_type = add_effect_type,
                 marker_effect = marker_effect,
                 SV_effect = SV_effect,
                 architecture = architecture,
                 diff_env_mean = diff_env_mean,
                 gxe = gxe,
                 diff_dist_gxe = diff_dist_gxe,
                 gen_cor_matrix = gen_cor_matrix,
                 res_cor_matrix = res_cor_matrix,
                 seed = seed_pop,
                 # ouput
                 QTN_variance = QTN_variance,
                 out_folder = out_folder_pop)

}




#### debug ----

# infile_name <- "data/usda_rils.all_markers.adjusted-n-markers.hmp.txt"
# sv_list <- "data/test_SVs_IDs.txt"
# # causal_variant <- "SNP"
# causal_variant <- "SV"
# # SNP_SV_ratio <- 0.5
# pops <- 10
# reps <- 3
# envs <- 5
# # envs <- 1
# h2 <- 0.2
# model <- "A"
# add_QTN_num <- 100
# add_effect_type <- "equal"
# # add_effect_type <- "geometric"
# marker_effect <- 0.1
# # SV_effect <- 0.2
# SV_effect <- NULL
# architecture <- "pleiotropic"
# seed <- 2020
# set.seed(seed); res_cor_matrix <- apply(randcorr(envs), MARGIN = c(1,2), function(x) as.numeric(as.character(x)))
#
# # gen_cor_matrix <- matrix(data = c(1, 0, 0, 0, 0,
# #                                   0, 1, 0, 0, 0,
# #                                   0, 0, 1, 0, 0,
# #                                   0, 0, 0, 1, 0,
# #                                   0, 0, 0, 0, 1),
# #                          nrow = 5, ncol = 5)
# # gen_cor_matrix <- matrix(data = c(1.0000000, 0.22824442, 0.4448667, 0.30036742, 0.5593773,
# #                                   0.2282444, 1.00000000, 0.2018070, 0.08436789, 0.4605805,
# #                                   0.4448667, 0.20180704, 1.0000000, 0.42500602, 0.2990179,
# #                                   0.3003674, 0.08436789, 0.4250060, 1.00000000, 0.4699402,
# #                                   0.5593773, 0.46058052, 0.2990179, 0.46994025, 1.0000000),
# #                          nrow = 5, ncol = 5)
# # gen_cor_matrix <- matrix(data = c(1.0000000, 0.8734700, 0.8893981, 0.5617946, 0.8260930,
# #                                   0.8734700, 1.0000000, 0.9116455, 0.6904279, 0.9076674,
# #                                   0.8893981, 0.9116455, 1.0000000, 0.6772922, 0.9060776,
# #                                   0.5617946, 0.6904279, 0.6772922, 1.0000000, 0.7586979,
# #                                   0.8260930, 0.9076674, 0.9060776, 0.7586979, 1.0000000),
# #                          nrow = 5, ncol = 5)
# # gen_cor_matrix <- matrix(data = c(1, 1, 1, 1, 1,
# #                                   1, 1, 1, 1, 1,
# #                                   1, 1, 1, 1, 1,
# #                                   1, 1, 1, 1, 1,
# #                                   1, 1, 1, 1, 1),
# #                          nrow = 5, ncol = 5)
#
# # mean_gen_cor <- round(mean(gen_cor_matrix[upper.tri(gen_cor_matrix)]), digits = 2)
#
# QTN_variance <- TRUE
#
#
# # out_folder <- paste0("analysis/trait_sim_mult-env/additive_model/test/correct_reps/",
# #                      add_QTN_num, "-QTNs_from_", causal_variant, "/",
# #                      h2, "-heritability/",
# #                      "effect0.1/",
# #                      "mean-gen-cor_", mean_gen_cor)
# # out_folder <- paste0("analysis/trait_sim_mult-env/additive_model/test/correct_reps/",
# #                      add_QTN_num, "-QTNs_from_", causal_variant, "/",
# #                      h2, "-heritability/",
# #                      "effect0.1-0.4-0.7/",
# #                      "mean-gen-cor_", mean_gen_cor)
#
#
# gxe <- TRUE
# gen_cor_matrix <- NULL
# diff_dist_gxe <- TRUE
#
# out_folder <- paste0("analysis/trait_sim_mult-env/additive_model/test/correct_reps_gxe/",
#                      add_QTN_num, "-QTNs_from_", causal_variant, "/",
#                      h2, "-heritability/",
#                      "effect0.1/",
#                      "with-gxe_new-rules/",
#                      "norm-dist_mean-0_sd-0.3")
#
# # out_folder <- paste0("analysis/trait_sim_mult-env/additive_model/test/correct_reps_gxe/",
# #                      add_QTN_num, "-QTNs_from_", causal_variant, "/",
# #                      h2, "-heritability/",
# #                      "effect0.1/",
# #                      "no-gxe/",
# #                      "norm-dist_mean-0_sd-0.3")
#
# # out_folder <- paste0("analysis/trait_sim/additive_model/test/correct_reps/",
# #                      add_QTN_num, "-QTNs_from_", causal_variant, "/",
# #                      h2, "-heritability/",
# #                      "effect0.1/")
