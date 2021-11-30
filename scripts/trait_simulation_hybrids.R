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
  --help                      show this helpful message
  --causal-variant=VALUE      marker type controlling trait variation ('SNP', 'SV' or 'both')
  --SNP-SV-ratio=VALUE        if both SNPs and SVs can be the causal variants, this option controls
                              the SNP/SV ratio to be sampled (default: 0.5)
  --pops=VALUE                number of populations to be simulated with different QTNs (default: 10)
  --reps=VALUE                number of reps (default: 3)
  --envs=VALUE                number of environments to simulate (default: 1)
  --h2=VALUE                  heritability (default: 0.5; comma-separated list of values also allowed)
  --impute-effect=VALUE       the marker effect ('Add' or 'Dom') when imputing missing data at the hapmap
                              numericalization step (default: 'Add')
  --impute-type=VALUE         the marker type ('Major', 'Middle' or 'Minor') when imputing missing data at 
                              the hapmap numericalization step (default: 'Major')
  --model=VALUE               additive, dominant or epistatic genetic model ('A', 'D' or 'E'), or any
                              combination of these models ('AE', 'DE' or 'ADE')
  --add-QTN-num=VALUE         number of additive loci controlling the trait (default: 10)
  --dom-QTN-num=VALUE         number of additive loci controlling the trait (default: 10)
  --effect-type=VALUE         type of effect among QTNs. If 'equal' is provided (default), effects across
                              all loci will be '--marker-effect'. If 'geometric' is provided, a geometric
                              series of effects starting at '--marker-effect' is applied to the loci
  --marker-effect-add=VALUE   size of additive marker effect (default: 0.1)
  --marker-effect-dom=VALUE   size of dominant marker effect (default: 0.1)
  --SV-effect-add=VALUE       if both SNPs and SVs can be the causal variants, this option allows
                              providing a different effect size for additive SVs (default: NULL; same as
                              '--marker-effect-add')
  --SV-effect-dom=VALUE       if both SNPs and SVs can be the causal variants, this option allows
                              providing a different effect size for dominant SVs (default: NULL; same as
                              '--marker-effect-dom')
  --architecture=VALUE        genetic architecture to be simulated ('pleiotropic' for traits being
                              controlled by the same QTNs, 'partially' for traits being controlled by
                              pleiotropic and trait-specific QTNs, 'LD' for traits being exclusively
                              controlled by different QTNs)
  --gen-cor-matrix=VALUE      if option is not provided, no genetic correlation matrix will be used.
                              If '/path/to/file', then the file should contain a genetic correlation
                              matrix among environments to be simulated (no column or row names). The
                              lower the genetic correlation among environments, higher the chance of
                              having GxE interactions. If '1', a perfect positive correlation matrix
                              (i.e. all environments have correlation = 1) will be generated with
                              'number of rows' = 'number of envs'
  --res-cor-matrix=VALUE      if option is not provided, a random residual correlation matrix will be
                              generated with 'number of rows' = 'number of envs'.  If '/path/to/file',
                              then the file should contain a residual correlation matrix among
                              environments to be simulated (no column or row names)
  --gxe                       add this option to simulate GxE by adding random effects (drawn from a
                              normal distribution) to the marker effects in each environment. This
                              option also overrides any genetic correlation matrix provided
  --diff-dist-gxe             add this option to draw random GxE effects from different normal
                              distributions depending on marker type (SNP or SV). The normal distribution
                              used to sample effects will have mean 0 for both SNPs and SVs, but the
                              standard deviation for SVs will be multiplied by the value of '--SV-effect'
                              option.
  --diff-env-mean             add this option to force simulated traits to have different means in
                              each environment
  --QTN-variance              add this option to write files with percent variance explained per QTN
  --seed=VALUE                value for set.seed (default: NULL; random number is selected)

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
                           model = c("A", "D", "E", "AD", "AE", "DE", "ADE"),
                           add_QTN_num = 3,
                           dom_QTN_num = 3,
                           effect_type = c("geometric", "equal"),
                           marker_effect_add = 0.1,
                           marker_effect_dom = 0.1,
                           SV_effect_add = NULL,
                           SV_effect_dom = NULL,
                           architecture = "pleiotropic",
                           diff_env_mean = FALSE,
                           gxe = FALSE,
                           diff_dist_gxe = FALSE,
                           gen_cor_matrix = NULL,
                           res_cor_matrix = NULL,
                           seed = 2020,
                           # ouput
                           QTN_variance = QTN_variance,
                           out_folder = getwd(),
                           # numericalization
                           SNP_effect = "Add",
                           SNP_impute = "Major") {

  # make sure list of SVs is provided
  if (is.null(list_of_SV_IDs)) stop("No list of SV IDs provided!")

  # separate data into SNPs and SVs
  geno_data_SVs <- geno_data[which(geno_data[, 1] %in% list_of_SV_IDs), ]
  geno_data_SNPs <- geno_data[which(!geno_data[, 1] %in% list_of_SV_IDs), ]

  # filter data based on sorce of which marker type is the causative variant
  if (source_trait_variation == "SV") geno2trait_sim <- geno_data_SVs
  if (source_trait_variation == "SNP") geno2trait_sim <- geno_data_SNPs
  if (source_trait_variation == "both") {

    splitCausalVars <- function(geno_data_SVs, geno_data_SNPs, QTN_num, SNP_SV_ratio, seed) {

      # get how many additive SNPs and SVs will be selected
      n_SNPs_to_sample <- ceiling(QTN_num * SNP_SV_ratio)
      n_SVs_to_sample <- ceiling(QTN_num * (1 - SNP_SV_ratio))
      # adjust numbers if they don't sum to total number QTNs (due to rounding above)
      if (n_SNPs_to_sample + n_SVs_to_sample != QTN_num) {
        n_markers_to_fix <- QTN_num - (n_SNPs_to_sample + n_SVs_to_sample)
        n_SVs_to_sample <- n_SVs_to_sample + n_markers_to_fix
      }
      # sample dataset
      set.seed(seed)
      SVs_to_keep <- sort(sample(1:NROW(geno_data_SVs), size = n_SVs_to_sample, replace = FALSE))
      set.seed(seed)
      SNPs_to_keep <- sort(sample(1:NROW(geno_data_SNPs), size = n_SNPs_to_sample, replace = FALSE))
      # filter datasets for each marker type separately
      geno_data_SVs <- geno_data_SVs[SVs_to_keep, ]
      geno_data_SNPs <- geno_data_SNPs[SNPs_to_keep, ]

      return(list(SVs = geno_data_SVs, SNPs = geno_data_SNPs))

    }

    # get additive SNPs and SVs
    geno_data_add <- splitCausalVars(geno_data_SVs, geno_data_SNPs, add_QTN_num, SNP_SV_ratio, seed)
    causal_vars_add <- c(geno_data_add$SVs[, 1], geno_data_add$SNPs[, 1])
    # remove additive QTNs from dataset
    geno_data_SVs <- geno_data_SVs[which(!geno_data_SVs[, 1] %in% geno_data_add$SVs[, 1]), ]
    geno_data_SNPs <- geno_data_SNPs[which(!geno_data_SNPs[, 1] %in% geno_data_add$SNPs[, 1]), ]
    # get dominant SNPs and SVs
    geno_data_dom <- splitCausalVars(geno_data_SVs, geno_data_SNPs, dom_QTN_num, SNP_SV_ratio, seed + seed)
    causal_vars_dom <- c(geno_data_dom$SVs[, 1], geno_data_dom$SNPs[, 1])

    # combine and order dataset
    geno2trait_sim <- rbind(do.call(rbind, geno_data_add), do.call(rbind, geno_data_dom))
    geno2trait_sim <- geno2trait_sim[order(geno2trait_sim$chr, geno2trait_sim$pos), ]

  }

  # select QTNs
  set.seed(seed)
  marker_qtns <- sort(sample(1:NROW(geno2trait_sim), size = add_QTN_num + dom_QTN_num, replace = FALSE))
  geno2trait_sim <- geno2trait_sim[marker_qtns, ]

  # check which type of effect will be used
  if (effect_type == "geometric") {

    # create a geometric series where first QTN effect is 'marker_effect_add"
    add_effect_value <- createGeomSeries(qtn_number = add_QTN_num, largest_effect = marker_effect_add)
    dom_effect_value <- createGeomSeries(qtn_number = dom_QTN_num, largest_effect = marker_effect_dom)

  }

  if (effect_type == "equal") {

    # make sure all QTNs have the same 'marker_effect_add' value
    add_effect_value <- rep(marker_effect_add, times = add_QTN_num)
    dom_effect_value <- rep(marker_effect_dom, times = dom_QTN_num)
    # transform to list for compatibility to simplePHENOTYPES
    add_effect_value <- list(add_effect_value)
    dom_effect_value <- list(dom_effect_value)

    # adjust additive effect size if SVs have different effects than SNPs
    if (!is.null(SV_effect_add)) {
      # find which qtns are SVs
      if (source_trait_variation == "both") {
        sv_qtns_add <- which(causal_vars_add %in% list_of_SV_IDs)
      } else {
        sv_qtns_add <- which(geno2trait_sim[, 1] %in% list_of_SV_IDs)
      }
      # modify effect size of svs
      if (length(sv_qtns_add) > 0) add_effect_value[[1]][sv_qtns_add] <- SV_effect_add
    }

    # adjust dominant effect size if SVs have different effects than SNPs
    if (!is.null(SV_effect_dom)) {
      # find which qtns are SVs
      if (source_trait_variation == "both") {
        sv_qtns_dom <- which(causal_vars_dom %in% list_of_SV_IDs)
      } else {
        sv_qtns_dom <- which(geno2trait_sim[, 1] %in% list_of_SV_IDs)
      }
      # modify effect size of svs
      if (length(sv_qtns_dom) > 0) dom_effect_value[[1]][sv_qtns_dom] <- SV_effect_dom
    }

  }

  # adjust some parameters if multiple traits
  if (envs > 1) {

    # make sure architecture option is provided if number of traits is bigger than 1
    if(is.null(architecture)) stop("Provide genetic architecture (pleiotropic, partially, or LD) when more than 1 trait is simulated")
    # set the same additive effects for all traits
    add_effect_value <- rep(add_effect_value, envs)
    dom_effect_value <- rep(dom_effect_value, envs)

    if (diff_env_mean) {

      # scale effect sizes among environments so they have different means
      for (env in 1:envs) {
        # use the absolute mean correlation among environments as the scaling factor
        mean_cor_env <- abs(mean(res_cor_matrix[env, -env]))
        add_effect_value[[env]] <- round(add_effect_value[[env]] * mean_cor_env, digits = 4)
        dom_effect_value[[env]] <- round(dom_effect_value[[env]] * mean_cor_env, digits = 4)
      }

    }

    # set the same heritability for all traits
    h2 <- rep(h2, envs)

    # if GxE, change effect sizes for each environment
    if (gxe) {

      addGxEeffects <- function(envs, QTN_num, effect_value, marker_effect, sv_qtns,
                                SV_effect, diff_dist_gxe = FALSE, seed) {

        # select QTNs to have constant effects across envs -- 10% QTNs without GxE
        n_constant_qtns <- ceiling(0.1 * QTN_num)
        constant_seed <- seed + seed
        set.seed(constant_seed)
        constant_qtns <- sort(sample(1:QTN_num, size = n_constant_qtns, replace = FALSE))

        # select QTNs that can change signs across environments -- 10% QTNs can change signs
        n_diff_sign_qtns <- ceiling(0.1 * QTN_num)
        diff_sign_seed <- seed * envs
        set.seed(diff_sign_seed)
        diff_sign_qtns <- sort(sample((1:QTN_num)[-constant_qtns], size = n_diff_sign_qtns, replace = FALSE))

        for (env in 1:length(effect_value)) {

          # select QTNs to have no effect in this env -- 10% QTNs with no effect at all
          n_zero_qtns <- ceiling(0.1 * QTN_num)
          zero_seed <- (seed + seed) * env
          set.seed(zero_seed)
          zero_qtns <- sort(sample((1:QTN_num)[-constant_qtns], size = n_zero_qtns, replace = FALSE))

          # draw specific effects from a normal distribution
          gxe_seed <- seed + env
          set.seed(gxe_seed)
          gxe_effect <- rnorm(n = QTN_num, mean = 0, sd = 0.3)
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
          effect_value[[env]][zero_qtns] <- 0
          effect_value[[env]][-c(zero_qtns, constant_qtns)] <- effect_value[[env]][-c(zero_qtns, constant_qtns)] + gxe_effect
          effect_value[[env]] <- sapply(effect_value[[env]], function(i) {
            if (i > 1) i = 1
            if (i < -1) i = -1
            return(i)
          })

        }

        # make sure only selected QTNs are allowed to change signs
        for (env in 1:length(effect_value)) {

          effect_value[[env]][-diff_sign_qtns] <- sapply(effect_value[[env]][-diff_sign_qtns], function(qtn, effect) {
            # change signs according to initial 'marker_effect' provided
            if (effect > 0) if (qtn < 0) qtn = qtn * (-1)
            if (effect < 0) if (qtn > 0) qtn = qtn * (-1)
            return(qtn)
          }, effect = marker_effect)

        }

        return(effect_value)

      }

      if (source_trait_variation != "both") {
        sv_qtns_add <- NULL
        sv_qtns_dom <- NULL
      }

      add_effect_value <- addGxEeffects(envs, add_QTN_num, add_effect_value, marker_effect_add, sv_qtns_add,
                                        SV_effect_add, diff_dist_gxe = diff_dist_gxe, seed)
      dom_effect_value <- addGxEeffects(envs, dom_QTN_num, dom_effect_value, marker_effect_dom, sv_qtns_dom,
                                        SV_effect_dom, diff_dist_gxe = diff_dist_gxe, seed + seed)

      # make sure to disable the use of any genetic correlation among environments
      gen_cor_matrix <- NULL

    }

  }

  cat("--- Simulating traits with ", add_QTN_num, " add QTNs, ", dom_QTN_num, " dom QTNs, and ",
      h2[1], " heritability (seed number: ", seed, ") ---\n", sep = "")

  if (source_trait_variation == "both") {
    QTN_list <- list(add = list(causal_vars_add), dom = list(causal_vars_dom))
  } else {
    QTN_list <- NULL
  }

  # simulate trait
  create_phenotypes(geno_obj = geno2trait_sim,
                    QTN_list = QTN_list,
                    rep = reps,
                    # architecture
                    ntraits = envs,
                    h2 = h2,
                    model = model,
                    add_QTN_num = add_QTN_num,
                    add_effect = add_effect_value,
                    dom_QTN_num = dom_QTN_num,
                    dom_effect = dom_effect_value,
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
                    SNP_effect = SNP_effect,
                    SNP_impute = SNP_impute)

}



#### command line options ----

# set default
causal_variant <- "SNP"
SNP_SV_ratio <- NULL
pops <- "10"
reps <- "3"
envs <- "1"
h2 <- "0.5"
model <- "AD"
impute_effect <- "Add"
impute_type <- "Major"
add_QTN_num <- "10"
dom_QTN_num <- "10"
effect_type <- "equal"
marker_effect_add <- "0.1"
marker_effect_dom <- "0.1"
SV_effect_add <- NULL
SV_effect_dom <- NULL
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
                        "--impute-effect", "--impute-type", "--model", "--add-QTN-num", "--dom-QTN-num",
                        "--effect-type", "--marker-effect-add", "--marker-effect-dom", "--SV-effect-add",
                        "--SV-effect-dom", "--architecture", "--gen-cor-matrix", "--res-cor-matrix",
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

if (!model %in% c("A", "D", "E", "AD", "AE", "DE", "ADE")) {
  stop("Optional argument '--model' should be 'A', 'D', 'E', 'AE', 'DE' or 'ADE'")
}

if (suppressWarnings(!any(is.na(as.integer(add_QTN_num))))) {
  add_QTN_num <- as.integer(add_QTN_num)
} else {
  stop("Optional argument '--add-QTN-num' should be an integer")
}

if (suppressWarnings(!any(is.na(as.integer(dom_QTN_num))))) {
  dom_QTN_num <- as.integer(dom_QTN_num)
} else {
  stop("Optional argument '--dom-QTN-num' should be an integer")
}

if (!effect_type %in% c("equal", "geometric")) {
  stop("Optional argument '--effect-type' should be a either 'equal' or 'geometric'")
}

if (suppressWarnings(!is.na(as.numeric(marker_effect_add)))) {
  marker_effect_add <- as.numeric(marker_effect_add)
} else {
  stop("Optional argument '--marker-effect-add' should be an integer")
}

if (suppressWarnings(!is.na(as.numeric(marker_effect_dom)))) {
  marker_effect_dom <- as.numeric(marker_effect_dom)
} else {
  stop("Optional argument '--marker-effect-dom' should be an integer")
}

if (!is.null(SV_effect_add)) {
  if (suppressWarnings(!is.na(as.numeric(SV_effect_add)))) {
    SV_effect_add <- as.numeric(SV_effect_add)
  } else {
    stop("Optional argument '--SV-effect-add' should be a number")
  }
}

if (!is.null(SV_effect_dom)) {
  if (suppressWarnings(!is.na(as.numeric(SV_effect_dom)))) {
    SV_effect_dom <- as.numeric(SV_effect_dom)
  } else {
    stop("Optional argument '--SV-effect-dom' should be a number")
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
cat("Number of additive QTNs per environment: ", add_QTN_num, "\n", sep = "")
cat("Number of dominant QTNs per environment: ", dom_QTN_num, "\n", sep = "")
cat("Causal variant: ", ifelse(causal_variant == "both", paste0("SNPs and SVs (SNP/SV ratio: ", SNP_SV_ratio, ")"), causal_variant), "\n", sep = "")
if (causal_variant != "both") {
  cat("Additive QTNs effect size: ", marker_effect_add, ifelse(effect_type == "geometric", " (geometric series)\n", "\n"), sep = "")
  cat("Dominant QTNs effect size: ", marker_effect_dom, ifelse(effect_type == "geometric", " (geometric series)\n", "\n"), sep = "")
} else {
  cat("Additive QTNs effect size: ", marker_effect_add, " (SNPs) and ", SV_effect_add, " (SVs)\n", sep = "")
  cat("Dominant QTNs effect size: ", marker_effect_dom, " (SNPs) and ", SV_effect_dom, " (SVs)\n", sep = "")
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
                 dom_QTN_num = dom_QTN_num,
                 effect_type = effect_type,
                 marker_effect_add = marker_effect_add,
                 marker_effect_dom = marker_effect_dom,
                 SV_effect_add = SV_effect_add,
                 SV_effect_dom = SV_effect_dom,
                 architecture = architecture,
                 diff_env_mean = diff_env_mean,
                 gxe = gxe,
                 diff_dist_gxe = diff_dist_gxe,
                 gen_cor_matrix = gen_cor_matrix,
                 res_cor_matrix = res_cor_matrix,
                 seed = seed_pop,
                 # ouput
                 QTN_variance = QTN_variance,
                 out_folder = out_folder_pop,
                 # numericalization
                 SNP_effect = impute_effect,
                 SNP_impute = impute_type)
  
}




#### debug ----

# infile_name <- "data/usda_hybrids_projected-SVs-SNPs.chr10.poly.low-missing.hmp.txt"
# sv_list <- "data/SVs_IDs_poly.txt"
# # causal_variant <- "SNP"
# causal_variant <- "both"
# SNP_SV_ratio <- 0.5
# impute_effect <- "Dom"
# impute_type <- "Middle"
# pops <- "2"
# reps <- "3"
# envs <- "5"
# # envs <- 1
# h2 <- "0.3"
# effect_type <- "equal"
# # effect_type <- "geometric"
# marker_effect_add <- "0.1"
# marker_effect_dom <- "0.1"
# # SV_effect_add <- 0.1
# # SV_effect_add <- 0.5
# # SV_effect_add <- NULL
# # SV_effect_dom <- 0.1
# # SV_effect_dom <- 0.5
# # SV_effect_dom <- NULL
# architecture <- "pleiotropic"
# 
# {
#   model <- "AD"
#   add_QTN_num <- 5
#   dom_QTN_num <- 5
#   
#   out_folder <- "analysis/trait_sim_hybrids/multi_env/with_gxe/AD_model/5-add-QTNs_5-dom-QTNs_from_both/0.3-heritability"
#   
# }
# 
# # {
# #   model <- "A"
# #   add_QTN_num <- "10"
# #   dom_QTN_num <- "0"
# #   
# #   out_folder <- "analysis/trait_sim_hybrids/multi_env/with_gxe/A_model/10-add-QTNs_0-dom-QTNs_from_SNP/0.3-heritability"
# # }
# 
# # out_folder <- paste0("tests/trait_sim_hybrids/",
# #                      add_QTN_num, "-add-QTNs_from_", causal_variant, "/",
# #                      dom_QTN_num, "-dom-QTNs_from_", causal_variant, "/",
# #                      h2, "-heritability/")
# 
# # out_folder <- paste0("tests/trait_sim_hybrids/",
# #                      add_QTN_num, "-add-QTNs_from_both/",
# #                      dom_QTN_num, "-dom-QTNs_from_both/",
# #                      "SNP-SV-ratio_", SNP_SV_ratio, "/effects_SNP-",
# #                      marker_effect_add, "_SV-", SV_effect_add, "/",
# #                      h2, "-heritability/")
# 
# 
# seed <- "75486"
# # set.seed(seed); res_cor_matrix <- apply(randcorr(envs), MARGIN = c(1,2), function(x) as.numeric(as.character(x)))
# res_cor_matrix <- "data/usda_envs_cor_matrix.txt"
# 
# QTN_variance <- TRUE
# gxe <- TRUE
# gen_cor_matrix <- NULL
# diff_dist_gxe <- FALSE
# # diff_dist_gxe <- TRUE
