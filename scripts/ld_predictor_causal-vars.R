library(data.table)
library(doParallel)

usage <- function() {
  cat("
description: get LD between predictors and causal variants.

usage: Rscript ld_predictor_causal-vars.R [predictors_file] [causal_vars_file] [ld_file] [outfile]

positional arguments:
  predictors_file       hapmap file with markers used as predictors in genomic prediction models
  causal_vars_file      file with markers selected by simplePHENOTYPES to be the causal variants
                        of a simulated trait
  ld_file               file with LD information for all markers of this population (should provide
                        only one chromosome file, but the other chromosome files should be in the
                        same folder)
  outfile               name of output file

optional argument:
  --help                show this helpful message

"
  )
}



#### command line options ----

args <- commandArgs(trailingOnly = TRUE)
if ("--help" %in% args) usage() & q(save = "no")

# assert to have the correct optional arguments
if (length(args) != 4) stop(usage(), "missing positional argument(s)")

# get positional arguments
predictors_file <- args[1]
causal_vars_file <- args[2]
ld_file <- args[3]
outfile <- args[4]



#### get ld between predictors and causal variants ----

# load predictors
pred <- fread(predictors_file, header = TRUE, data.table = FALSE)
pred <- pred[, c(1, 3, 4)]
colnames(pred) <- c("marker", "chrom", "pos")

# load causal variants
causal_vars <- fread(causal_vars_file, header = TRUE, data.table = FALSE)
causal_vars <- causal_vars[, c("snp", "chr", "pos")]
colnames(causal_vars) <- c("marker", "chrom", "pos")

# only look at chromosomes that have at least one predictor and one causal var
chrs_to_analyze <- sort(intersect(unique(pred[, "chrom"]), unique(causal_vars[, "chrom"])))

registerDoParallel(cores = length(chrs_to_analyze))
ld_results <- foreach(chr = chrs_to_analyze, .combine = rbind) %dopar% {

  # subset by chr
  pred_chr <- subset(pred, chrom == chr)
  causal_vars_chr <- subset(causal_vars, chrom == chr)

  # change ld file name for respective chr
  ld_file <- gsub("chr[0-9]+", paste0("chr", chr), ld_file, perl = TRUE)
  # load ld file
  ld <- fread(ld_file, header = TRUE, data.table = FALSE)

  # create empty df to store results of chromosome
  ld_results_chr <- data.frame()

  # get ld between predictors and causal variants
  for (i in 1:NROW(pred_chr)) {

    # get predictor info
    marker_predictor <- pred_chr[i, "marker"]
    pos_predictor <- pred_chr[i, "pos"]
    # keep only predictors in ld file
    ld_pred <- subset(ld, SNP_A == marker_predictor | SNP_B == marker_predictor)

    for (j in 1:NROW(causal_vars_chr)) {

      # get causal variants info
      marker_causal <- causal_vars_chr[j, "marker"]
      pos_causal <- causal_vars_chr[j, "pos"]
      # keep only predictor and causal variant in the ld file
      ld_pred_causal <- subset(ld_pred, SNP_A == marker_causal | SNP_B == marker_causal)
      # get their ld
      r2 <- ifelse(test = NROW(ld_pred_causal) > 0, yes = ld_pred_causal$R2, no = NA)

      # get pred and causal var MAFs
      if (NROW(ld_pred_causal) > 0) {
        maf_pred <- ifelse(test = ld_pred_causal$SNP_A == marker_predictor,
                           yes = ld_pred_causal$MAF_A, no = ld_pred_causal$MAF_B)
        maf_qtl <- ifelse(test = ld_pred_causal$SNP_A == marker_causal,
                          yes = ld_pred_causal$MAF_A, no = ld_pred_causal$MAF_B)
      } else {
        maf_pred <- NA
        maf_qtl <- NA
      }

      # append results to df
      ld_results_chr <- rbind(ld_results_chr,
                              data.frame(chr = chr, predictor = marker_predictor,
                                         causal_var = marker_causal, r2 = r2,
                                         maf_pred = maf_pred, maf_qtl = maf_qtl,
                                         bp_distance = pos_predictor - pos_causal))

    }
  }

  return(ld_results_chr)

}
stopImplicitCluster()

# make sure df is ordered by chromosome
ld_results <- ld_results[order(ld_results$chr), ]

# write results
fwrite(ld_results, outfile, sep = "\t", quote = FALSE, na = NA, row.names = FALSE)



#### debug ----

# predictors_file <- "analysis/trait_sim/datasets/iter17/usda_rils.all_markers.adjusted-n-markers.hmp.txt"
# causal_vars_file <- "analysis/trait_sim/multi_env/with_gxe/additive_model/equal_effects/10-QTNs_from_SNP/0.3-heritability/pop2/Additive_QTNs.txt"
# ld_file <- "analysis/trait_sim/multi_env/with_gxe/additive_model/equal_effects/10-QTNs_from_SNP/ld_causative-vars_predictors/pop2/all_markers.pred-iter17.causal-pop2.chr10.ld.gz"
# outfile <- "analysis/trait_sim/multi_env/with_gxe/additive_model/equal_effects/10-QTNs_from_SNP/ld_causative-vars_predictors/pop2/ld_summary.all_markers.pred-iter17.causal-pop2.txt"
