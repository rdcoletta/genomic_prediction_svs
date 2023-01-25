library(data.table)
library(tidyr)
library(dplyr)
library(ggplot2)
library(lme4)
library(lmerTest)
library(MuMIn)

usage <- function() {
  cat("
description: simulate trait based on user-defined genetic architecture.

usage: Rscript anova_sim_traits.R [dir_sim_traits] [...]

positional arguments:
  dir_sim_traits          path to folder containing simulated traits from simplePHENOTYPES

optional argument:
  --help                  show this helpful message

"
  )
}




#### command line options ----

args <- commandArgs(trailingOnly = TRUE)
if ("--help" %in% args) usage() & q(save = "no")

# assert to have the correct optional arguments
if (length(args) != 1) stop(usage(), "missing positional argument(s)")

# get positional arguments
sim_folder <- args[1]
# sim_folder <- "analysis/trait_sim_mult-env/additive_model/test/script_test/100-QTNs_from_SNP/0.2-heritability/marker-effect0.1/mean-gen-cor_1/pop1"
# sim_folder <- "analysis/trait_sim_mult-env/additive_model/test/correct_reps_gxe/100-QTNs_from_SNP/0.2-heritability/effect0.1/no-gxe/norm-dist_mean-0_sd-0.3/pop1"



#### QC simulated tratis ----

# get file names
file_pheno <- list.files(sim_folder, pattern = "Simulated_Data", full.names = TRUE)
# file_qtns <- list.files(sim_folder, pattern = "geno_info", full.names = TRUE)
# file_log <- list.files(sim_folder, pattern = "Log", full.names = TRUE)

# open files
pheno_pop <- fread(file_pheno, header = TRUE, data.table = FALSE)
# info_pop <- fread(file_qtns, header = TRUE, data.table = FALSE)
# log_pop <- scan(file_log, what = "character", sep = "\n")

# adjust column names
colnames(pheno_pop)[1] <- "genotype"
colnames(pheno_pop)[NCOL(pheno_pop)] <- "rep"
colnames(pheno_pop)[2:(NCOL(pheno_pop) - 1)] <- paste0("env", 1:length(2:(NCOL(pheno_pop) - 1)))
# transform df into long format
pheno_pop <- pivot_longer(pheno_pop, !c(genotype, rep), names_to = "environment", values_to = "trait_value") %>%
  mutate(rep = paste0("rep", rep)) %>% 
  relocate(rep, .after = environment) %>%
  arrange(environment, genotype)

# write summary
out_summary <- paste0(sim_folder, "/simulated_traits_summary.txt")
fwrite(x = pheno_pop, file = out_summary, quote = FALSE, sep = "\t", na = NA, row.names = FALSE)

# plot summary
plot_summary <- ggplot(pheno_pop, aes(x = environment, y = trait_value)) +
  geom_violin(color = "gray70", fill = "gray70") +
  stat_summary(fun.data = "mean_cl_normal", color = "black", show.legend = FALSE) 

out_plot_summary <- paste0(sim_folder, "/simulated_traits_distribution-per-env.pdf")
ggsave(filename = out_plot_summary, plot = plot_summary, device = "pdf")

# summarize individual performance per env
geno_mean_per_env <- pheno_pop %>% 
  group_by(genotype, environment) %>% 
  summarize(mean = mean(trait_value))

# reorder genotypes based on mean of env1
order_geno_env1 <- geno_mean_per_env[which(geno_mean_per_env$environment == "env1"), ]
order_geno_env1 <- order_geno_env1[order(order_geno_env1$mean), ]
order_geno_env1 <- pull(order_geno_env1, genotype)
geno_mean_per_env$genotype <- factor(geno_mean_per_env$genotype, levels = order_geno_env1)

# plot individual performance per env
plot_ind_perf <- ggplot(geno_mean_per_env, aes(x = genotype, y = mean)) +
  facet_grid(~environment) +
  geom_bar(stat = "identity") +
  theme(axis.text.y = element_text(size = 5)) +
  coord_flip()

out_plot_ind_perf <- paste0(sim_folder, "/simulated_traits_ind-performance-per-env.pdf")
ggsave(filename = out_plot_ind_perf, plot = plot_ind_perf, device = "pdf", height = 20)

# fit mixed linear model
# fixed effects -- genotypes
# random effects -- reps within envs, envs, GxE)
trait_lmer <- lmer(trait_value ~ genotype + (1 | environment/rep) + (1 | genotype:environment),
                   data = pheno_pop, REML = TRUE)
trait_lmer_summary <- summary(trait_lmer)

# calculate percent variation explained
pve_fixed <- as.numeric(r.squaredGLMM(trait_lmer)[1, "R2m"])
pve_random <- as.numeric(r.squaredGLMM(trait_lmer)[1, "R2c"] - pve_fixed)
pve_residual <- 1 - pve_fixed - pve_random

# extract variation from random effects
var_random <- data.frame(trait_lmer_summary$varcor) %>% 
  filter(grp != "Residual") %>% 
  select(source = grp, variance = vcov) %>% 
  mutate(pve = (variance / sum(variance)) * pve_random) %>% 
  select(source, pve)

# run anova on fixed and random effects separately
trait_lmer_anova_fixed <- anova(trait_lmer)
trait_lmer_anova_random <- ranova(trait_lmer)

# format anova results and add percent variance explained
pval_fixed <- cbind(source = rownames(trait_lmer_anova_fixed),
                    data.frame(trait_lmer_anova_fixed, row.names = NULL, check.names = FALSE)) %>% 
  select(source, pval = `Pr(>F)`) %>% 
  mutate(pve = pve_fixed)

pval_random <- cbind(source = rownames(trait_lmer_anova_random),
                     data.frame(trait_lmer_anova_random, row.names = NULL, check.names = FALSE)) %>% 
  select(source, pval = `Pr(>Chisq)`) %>% 
  mutate(source = gsub("(1 | ", "", source, fixed = TRUE),
         source = gsub(")", "", source, fixed = TRUE)) %>% 
  filter(source != "<none>") %>% 
  full_join(var_random, by = "source") 

mm_pve_pval <- rbind(pval_fixed, pval_random) %>% 
  mutate(pve = round(pve, digits = 3),
         signif = case_when(pval > 0.1 ~ "NS",
                            pval > 0.05 & pval <= 0.1 ~ ".",
                            pval > 0.01 & pval <= 0.05 ~ "*",
                            pval <= 0.01 ~ "***")) %>% 
  select(source, pve, pval, signif)

# add residual info into anova results
mm_pve_pval <- rbind(mm_pve_pval, data.frame(source = "residual",
                                             pve = round(pve_residual, digits = 3),
                                             pval = NA,
                                             signif = NA))

# write anova results
out_pve <- paste0(sim_folder, "/simulated_traits_pve.txt")
fwrite(x = mm_pve_pval, file = out_pve, quote = FALSE, sep = "\t", na = NA, row.names = FALSE)

# also print full anova results
options(max.print=2000)
sink(paste0(sim_folder, "/simulated_traits_anova.txt"))
print(trait_lmer_summary)
cat("\n---- ANOVA fixed effects ----\n")
print(trait_lmer_anova_fixed)
cat("\n---- ANOVA random effects ----\n")
print(trait_lmer_anova_random)
sink()


# More on calculating PVE for mixed models:
#   https://stats.stackexchange.com/questions/7240/proportion-of-explained-variance-in-a-mixed-effects-model
#   https://ecologyforacrowdedplanet.wordpress.com/2013/08/27/r-squared-in-mixed-models-the-easy-way/
#   https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/j.2041-210x.2012.00261.x
#   Marginal R2 is concerned with variance explained by fixed factors
#   Conditional R2 is concerned with variance explained by both fixed and random factors. 
