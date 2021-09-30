library(data.table)
library(ggplot2)
library(dplyr)

folder_base <- "analysis/trait_sim/multi_env"

gxe_results <- data.frame(gxe = as.character(),
                          h2 = as.numeric(),
                          qtn = as.numeric(),
                          var = as.character(),
                          ratio = as.numeric(),
                          sv_effect = as.numeric(),
                          diff_dist = as.logical(),
                          pop = as.numeric(),
                          gxe_pve = as.numeric(),
                          gxe_pval = as.numeric())



for (gxe in c("no", "with")) {
  for (qtn in c(10, 100)) {
    for (var in c("SNP", "SV", "both")) {
      for (h2 in c(0.3, 0.7)) {
        for (pop in 1:20) {

          # begin var both
          if (var == "both") {

            for (ratio in c(0.5, 0.8)) {
              for (sv_effect in c(0.1, 0.2, 0.5)) {

                if (gxe == "with") {

                  if (sv_effect %in% c(0.2, 0.5)) {
                    for (diff_dist in c(FALSE, TRUE)) {
                      if (!diff_dist) {
                        # with gxe, same distribution
                        folder_traits <- paste0(folder_base, "/", gxe, "_gxe/additive_model/equal_effects/", qtn,
                                                "-QTNs_from_", var, "/SNP-SV-ratio_", ratio, "/effects_SNP-0.1_SV-",
                                                sv_effect, "/", h2, "-heritability/pop", pop)

                        trait_pve <- fread(paste0(folder_traits, "/simulated_traits_pve.txt"),
                                           header = TRUE, data.table = FALSE)

                        pve <- trait_pve[trait_pve$source == "genotype:environment", "pve"]
                        pval <- trait_pve[trait_pve$source == "genotype:environment", "pval"]

                        gxe_results <- rbind(gxe_results,
                                             data.frame(gxe = gxe,
                                                        h2 = h2,
                                                        qtn = qtn,
                                                        var = var,
                                                        ratio = ratio,
                                                        sv_effect = sv_effect,
                                                        diff_dist = diff_dist,
                                                        pop = pop,
                                                        gxe_pve = pve,
                                                        gxe_pval = pval))

                      } else {
                        # with gxe, different distribution
                        folder_traits <- paste0(folder_base, "/", gxe, "_gxe/additive_model/equal_effects/", qtn,
                                                "-QTNs_from_", var, "/SNP-SV-ratio_", ratio, "/effects_SNP-0.1_SV-",
                                                sv_effect, "_diff-dist-gxe/", h2, "-heritability/pop", pop)

                        trait_pve <- fread(paste0(folder_traits, "/simulated_traits_pve.txt"),
                                           header = TRUE, data.table = FALSE)

                        pve <- trait_pve[trait_pve$source == "genotype:environment", "pve"]
                        pval <- trait_pve[trait_pve$source == "genotype:environment", "pval"]

                        gxe_results <- rbind(gxe_results,
                                             data.frame(gxe = gxe,
                                                        h2 = h2,
                                                        qtn = qtn,
                                                        var = var,
                                                        ratio = ratio,
                                                        sv_effect = sv_effect,
                                                        diff_dist = diff_dist,
                                                        pop = pop,
                                                        gxe_pve = pve,
                                                        gxe_pval = pval))
                      }
                    }
                  } else {
                    # with gxe, same distribution
                    folder_traits <- paste0(folder_base, "/", gxe, "_gxe/additive_model/equal_effects/", qtn,
                                            "-QTNs_from_", var, "/SNP-SV-ratio_", ratio, "/effects_SNP-0.1_SV-",
                                            sv_effect, "/", h2, "-heritability/pop", pop)

                    trait_pve <- fread(paste0(folder_traits, "/simulated_traits_pve.txt"),
                                       header = TRUE, data.table = FALSE)

                    pve <- trait_pve[trait_pve$source == "genotype:environment", "pve"]
                    pval <- trait_pve[trait_pve$source == "genotype:environment", "pval"]

                    gxe_results <- rbind(gxe_results,
                                         data.frame(gxe = gxe,
                                                    h2 = h2,
                                                    qtn = qtn,
                                                    var = var,
                                                    ratio = ratio,
                                                    sv_effect = sv_effect,
                                                    diff_dist = FALSE,
                                                    pop = pop,
                                                    gxe_pve = pve,
                                                    gxe_pval = pval))
                  }

                } else {
                  # no gxe
                  folder_traits <- paste0(folder_base, "/", gxe, "_gxe/additive_model/equal_effects/", qtn,
                                          "-QTNs_from_", var, "/SNP-SV-ratio_", ratio, "/effects_SNP-0.1_SV-",
                                          sv_effect, "/", h2, "-heritability/pop", pop)

                  trait_pve <- fread(paste0(folder_traits, "/simulated_traits_pve.txt"),
                                     header = TRUE, data.table = FALSE)

                  pve <- trait_pve[trait_pve$source == "genotype:environment", "pve"]
                  pval <- trait_pve[trait_pve$source == "genotype:environment", "pval"]

                  gxe_results <- rbind(gxe_results,
                                       data.frame(gxe = gxe,
                                                  h2 = h2,
                                                  qtn = qtn,
                                                  var = var,
                                                  ratio = ratio,
                                                  sv_effect = sv_effect,
                                                  diff_dist = FALSE,
                                                  pop = pop,
                                                  gxe_pve = pve,
                                                  gxe_pval = pval))
                }

              }
            }
            # end of var both
          } else {
            # begin var snp or sv
            folder_traits <- paste0(folder_base, "/", gxe, "_gxe/additive_model/equal_effects/", qtn,
                                    "-QTNs_from_", var, "/", h2, "-heritability/pop", pop)

            trait_pve <- fread(paste0(folder_traits, "/simulated_traits_pve.txt"),
                            header = TRUE, data.table = FALSE)

            pve <- trait_pve[trait_pve$source == "genotype:environment", "pve"]
            pval <- trait_pve[trait_pve$source == "genotype:environment", "pval"]

            gxe_results <- rbind(gxe_results,
                                 data.frame(gxe = gxe,
                                            h2 = h2,
                                            qtn = qtn,
                                            var = var,
                                            ratio = NA,
                                            sv_effect = NA,
                                            diff_dist = FALSE,
                                            pop = pop,
                                            gxe_pve = pve,
                                            gxe_pval = pval))
          }

        }
      }
    }
  }
}

# gxe_pve_summary <- gxe_results %>%
#   group_by(gxe, h2, qtn, var) %>%
#   summarize(gxe_pve_mean = mean(gxe_pve),
#             gxe_pve_se = sd(gxe_pve)/sqrt(n()),
#             gxe_pval_mean = mean(gxe_pval),
#             gxe_pval_se = sd(gxe_pval)/sqrt(n()))

gxe_pve_summary <- gxe_results %>%
  group_by(gxe) %>%
  summarize(gxe_pve_mean = mean(gxe_pve),
            gxe_pve_se = sd(gxe_pve)/sqrt(n()),
            gxe_pval_mean = mean(gxe_pval),
            gxe_pval_se = sd(gxe_pval)/sqrt(n()))
# no gxe mean (se) = 0.04 (+-0.001)
# with gxe mean (se) = 0.18 (+-0.003)

fwrite(gxe_pve_summary, file = paste0(folder_base, "/gxe_pve_summary.txt"),
       sep = "\t", quote = FALSE, na = NA, row.names = FALSE)

plot_summary <- ggplot(gxe_pve_summary, aes(x = gxe, y = gxe_pve_mean)) +
  geom_col(fill = "firebrick") +
  geom_errorbar(aes(ymin = gxe_pve_mean - gxe_pve_se, ymax = gxe_pve_mean + gxe_pve_se), width = 0.2) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "GxE",
       y = "Mean PVE by GxE")

ggsave(plot_summary, filename = paste0(folder_base, "/gxe_pve_summary.pdf"), device = "pdf")
